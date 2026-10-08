#!/usr/bin/env Rscript
#
# Non-split (read-through) counting for FRASER2 driven by a FIXED splice-site
# universe built from STAR SJ.out.tab, using FRASER's own internals + Rsubread.
#
# Replaces the standalone-featureCounts route (sjtab_to_fraser_saf.py +
# CLI featureCounts): FRASER::: does the counting, so there is no SAF window
# math, no -s/-p/--countReadPairs flag archaeology, and no version skew.
#
# Why a FIXED universe: countNonSplicedReads' per-sample HDF5 cache is validated
# by ROW COUNT ALONE (countRNAseqData.R: `length(spliceSiteCoords) == nrow(cache)`).
# Derive the universe from the current cohort's junctions and it changes every
# time a sample is added -> every sample recounts from BAM. Freeze it and the
# cache survives cohort growth.
#
# THE INVARIANT THAT MAKES THIS SAFE
#   Junction<->splice-site joins (calculateIntronNonsplitSum) are on the INTEGER
#   spliceSiteID, which annotateSpliceSite re-derives from scratch on whatever
#   set it is handed. The same coordinate gets a different integer in a
#   different universe, and the join is merge(all.x=TRUE) -> a stale ID does not
#   error, it silently attaches counts to the WRONG splice site.
#   Therefore: the universe is frozen as COORDINATES ONLY. IDs are re-derived
#   every run from the current junction set, and cached counts are matched back
#   in by coordinate. Never persist spliceSiteID/startID/endID.
#
# Sites present in the cohort but absent from the frozen universe get 0, applied
# uniformly to every sample. That biases jaccard UP (denominator shrinks), so it
# can only ever MISS intron retention -- never manufacture it. The asymmetric
# alternative (counting a new site for new samples only, leaving old ones at 0)
# fabricates exactly the outlier signal FRASER hunts for. Do not do it.

suppressPackageStartupMessages({
    library(FRASER)
    library(GenomicRanges)
    library(GenomicAlignments)   # readGAlignments / summarizeJunctions (strand recovery)
    library(Rsamtools)           # ScanBamParam
    library(SummarizedExperiment)
    library(BiocParallel)        # bplapply (parallel per-sample strand recovery)
    library(data.table)
})

STAR_STRAND <- c("*", "+", "-")   # SJ.out.tab col4: 0=undefined, 1=+, 2=-

#' Union of junctions across SJ.out.tab files, filtered like FRASER's
#' minExpressionInOneSample gate.
#'
#' min_count here is the DISCOVERY threshold (fraser_cohort.R --sj-min-count,
#' default 10), which is deliberately lower than the fit's testing threshold
#' (--min-count, 20). It decides what gets counted and can therefore be seen at
#' all; it does not decide what gets a p value.
#'
#' SJ.out.tab cols 2-3 are the 1-based first/last intron base, which is exactly
#' FRASER's junction GRanges -- no offset conversion.
#'
#' count_col="unique" (col7) is the default, and it REQUIRES
#'     scanBamParam(fds) <- ScanBamParam(mapqFilter = 255)
#' STAR gives unique alignments MAPQ 255, so that filter selects exactly STAR's
#' col7 set. Measured on a real sample (chr22, 3.4M PE fragments): FRASER's own
#' BAM counting = 2,255,122 vs STAR col7 = 2,255,120, with 3557/3557 junctions
#' exact (r = 1.000000) among those FRASER actually tests. Bit-identical.
#'
#' Do NOT use "both" (col7+col8) with the default mapqFilter=0. That pairing
#' looks right on paper -- FRASER's default keeps multimappers -- but FRASER's
#' multimapper handling diverges from STAR's col8 (readGAlignmentPairs dumps
#' ambiguously-paired alignments), so it only matches 86% of junctions and blows
#' up at repeat loci (IGL 22:22905026-22922593 measured 213 vs 2028). The
#' aggregate correlation stays at 0.9996 and HIDES this -- don't trust it.
read_sj_union <- function(sj_files, min_count = 20,
                          count_col = c("unique", "both", "multi")) {
    count_col <- match.arg(count_col)
    cols <- c("seqnames", "start", "end", "strand_code", "motif", "annot",
              "n_unique", "n_multi", "overhang")

    per_file <- lapply(sj_files, function(f) {
        if (file.size(f) == 0) {          # STAR produced no junctions (empty sample)
            warning("empty SJ.out.tab, contributes no junctions: ", f)
            return(data.table(seqnames = character(), start = integer(),
                              end = integer(), strand_code = integer()))
        }
        dt <- fread(f, header = FALSE, col.names = cols, showProgress = FALSE)
        dt[, k := switch(count_col,
                         both   = n_unique + n_multi,
                         unique = n_unique,
                         multi  = n_multi)]
        unique(dt[k >= min_count, .(seqnames, start, end, strand_code)])
    })

    jx <- unique(rbindlist(per_file))
    message(date(), ": ", length(sj_files), " SJ.out.tab file(s) -> ",
            nrow(jx), " junctions passing >=", min_count, " (", count_col,
            ") in >=1 sample")

    gr <- GRanges(seqnames = jx$seqnames,
                  ranges   = IRanges(start = jx$start, end = jx$end),
                  strand   = STAR_STRAND[jx$strand_code + 1L])
    sort(gr)
}

#' Recover the +/- strand of STAR col4=0 ("undefined motif") junctions from BAM
#' read orientation, reproducing FRASER's own strand assignment.
#'
#' STAR strands junctions by SPLICE MOTIF, so non-canonical junctions get col4=0
#' -> '*'. FRASER instead strands EVERY junction by READ ORIENTATION. A '*'
#' mixed into a +/- set is a wildcard that breaks annotateSpliceSite's
#' one-hit-per-side invariant -> "mix of stranded and unstranded (+,-,*)" error.
#' This recovers the read-based strand so the set stays uniformly +/-.
#'
#' METHOD (rewritten 2026-09-18; the previous one is described under WHAT THIS
#' REPLACED). For each candidate junction, read the alignments overlapping its
#' DONOR BASE -- the last exonic base before the intron, start-1 -- with
#' readGAlignments(what="flag"); find the reads that actually splice across this
#' exact junction via summarizeJunctions(with.revmap=TRUE); derive each such
#' read's transcript strand from its SAM flag bits. Under strandMode 2 the
#' SECOND mate's own strand IS the transcript strand and the first mate is
#' flipped (0x80 = second in pair, 0x10 = reverse); strandMode 1 is the mirror
#' image. Majority vote, tie -> "+", purity = max(n+,n-)/(n+ + n-).
#'
#' WHAT THIS REPLACED, AND WHY IT HAD TO GO. The old implementation built read
#' PAIRS (readGAlignmentPairs(strandMode=)) over the whole reduced intron body.
#' PAIRING, not read volume, is the superlinear cost: measured, 968k reads
#' finished in 181 s while 497k TIMED OUT at 530-1187 s (both timeouts EEF1A1 --
#' huge depth in a tiny span, with `ambiguous pairing ... dumped` climbing on
#' exactly those). To stop that hanging a genome-wide run, a --max-region-reads
#' gate counted reads over the intron body with countBam and SKIPPED anything
#' above it. But that gate measures EXPRESSION, not strand ambiguity, so it
#' systematically discarded col4=0 junctions in highly expressed genes -- 83.7%
#' of candidates at the 50000 the slurm pinned -- and it skipped them BEFORE
#' reading a single alignment, so they were never tested at all. Downstream,
#' that silence was indistinguishable from "this sample does not carry the
#' junction", which is how junctions carried by all 50 controls (GNAS, MAZ,
#' IFI44L) were reported as private to a case.
#'
#' Removing pair construction removes the cost, so the gate can go. Benchmarked
#' over ALL 853 candidates of a real add (job 30524195): 853/853 resolved, 850
#' at frac >= 0.9, 4.4 min wall on 15 workers, per-junction median 0.91 s,
#' p99 7.70 s, max 9.42 s -- against which the OLD code spent 11.5 min on the
#' countBam gate alone, merely deciding what to skip. It reproduces every
#' junction the paired method resolved (vote_frac agreeing to 4 decimals, same
#' strand call) and GAINS votes, because both mates vote independently instead
#' of one vote per fragment: purity is unaffected (mates agree), only the
#' denominator grows, by up to 2x.
#'
#' NOTHING IS SILENTLY SKIPPED ANY MORE. Every candidate is attempted.
#' `budget_secs` bounds each junction so a pathological locus still cannot hang
#' a run, but a junction over budget is RECORDED as unresolved_timeout instead
#' of vanishing. With `attempts_out` set, every candidate and its outcome is
#' written there -- which is what lets a later `add` tell "tested here and
#' genuinely absent" from "never tested here", the distinction whose absence
#' produced those false private calls.
#'
#' Each junction is resolved from the sample with the most reads over it (the
#' most vote evidence). That choice is free: it comes from SJ.out.tab column 7,
#' not from the BAM. If that sample yields no votes, the next-best untried
#' sample is tried, up to `max_rounds` in total -- the old code was single-shot,
#' so one uncooperative sample lost the junction for the whole cohort. A
#' junction that DOES return votes is never re-asked elsewhere: low purity is a
#' real answer about the reads, and retrying it in another sample would just be
#' shopping for a sample that agrees.
#'
#' @param min_vote_frac  optional purity gate: a junction's vote must be at
#'   least this one-sided or it is treated as unrecovered, same as no read
#'   orientation at all. NULL (default, what the frozen-universe caller below
#'   uses) keeps majority-vote-only behaviour.
#' @param return_metadata  attach vote_frac/carrier_sampleID/carrier_count as
#'   mcols on the returned GRanges. FALSE (default) returns a bare GRanges, so
#'   `c(strand_map, recover_star_strand(...))` below (which combines against an
#'   mcols-less GRanges) is unaffected. The private-novel recovery path sets
#'   this TRUE because it needs to report WHICH sample and read count.
#' @param budget_secs  per-junction elapsed budget. Insurance only -- the
#'   benchmark's slowest junction was 9.42 s against this 60 -- but it is what
#'   guarantees no single locus can hang a run now that nothing is pre-skipped.
#' @param max_rounds  how many samples to try per junction before giving up.
#' @param attempts_out  optional TSV (.gz honoured) recording one row per
#'   candidate junction with its final outcome. `build` writes this into the
#'   cohort bundle; it is the record `add` reads to judge whether a control's
#'   0 at some coordinate means "absent" or "never tested".
#'
#' Returns a GRanges of ONLY the successfully recovered junctions (all +/-) --
#' the strand map to persist and apply to the universe and the split count files.

# Diagnostics computed deep inside the counting path that the CALLER wants to
# record (fraser_cohort.R writes them to runs.tsv). Kept in an environment
# rather than threaded back through three function signatures, because they are
# reporting only -- nothing here changes what gets counted. SCALARS ONLY: the
# per-junction attempt record is a table and goes to `attempts_out`, never here.
.fraser_run_stats <- new.env(parent = emptyenv())

# The strand-recovery method, stamped into run_params.tsv and pinned as a
# CRITICAL parameter. Not a knob -- there is exactly one method and no CLI flag
# to change it. It exists so that a cohort frozen with the pre-2026-09-18 gated
# pair-based recovery cannot be extended by an `add` using this one: that
# asymmetry (case tested thoroughly, controls never tested) is precisely what
# manufactured false private junctions, and it gets WORSE, not better, when only
# one side is upgraded.
STRAND_METHOD <- "donor_flag_v2"

#' Transcript strand of each read from its SAM flag, under FRASER's strandMode.
#' Orphan / single-end reads carry no mate convention, so they do not vote.
tstrand_from_flag <- function(flag, mode) {
    rev    <- bitwAnd(flag, 0x10L) != 0L
    sec    <- bitwAnd(flag, 0x80L) != 0L
    paired <- bitwAnd(flag, 0x1L)  != 0L
    s <- ifelse(sec, ifelse(rev, "-", "+"), ifelse(rev, "+", "-"))
    if (mode == 1L) s <- ifelse(s == "+", "-", "+")
    s[!paired] <- NA_character_
    s
}

recover_star_strand <- function(star_gr, ss, strandMode, mapq = 255L,
                                count_col = c("unique", "both", "multi"),
                                BPPARAM = bpparam(),
                                min_vote_frac = NULL,
                                return_metadata = FALSE,
                                budget_secs = 60,
                                max_rounds = 3L,
                                attempts_out = NULL) {
    count_col <- match.arg(count_col)
    if (length(star_gr) == 0) return(star_gr)
    if (!strandMode %in% c(1L, 2L)) {   # unstranded: no read orientation to use
        warning("recover_star_strand: strandMode=", strandMode,
                " is not stranded (1/2); cannot recover col4=0 strand from reads.")
        return(star_gr[0])
    }
    stopifnot(all(c("sampleID", "bamFile", "sjFile") %in% names(ss)))
    cols <- c("seqnames", "start", "end", "strand_code", "motif", "annot",
              "n_unique", "n_multi", "overhang")

    # Per-junction read count in EVERY sample, from SJ.out.tab (text only, no
    # BAM). The old code kept just the running best; the full matrix is what
    # makes the next-best fallthrough possible, and it is small -- one integer
    # per (candidate junction x sample).
    tk <- paste(seqnames(star_gr), start(star_gr), end(star_gr), sep = ":")
    cnt_mat <- matrix(0L, nrow = length(star_gr), ncol = nrow(ss))
    for (i in seq_len(nrow(ss))) {
        if (file.size(ss$sjFile[i]) == 0) next
        dt <- fread(ss$sjFile[i], header = FALSE, col.names = cols, showProgress = FALSE)
        dt[, cnt := switch(count_col, both = n_unique + n_multi,
                           unique = n_unique, multi = n_multi)]
        dt <- dt[strand_code == 0 & cnt > 0]
        if (!nrow(dt)) next
        m  <- match(tk, paste(dt$seqnames, dt$start, dt$end, sep = ":"))
        hit <- !is.na(m)
        cnt_mat[hit, i] <- dt$cnt[m[hit]]
    }

    # Resolve one sample's assigned junctions. Parallel across samples (each
    # worker reads only its own BAM, disjoint file handles, fork-safe like
    # count_nonsplit); serial per junction inside, which is what keeps memory
    # bounded -- a donor base's reads, not a whole reduced intron region's.
    vote_at <- function(bam, j, d) {
        ga <- readGAlignments(bam, param = ScanBamParam(which = d, mapqFilter = mapq,
                                                        what = "flag"))
        if (!length(ga)) return(list(v = 0L, s = NA_character_, f = NA_real_))
        jc <- summarizeJunctions(ga, with.revmap = TRUE)
        if (!length(jc)) return(list(v = 0L, s = NA_character_, f = NA_real_))
        h <- findOverlaps(j, jc, type = "equal", ignore.strand = TRUE)
        if (!length(h)) return(list(v = 0L, s = NA_character_, f = NA_real_))
        ridx <- mcols(jc)[["revmap"]][[subjectHits(h)[1]]]
        ts   <- tstrand_from_flag(mcols(ga)$flag[ridx], strandMode)
        ts   <- ts[!is.na(ts)]
        np <- sum(ts == "+"); nm <- sum(ts == "-"); nv <- np + nm
        list(v = nv,
             s = if (!nv) NA_character_ else if (np >= nm) "+" else "-",
             f = if (nv) max(np, nm) / nv else NA_real_)
    }

    resolve_chunk <- function(i, idx) {
        bam <- ss$bamFile[i]; sid <- ss$sampleID[i]
        out <- lapply(idx, function(p) {
            j  <- star_gr[p]
            db <- max(start(j) - 1L, 1L)          # last exonic base before the intron
            # as.character() is load-bearing: ScanBamParam(which=) validates
            # EVERY seqlevel of the GRanges against the BAM header, not just the
            # ones it uses. Passing seqnames(j) carries the whole candidate
            # set's seqlevels along (ERCC/SIRV spike-ins, unplaced contigs), so
            # a single contig missing from one sample's header would fail every
            # junction assigned to that sample. Building from a bare character
            # gives a one-seqlevel range.
            d  <- GRanges(as.character(seqnames(j)), IRanges(db, db))
            t0 <- Sys.time(); timed_out <- FALSE; emsg <- ""
            r <- tryCatch({
                setTimeLimit(elapsed = budget_secs, transient = TRUE)
                on.exit(setTimeLimit(), add = TRUE)
                vote_at(bam, j, d)
            }, error = function(e) {
                # Do NOT call every failure a timeout. A BAM whose header lacks
                # this junction's seqlevel, an unindexed file, a truncated read
                # -- all land here, and reporting them as "over budget" hides a
                # real fault behind a knob, which is how the old gate's damage
                # stayed invisible for so long.
                emsg <<- conditionMessage(e)
                timed_out <<- grepl("elapsed time limit|CPU time limit", emsg)
                NULL
            })
            setTimeLimit()
            if (is.null(r)) r <- list(v = 0L, s = NA_character_, f = NA_real_)
            data.table(pos = p, votes = r$v, strand = r$s, vote_frac = r$f,
                       timed_out = timed_out,
                       error = if (timed_out) "" else emsg,
                       sampleID = sid, count = cnt_mat[p, i],
                       secs = as.numeric(Sys.time() - t0, "secs"))
        })
        rbindlist(out)
    }

    has_reads <- rowSums(cnt_mat > 0L) > 0L
    pending   <- which(has_reads)
    tried     <- matrix(FALSE, nrow(cnt_mat), ncol(cnt_mat))
    got       <- vector("list", 0L)
    message(date(), ": recovering col4=0 strand for ", length(pending),
            " junction(s) with reads in >=1 sample [method=", STRAND_METHOD,
            ", budget=", budget_secs, "s/junction, up to ", max_rounds, " sample(s) each]")

    for (round in seq_len(max_rounds)) {
        if (!length(pending)) break
        avail <- cnt_mat
        avail[tried] <- 0L
        sub  <- avail[pending, , drop = FALSE]
        pick <- max.col(sub, ties.method = "first")
        pcnt <- sub[cbind(seq_along(pending), pick)]
        ok   <- pcnt > 0L
        if (!any(ok)) break
        jsel <- pending[ok]; psel <- pick[ok]
        tried[cbind(jsel, psel)] <- TRUE
        groups <- split(jsel, psel)
        if (round > 1L)
            message(date(), ": retry round ", round, " -- ", length(jsel),
                    " junction(s) with no votes yet, across ", length(groups),
                    " further sample(s)")
        res <- rbindlist(bplapply(names(groups), function(gi) {
            i   <- as.integer(gi)
            idx <- groups[[gi]]
            r   <- resolve_chunk(i, idx)
            nerr <- sum(nzchar(r$error))
            message(date(), ": [", ss$sampleID[i], "] ", length(idx),
                    " col4=0 junction(s); ", sum(r$votes > 0L), " resolved, ",
                    sum(r$timed_out), " over budget",
                    if (nerr) paste0(", ", nerr, " errored (", r$error[nzchar(r$error)][1], ")") else "")
            r
        }, BPPARAM = BPPARAM))
        got[[length(got) + 1L]] <- res
        pending <- res$pos[res$votes == 0L]
    }

    res <- if (length(got)) rbindlist(got) else
        data.table(pos = integer(0), votes = integer(0), strand = character(0),
                   vote_frac = numeric(0), timed_out = logical(0),
                   error = character(0), sampleID = character(0),
                   count = integer(0), secs = numeric(0))
    # One row per junction: the attempt that produced votes, else the last tried.
    setorder(res, pos, -votes)
    res_final <- res[, .SD[1L], by = pos]
    n_tries   <- res[, .(n_samples_tried = .N), by = pos]
    res_final <- merge(res_final, n_tries, by = "pos")

    strand_out <- rep(NA_character_, length(star_gr))
    frac_out   <- rep(NA_real_,      length(star_gr))
    sid_out    <- rep(NA_character_, length(star_gr))
    cnt_out    <- rep(NA_integer_,   length(star_gr))
    strand_out[res_final$pos] <- res_final$strand
    frac_out[res_final$pos]   <- res_final$vote_frac
    sid_out[res_final$pos]    <- res_final$sampleID
    cnt_out[res_final$pos]    <- res_final$count

    low_purity <- rep(FALSE, length(star_gr))
    if (!is.null(min_vote_frac)) {
        weak <- !is.na(frac_out) & frac_out < min_vote_frac
        if (any(weak)) {
            message(date(), ": dropping ", sum(weak), " junction(s) with vote ",
                    "purity < ", min_vote_frac, " (ambiguous strand evidence)")
            strand_out[weak] <- NA_character_
            low_purity[weak] <- TRUE
        }
    }

    # The attempt record. Every candidate appears exactly once, including the
    # ones with no reads anywhere and the ones that timed out -- "not in this
    # file" must never again be the only trace of a junction that was skipped.
    if (!is.null(attempts_out)) {
        # not_attempted_no_sj_reads: no sample had a single col4=0 read here, so
        #   there was nothing to read a BAM for. Distinct from unresolved_no_votes,
        #   where reads WERE reported by STAR but none span the junction in the BAM.
        outcome <- rep("not_attempted_no_sj_reads", length(star_gr))
        outcome[res_final$pos] <- ifelse(res_final$votes > 0L, "resolved",
                                  ifelse(res_final$timed_out, "unresolved_timeout",
                                  ifelse(nzchar(res_final$error), "unresolved_error",
                                         "unresolved_no_votes")))
        outcome[low_purity] <- "unresolved_low_purity"
        att <- data.table(
            seqnames = as.character(seqnames(star_gr)),
            start = start(star_gr), end = end(star_gr), width = width(star_gr),
            carrier_sampleID = sid_out, carrier_count = cnt_out,
            votes = NA_integer_, vote_frac = frac_out,
            strand = ifelse(is.na(strand_out), "*", strand_out),
            outcome = outcome, error = "", n_samples_tried = 0L, secs = NA_real_,
            max_sj_count = as.integer(apply(cnt_mat, 1L, max)),
            n_samples_with_reads = as.integer(rowSums(cnt_mat > 0L)))
        att$votes[res_final$pos]           <- res_final$votes
        att$error[res_final$pos]           <- res_final$error
        att$n_samples_tried[res_final$pos] <- res_final$n_samples_tried
        att$secs[res_final$pos]            <- round(res_final$secs, 3)
        # NB vote_frac is reported as MEASURED even where the purity gate
        # rejected it -- the outcome column says it was rejected; blanking the
        # number would hide why.
        fwrite(att, attempts_out, sep = "\t")
        message(date(), ": recorded ", nrow(att), " strand-recovery attempt(s) -> ",
                attempts_out)
    }

    out <- star_gr
    strand(out) <- ifelse(is.na(strand_out), "*", strand_out)
    if (return_metadata) {
        mcols(out)$vote_frac        <- frac_out
        mcols(out)$carrier_sampleID <- sid_out
        mcols(out)$carrier_count    <- cnt_out
    }
    n_rec <- sum(strand_out %in% c("+", "-"))
    .fraser_run_stats$strand_attempted <- length(star_gr)
    .fraser_run_stats$strand_resolved  <- n_rec
    message(date(), ": recovered strand for ", n_rec, "/", length(star_gr),
            " col4=0 junction(s) from BAM read orientation")
    out[as.character(strand(out)) != "*"]     # keep only successfully recovered
}

#' Apply a saved strand map (recovered col4=0 junctions, +/-) to a junction set:
#' override the '*' strand of matching coordinates, then DROP any junction whose
#' strand is still '*' (unrecoverable). Match is by coordinate only.
apply_strand_map <- function(gr, strand_map) {
    is_star <- as.character(strand(gr)) == "*"
    if (any(is_star) && length(strand_map) > 0) {
        h  <- findOverlaps(gr[is_star], strand_map, type = "equal", ignore.strand = TRUE)
        st <- as.character(strand(gr))
        si <- which(is_star)
        st[si[queryHits(h)]] <- as.character(strand(strand_map))[subjectHits(h)]
        strand(gr) <- st
    }
    kept    <- gr[as.character(strand(gr)) != "*"]
    dropped <- length(gr) - length(kept)
    if (dropped > 0)
        message(date(), ": dropped ", dropped,
                " junction(s) with unrecoverable strand (col4=0, no read support)")
    kept
}

#' SJ.out.tab -> per-sample split-count TSVs for getSplitReadCountsForAllSamples(countFiles=).
#'
#' This is the SPLIT side: it replaces FRASER's BAM scan entirely (K, and both
#' N_psi5/N_psi3, are reconstructable from SJ.out.tab; FRASER derives the N's
#' from K itself).
#'
#' DO NOT pre-filter these files. FRASER unions the per-sample files with 0-fill
#' and only then applies minExpressionInOneSample. Dropping a junction here
#' because it is low in THIS sample makes it 0 rather than its true count, which
#' corrupts K for every junction that is low in one sample and high in another.
#' The min_count filter belongs on the universe, not on the split counts.
#'
#' strand_map: recovered strand for col4=0 junctions (see recover_star_strand).
#' getSplitReadCountsForAllSamples runs annotateSpliceSite on the UNION of these
#' files, so a '*' here crashes it exactly like in the universe. col4=0 rows are
#' therefore rewritten to their recovered +/- strand; any col4=0 row NOT in the
#' map is dropped (unrecoverable, or sub-threshold in every sample so it never
#' reaches the universe -- either way it would be filtered, so dropping is safe).
#' With strand_map=NULL the behaviour is the legacy pass-through (col4 strand).
write_split_countfiles <- function(sj_files, sample_ids, outdir,
                                   count_col = c("unique", "both", "multi"),
                                   strand_map = NULL) {
    count_col <- match.arg(count_col)
    stopifnot(length(sj_files) == length(sample_ids))
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    cols <- c("seqnames", "start", "end", "strand_code", "motif", "annot",
              "n_unique", "n_multi", "overhang")

    mapdt <- NULL
    if (!is.null(strand_map) && length(strand_map) > 0) {
        mapdt <- data.table(seqnames = as.character(seqnames(strand_map)),
                            start = start(strand_map), end = end(strand_map),
                            map_strand = as.character(strand(strand_map)))
    }

    vapply(seq_along(sj_files), function(i) {
        f <- file.path(outdir, paste0(sample_ids[i], ".splitCounts.tsv"))
        if (file.size(sj_files[i]) == 0) {          # empty sample -> empty countfile
            fwrite(data.table(seqnames = character(), start = integer(),
                              end = integer(), strand = character(),
                              count = integer()), f, sep = "\t")
            return(f)
        }
        dt <- fread(sj_files[i], header = FALSE, col.names = cols, showProgress = FALSE)
        dt[, count := switch(count_col,
                             both   = n_unique + n_multi,
                             unique = n_unique,
                             multi  = n_multi)]
        dt <- dt[count > 0]
        dt[, seqnames := as.character(seqnames)]   # fread parses "22" as integer
        dt[, strand := STAR_STRAND[strand_code + 1L]]
        if (!is.null(mapdt)) {
            # canonical rows pass through; col4=0 rows take the recovered strand,
            # and col4=0 rows absent from the map are dropped.
            nonstar <- dt[strand_code != 0, .(seqnames, start, end, strand, count)]
            star    <- dt[strand_code == 0, .(seqnames, start, end, count)]
            if (nrow(star)) {
                star <- merge(star, mapdt, by = c("seqnames", "start", "end"),
                              all.x = TRUE)
                star <- star[!is.na(map_strand),
                             .(seqnames, start, end, strand = map_strand, count)]
            } else {
                star <- nonstar[0]
            }
            dt <- rbind(nonstar, star)
        } else {
            dt <- dt[, .(seqnames, start, end, strand, count)]
        }
        fwrite(dt, f, sep = "\t")
        f
    }, character(1))
}

#' Freeze the splice-site universe. Coordinates only -- mcols are dropped on
#' purpose so a stale spliceSiteID can never be resurrected from this file.
build_frozen_universe <- function(jx, path) {
    universe <- FRASER:::extractSpliceSiteCoordinates(
        FRASER:::annotateSpliceSite(jx))
    universe <- granges(sort(universe))   # sort() is a no-op here; assert intent
    saveRDS(universe, path)
    message(date(), ": froze ", length(universe), " splice sites -> ", path)
    universe
}

load_frozen_universe <- function(path) {
    universe <- readRDS(path)
    stopifnot(is(universe, "GRanges"), is.null(universe$spliceSiteID))
    # SORT ORDER IS LOAD-BEARING, and nothing downstream would notice it broken.
    # countNonSplicedReads assigns featureCounts' output to spliceSiteCoords in
    # the order it was handed, then does `spliceSiteCoords <- sort(...)` before
    # caching -- so the cached vector is in SORTED order of whatever we passed.
    # assemble_nonsplit_se re-keys it by position against this object exactly as
    # loaded. Both are sorted today (extractSpliceSiteCoordinates returns
    # unique(sort(.)) and build_frozen_universe re-sorts), so this holds -- but
    # an unsorted universe.RDS of the right length would silently permute every
    # non-split count onto the wrong splice site, with no error anywhere. That
    # is the same silent-corruption class as the spliceSiteID trap above, so it
    # gets the same treatment: assert, do not assume.
    if (is.unsorted(universe))
        stop("frozen universe is not sorted: ", path,
             "\n  Non-split counts are matched to it BY POSITION against a cache ",
             "FRASER writes in sorted order, so an unsorted universe silently ",
             "assigns every count to the wrong splice site. Rebuild it with ",
             "build_frozen_universe(), which sorts.")
    universe
}

#' Count non-split reads per sample against the frozen universe.
#'
#' countNonSplicedReads only reads COORDINATES off spliceSiteCoords (to build the
#' SAF it hands Rsubread), so passing an ID-free universe is fine. It writes
#' nonSplicedCounts-<id>.h5 and returns counts in sorted-coordinate order.
#'
#' recount=FALSE + a frozen universe => a sample already counted is never
#' re-read from BAM, no matter how much the cohort grows.
count_nonsplit <- function(fds, universe, minAnchor = 5, recount = FALSE,
                           NcpuPerSample = 1, BPPARAM = bpparam()) {
    counts <- bplapply(samples(fds), FUN = FRASER:::countNonSplicedReads,
                       fds              = fds,
                       splitCountRanges = NULL,   # ignored; universe is supplied
                       spliceSiteCoords = universe,
                       minAnchor        = minAnchor,
                       recount          = recount,
                       NcpuPerSample    = NcpuPerSample,
                       BPPARAM          = BPPARAM)
    names(counts) <- samples(fds)
    counts
}

#' Assemble the rawCountsSS SummarizedExperiment for THIS run.
#'
#' splitCounts must be the SE from getSplitReadCountsForAllSamples(), whose
#' rowRanges were annotated in ONE shot -- its startID/endID are the integers the
#' downstream join uses. The run universe is derived from a SUBSET of that same
#' annotated object, so it inherits those IDs rather than re-deriving (and
#' shifting) them.
#'
#' The filter here mirrors countRNAData() (filter=TRUE); FRASER builds its
#' splice-site table from the FILTERED junctions, not from every junction in
#' splitCounts.
#'
#' NOTE the value fraser_cohort.R passes for minExpressionInOneSample here is
#' --sj-min-count (the DISCOVERY threshold that froze the universe), NOT
#' --min-count (FRASER's testing threshold). The two were one number until
#' 2026-09-01. Pass the testing threshold here instead and every junction
#' between the two thresholds keeps its split counts but loses its splice sites:
#' still in the junction set, still in the private-novel candidates, non-split
#' term 0. Whatever value froze the universe is the value that belongs here.
assemble_nonsplit_se <- function(fds, splitCounts, universe, countList,
                                 filter = TRUE, minExpressionInOneSample = 20) {
    splitCountRanges <- FRASER:::extractSplitCountRanges(
        splitCounts = splitCounts, filter = filter,
        minExpressionInOneSample = minExpressionInOneSample)
    run_ss <- FRASER:::extractSpliceSiteCoordinates(splitCountRanges)

    # Cached rows are in frozen-universe order; re-key them onto this run's
    # sites by COORDINATE. Equality of (seqnames, start, end, strand).
    hits <- findOverlaps(run_ss, universe, type = "equal")
    stopifnot(!any(duplicated(queryHits(hits))))

    cached <- do.call(cbind, lapply(countList, as.matrix))
    # Positional re-keying again: row i of `cached` IS site i of `universe`.
    # FRASER validates its per-sample cache by row count alone, so a cache built
    # against a different universe of the SAME length is accepted and lands here
    # silently. A length mismatch is the one part of that we can still catch.
    if (nrow(cached) != length(universe))
        stop("non-split cache has ", nrow(cached), " row(s) but the frozen ",
             "universe has ", length(universe), " site(s). The cached counts ",
             "were produced against a DIFFERENT universe; delete ",
             "cache/nonSplicedCounts/ and recount, or restore the universe.RDS ",
             "they were made with.")
    counts <- matrix(0L, nrow = length(run_ss), ncol = length(countList),
                     dimnames = list(NULL, names(countList)))
    counts[queryHits(hits), ] <- cached[subjectHits(hits), , drop = FALSE]

    n_missing <- length(run_ss) - length(hits)
    .fraser_run_stats$sites_zeroed <- n_missing
    .fraser_run_stats$sites_total  <- length(run_ss)
    if (n_missing > 0) {
        warning(n_missing, " of ", length(run_ss), " splice sites are absent ",
                "from the frozen universe and are set to 0 for ALL samples ",
                "(uniform -> biases jaccard up, cannot create false outliers). ",
                "Rebuild + recount the universe if this grows large.")
    }

    h5 <- FRASER:::saveAsHDF5(fds, "rawCountsSS", counts[, samples(fds), drop = FALSE])
    colnames(h5) <- samples(fds)
    SummarizedExperiment(colData   = colData(fds),
                         assays    = list(rawCountsSS = h5),
                         rowRanges = run_ss)
}

#' Full path: split counts (from your SJ.out.tab-derived countFiles) + non-split
#' counts (frozen universe) -> a countRNAData()-equivalent fds, zero BAM reads
#' for any sample already cached.
#'
#' countRNAData() does NOT forward countFiles, which is why this is assembled by
#' hand rather than by calling it.
build_fds_counts <- function(fds, countFiles, universe_path,
                             minAnchor = 5, recount = FALSE,
                             minExpressionInOneSample = 20,
                             count_col = "unique",
                             NcpuPerSample = 1, BPPARAM = bpparam()) {
    # The split counts you injected and the non-split counts FRASER is about to
    # take from the BAM must agree on multimappers, or the jaccard numerator and
    # denominator are drawn from different read sets.
    mapq <- bamMapqFilter(scanBamParam(fds))
    if (identical(count_col, "unique") && !identical(as.integer(mapq), 255L)) {
        stop("count_col='unique' (SJ.out.tab col7) requires ",
             "scanBamParam(fds) <- ScanBamParam(mapqFilter = 255), but ",
             "mapqFilter is ", mapq, ". STAR gives unique alignments MAPQ 255; ",
             "without this filter FRASER counts non-split reads from the ",
             "multimapper-inclusive BAM while your split counts are unique-only.")
    }
    if (!identical(count_col, "unique") && identical(as.integer(mapq), 255L)) {
        warning("count_col='", count_col, "' includes multimappers but ",
                "mapqFilter=255 excludes them on the non-split side.")
    }

    splitCounts <- getSplitReadCountsForAllSamples(
        fds = fds, countFiles = countFiles, recount = recount, BPPARAM = BPPARAM)

    universe <- load_frozen_universe(universe_path)
    countList <- count_nonsplit(fds, universe, minAnchor = minAnchor,
                                recount = recount, NcpuPerSample = NcpuPerSample,
                                BPPARAM = BPPARAM)
    nonSplitCounts <- assemble_nonsplit_se(
        fds, splitCounts, universe, countList,
        minExpressionInOneSample = minExpressionInOneSample)

    # validObject() runs validateNonSplicedReadsSanity here. Do NOT read more
    # into it than it does: it is
    #   nrow(nonSplicedReads) == length(intersect(spliceSiteID, c(startID, endID)))
    # i.e. a set-MEMBERSHIP count. run_ss derives its IDs fresh from the current
    # junction set, so they are always valid members and this always passes --
    # including when the count VALUES have been joined to the wrong rows, which
    # is the actual risk in assemble_nonsplit_se. It catches a merged/mismatched
    # ID space; it cannot catch a permuted universe. The sortedness and row-count
    # assertions above are what guard that, and they are not redundant with this.
    addCountsToFraserDataSet(fds, splitCounts, nonSplitCounts)
}
