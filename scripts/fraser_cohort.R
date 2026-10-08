#!/usr/bin/env Rscript
#
# FRASER2 cohort runner on STAR SJ.out.tab + a frozen splice-site universe.
#
#   build : first run. Freezes the splice-site universe, counts every sample.
#   add   : add sample(s). Reuses the frozen universe -> already-counted samples
#           are NEVER re-read from BAM. One featureCounts pass for the new sample.
#   fit   : PSI -> filter -> q -> fit -> results. No counting, no BAM access.
#
# Usage
#   Rscript fraser_cohort.R build --samplesheet cohort.tsv --workdir WD
#   Rscript fraser_cohort.R add   --samplesheet cohort.tsv --workdir WD
#   Rscript fraser_cohort.R fit   --workdir WD
#
#   `add` will also take one sample inline, in place of a one-line samplesheet:
#   Rscript fraser_cohort.R add --sample-id ID --bam F --sj F --workdir WD
#
#   `fit` takes NO sample input -- it reads the counts object build/add wrote.
#   --samplesheet is accepted there for backwards compatibility and ignored.
#
# The samplesheet for `add` is the FULL cohort (old rows + the new one), not just
# the new sample -- FRASER fits across the whole cohort. `add` is `build` minus
# the universe rebuild; that single difference is what preserves the caches.
#
# ---------------------------------------------------------------------------
# FROZEN CONTROL COHORT: --cohort-dir (read-only) + --workdir (writable)
# ---------------------------------------------------------------------------
#   Rscript fraser_cohort.R add --cohort-dir COHORT --workdir WD \
#                               --sample-id NEWID --bam F --sj F
#   Rscript fraser_cohort.R fit --cohort-dir COHORT --workdir WD \
#                               --results-samples NEWID --out-dir fraser_out
#
#   (--samplesheet just_the_new_sample.tsv still works in place of the three
#   inline flags on `add`, and is simply ignored on `fit`.)
#
# For the one-case-against-a-fixed-panel workflow (one Snakemake run per sample,
# many samples in flight at once). COHORT is NEVER written to, so any number of
# samples can run against it concurrently, and every case is fitted against a
# cohort that does not move.
#
# What --cohort-dir changes:
#   * universe.RDS, run_params.tsv, strand_map.RDS and splitCountFiles/ are READ
#     from COHORT. Everything produced is written under WD.
#   * THE SAMPLESHEET IS THE NEW SAMPLE(S) ONLY. The roster of the frozen cohort
#     comes off its own saved counts object, which is the only authoritative
#     statement of what the cohort contains -- restating 20+ frozen paths in
#     every per-sample samplesheet is how they drift apart. Rows naming a sample
#     already in the cohort are dropped with a message, so passing a full sheet
#     still works.
#   * Only the new sample's SJ.out.tab is read and only its split-count file is
#     written; the frozen ones are used where they lie. Without this the `add`
#     re-read every control's SJ.out.tab and rewrote all of their count files --
#     ~12 of the 13 minutes of the 2026-08-28 run, all of it discarded.
#   * The frozen per-sample non-split caches are SYMLINKED into WD (see
#     link_cohort_cache): FRASER addresses those caches through workingDir(fds)
#     alone, so they have to be reachable under the writable workdir. Links are
#     bytes; the new sample's h5 is a real file created beside them, in WD.
#   * strand_map.RDS is FROZEN, not extended. Recovering a new col4=0 strand for
#     the case alone would hand it real reads at a junction where every control
#     keeps the 0 the build-time drop gave it -- a control-carried junction then
#     reports as private to the case. --max-star-span is inert on this path as a
#     result. See the FREEZE_STRANDS block.
#
# Without --cohort-dir every path resolves to WD exactly as before.
#
# SAVED OBJECTS
#   savedObjects/<name>       counts only, written by build/add. The fit never
#                             writes here, so every fit starts from the same
#                             full junction set.
#   savedObjects/<name>_fit   filtered + annotated + fitted, written by fit.
#   Keeping these apart is what makes `fit` re-runnable. When fit saved over
#   savedObjects/<name> it consumed its own input: the next fit loaded the
#   already-filtered set, and the private-novel table -- built from the
#   junctions the filter rejects -- silently came out empty.
#
# OUTPUTS (fit). With --results-samples ID these move to samples/<ID>/<ID>_*
# and the cohort-wide copies are not written at all; without it they are
# <workdir>/<name>_* exactly as before. --out-dir DIR redirects every table to
# DIR with NO per-sample subdirectory -- the file names already carry the sample
# ID, so one flat directory holds the deliverables for any number of samples.
#   <name>_results_all.tsv.gz  every junction  (all=TRUE)   -> parse_fraser.R
#   <name>_results.tsv         outliers only   (all=FALSE)
#   <name>_results_gene.tsv    outliers aggregated per gene, carrying
#                              pValueGene/padjustGene. Gene-level stats are
#                              computed inside FRASER() and ONLY exist if the
#                              annotation ran first -- see the fit block.
#   <name>_results_gene_all.tsv.gz   the all=TRUE gene table. OPT-IN
#                              (--gene-all-table): heaviest call in the fit,
#                              OOM-killed job 30219333 at 64G.
#   <name>_private_novel.tsv   junctions FRASER COULD NOT TEST that are carried
#                              by exactly one sample -- no p-value, by design.
#                              See write_private_novel_table() below for why this
#                              class is invisible in the results tables and why it
#                              matters most in a proband-vs-frozen-cohort `add`.
#                              Disable with --no-private-table. With
#                              --results-samples a per-sample copy is also
#                              written to samples/<ID>/ -- a filter of this
#                              table, since a private junction has one carrier.
#                              Under --out-dir this cohort-wide table is named
#                              <ID>_<name>_private_novel.tsv: many samples share
#                              one --out-dir and its content differs per case, so
#                              a bare <name>_ prefix had them overwrite each
#                              other. See cohort_prefix().
#   runs.tsv                   append-only record of every build/add/fit in this
#                              workdir: sample count, junctions in/passed, q,
#                              filter quantile, counting parameters, universe
#                              drift. Every add REFITS the whole cohort, so a
#                              sample's p values depend on which cohort it was
#                              called against; this is what ties a released
#                              report back to the run that produced it. A copy
#                              lands in each samples/<ID>/ as <ID>_fit_info.tsv.
#
# PARAMETER PINNING
#   `build` records count_col/mapq/strand/min_count/sj_min_count (+ the two
#   col4=0 limits) in
#   <workdir>/run_params.tsv; `add` refuses to run if the critical four differ.
#   The non-split cache is validated by ROW COUNT ALONE, so without this an
#   `add --mapq 0` would silently mix read sets across samples. See
#   check_run_params() for the full argument.
#
# Samplesheet: TSV, no header, positional columns sampleID, bamFile, sjFile
#              (data assumed paired-end; a header row is tolerated if present)
#              Required by `build`. On `add`, interchangeable with the inline
#              --sample-id/--bam/--sj trio. On `fit`, not needed at all.
#
# WHY THIS IS FAST
#   FRASER's counting is 2 full BAM passes/sample: split reads + non-split reads.
#   - Split: STAR already computed it. SJ.out.tab col7 == FRASER's own BAM counts
#     EXACTLY, given mapqFilter=255 (measured: 3557/3557 junctions, r=1.000000).
#     -> that pass drops to zero.
#   - Non-split: no STAR equivalent (STAR never reports reads that overlap a
#     splice site without splicing), so 1 BAM pass is the floor. But FRASER
#     derives its splice-site universe from the FILTERED COHORT junctions, which
#     grow every time you add a sample; the per-sample cache is validated by row
#     count alone, so a growing universe silently invalidates EVERY sample's
#     cache. Freezing the universe is what stops that.
#   Adding 1 sample to 100: ~101 BAM passes -> 1.
#
# TWO THRESHOLDS, AND THEY ARE NOT THE SAME NUMBER
#   --sj-min-count (default 10) DISCOVERY. A junction enters the frozen
#       splice-site universe -- and contributes its two splice sites to the
#       per-run site set -- if it carries >= this many reads in >= 1 sample.
#       Everything at or above it is COUNTED, appears in K, gets a real
#       non-split denominator, and is eligible for the private-novel table.
#   --min-count (default 20)    TESTING. FRASER's minExpressionInOneSample, the
#       fit's expression gate. This is what decides which junctions get a
#       p-value, so it also fixes the multiple-testing burden. Left at DROP's
#       20 on purpose.
#   These were ONE flag until 2026-09-01, which meant "see more junctions" and
#   "test more junctions" could not be asked for separately. They can now.
#   The ordering constraint is --sj-min-count <= --min-count: the per-run site
#   set is re-derived from --sj-min-count on every run
#   (assemble_nonsplit_se -> extractSplitCountRanges), so a junction above the
#   testing threshold but below the discovery one would be tested with a
#   non-split denominator of 0. The script warns if you invert them.
#   --sj-min-count is PINNED (it defines the universe); --min-count is pinned
#   too, because it is what a released p-value was computed under.
#
#   WHY 10 AND NOT 5 OR 20. Measured on the 22-sample LCL control cohort
#   (2026-09-01, leave-one-out over 6 samples treated as probands). Lowering the
#   discovery threshold does NOT change any junction already >= 20 in a control:
#   calculateIntronNonsplitSum reads a junction's OWN two sites, and those are in
#   the 20-universe by construction. The entire gain is one class -- a junction
#   carried WEAKLY by the controls (below the threshold) and STRONGLY by the
#   case. At 20 its sites never entered the frozen universe, so its non-split
#   term is 0 for everyone, its jaccard is biased UP, and real intron retention
#   there is invisible (a false negative; the bias is uniform, so it can never
#   fabricate an outlier).
#     threshold   splice sites   col4=0 to recover   median junctions rescued
#        20         357,918            2,340                  --
#        10         453,537 (1.27x)    5,618  (2.4x)         491
#         5         628,586 (1.76x)   14,783  (6.3x)         612
#         3         873,344 (2.44x)   33,319 (14.2x)         660
#   Marginal cost per extra rescued junction: 195 splice sites for 20->10, 1,447
#   for 10->5, 5,141 for 5->3. 10 buys ~84% of the available correction (77-90%
#   across the 6 cases) for ~27% of the incremental cost, so that is the default.
#   Splice sites are the cost driver: they become featureCounts SAF rows, so
#   per-sample non-split counting time and .h5 cache size both scale with them.
#   The residual is irreducible by this knob -- even at 3, every case still has
#   200-800 testable junctions with a site off-universe, because the controls
#   carry ZERO reads there. That class is what <name>_private_novel.tsv is for.
#
# MAPQ IS NOT OPTIONAL
#   --mapq 255 keeps both sides unique-only, matching STAR col7. With the default
#   mapqFilter=0 you would inject unique-only split counts while FRASER counts
#   non-split reads from the multimapper-inclusive BAM -- numerator and
#   denominator drawn from different read sets. build_fds_counts() hard-errors.

suppressPackageStartupMessages({
    library(FRASER); library(GenomicRanges); library(SummarizedExperiment)
    library(data.table); library(Rsamtools); library(BiocParallel); library(HDF5Array)
})

self <- tryCatch(normalizePath(sub("^--file=", "",
            grep("^--file=", commandArgs(FALSE), value = TRUE)[1])),
         error = function(e) NA)
source(file.path(dirname(self), "fraser_nonsplit_from_sjtab.R"))
source(file.path(dirname(self), "fraser_gene_annotation.R"))

# ---------------------------------------------------------------- args -------
args <- commandArgs(TRUE)
mode <- args[1]
if (!mode %in% c("build", "add", "fit")) {
    stop("usage: fraser_cohort.R <build|add|fit> --samplesheet F --workdir D ",
         "[--name N] [--strand 0|1|2] [--mapq 255] [--threads T] ",
         "[--min-count 20] [--sj-min-count 10] [--count-col unique|both] ",
         "[--delta-jaccard 0.1] ",
         "[--ncpu-per-sample auto|N] [--max-star-span 100000] [--strand-budget-secs 60] ",
         "[--gtf F | --txdb F] [--gene-map F] [--rebuild-txdb]   # fit: annotate genes (default txdb auto-found)\n",
         "  add only: --sample-id ID --bam F --sj F   one sample, in place of ",
         "--samplesheet (all three required together)\n",
         "  fit only: --samplesheet is OPTIONAL (the fit reads no BAM and no ",
         "SJ.out.tab)\n",
         "  add/fit: [--cohort-dir D]  frozen, read-only control cohort; the ",
         "samplesheet is then the NEW sample(s) only\n",
         "  fit also: [--filter-quantile 0.95] [--quantile-min-expression 10] ",
         "[--private-min-count N] ",
         "[--private-max-other 2] [--no-private-table] [--gene-table] [--gene-all-table]\n",
         "  add also (against --cohort-dir): [--strand-vote-min-frac 0.9]  purity ",
         "gate for recovering the new sample's private col4=0 junctions into ",
         "the private-novel table (never into strand_map/universe)\n",
         "            [--private-into-counts]  also fold those junctions into ",
         "THIS sample's split counts, so they get a real jaccard and count ",
         "toward neighbouring junctions' N (still never tested: one carrier ",
         "fails the expression filter)\n",
         "            [--results-samples ID[,ID]]  report only these samples ",
         "-> <workdir>/samples/<ID>/\n",
         "            [--out-dir D]                write every results table to ",
         "D, flat (no samples/<ID>/ subdirectory)\n",
         "            [--allow-filtered-input]     fit an object an older fit ",
         "already filtered (private-novel table will be empty)\n",
         "            [--q, --q-search]  --q-search writes ",
         "<prefix>_q_search.tsv/.png (the injected-outlier AROC-vs-q grid, ",
         "i.e. estimateBestQ(useOHT=FALSE) -- what used to be ",
         "optimHyperParams()) and, when --q was not given, ",
         "<prefix>_q_oht_spectrum.png (the OHT singular-value plot). When --q ",
         "was not given, the grid's best-AROC q REPLACES the OHT estimate as ",
         "the q the fit uses; an explicit --q always overrides both. Slower: ",
         "up to ~15 extra PCA fits over the candidate q range.")
}
# LAST occurrence wins. match() returned the FIRST, which silently discarded
# any flag a wrapper had already hardcoded: run_fraser_cohort.slurm writes its
# own --min-count/--sj-min-count BEFORE "${EXTRA[@]}", so FRASER_EXTRA="--sj-min-count 10"
# parsed cleanly, was echoed in the log, and then had no effect at all.
opt <- function(flag, default = NULL) {
    i <- which(args == flag); if (!length(i)) return(default)
    args[max(i) + 1L]
}
SHEET   <- opt("--samplesheet"); WD <- opt("--workdir")
# One sample given inline, in place of a one-line samplesheet (`add` only). The
# per-sample workflow builds a tsv holding exactly one row whose fields the
# caller already has in hand; a Snakemake rule has them as {wildcards.sampleID},
# {input.bam}, {input.sj}. See the SAMPLE INPUT block below for the rules.
SID     <- opt("--sample-id", NULL)
BAM     <- opt("--bam", NULL)
SJ      <- opt("--sj", NULL)
NAME    <- opt("--name", "cohort")
STRAND  <- as.integer(opt("--strand", "2"))     # 2 = reverse (dUTP). SET THIS.
MAPQ    <- as.integer(opt("--mapq", "255"))
THREADS <- as.integer(opt("--threads", "4"))
MINCNT  <- as.integer(opt("--min-count", "20"))    # TESTING threshold (FRASER)
# DISCOVERY threshold: which junctions are counted at all. See the header block
# "TWO THRESHOLDS". Feeds read_sj_union (what freezes the universe) and
# assemble_nonsplit_se's extractSplitCountRanges (the per-run site set) -- those
# two MUST use the same number or the sites in the gap are zeroed for everyone.
# It does NOT feed filterExpressionAndVariability: --min-count still does.
SJMIN   <- as.integer(opt("--sj-min-count", "10"))
if (SJMIN > MINCNT)
    warning("--sj-min-count (", SJMIN, ") is above --min-count (", MINCNT,
            "). Junctions between the two are TESTED but their splice sites ",
            "never entered the universe, so their non-split denominator is 0 ",
            "and their jaccard is biased up. Set --sj-min-count <= --min-count.",
            immediate. = TRUE)
CCOL    <- opt("--count-col", "unique")
DPSI    <- as.numeric(opt("--delta-jaccard", "0.1"))
# fit filter: rowQuantiles(N_jaccard, FQUANT) >= 10. We default to 0.95 to MATCH
# DROP (quantileForFiltering, AberrantSplicing.py:47), not to FRASER2's own 0.75
# default -- higher is MORE PERMISSIVE, so 0.75 was discarding junctions a DROP
# run keeps. Measured for n=40: 0.95 needs 3/40 samples with N>=10, 0.75 needs
# 11/40. Pass --filter-quantile 0.75 for stock FRASER2 behaviour.
FQUANT  <- as.numeric(opt("--filter-quantile", "0.95"))
# The OTHER half of the expression filter: rowQuantiles(N_jaccard, FQUANT) >=
# quantileMinExpression. FRASER hardcodes 10 and we never passed it; this
# exposes it. 0 disables that half entirely (quantileValueN >= 0 is always
# true), leaving the gate as rowMaxs(K) >= --min-count alone.
#
# ⚠ WHAT THIS GUARD ACTUALLY DOES. The quantile term tests the junction's
# SPLICE SITES across the cohort, not the junction. A junction carried by ONE
# sample, whose sites no other sample covers, has N=0 everywhere else, so the
# quantile lands among the zeros and it fails at EVERY FQUANT -- even at
# K=1e6. That is deliberate: such a junction has no cohort distribution to be
# an outlier AGAINST, and the private-novel table exists to report it instead.
# Setting this to 0 lets that whole class into the fit, where it WILL tend to
# look significant, because the model sees ~n-1 zeros and one large value.
# Some of those zeros are real absence; others are structural (a site never
# counted, a col4=0 junction whose strand was never resolved). The model cannot
# tell them apart -- which is the asymmetry fraser_nonsplit_from_sjtab.R's
# header forbids. Use for experiments, read the results knowing this, and do
# not make it the default.
QMINEXPR <- as.numeric(opt("--quantile-min-expression", "10"))
QOPT    <- opt("--q", NULL)                      # override estimateBestQ
# Runs estimateBestQ(useOHT=FALSE) -- the injected-outlier grid search
# FRASER2 otherwise skips (see the usage block). When --q was not given, its
# winner (max AROC) REPLACES the fast OHT estimate as the q the fit uses.
QSEARCH <- "--q-search" %in% args
NCPS    <- opt("--ncpu-per-sample", "auto")      # featureCounts nthreads/sample
MAXSPAN <- as.integer(opt("--max-star-span", "100000")) # drop col4=0 junctions wider than this bp (0=off)
# Per-junction elapsed budget for col4=0 strand recovery. This REPLACES
# --max-region-reads, which was removed 2026-09-18 along with the pair-based
# recovery whose cost it existed to contain. The old flag pre-skipped any
# junction whose intron body carried more than N reads -- a measure of
# EXPRESSION, not of strand ambiguity -- and so discarded 83.7% of candidates
# without reading one alignment, which is how junctions carried by every
# control were reported as private to a case. Nothing is pre-skipped now; this
# budget is a hang guard only (the benchmark's slowest junction: 9.42 s), and a
# junction that exceeds it is RECORDED as unresolved_timeout, not silently
# dropped. See recover_star_strand() in fraser_nonsplit_from_sjtab.R.
STRANDBUDGET <- as.numeric(opt("--strand-budget-secs", "60"))
GTF     <- opt("--gtf", NULL)                    # fit: gene-annotate before results()
# --txdb points at a txdb.db built ONCE by build_fraser_txdb.R (beside the GTF).
# Prefer it: makeTxDbFromGFF on a full Ensembl GTF is slow and memory-hungry and
# nothing about it depends on the cohort. With --gtf alone the cache is still
# looked up (and written) beside the GTF, so the second fit is just as fast.
TXDB    <- opt("--txdb", NULL)                   # prebuilt txdb.db; no --gtf needed
GMAP    <- opt("--gene-map", NULL)               # its gene-name mapping tsv (default: sibling)
REBUILD <- "--rebuild-txdb" %in% args            # force a txdb/mapping rebuild

# Default annotation, so the ordinary fit needs NO annotation flag at all. The
# txdb + gene-name mapping are built once per annotation release
# (build_fraser_txdb.R) and shared by every workdir, so there is one obvious
# location -- same convention as parse_fraser.R's blacklist/exon tables: a
# hardcoded default, an env override, and a warning (never a hard error) when
# nothing is found. Not having this is how a fit ends up rebuilding a 266 MB
# txdb into its own workdir, or silently producing results with no hgncSymbol.
DEFAULT_TXDB <- paste0("/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/",
                       "annotation_files/Homo_sapiens.GRCh38.115.txdb.db")

# The mapping is NOT optional: annotate_fds_genes() needs it for the
# gene_id -> gene_name lookup, so a txdb without its sibling mapping is not a
# usable candidate. Mirrors build_txdb_and_map()'s own derivation, including the
# fallback for a DROP-style plain "txdb.db".
map_for_txdb <- function(f) {
    m <- sub("\\.txdb\\.db$", ".gene_name_mapping.tsv", f)
    if (identical(m, f)) m <- file.path(dirname(f), "gene_name_mapping.tsv")
    m
}

#' Locate a prebuilt txdb when neither --txdb nor --gtf was given.
#' $FRASER_TXDB > the shared annotation_files cache > a txdb.db an older run
#' left in this workdir. Only a candidate whose mapping is also present counts.
find_txdb <- function(wd) {
    cands <- character(0)
    if (nzchar(Sys.getenv("FRASER_TXDB"))) cands <- c(cands, Sys.getenv("FRASER_TXDB"))
    cands <- path.expand(c(cands, DEFAULT_TXDB, file.path(wd, "txdb.db")))
    ok <- cands[file.exists(cands) & file.exists(vapply(cands, map_for_txdb, ""))]
    if (length(ok)) ok[1] else NA_character_
}
PRIVMIN <- as.integer(opt("--private-min-count", MINCNT))  # fit: private-novel reporting
# max reads allowed in any OTHER sample before a junction stops counting as
# private. At the old default of 0 a single stray read anywhere in the cohort
# erased a junction from the table entirely -- and because 0 also forces
# max_other_sample_count and n_other_samples_nonzero to be 0 on every surviving
# row, the two columns that exist to let you judge how private a junction really
# is were constants. 2 keeps the table a triage list rather than a pre-filtered
# verdict; the columns now carry the evidence and you sort on them.
PRIVOTH <- as.integer(opt("--private-max-other", "2"))
PRIVTAB <- !("--no-private-table" %in% args)
# fit: also write the aggregate/collapsed GENE table (_results_gene.tsv). OFF by
# default -- aggregate=TRUE still has to join every junction to its gene and
# collapse across the whole cohort to get correct Holm/BY-corrected gene-level
# p-values BEFORE any deltaPsiCutoff filtering is applied, so it's expensive
# regardless of how few rows survive the cutoff. This OOM-killed jobs 30502778
# and 30504802/3/4 (Blood/Fibroblasts/LCL, Sept 2026) even at all=FALSE. Give
# this step more memory, then rerun results-extraction only with --gene-table.
GENETAB <- "--gene-table" %in% args
# add (against --cohort-dir): purity gate for recovering a new sample's own
# col4=0 junctions into the private-novel table (see FREEZE_STRANDS below).
# 0.9 means the winning strand needs >=90% of the supporting read-pairs --
# below that the evidence is too ambiguous to call it resolved, same as no
# orientation info at all.
STRANDVOTE <- as.numeric(opt("--strand-vote-min-frac", "0.9"))
# Fold the new sample's recovered private col4=0 junctions into ITS OWN split
# count file, so FRASER computes a real jaccard for them and their reads enter
# N_psi5/N_psi3 for neighbouring junctions that ARE tested. Off by default: it
# changes the denominator of tested junctions in the case and not in the
# controls, which is only correct if the controls genuinely lack those reads.
# See the --max-star-span gate in the FREEZE_STRANDS block, which is what makes
# that true. The frozen strand_map/universe/control caches stay untouched either
# way -- only this sample's countfile changes.
PRIVCOUNTS <- "--private-into-counts" %in% args
# fit: also write the all=TRUE GENE table. OFF by default -- it is the single
# heaviest call in the fit (dense gene x sample p-value matrices, plus a per-gene
# sort inside every forked worker) and it OOM-killed job 30219333 at 64G. Needs
# --gene-table too since it builds on the same aggregate=TRUE, collapse=FALSE
# results() call.
GENEALL <- "--gene-all-table" %in% args

# fit: report only these samples instead of the whole cohort. Default (NULL) is
# every sample, i.e. exactly the old behaviour.
#
# This is purely a row SELECTION: every p value and FDR is fixed by FRASER()
# before results() is ever called, and both the junction-level FDR
# (pvalsNzscore.R:270) and the gene-level one (:353) are `apply(pvals, 2,
# p.adjust)` -- column-wise, i.e. WITHIN a sample. So one sample's table is
# identical to that sample's slice of the cohort table. What it buys is size and
# memory: the all=TRUE junction pass drops from nSamples x nJunctions rows
# (5,052,663 / 170 MB for the 23-sample LCL cohort) to nJunctions, and the
# aggregate=TRUE gene pass -- which OOM-killed job 30219333 at 64G by
# materialising dense gene x sample matrices -- shrinks by the same factor.
FOCUS <- opt("--results-samples", NULL)
if (!is.null(FOCUS)) {
    FOCUS <- trimws(strsplit(FOCUS, ",", fixed = TRUE)[[1]])
    FOCUS <- FOCUS[nzchar(FOCUS)]
    if (!length(FOCUS)) stop("--results-samples was given but lists no sample")
}
# Escape hatch for the guard that refuses to fit an already-filtered object.
ALLOWFILT <- "--allow-filtered-input" %in% args

# Frozen, read-only control cohort. Everything the run READS that belongs to the
# cohort is resolved here; everything it WRITES goes to --workdir. That is what
# lets many samples run against one cohort at once, and what keeps a case out of
# the panel it was measured against.
COHORT <- opt("--cohort-dir", NULL)
# fit deliverables go here, flat. Default NULL = the historical layout
# (<workdir>/<name>_* , or <workdir>/samples/<ID>/<ID>_* under --results-samples).
OUTDIR <- opt("--out-dir", NULL)

# ---------------------------------------------------------- sample input -----
# Two ways to name the samples this run operates on:
#   --samplesheet F                     a TSV, any number of rows (all modes)
#   --sample-id ID --bam F --sj F       one sample, inline (`add` only)
# They are mutually exclusive: two statements of what to add are two things to
# keep in step. The inline form exists because the per-sample frozen-cohort
# workflow otherwise writes a one-line tsv whose three fields the caller already
# holds -- a Snakemake rule has them as {wildcards.sampleID}/{input.bam}/{input.sj}.
#
# --sample-id is REQUIRED rather than derived from the BAM basename. It is the
# primary key for the colData row, <ID>.splitCounts.tsv, nonSplicedCounts-<ID>.h5
# and the fit's --results-samples/output prefix, and `X.mdup.bam` would derive
# `X.mdup`, which nothing notices until the fit cannot find the sample. Same
# reasoning as parse_fraser.R's required -s.
INLINE <- !is.null(BAM) || !is.null(SJ) || !is.null(SID)

stopifnot(!is.null(WD))

if (INLINE) {
    if (mode != "add")
        stop("--sample-id/--bam/--sj are for `add` only (this is `", mode,
             "`).\n  build: use --samplesheet, it takes the whole cohort.",
             "\n  fit:   needs no sample input at all -- --samplesheet is ",
             "optional and --bam/--sj are never read.")
    if (!is.null(SHEET))
        stop("give EITHER --samplesheet OR --sample-id/--bam/--sj, not both.\n",
             "  Two statements of what to add is one too many; drop whichever ",
             "is not the one you meant.")
    absent <- c("--sample-id", "--bam", "--sj")[c(is.null(SID), is.null(BAM),
                                                  is.null(SJ))]
    if (length(absent))
        stop("the inline form needs all three of --sample-id, --bam and --sj; ",
             "missing: ", paste(absent, collapse = ", "), ".\n",
             "  --sample-id is not derived from the BAM name on purpose: it ",
             "names the split-count file, the non-split cache and the fit's ",
             "output prefix, so a guess only surfaces as a failure at fit time.")
    # A trailing flag with no value makes opt() return NA, not NULL, so the
    # is.null checks above pass and the NA travels all the way into colData.
    blank <- c("--sample-id", "--bam", "--sj")[c(is.na(SID) || !nzchar(SID),
                                                 is.na(BAM) || !nzchar(BAM),
                                                 is.na(SJ)  || !nzchar(SJ))]
    if (length(blank))
        stop("no value given for: ", paste(blank, collapse = ", "),
             " -- check the variables in your batch loop")
} else if (is.null(SHEET) && mode != "fit") {
    # `fit` reads neither a BAM nor an SJ.out.tab: it loads the counts object
    # written by build/add and works from there. The samplesheet it used to
    # demand was fread, existence-checked, and then never referenced again, so
    # the per-sample workflow was writing a tsv purely to satisfy this check.
    stop("`", mode, "` needs --samplesheet F",
         if (mode == "add") "  (or --sample-id ID --bam F --sj F)" else "")
}

if (!is.null(COHORT)) {
    if (mode == "build")
        stop("--cohort-dir is for `add`/`fit` against an existing frozen cohort. ",
             "`build` creates one: run it with --workdir <the cohort dir> alone.")
    COHORT <- normalizePath(COHORT, mustWork = FALSE)
    if (!dir.exists(COHORT)) stop("--cohort-dir does not exist: ", COHORT)
    # Same directory for both would silently reintroduce exactly what the split
    # exists to prevent: the case written into the panel it is measured against.
    if (identical(COHORT, normalizePath(WD, mustWork = FALSE)))
        stop("--cohort-dir and --workdir are the same directory (", COHORT,
             ").\n  The cohort must stay read-only; give --workdir a path of ",
             "its own.")
}

#' Resolve a cohort INPUT: the frozen cohort when there is one and it has the
#' file, otherwise the workdir. The fallback is what keeps every pre-existing
#' single-directory invocation working unchanged.
cohort_in <- function(rel) {
    if (!is.null(COHORT)) {
        p <- file.path(COHORT, rel)
        if (file.exists(p)) return(p)
    }
    file.path(WD, rel)
}

UNIVERSE   <- cohort_in("universe.RDS")
PARAMS     <- cohort_in("run_params.tsv")      # read
PARAMS_OUT <- file.path(WD, "run_params.tsv")  # write -- never into the cohort
SCF_IN     <- cohort_in("splitCountFiles")     # frozen per-sample split counts
SCF        <- file.path(WD, "splitCountFiles") # only the new sample's is written
RUNS       <- file.path(WD, "runs.tsv")
# One row per col4=0 junction the build TRIED to strand-resolve, with its
# outcome. Part of the bundle, read (never written) by a later `add`: it is what
# distinguishes "this control was tested here and genuinely lacks the junction"
# from "this control was never tested here", which a bare 0 in a split count
# file cannot express. Before it existed, the second case was read as the first,
# and 17% of the strand-recovered private calls in the 2026-09-17 add were
# junctions their controls actually carried.
ATTEMPTS    <- cohort_in("strand_attempts.tsv.gz")      # read
ATTEMPTS_OUT <- file.path(WD, "strand_attempts.tsv.gz") # write (build only)

# TWO saved objects, and the split is load-bearing.
#   savedObjects/<name>      counts only, written by build/add. The fit NEVER
#                            writes here, so it stays the full ~1.7M-junction
#                            set and every fit starts from the same input.
#   savedObjects/<name>_fit  filtered + annotated + fitted, written by fit.
# Before this, fit saved its filtered object back over savedObjects/<name>, so
# the fit consumed its own input: the second fit on a workdir loaded ~220k
# junctions instead of ~1.7M and the private-novel table -- which is built from
# the junctions the filter REJECTS -- came out empty (job 30219333 overwrote
# 1,395 rows with a header). The rename happens BEFORE the fit starts, not at
# the final save, because FRASER() writes assays to disk as it goes and those
# writes are what truncated rawCountsJ.h5 to the filtered rows.
FITNAME  <- paste0(NAME, "_fit")
dir.create(WD, showWarnings = FALSE, recursive = TRUE)

# DROP's speed settings: FRASER standalone defaults to SerialParam (single core),
# which is the single biggest slowdown. HDF5 always-on avoids lock contention.
register(MulticoreParam(THREADS)); setAutoBPPARAM(MulticoreParam(THREADS))
options("FRASER.maxSamplesNoHDF5" = 0, "FRASER.maxJunctionsNoHDF5" = -1)
h5disableFileLocking()
BPP <- MulticoreParam(THREADS)

# ---------------------------------------------------------- samplesheet ------
# Both input forms converge here on ONE shape -- a data.table of sampleID,
# bamFile, sjFile -- so nothing downstream has to know which was used. A
# single-row ss is the ordinary case under --cohort-dir, not a special one.
#
# Samplesheet: no header required, columns are positional. A header row, if
# present, is detected and dropped so both forms work.
ss <- if (!is.null(SHEET)) {
    x <- fread(SHEET, header = FALSE, colClasses = "character")
    if (ncol(x) < 3)
        stop("samplesheet needs >=3 tab-separated columns: sampleID, bamFile, sjFile")
    setnames(x, 1:3, c("sampleID", "bamFile", "sjFile"))
    if (x$sampleID[1] == "sampleID" && x$bamFile[1] == "bamFile")
        x <- x[-1]                               # tolerate an optional header row
    x
} else if (INLINE) {
    data.table(sampleID = SID, bamFile = BAM, sjFile = SJ)
} else {
    # fit with no --samplesheet: no sample input is needed. An empty table with
    # the right columns keeps the roster/path-check code below unconditional.
    data.table(sampleID = character(), bamFile = character(), sjFile = character())
}
ss[, sampleID := as.character(sampleID)]
if (anyDuplicated(ss$sampleID)) stop("duplicate sampleIDs in samplesheet")

# --------------------------------------------------------- frozen roster -----
# The frozen cohort states its own membership, in its own saved counts object.
# That is the only place the roster and the cached counts cannot disagree: a
# samplesheet restating 20+ frozen paths in every per-sample run is a second
# copy of the truth, and the two drift. Read it off the object instead, and let
# the samplesheet carry only what this run is actually adding.
CTL <- data.table(sampleID = character(), bamFile = character())
if (!is.null(COHORT) && mode == "add") {
    cf_obj <- file.path(COHORT, "savedObjects", NAME, "fds-object.RDS")
    if (!file.exists(cf_obj))
        stop("--cohort-dir has no counts object for --name ", NAME, ":\n  ",
             cf_obj, "\n  A frozen cohort is a workdir that `build` (and any ",
             "later `add`) has already written. Check --name.")
    CTL <- as.data.table(as.data.frame(
        colData(loadFraserDataSet(COHORT, NAME))))[, .(sampleID, bamFile)]
    CTL[, sampleID := as.character(sampleID)]
    CTL[, bamFile  := as.character(bamFile)]

    # A full sheet (cohort rows + the new one) still works: rows naming a sample
    # the cohort already holds are dropped rather than rejected, so the same
    # sheet can be used with or without --cohort-dir.
    dup <- ss$sampleID %in% CTL$sampleID
    if (any(dup)) {
        message("samplesheet: dropping ", sum(dup), " row(s) already in the ",
                "frozen cohort (", paste(ss$sampleID[dup], collapse = ", "), ")")
        ss <- ss[!dup]
    }
    if (!nrow(ss))
        stop("nothing to add: every samplesheet row is already in the frozen ",
             "cohort at ", COHORT)
    message("frozen cohort: ", nrow(CTL), " sample(s); adding ", nrow(ss), ": ",
            paste(ss$sampleID, collapse = ", "))
}

# Only the rows this run actually reads are path-checked. The frozen cohort's
# BAMs are never opened -- its non-split counts are cached and its split counts
# are already written -- so a control whose BAM has since moved must not fail
# the run.
missing <- ss[!file.exists(bamFile) | !file.exists(sjFile)]
if (nrow(missing)) stop("missing files for: ", paste(missing$sampleID, collapse = ", "))

# ------------------------------------------------------------ run record -----
# Append-only provenance for the workdir. Every add REFITS the whole cohort, so
# a sample's p values depend on which cohort it was called against -- a released
# report is only interpretable next to the run that produced it. Kept as a flat
# TSV of all-character columns so appending can never misalign, and so a run
# from an older version of this script still reads.
RUN_COLS <- c("timestamp", "jobid", "mode", "name", "n_samples", "junctions_in",
              "junctions_passed", "sites_zeroed", "sites_total", "q",
              "filter_quantile", "min_count", "sj_min_count", "count_col",
              "mapq", "strand",
              "focus", "fraser_version", "quantile_min_expression")

run_record <- function(...) {
    v <- list(...)
    rec <- lapply(RUN_COLS, function(k) {
        x <- v[[k]]
        if (is.null(x) || !length(x)) NA_character_ else paste(as.character(x), collapse = ",")
    })
    names(rec) <- RUN_COLS
    as.data.table(rec)
}

write_run_record <- function(rec) {
    fresh <- !file.exists(RUNS)
    if (!fresh) {
        # A runs.tsv written by an older version of this script has fewer
        # columns. Appending a wider row than its header would silently shift
        # every field of every later read, so conform to the header that is
        # already there. NB a column added since that file was created is then
        # NOT recorded -- use a fresh workdir when the new field matters.
        hdr <- names(fread(RUNS, nrows = 0L))
        if (!identical(hdr, names(rec))) {
            for (k in setdiff(hdr, names(rec))) rec[[k]] <- NA_character_
            rec <- rec[, hdr, with = FALSE]
        }
    }
    fwrite(rec, RUNS, sep = "\t", append = !fresh, col.names = fresh)
    message(date(), ": run recorded -> ", RUNS)
    invisible(rec)
}

# Extra fields are passed straight through, NOT patched in afterwards with
# `:=`: inside a data.table `:=` the RHS is evaluated in the table's scope, so a
# bare `q` resolves to the (still NA) COLUMN q instead of the variable, and the
# value is silently lost.
base_record <- function(mode, n_samples, ...) run_record(
    ...,
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    jobid = Sys.getenv("SLURM_JOB_ID", NA_character_),
    mode = mode, name = NAME, n_samples = n_samples,
    min_count = MINCNT, sj_min_count = SJMIN, count_col = CCOL,
    mapq = MAPQ, strand = STRAND,
    fraser_version = as.character(packageVersion("FRASER")))

# The cohort this run operates on: the frozen roster (empty without
# --cohort-dir) followed by the samplesheet rows. Order is load-bearing --
# getSplitReadCountsForAllSamples consumes countFiles POSITIONALLY
# (`names(countList) <- samples(fds)`, countRNAseqData.R), so the count-file
# vector is ordered against samples(fds) rather than against the samplesheet.
ROSTER <- rbind(CTL[, .(sampleID, bamFile)], ss[, .(sampleID, bamFile)])

make_fds <- function(roster = ROSTER) {
    fds <- FraserDataSet(
        colData    = data.table(sampleID = roster$sampleID,
                                bamFile  = roster$bamFile,
                                pairedEnd = TRUE),   # data assumed paired-end
        workingDir = WD, name = NAME)
    strandSpecific(fds) <- STRAND
    scanBamParam(fds)   <- ScanBamParam(mapqFilter = MAPQ)
    fds
}

#' Make the frozen cohort's per-sample non-split caches reachable from the
#' writable workdir.
#'
#' FRASER resolves those caches through workingDir(fds) alone
#' (getNonSplicedCountCacheFolder -> <workingDir>/cache/nonSplicedCounts/<name>),
#' so with the cohort read-only and workingDir pointed at WD it would find
#' nothing and recount EVERY control from BAM -- the one cost this whole design
#' exists to avoid. Symlinking the individual FILES into a REAL directory keeps
#' the reads resolving while the new sample's h5 is created beside them, inside
#' WD. Linking the directory instead would send that write into the cohort.
#'
#' Only nonSplicedCounts is linked. otherCounts is derived from the junction
#' set, which changes the moment a sample is added, so a shared copy would be
#' stale by construction; the fit recomputes it (~5 min) instead.
link_cohort_cache <- function() {
    if (is.null(COHORT)) return(invisible(0L))
    src <- file.path(COHORT, "cache", "nonSplicedCounts", NAME)
    dst <- file.path(WD,     "cache", "nonSplicedCounts", NAME)
    if (!dir.exists(src)) {
        warning("frozen cohort has no non-split cache at ", src,
                " -- every control sample will be recounted from BAM ",
                "(hours, and it needs their BAMs to still exist).",
                immediate. = TRUE)
        return(invisible(0L))
    }
    dir.create(dst, recursive = TRUE, showWarnings = FALSE)
    fs <- list.files(src, pattern = "\\.h5$", full.names = FALSE)
    n  <- 0L
    for (f in fs) {
        to <- file.path(dst, f)
        # file.exists() follows links, so an existing link or a real h5 this
        # workdir already wrote is left alone either way.
        if (file.exists(to)) next
        if (file.symlink(normalizePath(file.path(src, f)), to)) n <- n + 1L
    }
    message("linked ", n, " frozen non-split cache file(s) (", length(fs),
            " in the cohort) -> ", dst)
    invisible(n)
}

# ---------------------------------------------------------- run params -------
# THE ONE SILENT-CORRUPTION PATH THE MAPQ GUARD IN build_fds_counts DOES NOT
# COVER. That guard checks count_col against mapq WITHIN a run. It cannot see
# the CACHE, and the non-split cache is validated by row count alone
# (FRASER countRNAseqData.R:898 `length(spliceSiteCoords) == nrow(cache)`) --
# the mapq and strand that produced those counts are never recorded or checked.
#
# So `add --mapq 0` reuses cohort counts made at MAPQ 255 while counting the new
# sample multimapper-inclusive: different read sets in one rawCountsSS matrix,
# no error, no warning. That is exactly the per-sample asymmetry the frozen
# universe exists to prevent, and it lands on the sample you just added. Same
# trap for --strand, which becomes featureCounts -s and is equally unverified.
#
# The split side self-heals (write_split_countfiles rewrites EVERY sample from
# SJ.out.tab each run), so only the cached non-split side needs this guard --
# but sj_min_count and count_col also define the frozen universe, so they are
# pinned too. sj_min_count doubly so: the per-run splice-site set is re-derived
# with it on EVERY run (assemble_nonsplit_se), so a value that disagrees with
# the build's turns every junction in the gap into a 0-denominator site.
# min_count no longer touches the cache at all -- since 2026-09-01 it is the
# fit's testing threshold only -- but it stays pinned, because it is what a
# released p value was computed under.
# max_star_span only changes which col4=0 junctions are dropped, uniformly
# across samples, so a mismatch is a warning not an error.
#
# strand_method IS critical, and it is not a knob -- there is one method and no
# flag to change it (STRAND_METHOD in fraser_nonsplit_from_sjtab.R). It is
# pinned so that a cohort frozen with the pre-2026-09-18 gated pair-based
# recovery cannot be extended by an `add` using the current one. That
# combination is worse than either alone: the case would be tested thoroughly
# while every control keeps the structural zero its own skipped build gave it,
# so MORE junctions surface as spuriously private, concentrated in exactly the
# high-expression loci where the controls were skipped. The only correct
# response to a mismatch is to rebuild the cohort.
PARAM_CRITICAL <- c("count_col", "mapq", "strand", "min_count", "sj_min_count",
                    "strand_method")
run_params <- function() data.table(
    param = c("count_col", "mapq", "strand", "min_count", "sj_min_count",
              "strand_method", "max_star_span"),
    value = c(CCOL, as.character(MAPQ), as.character(STRAND),
              as.character(MINCNT), as.character(SJMIN), STRAND_METHOD,
              as.character(MAXSPAN)))

write_run_params <- function() {
    fwrite(run_params(), PARAMS_OUT, sep = "\t")
    message(date(), ": recorded counting parameters -> ", PARAMS_OUT)
}

#' Compare this run's parameters against the ones the cohort was pinned with.
#'
#' @param only     restrict the comparison to these params. `fit` passes
#'   "min_count", because it is the ONLY pinned value a fit consumes: since the
#'   discovery/testing split, --min-count feeds filterExpressionAndVariability
#'   and --private-min-count and nothing else. Every other row is a COUNTING
#'   parameter the fit never touches, and comparing them would fail a perfectly
#'   valid fit -- e.g. a legacy cohort backfills sj_min_count=20 against today's
#'   default 10.
#' @param critical which of `only` are hard errors rather than warnings.
#' @param record   write run_params.tsv when the cohort has none. `add` does
#'   (so later runs ARE checked); `fit` must not -- it would stamp the workdir
#'   with counting parameters no counting run ever used.
check_run_params <- function(only = NULL, critical = PARAM_CRITICAL,
                             record = TRUE) {
    want <- run_params()
    if (!is.null(only)) want <- want[param %in% only]
    if (!file.exists(PARAMS)) {
        if (!record) {
            warning("no ", basename(PARAMS), " beside the counts object, so ",
                    "--min-count ", MINCNT, " cannot be checked against the ",
                    "value the cohort was built with. If they differ, these ",
                    "results are not comparable to that cohort's earlier ones.",
                    immediate. = TRUE)
            return(invisible(NULL))
        }
        warning("no ", basename(PARAMS), " in the workdir -- this cohort was ",
                "built before parameter pinning existed, so the counting ",
                "parameters of the cached samples CANNOT be verified. Confirm ",
                "by hand that build used --count-col ", CCOL, " --mapq ", MAPQ,
                " --strand ", STRAND, " --min-count ", MINCNT,
                " --sj-min-count ", SJMIN,
                ". A cohort this old predates the discovery/testing split, so ",
                "it was discovered at its --min-count: --sj-min-count must be ",
                "that same value here. Recording the current values so later ",
                "`add` runs are checked.",
                immediate. = TRUE)
        write_run_params()
        return(invisible(NULL))
    }
    old <- fread(PARAMS, colClasses = "character")
    # A cohort built before the discovery/testing split (2026-09-01) has a single
    # `min_count` row that played both roles, so its DISCOVERY threshold is that
    # same value. Backfill it, rather than reporting `build=NA`, which reads like
    # a corrupt params file instead of what it is: a cohort discovered at 20.
    if (nrow(old) && !("sj_min_count" %in% old$param) &&
        "min_count" %in% old$param)
        old <- rbind(old, data.table(param = "sj_min_count",
                                     value = old[param == "min_count", value[1]]))
    # Params this version no longer pins -- max_region_reads, retired 2026-09-18
    # with the gated pair-based strand recovery whose cost it existed to contain.
    # A bundle written by an older build still carries the row; it has no current
    # counterpart, so comparing it would warn on every legacy cohort about a flag
    # that no longer exists. A param that is MISSING from the bundle but present
    # here is NOT dropped: that is the strand_method case, and it must be caught.
    retired <- setdiff(old$param, run_params()$param)
    if (length(retired)) old <- old[!param %in% retired]
    if (!is.null(only)) old <- old[param %in% only]
    cmp <- merge(old, want, by = "param", all = TRUE,
                 suffixes = c("_build", "_now"))
    bad <- cmp[is.na(value_build) | is.na(value_now) | value_build != value_now]
    if (!nrow(bad)) return(invisible(NULL))
    fmt <- bad[, paste0("  ", param, ": build=", value_build, "  now=", value_now)]
    crit <- bad[param %in% critical]
    # min_count (only=="min_count", the fit-only call below) is NOT a counting
    # parameter -- build/add never read MINCNT, so nothing about the universe,
    # caches or split counts depends on it (see the @param only comment above).
    # What's pinned here is purely a cross-fit convention: every fit against
    # this cohort should use the same testing threshold, so p values/FDR stay
    # comparable case to case. A mismatch is therefore a warning, never a hard
    # stop -- the fit itself is perfectly valid at whatever --min-count was
    # actually passed, just not directly comparable to the cohort's other fits.
    if (nrow(crit) && !is.null(only)) {
        warning("this fit uses different thresholds than the cohort was ",
             "pinned with:\n", paste(fmt, collapse = "\n"),
             "\n  --min-count is the fit's expression gate ",
             "(filterExpressionAndVariability) and the default for ",
             "--private-min-count, so it changes which junctions are tested, ",
             "the multiple-testing burden, and therefore every p value and FDR ",
             "in the output. Proceeding with --min-count ", MINCNT,
             " as given; results are not directly comparable to this cohort's ",
             "other fits at the pinned value.", immediate. = TRUE)
        return(invisible(NULL))
    }
    if ("strand_method" %in% crit$param) {
        stop("this cohort was frozen with a DIFFERENT col4=0 strand-recovery ",
             "method:\n", paste(fmt, collapse = "\n"),
             "\n  build=NA means the bundle predates 2026-09-18, when the ",
             "--max-region-reads gate was removed: that build SKIPPED 83.7% of ",
             "its col4=0 candidates without reading an alignment, so its ",
             "controls carry structural zeros at coordinates they were never ",
             "tested at. Adding a sample with the current method would test the ",
             "CASE thoroughly against controls that were never asked, which ",
             "surfaces junctions the controls do carry as private to the case -- ",
             "the very failure this method was rewritten to fix, made worse. ",
             "There is no flag to downgrade, and there should not be: rebuild ",
             "the cohort (delete universe.RDS and re-run build).")
    }
    if (nrow(crit)) {
        stop("counting parameters differ from the frozen cohort:\n",
             paste(fmt, collapse = "\n"),
             "\n  The cached non-split counts were produced with the build ",
             "values and are NOT re-validated against them, so proceeding ",
             "would mix read sets across samples and manufacture outliers in ",
             "the sample you are adding. Re-run with the build values, or ",
             "delete universe.RDS + the cache and rebuild the whole cohort.")
    }
    warning("non-critical parameters differ from the frozen cohort:\n",
            paste(fmt, collapse = "\n"),
            "\n  These only steer col4=0 strand recovery, which an `add` ",
            "against --cohort-dir does not run at all (the strand map is frozen ",
            "with the cohort), and which a plain `add` applies to every sample ",
            "at once. Either way the cache stays valid and no sample is treated ",
            "differently from the others.",
            immediate. = TRUE)
    invisible(NULL)
}

# ------------------------------------------------------------- counting ------
if (mode %in% c("build", "add")) {
    fds <- make_fds()
    STRANDMAP     <- cohort_in("strand_map.RDS")           # read
    STRANDMAP_OUT <- file.path(WD, "strand_map.RDS")       # write

    if (mode == "build") {
        if (file.exists(UNIVERSE))
            stop("universe already exists: ", UNIVERSE, "\n  Use `add` to extend ",
                 "the cohort. Rebuilding forces a full non-split recount of EVERY ",
                 "sample. If that is really what you want, delete it explicitly.")
        message("=== BUILD: freezing splice-site universe from ", nrow(ss), " samples ===")
    } else {
        if (!file.exists(UNIVERSE)) stop("no universe at ", UNIVERSE, " -- run `build` first")
        message("=== ADD: reusing frozen universe (cached samples will not be recounted) ===")
        # Before any BAM work: the cached counts must have been made the same way.
        check_run_params()
    }

    # THE STRAND MAP IS FROZEN WITH THE COHORT.
    #
    # STAR strands junctions by MOTIF, so non-canonical junctions come out col4=0
    # -> '*'. A '*' mixed with +/- crashes annotateSpliceSite (both in the
    # universe and in the split count files). The map records the read-based
    # strand FRASER would assign, recovered from BAM read orientation; a col4=0
    # junction absent from it is DROPPED from every split count file.
    #
    # Against a frozen cohort we must NOT extend it. The controls' split count
    # files were written once, with the map as it stood at `build`, and they are
    # never rewritten. Recovering a NEW col4=0 junction here would put it in the
    # case's count file with real reads while every control keeps the 0 that the
    # build-time drop gave it -- so a junction the controls genuinely carry
    # surfaces as private to the case, with no p-value to moderate it, straight
    # into <ID>_private_novel.tsv. That is precisely the asymmetry
    # fraser_nonsplit_from_sjtab.R's header forbids ("counting a new site for new
    # samples only, leaving old ones at 0 ... fabricates exactly the outlier
    # signal FRASER hunts for"). It was reachable in practice, and measured: of
    # the 138 junctions the 2026-09-17 add recovered this way, 24 (17%) were
    # carried by controls that had simply never been tested at that coordinate.
    # The 2026-09-18 recovery rewrite removes the skipping that caused it, and
    # strand_attempts.tsv.gz records what was tested, but the freeze stays --
    # the case and the controls must be read the same way, by the same build.
    #
    # Freezing costs only novel NON-CANONICAL junctions in the case: col4=0 is a
    # non-canonical motif, and a real cryptic splice site is GT-AG, i.e. col4=1/2,
    # which never touches this path. To test a sample against a cohort that DOES
    # know its col4=0 junctions, build a cohort that includes it.
    #
    # Without --cohort-dir nothing changes: `add` rewrites EVERY sample's split
    # count file from SJ.out.tab on the same run (see write_split_countfiles
    # below), so an extended map is applied to all of them at once and stays
    # uniform. That is also why max_star_span remains non-critical -- there it
    # moves every sample together, and here it does nothing at all.
    # Recovered private col4=0 junctions kept for --private-into-counts. NULL on
    # every other path, so write_split_countfiles below sees the frozen map
    # unchanged and the legacy behaviour is byte-identical.
    priv_extra <- NULL
    FREEZE_STRANDS <- !is.null(COHORT) && mode == "add"
    if (FREEZE_STRANDS) {
        if (!file.exists(STRANDMAP))
            stop("--cohort-dir has no strand_map.RDS:\n  ", STRANDMAP,
                 "\n  A cohort bundle must carry it; without it the case's ",
                 "col4=0 junctions cannot be treated the same way the frozen ",
                 "controls' were.")
        strand_map <- readRDS(STRANDMAP)
        message("strand map: frozen with the cohort, ", length(strand_map),
                " col4=0 junction(s); NOT extended for the new sample ",
                "(--max-star-span is inert on this path)")
        saveRDS(strand_map, STRANDMAP_OUT)

        # ---- sample-private strand recovery, for private_novel.tsv only ------
        #
        # The map above is deliberately NOT extended -- see the block just
        # above this one: recovering a new col4=0 junction into the SHARED map
        # would create exactly the case-vs-control asymmetry that block
        # describes (the case gets real reads, every control keeps the 0 its
        # own build-time drop gave it, straight into an untested "outlier").
        #
        # But dropping it outright throws away real, well-supported, cleanly
        # resolvable evidence that <name>_private_novel.tsv exists to catch
        # (see "private novel junctions" below: "the most interesting class
        # there is: a genuinely novel junction created by a variant"). Today
        # it never reaches that table, because the drop above happens
        # upstream of it -- before write_split_countfiles even runs, so the
        # junction never enters fds at all.
        #
        # This recovers those candidates from the NEW sample's own BAM only
        # (never a control's), gated on:
        #   - read count >= --private-min-count (PRIVMIN) -- cheap, applied
        #     before any BAM is touched
        #   - strand-vote purity >= --strand-vote-min-frac (STRANDVOTE) -- the
        #     frozen map's own recovery is majority-vote-only (any margin
        #     wins), but "private and untested" should mean unambiguous, not
        #     merely tie-broken
        #   - support in every OTHER sample (this add batch's siblings, by
        #     their own SJ.out.tab; the frozen controls, by their already-
        #     written split count files -- a read, never a recount) capped at
        #     --private-max-other (PRIVOTH), same bar every other row in that
        #     table already clears
        # Nothing here touches strand_map/STRANDMAP_OUT, universe.RDS, or
        # fds -- the tested pipeline and every control's cache stay exactly as
        # frozen as before. Survivors are written to a side file that `fit`
        # reads and folds into the private-novel table with no pValue/padjust,
        # exactly like every other row there.
        #
        # KNOWN BLIND SPOT, shared with every other row in that table: a
        # control whose OWN reads at this coordinate never resolved a strand
        # at ITS build time would itself show 0 here even if it truly carries
        # support -- this is not a new gap, write_private_novel_table's
        # max_other already has it for the both_sites_novel class.
        priv_out <- file.path(WD, "strand_recovered_private.tsv")
        if (file.exists(priv_out)) file.remove(priv_out)
        if (!PRIVTAB) {
            message("--no-private-table: skipping sample-private strand recovery")
        } else {
        priv_jx   <- read_sj_union(ss$sjFile, min_count = PRIVMIN, count_col = CCOL)
        priv_star <- priv_jx[as.character(strand(priv_jx)) == "*"]
        # ⭐ THE SAME --max-star-span GATE THE BUILD APPLIES. Without it this
        # path recovers long-span col4=0 junctions that the CONTROLS could
        # never have had: the build drops them from jx BEFORE strand recovery
        # (see the !FREEZE_STRANDS branch), so they never reach the control
        # strand map or any control's split count file, and their 0 there is
        # structural -- exactly the GNAS failure mode, one filter upstream.
        # MEASURED on the 2026-09-20 add of DSK_260_03 against the 50-sample
        # 09182026 LCL bundle, by scanning all 50 controls' RAW SJ.out.tab
        # (col4=0 rows included): of 476 recovered candidates, ALL 213 refuted
        # by real control support had span > 100 kb (median 238,936 bp, median
        # 698 control reads, raw_n_other_resolved == 0 throughout), while of the
        # 235 candidates at span <= 100 kb, ZERO were refuted. So this one line
        # removes 213/213 of the false positives and costs 23 long-span
        # confirmed-private junctions, which the July span analysis found are
        # mostly chimeric/misalignment anyway. (UBA5 3:132670285-132670967, the
        # junction this whole path exists to catch, is 683 bp.)
        if (MAXSPAN > 0 && length(priv_star)) {
            wide <- width(priv_star) > MAXSPAN
            if (any(wide)) {
                message(date(), ": dropping ", sum(wide), " col4=0 ('*') ",
                        "private candidate(s) spanning > ", MAXSPAN, " bp -- the ",
                        "controls' build dropped these before strand recovery, so ",
                        "their absence there is structural, not evidence of privacy")
                priv_star <- priv_star[!wide]
            }
        }
        if (length(priv_star) && length(strand_map)) {
            known <- !is.na(findOverlaps(priv_star, strand_map, type = "equal",
                                         ignore.strand = TRUE, select = "first"))
            priv_star <- priv_star[!known]
        }
        if (!length(priv_star)) {
            message("no sample-private col4=0 junctions >= --private-min-count (",
                    PRIVMIN, ") outside the frozen strand map")
        } else {
            message(date(), ": attempting strand recovery for ", length(priv_star),
                    " sample-private col4=0 junction(s) >= --private-min-count (",
                    PRIVMIN, ") absent from the frozen strand map (private-novel ",
                    "reporting only; the frozen map/universe are untouched)")
            rec <- recover_star_strand(priv_star, ss, strandMode = STRAND,
                                       mapq = MAPQ, count_col = CCOL,
                                       BPPARAM = BPP,
                                       budget_secs = STRANDBUDGET,
                                       min_vote_frac = STRANDVOTE,
                                       return_metadata = TRUE)
            if (!length(rec)) {
                message("none recovered with vote purity >= ", STRANDVOTE)
            } else {
                # Other-sample support lookup: two file formats, because
                # siblings in this add batch are still raw SJ.out.tab (9 STAR
                # columns) while the frozen controls' split count files are
                # already the reduced 5-column form write_split_countfiles
                # produces (seqnames/start/end/strand/count).
                lookup_sjtab <- function(f, keys, count_col) {
                    if (!file.exists(f) || file.size(f) == 0) return(rep(0L, length(keys)))
                    dt <- fread(f, header = FALSE,
                               col.names = c("seqnames", "start", "end", "strand_code",
                                            "motif", "annot", "n_unique", "n_multi",
                                            "overhang"),
                               showProgress = FALSE)
                    dt[, cnt := switch(count_col, both = n_unique + n_multi,
                                      unique = n_unique, multi = n_multi)]
                    m <- match(keys, paste(dt$seqnames, dt$start, dt$end, sep = ":"))
                    ifelse(is.na(m), 0L, dt$cnt[m])
                }
                lookup_splitcounts <- function(f, keys) {
                    if (!file.exists(f) || file.size(f) == 0) return(rep(0L, length(keys)))
                    dt <- fread(f, header = TRUE, showProgress = FALSE)
                    if (!nrow(dt)) return(rep(0L, length(keys)))
                    m <- match(keys, paste(dt$seqnames, dt$start, dt$end, sep = ":"))
                    ifelse(is.na(m), 0L, dt$count[m])
                }
                # THE BLIND SPOT THIS CLOSES. A frozen control's split count
                # file cannot contain a col4=0 junction its own build failed to
                # strand-resolve -- write_split_countfiles drops exactly those
                # rows. So lookup_splitcounts returns a STRUCTURAL 0 for them,
                # indistinguishable from "this control does not carry the
                # junction", and the candidate is reported private. Measured on
                # the 2026-09-17 add: 24 of 138 strand-recovered rows (17%) were
                # refuted this way, every one with raw_n_other_resolved == 0 --
                # the worst, GNAS 20:58895685-58903527, is carried by ALL 50
                # controls at up to 12,146 reads and would have been released as
                # a private novel junction in a disease gene.
                #
                # strand_attempts.tsv.gz is the build's own record of what it
                # tried, and it carries max_sj_count: the highest col4=0 read
                # count for that junction across the cohort's samples, whether
                # or not a strand was ever resolved. That is the number the
                # privacy check always needed and could not see.
                #
                # RESIDUAL GAP, deliberately left and reported rather than
                # hidden: the record covers junctions that reached the build's
                # --sj-min-count in at least one sample. A junction no control
                # reached that threshold on is absent from it, so controls
                # carrying 1-9 reads there are still invisible (8 of the 138
                # were in that 3-9 read band). Those rows get
                # control_evidence="not_in_build_universe" -- weaker evidence of
                # privacy than "tested_absent", not a silent equivalent.
                cand_key  <- paste(seqnames(rec), start(rec), end(rec), sep = ":")
                other_max <- integer(length(rec))
                sib <- setdiff(ss$sampleID, mcols(rec)$carrier_sampleID)
                for (f in ss$sjFile[match(sib, ss$sampleID)])
                    other_max <- pmax(other_max, lookup_sjtab(f, cand_key, CCOL))
                if (nrow(CTL))
                    for (f in file.path(SCF_IN, paste0(CTL$sampleID, ".splitCounts.tsv")))
                        other_max <- pmax(other_max, lookup_splitcounts(f, cand_key))

                ctl_raw  <- rep(NA_integer_,   length(rec))  # build-time col4=0 support
                ctl_out  <- rep(NA_character_, length(rec))  # what the build made of it
                if (file.exists(ATTEMPTS)) {
                    att <- fread(ATTEMPTS, showProgress = FALSE,
                                 colClasses = list(character = "seqnames"))
                    ai  <- match(cand_key, paste(att$seqnames, att$start, att$end, sep = ":"))
                    ctl_raw <- att$max_sj_count[ai]
                    ctl_out <- att$outcome[ai]
                    message(sprintf("control strand-attempt record: %d/%d candidate(s) found in %s",
                                    sum(!is.na(ai)), length(ai), basename(ATTEMPTS)))
                } else {
                    warning("no strand_attempts.tsv.gz in the cohort bundle:\n  ", ATTEMPTS,
                            "\n  This bundle was built before the record existed (or by a ",
                            "build that resolved no col4=0 junctions). The controls' own ",
                            "unresolved col4=0 reads are therefore invisible here, so a ",
                            "junction they carry can still be reported as private to this ",
                            "sample -- 17% of such rows were, when this was measured. The ",
                            "column control_evidence will read 'no_record' for every row. ",
                            "Rebuild the cohort to close it.",
                            immediate. = TRUE)
                    ctl_out <- rep("no_record", length(rec))
                }
                # No control had a single col4=0 read here -> the 0 downstream is
                # real absence. Reads existed but the build never got a verdict
                # (timeout, or a BAM-level error) -> the 0 is structural and says
                # nothing about carriage, so the candidate is held back rather
                # than released as private. Reads existed and were read but do
                # not span (unresolved_no_votes) -> max_sj_count still judges it.
                evidence <- ifelse(is.na(ctl_out), "not_in_build_universe",
                            ifelse(ctl_out == "no_record", "no_record",
                            ifelse(ctl_out == "not_attempted_no_sj_reads", "tested_absent",
                            ifelse(ctl_out %in% c("unresolved_timeout", "unresolved_error"),
                                   "untested", "tested_present"))))
                ctl_raw[is.na(ctl_raw)] <- 0L
                other_max <- pmax(other_max, ctl_raw)

                keep <- other_max <= PRIVOTH & evidence != "untested"
                n_ref <- sum(other_max > PRIVOTH)
                n_unk <- sum(evidence == "untested" & other_max <= PRIVOTH)
                message(sprintf("strand-recovered %d/%d candidate(s); %d pass --private-max-other (%d)%s",
                                length(rec), length(priv_star), sum(keep), PRIVOTH,
                                if (n_ref + n_unk)
                                    sprintf(" [%d refuted by control support, %d dropped as never tested in the controls]",
                                            n_ref, n_unk) else ""))
                if (any(keep)) {
                    rec <- rec[keep]
                    priv_dt <- data.table(
                        seqnames = as.character(seqnames(rec)), start = start(rec),
                        end = end(rec), width = width(rec),
                        strand = as.character(strand(rec)),
                        sampleID = mcols(rec)$carrier_sampleID,
                        count = mcols(rec)$carrier_count,
                        vote_frac = mcols(rec)$vote_frac,
                        max_other_sample_count = other_max[keep],
                        control_evidence = evidence[keep])
                    fwrite(priv_dt, priv_out, sep = "\t")
                    message("strand-recovered private junctions: ", nrow(priv_dt),
                            "  ->  ", priv_out)
                    # granges() strips the recovery mcols (vote_frac,
                    # carrier_*). Not for safety -- c() does NOT error on an
                    # mcols mismatch, it silently fills the bare operand's rows
                    # with NA (checked) -- but so the map stays the plain
                    # coordinate+strand object every other caller expects,
                    # rather than growing a mostly-NA column. write_split_
                    # countfiles reads seqnames/start/end/strand only.
                    if (PRIVCOUNTS) priv_extra <- granges(rec)
                }
            }
        }
        }
    } else {

    # DISCOVERY threshold (--sj-min-count), not the fit's testing threshold.
    jx   <- read_sj_union(ss$sjFile, min_count = SJMIN, count_col = CCOL)

    # Drop col4=0 ('*') junctions with an intron span > --max-star-span. STAR emits
    # col4=0 for non-canonical splice motifs; genome-wide these are dominated by
    # huge-span (up to ~900kb) junctions clustered in the chr6 MHC -- almost all
    # artifact, never testable, and the strand-recovery hotspot. Canonical (+/-)
    # junctions are never touched, no matter how long the intron (real large genes).
    if (MAXSPAN > 0) {
        drop <- as.character(strand(jx)) == "*" & width(jx) > MAXSPAN
        if (any(drop)) {
            message(date(), ": dropping ", sum(drop), " col4=0 ('*') junction(s) ",
                    "spanning > ", MAXSPAN, " bp (non-canonical, likely artifact)")
            jx <- jx[!drop]
        }
    }
    star <- jx[as.character(strand(jx)) == "*"]
    strand_map <- if (file.exists(STRANDMAP)) readRDS(STRANDMAP) else jx[0]
    if (length(star)) {
        if (length(strand_map)) {                 # only recover NEW '*' junctions
            known <- !is.na(findOverlaps(star, strand_map, type = "equal",
                                         ignore.strand = TRUE, select = "first"))
            new_star <- star[!known]
        } else new_star <- star
        if (length(new_star))
            strand_map <- c(strand_map,
                            recover_star_strand(new_star, ss, strandMode = STRAND,
                                                mapq = MAPQ, count_col = CCOL,
                                                BPPARAM = BPP,
                                                budget_secs = STRANDBUDGET,
                                                attempts_out = if (mode == "build")
                                                    ATTEMPTS_OUT else NULL))
    }
    saveRDS(strand_map, STRANDMAP_OUT)
    }   # end !FREEZE_STRANDS

    if (mode == "build") {
        build_frozen_universe(apply_strand_map(jx, strand_map), UNIVERSE)
        # Pin the parameters the universe and the caches were made with, so a
        # later `add` cannot silently count the new sample a different way.
        write_run_params()
    }

    # Split counts: pure text, no BAM. Only the samplesheet's own samples are
    # written. Against a frozen cohort the controls' count files were written
    # once, by the run that built it, from SJ.out.tab files and a strand map
    # that are themselves frozen -- rewriting them can only reproduce the same
    # bytes. (The 2026-08-28 add spent ~12 of its 13 minutes re-reading 23
    # control SJ.out.tab files and rewriting 248 MB of their count files, then
    # discarded all of it when the split-count cache short-circuited.)
    # --private-into-counts: extend the map used to WRITE this sample's count
    # file with its own recovered private junctions, without touching the map
    # that was saved to the workdir (already written, frozen, above). The
    # controls' count files are read from the bundle and never rewritten, so
    # these junctions enter the fds with K>0 in this sample and 0 elsewhere:
    # they gain a real jaccard and contribute to N_psi5/N_psi3 at their splice
    # sites, and the expression filter still rejects them (one carrier cannot
    # clear the quantile term), so they are reported, never tested.
    map_for_counts <- strand_map
    if (!is.null(priv_extra) && length(priv_extra)) {
        map_for_counts <- c(strand_map, priv_extra)
        message("--private-into-counts: ", length(priv_extra), " recovered private ",
                "col4=0 junction(s) folded into ", paste(ss$sampleID, collapse = ","),
                "'s split counts (frozen strand_map/universe untouched); their ",
                "reads now count toward neighbouring junctions' N in this sample")
    }
    cf_new <- write_split_countfiles(ss$sjFile, ss$sampleID, SCF,
                                     count_col = CCOL, strand_map = map_for_counts)
    cf_ctl <- character(0)
    if (nrow(CTL)) {
        cf_ctl <- file.path(SCF_IN, paste0(CTL$sampleID, ".splitCounts.tsv"))
        gone   <- !file.exists(cf_ctl)
        if (any(gone))
            stop("the frozen cohort is missing split count files for ",
                 sum(gone), " of its own sample(s):\n  ",
                 paste(head(CTL$sampleID[gone], 5), collapse = ", "),
                 "\n  A cohort bundle must carry splitCountFiles/ alongside ",
                 "savedObjects/, universe.RDS and cache/. Looked in ", SCF_IN)
    }
    # POSITIONAL: countFiles is zipped against samples(fds) by position, so this
    # must be ordered by the fds colData, not by the samplesheet.
    cf <- c(cf_ctl, cf_new)[match(samples(fds), c(CTL$sampleID, ss$sampleID))]
    stopifnot(!anyNA(cf), length(cf) == ncol(fds))

    # Frozen non-split caches have to be visible under the writable workdir
    # before any counting starts, or every control is recounted from BAM.
    link_cohort_cache()

    # Two parallelism axes over the same core budget (THREADS): W sample-workers
    # (bplapply) x C featureCounts threads/sample (NcpuPerSample). W*C must not
    # exceed THREADS or the cores oversubscribe. `add` (recount only the new
    # sample) and small cohorts leave workers idle unless C picks up the slack.
    # Against a frozen cohort exactly nrow(ss) samples are actually counted --
    # every other one returns from cache in milliseconds. Sizing `auto` by the
    # full cohort would hand the single new sample ONE featureCounts thread and
    # leave the other 14 idle for the whole BAM pass (12.7 min on a deep LCL
    # sample in job 29824253). Without --cohort-dir which samples are cached is
    # not knowable here, so that path keeps the old sizing exactly.
    # The worker pool still spans the WHOLE cohort: the same BPPARAM drives
    # mergeCounts over every sample's split counts, which is the expensive half
    # of a first `add` into a fresh workdir. Only the per-sample featureCounts
    # thread count follows ncount. The two do not oversubscribe, because a
    # worker whose sample is cached returns without spawning any threads:
    # real usage is ncount x ncps ~ THREADS in both the cached and cold cases.
    nsamp   <- ncol(fds)
    ncount  <- if (!is.null(COHORT)) max(1L, nrow(ss)) else nsamp
    if (identical(NCPS, "auto")) {
        ncps    <- max(1L, THREADS %/% ncount)
        workers <- min(THREADS, nsamp)
    } else {
        ncps    <- max(1L, as.integer(NCPS))
        workers <- max(1L, THREADS %/% ncps)
    }
    message(sprintf(paste0("counting layout: %d sample-worker(s) x %d featureCounts ",
                           "thread(s)/sample; %d sample(s) to count = %d core(s) of %d"),
                    workers, ncps, ncount, ncount * ncps, THREADS))
    COUNTBPP <- MulticoreParam(workers)

    # SJMIN, not MINCNT: this argument selects the junctions whose splice sites
    # form THIS run's site set, and that has to be the same threshold that froze
    # the universe, or every junction in the gap gets a 0 non-split denominator.
    fds <- build_fds_counts(fds, countFiles = cf, universe_path = UNIVERSE,
                            minExpressionInOneSample = SJMIN,
                            count_col = CCOL, NcpuPerSample = ncps,
                            BPPARAM = COUNTBPP)
    saveFraserDataSet(fds)

    message("\n=== counts ready ===")
    message("thresholds: discovery (--sj-min-count) ", SJMIN,
            " | testing (--min-count, applied at fit) ", MINCNT)
    message("samples  : ", ncol(fds))
    message("junctions: ", nrow(fds))
    message("sites    : ", nrow(nonSplicedReads(fds)))

    # Sites the cohort now uses that the frozen universe never counted are set to
    # 0 for EVERY sample. That is uniform, so it can only lose retention
    # sensitivity, never fabricate an outlier -- but it grows with each add, and
    # the only fix is a rebuild + full recount. Record it so the drift is
    # visible over time rather than living in one job's stderr.
    zeroed <- get0("sites_zeroed", envir = .fraser_run_stats)
    total  <- get0("sites_total",  envir = .fraser_run_stats)
    if (!is.null(zeroed) && !is.null(total) && total > 0) {
        frac <- 100 * zeroed / total
        message(sprintf("universe drift: %d of %d sites (%.2f%%) absent from the ",
                        zeroed, total, frac),
                "frozen universe -> 0 for all samples")
        if (frac >= 1)
            warning(sprintf("universe drift is %.2f%% of splice sites. Past ~1%% ",
                            frac), "consider rebuilding the universe (full recount).",
                    immediate. = TRUE)
    }
    write_run_record(base_record(mode, ncol(fds),
                                 junctions_in = nrow(fds),
                                 sites_zeroed = zeroed, sites_total = total))

    # The fit takes no sample input -- it reads the counts object this run just
    # wrote -- so the hint names only the workdir, whichever form got us here.
    message("\nnext: Rscript fraser_cohort.R fit --workdir ", WD,
            if (!is.null(COHORT)) paste0(" --cohort-dir ", COHORT) else "",
            if (mode == "add" && nrow(ss) == 1L)
                paste0(" --results-samples ", ss$sampleID) else "")
}

# ------------------------------------------- private novel junctions ---------
#' Junctions FRASER cannot test, carried by exactly one sample.
#'
#' A junction with reads in one sample and none anywhere else has no cohort
#' dispersion to test against, so the jaccard model cannot produce a meaningful
#' p-value and the expression filter drops it. That is statistically correct,
#' but in a proband-against-frozen-cohort run it silently discards the most
#' interesting class there is: a genuinely novel junction created by a variant.
#' Report those on evidence instead of on a p-value.
#'
#' Why they fail the filter, specifically: the gate is
#' `rowMaxs(K) >= minExpression AND rowQuantiles(N_jaccard, --filter-quantile) >= 10`.
#' The first gate is a max ACROSS samples, so it does not reject single-sample
#' junctions -- the second one does, and no quantile value can save a junction
#' whose sites have no cohort coverage (the quantile position lands among the
#' zeros, so the carrier's own N is never read). For a
#' junction nobody else carries, K = 0 in every other sample AND no other
#' junction uses its splice sites, so N collapses to the non-split read-through
#' term alone. If both sites are absent from the frozen universe their non-split
#' counts were set to 0 for EVERY sample (uniformly -- see
#' fraser_nonsplit_from_sjtab.R), N is 0 across the cohort, and the 75%-quantile
#' gate can never be met. The *_site_in_universe columns distinguish that case
#' from a site that WAS counted and genuinely has no read-through.
#'
#' This is NOT an outlier list: no pValue, no padjust, deliberately. Rank it on
#' read support, gene context and variant proximity downstream. Expect mapping
#' artifacts (IG/TR, MHC, blacklist regions) to dominate until parse_fraser.R's
#' artifact filters are applied.
#' The empty table, as its own function: a no-findings run must still write a
#' readable header (downstream fread() breaks on a 0-byte file), and the
#' per-sample split needs the same header when a sample carries nothing.
private_novel_empty <- function() data.table(
    seqnames = character(), start = integer(), end = integer(),
    width = integer(), strand = character(), sampleID = character(),
    count = integer(), jaccard = numeric(),
    max_other_sample_count = integer(), n_other_samples_nonzero = integer(),
    startSite_in_universe = logical(), endSite_in_universe = logical(),
    both_sites_novel = logical(), nonsplit_sample = integer(),
    nonsplit_cohort_max = integer(), maxCount = integer(),
    quantileValueN = numeric(), passedExpression = logical(),
    passedVariability = logical(), passed = logical(),
    pValue = numeric(), padjust = numeric(),
    hgnc_symbol = character(),
    strand_recovered_private = logical(), vote_frac = numeric(),
    control_evidence = character())

#' Build private-novel-style rows directly off a set of fds junction indices.
#' `sel` indexes rowRanges(fds,type="j"); `sidx` is, per row of `sel`, the
#' 1-based sample-column its value should be read from.
#'
#' Factored out 2026-09-21 so a junction can be described this way from TWO
#' entry points: the original `!passed` candidate scan below, and (new) the
#' extra_file lookup for --private-into-counts rows that DID enter the fds but
#' PASSED the filter -- a strand-recovered junction sharing one splice site
#' with cohort-active splicing gets real N there and becomes genuinely
#' testable (its own row in <ID>_results_all.tsv.gz), which the plain
#' "!passed" scan can never see, so those rows were falling through to the
#' NA-filled append path below despite having real jaccard/counts sitting in
#' the fds the whole time. Real jaccard/counts/N/filter-flags are pulled
#' identically either way; only provenance (strand_recovered_private/
#' vote_frac/control_evidence) differs by caller.
#'
#' `pValue`/`padjust` are always NA here -- this runs BEFORE FRASER() fits the
#' model (see write_private_novel_table's own doc comment), so no p-value
#' exists yet for anyone, `passed` or not. `enrich_private_novel_pvalues()`,
#' called once the fit is done, patches them in afterward for `passed=TRUE`
#' rows only (a `!passed` row has no model p-value, ever -- that is the
#' entire reason this table exists instead of just reading results_all).
rows_from_fds <- function(fds, sel, sidx, universe_path, anno = NULL) {
    empty_cols <- setdiff(names(private_novel_empty()),
                          c("strand_recovered_private", "vote_frac", "control_evidence"))
    if (!length(sel)) return(private_novel_empty()[, ..empty_cols])
    rows <- seq_along(sel)
    jr  <- rowRanges(fds, type = "j")[sel]
    key <- function(gr) paste(seqnames(gr), start(gr), end(gr), strand(gr), sep = ":")
    if (file.exists(universe_path)) {
        uk   <- key(load_frozen_universe(universe_path))
        in_l <- key(GRanges(seqnames(jr), IRanges(start(jr) - 1L, start(jr)),
                            strand(jr))) %in% uk
        in_r <- key(GRanges(seqnames(jr), IRanges(end(jr), end(jr) + 1L),
                            strand(jr))) %in% uk
    } else {
        in_l <- NA; in_r <- NA
    }
    jmc <- mcols(fds, type = "j")
    col <- function(n) if (n %in% colnames(jmc)) jmc[[n]][sel] else NA
    ns <- NULL
    if (all(c("startID", "endID") %in% colnames(jmc))) {
        ssr  <- nonSplicedReads(fds)
        ssc  <- as.matrix(assay(ssr, "rawCountsSS"))
        ssid <- mcols(ssr)$spliceSiteID
        pick <- function(id) {
            m  <- match(id, ssid)
            x  <- matrix(0L, nrow = length(id), ncol = ncol(ssc))
            ok <- !is.na(m)
            x[ok, ] <- ssc[m[ok], , drop = FALSE]
            x
        }
        ns <- pick(jmc$startID[sel]) + pick(jmc$endID[sel])
    }
    k      <- as.matrix(K(fds, type = "j"))[sel, , drop = FALSE]
    kother <- k; kother[cbind(rows, sidx)] <- -1L
    jac <- tryCatch(as.matrix(assay(fds, "jaccard"))[sel, , drop = FALSE],
                    error = function(e) NULL)
    hgnc <- if (!is.null(anno)) annotate_gene_symbol(jr, anno$txdb, anno$map)
            else rep(NA_character_, length(jr))
    data.table(
        seqnames = as.character(seqnames(jr)),
        start    = start(jr),
        end      = end(jr),
        width    = width(jr),
        strand   = as.character(strand(jr)),
        sampleID = samples(fds)[sidx],
        count    = k[cbind(rows, sidx)],
        jaccard  = if (is.null(jac)) NA_real_ else jac[cbind(rows, sidx)],
        max_other_sample_count  = matrixStats::rowMaxs(kother),
        n_other_samples_nonzero = rowSums(k > 0) - 1L,
        startSite_in_universe = in_l,
        endSite_in_universe   = in_r,
        both_sites_novel      = !in_l & !in_r,
        nonsplit_sample     = if (is.null(ns)) NA_integer_ else ns[cbind(rows, sidx)],
        nonsplit_cohort_max = if (is.null(ns)) NA_integer_ else matrixStats::rowMaxs(ns),
        maxCount          = col("maxCount"),
        quantileValueN    = col("quantileValueN"),
        passedExpression  = col("passedExpression"),
        passedVariability = col("passedVariability"),
        passed            = col("passed"),
        pValue = NA_real_, padjust = NA_real_,
        hgnc_symbol       = hgnc)
}

#' Gene symbol for arbitrary ranges, by overlap with a txdb -- the same lookup
#' annotateRangesWithTxDb (fraser_gene_annotation.R) does internally
#' (getAnnotationFeature is FRASER's own exported helper for exactly this),
#' but callable on any GRanges rather than only an fds's rowRanges.
#'
#' Needed here because write_private_novel_table runs BEFORE the fds itself
#' gets annotated (annotation is positioned right before FRASER(), see the
#' note at that block -- gene-level p-values need it there, not here), so
#' fds's own mcols have no hgnc_symbol yet at this point; and because the
#' strand-recovered private rows (FREEZE_STRANDS above) never enter fds at
#' all, so they could never inherit it from fds either way.
annotate_gene_symbol <- function(gr, txdb, map, keytype = "gene_id",
                                 feature = "gene_name", featureName = "hgnc_symbol") {
    if (!length(gr)) return(character(0))
    style <- seqlevelsStyle(gr)[1]
    suppressWarnings({
        seqlevelsStyle(map$seqnames) <- style
        seqlevelsStyle(txdb)         <- style
    })
    anno <- genes(txdb)
    tmp <- if (is.data.table(map)) {
        merge(x = as.data.table(anno)[, .(gene_id)], y = map,
             by.y = keytype, by.x = "gene_id", all.x = TRUE, sort = FALSE)[
             , c("gene_id", feature), with = FALSE]
    } else {
        as.data.table(AnnotationDbi::select(map, keys = mcols(anno)[, "gene_id"],
                                            columns = feature, keytype = keytype))
    }
    tmp[, uniqueID := .GRP, by = keytype]
    anno <- anno[tmp$uniqueID]
    mcols(anno)[[featureName]] <- tmp[[feature]]
    anno <- anno[!is.na(mcols(anno)[[featureName]])]
    anno <- anno[mcols(anno)[[featureName]] != ""]
    # Not exported by FRASER (checked: absent from getNamespaceExports("FRASER"))
    # -- it is the same internal helper annotateRangesWithTxDb itself calls.
    FRASER:::getAnnotationFeature(gr, featureName, anno)
}

#' @param anno        list(txdb=, map=) from build_txdb_and_map(), or NULL to
#'   skip gene annotation (hgnc_symbol stays NA).
#' @param extra_file  optional TSV of sample-private, strand-recovered col4=0
#'   junctions (written by the FREEZE_STRANDS block above, one add run
#'   earlier) -- seqnames/start/end/width/strand/sampleID/count/vote_frac/
#'   max_other_sample_count, already filtered to max_other <= max_other.
#'   Without --private-into-counts these never entered fds, so every
#'   fds-derived column is NA by construction. WITH it (2026-09-21) three
#'   outcomes are possible, checked in order: (1) the junction is a `!passed`
#'   base row already -- provenance is merged onto it; (2) it entered fds but
#'   PASSED the filter (one splice site shared with cohort-active splicing
#'   gives it real N) -- real jaccard/counts/N are pulled from the fds
#'   directly via rows_from_fds(), same as a base row, `passed` will read
#'   TRUE, and it has its own tested row in <ID>_results_all.tsv.gz too; (3)
#'   it never entered fds at all -- the old NA-filled append, `passed=NA`.
write_private_novel_table <- function(fds, passed, universe_path, out,
                                      min_count = 20, max_other = 0,
                                      anno = NULL, extra_file = NULL) {
    if (ncol(fds) < 2) {
        message("private-novel table needs >= 2 samples; skipped")
        return(invisible(NULL))
    }
    k    <- as.matrix(K(fds, type = "j"))
    kmax <- matrixStats::rowMaxs(k)      # matrixStats is a hard FRASER dependency

    # Narrow before doing anything wide: only filter-rejected, well-supported rows.
    cand <- which(!passed & kmax >= min_count)
    keep <- logical(0)
    if (length(cand)) {
        kc     <- k[cand, , drop = FALSE]
        top    <- max.col(kc, ties.method = "first")     # the sole carrier
        ko     <- kc; ko[cbind(seq_len(nrow(ko)), top)] <- -1L
        kother <- matrixStats::rowMaxs(ko)               # best support anywhere else
        keep   <- rowSums(kc >= min_count) == 1L & kother <= max_other
    }
    # dt starts as the (possibly empty) fds-derived table; the extra_file
    # block below can still add rows to it even when this part is empty --
    # zero base candidates must not drop a strand-recovered private junction
    # that DID clear its own gate.
    dt <- private_novel_empty()
    if (any(keep)) {
    sel <- cand[keep]; sidx <- top[keep]
    dt <- cbind(rows_from_fds(fds, sel, sidx, universe_path, anno), data.table(
        strand_recovered_private = FALSE,
        vote_frac = NA_real_,
        # A base row's junction is in the frozen strand map by definition, so no
        # control's read at that coordinate was ever dropped for want of a
        # strand: max_other_sample_count is a real count off the fds, not a
        # structural 0. Measured: raw_n_other_col4zero == 0 across all 2,082
        # base rows of the 2026-09-17 add. This class never had the blind spot.
        control_evidence = "counted_in_fds"))
    }

    # Sample-private col4=0 junctions recovered from THIS sample's own BAM
    # against a frozen --cohort-dir (see the FREEZE_STRANDS block, far above):
    # never counted into fds, so every fds-derived column here is NA by
    # construction, not a gap in this function.
    n_extra <- 0L
    if (!is.null(extra_file) && file.exists(extra_file) && file.size(extra_file) > 0) {
        ex <- fread(extra_file, showProgress = FALSE, colClasses = c(
            seqnames = "character", strand = "character", sampleID = "character"))
        # Under --private-into-counts these junctions ARE in the fds, so they
        # have already come through as ordinary base rows -- with a real
        # jaccard, real site-membership flags and real non-split numbers. The
        # side file then carries only provenance the fds cannot know
        # (strand_recovered_private, vote_frac, control_evidence), so it is
        # merged ONTO those rows instead of appended, or the junction would
        # appear twice with the richer copy indistinguishable from the poorer.
        # Keyed on junction + sample, since a private junction has exactly one
        # carrier. Without the flag nothing matches and this is a no-op.
        if (nrow(ex) && nrow(dt)) {
            exk <- paste(ex$seqnames, ex$start, ex$end, ex$sampleID, sep = ":")
            dtk <- paste(dt$seqnames, dt$start, dt$end, dt$sampleID, sep = ":")
            hit <- match(dtk, exk)
            if (any(!is.na(hit))) {
                i <- which(!is.na(hit)); j <- hit[i]
                dt[i, `:=`(strand_recovered_private = TRUE,
                           vote_frac = ex$vote_frac[j],
                           control_evidence = if (!is.null(ex$control_evidence))
                               as.character(ex$control_evidence[j]) else "unknown")]
                message("private-novel: ", length(i), " strand-recovered junction(s) ",
                        "already present as counted rows (--private-into-counts); ",
                        "merged provenance rather than duplicating")
                ex <- ex[-j]
            }
        }
        # --private-into-counts folded these into the sample's own countfile,
        # so a junction can be in fds WITHOUT being a `!passed` candidate: one
        # shared splice site is enough real N to clear the quantile term and
        # pass, at which point it's not "private" in the statistical sense any
        # more, it's a genuinely tested junction (own p value in
        # <ID>_results_all.tsv.gz) -- the match above only ever looks at
        # `!passed` rows, so this class fell through to the NA-append below
        # with a real jaccard sitting unused in the fds. Measured on the
        # 2026-09-21 DSK_260_03 --private-into-counts run: 85/235
        # strand-recovered junctions were exactly this case (incl. UBA5
        # itself), all with passed=TRUE.
        if (nrow(ex)) {
            full_jr  <- rowRanges(fds, type = "j")
            full_key <- paste(seqnames(full_jr), start(full_jr), end(full_jr), sep = ":")
            ex_key   <- paste(ex$seqnames, ex$start, ex$end, sep = ":")
            fds_row  <- match(ex_key, full_key)
            fds_col  <- match(ex$sampleID, samples(fds))
            found    <- !is.na(fds_row) & !is.na(fds_col)
            if (any(found)) {
                recov <- cbind(rows_from_fds(fds, fds_row[found], fds_col[found],
                                             universe_path, anno), data.table(
                    strand_recovered_private = TRUE,
                    vote_frac = ex$vote_frac[found],
                    control_evidence = if (!is.null(ex$control_evidence))
                        as.character(ex$control_evidence[found]) else "unknown"))
                dt <- rbind(dt, recov)
                message("private-novel: ", sum(found), " more strand-recovered ",
                        "junction(s) found in the fds under a filter outcome the ",
                        "!passed scan can't see (--private-into-counts); real ",
                        "jaccard recovered from the count matrix rather than ",
                        "reported NA -- ", sum(recov$passed, na.rm = TRUE),
                        " of them PASSED and have their own tested row in ",
                        "<ID>_results_all.tsv.gz")
                ex <- ex[!found]
            }
        }
        # Genuinely never entered fds -- e.g. dropped before assembly, or a
        # workdir without --private-into-counts. Every fds-derived column NA
        # by construction, not a gap in this function.
        if (nrow(ex)) {
            n_extra <- nrow(ex)
            ex_gr <- GRanges(ex$seqnames, IRanges(ex$start, ex$end), strand = ex$strand)
            ex_hgnc <- if (!is.null(anno)) annotate_gene_symbol(ex_gr, anno$txdb, anno$map)
                      else rep(NA_character_, nrow(ex))
            dt <- rbind(dt, data.table(
                seqnames = ex$seqnames, start = ex$start, end = ex$end,
                width = ex$width, strand = ex$strand, sampleID = ex$sampleID,
                count = ex$count, jaccard = NA_real_,
                max_other_sample_count = ex$max_other_sample_count,
                n_other_samples_nonzero = NA_integer_,
                startSite_in_universe = NA, endSite_in_universe = NA,
                both_sites_novel = NA,
                nonsplit_sample = NA_integer_, nonsplit_cohort_max = NA_integer_,
                maxCount = NA_integer_, quantileValueN = NA_real_,
                passedExpression = NA, passedVariability = NA, passed = NA,
                pValue = NA_real_, padjust = NA_real_,
                hgnc_symbol = ex_hgnc, strand_recovered_private = TRUE,
                vote_frac = ex$vote_frac,
                # Absent from a side file written before the control-evidence
                # check existed; "unknown" is the honest reading of that, not
                # "tested_absent".
                control_evidence = if (!is.null(ex$control_evidence))
                    as.character(ex$control_evidence) else "unknown"))
        }
    }

    setorder(dt, -count)
    fwrite(dt, out, sep = "\t")
    novel <- if (all(is.na(dt$both_sites_novel))) "unknown (no universe.RDS)"
             else as.character(sum(dt$both_sites_novel, na.rm = TRUE))
    # Counted off the column, not off n_extra: under --private-into-counts a
    # strand-recovered junction arrives as a COUNTED base row and is merged,
    # not appended, so n_extra alone would report 0 of them.
    n_sr      <- sum(dt$strand_recovered_private, na.rm = TRUE)
    n_sr_pass <- sum(dt$strand_recovered_private & dt$passed, na.rm = TRUE)
    message("private novel junctions: ", nrow(dt), " (",
            uniqueN(dt$sampleID), " sample(s); ", novel,
            " with both sites absent from the frozen universe; ", n_sr,
            " strand-recovered col4=0 -- ", n_sr_pass, " PASSED the fit filter ",
            "(own p value in <ID>_results_all.tsv.gz), ",
            n_sr - n_sr_pass - n_extra, " counted but filtered (real jaccard ",
            "here), ", n_extra, " never entered the count matrix)  ->  ", out)
    invisible(dt)
}

#' Retrofit pValue/padjust into an already-written private-novel file, for
#' rows that turned out to be genuinely tested (passed=TRUE).
#'
#' write_private_novel_table() necessarily runs BEFORE FRASER() fits the
#' model (it needs to see the junctions the filter is about to reject, and
#' losing it to a crash during the fit would be worse than a stale one), so
#' `passed=TRUE` rows are written with pValue/padjust = NA even though a real
#' beta-binomial test now exists for them (see rows_from_fds()'s doc comment).
#' Call this once per output file AFTER `fds <- FRASER(...)` -- `fds` must be
#' the fitted, filtered object (the one write_results_tables() itself uses),
#' so a coordinate match against its rowRanges lands on the SAME row/column
#' the fit actually used. Mirrors extract_gene_results.R's --recalc-padj
#' philosophy: compute once, patch cheaply afterward, no refit.
#'
#' `!passed` rows are untouched by design -- they have no model p-value, ever
#' (no dispersion to test a single carrier against), which is the entire
#' reason this table exists instead of pointing at results_all.
enrich_private_novel_pvalues <- function(fds, path) {
    if (!file.exists(path) || file.size(path) == 0) return(invisible(NULL))
    # Explicit colClasses off the schema, not fread's own guess: an all-NA
    # column (pValue/padjust, before this function ever runs) reads back as
    # `logical` by default, and assigning a real double into a logical
    # data.table column SILENTLY coerces it to TRUE/FALSE rather than erroring
    # -- caught by the unit test, not by inspection. Wrong on any all-NA
    # column, not just these two, so the whole schema is pinned here.
    dt <- fread(path, showProgress = FALSE,
               colClasses = vapply(private_novel_empty(), function(x) class(x)[1],
                                   character(1)))
    if (!nrow(dt) || !"passed" %in% names(dt) || !any(dt$passed, na.rm = TRUE))
        return(invisible(dt))
    sel  <- which(dt$passed)
    jr   <- rowRanges(fds, type = "j")
    fkey <- paste(seqnames(jr), start(jr), end(jr), sep = ":")
    dkey <- paste(dt$seqnames[sel], dt$start[sel], dt$end[sel], sep = ":")
    row  <- match(dkey, fkey)
    col  <- match(dt$sampleID[sel], samples(fds))
    ok   <- !is.na(row) & !is.na(col)
    if (any(ok)) {
        idx <- cbind(row[ok], col[ok])
        dt[sel[ok], `:=`(
            pValue  = assay(fds, "pvaluesBetaBinomial_jaccard")[idx],
            padjust = assay(fds, "padjBetaBinomial_jaccard")[idx])]
        fwrite(dt, path, sep = "\t")
    }
    message("private-novel: patched pValue/padjust for ", sum(ok), "/", length(sel),
            " passed=TRUE row(s) in ", basename(path),
            if (any(!ok)) paste0(" (", sum(!ok), " unmatched against the post-fit ",
                                 "fds -- kept NA, check for a re-annotation or ",
                                 "coordinate drift)") else "")
    invisible(dt)
}

#' Where one reporting scope's tables go.
#'
#' Default (no --out-dir) is the historical layout: the cohort tables beside the
#' saved objects in the workdir, a per-sample scope in its own samples/<ID>/.
#' With --out-dir everything lands in ONE flat directory instead -- every file
#' name is already prefixed with its scope (<ID>_results.tsv, <name>_*), so a
#' single directory can hold the deliverables for any number of samples without
#' collision, and nothing downstream has to walk a per-sample subdirectory to
#' find them.
res_dir <- function(scope) {
    if (!is.null(OUTDIR)) return(OUTDIR)
    if (is.null(scope)) WD else file.path(WD, "samples", scope)
}

#' File-name prefix for a COHORT-scope table.
#'
#' Per-sample tables carry their sample ID, so a flat --out-dir can hold any
#' number of them. Cohort-scope tables were named <name>_* -- `cohort_*` for
#' every run -- and the documented per-tissue workflow points many concurrent
#' samples at ONE --out-dir. Two cases running against the same frozen cohort
#' therefore wrote different content to the same path: last writer won, readers
#' mid-write saw a truncated file, and whatever survived was silently attributed
#' to the wrong case (its private-novel context is 22 controls + THAT case).
#' Under --out-dir the cohort tables are prefixed with the reporting scope, so
#' they stay distinct and stay attributable. The workdir layout is unchanged --
#' there each run owns its own directory already.
cohort_prefix <- function() {
    if (!is.null(OUTDIR) && !is.null(FOCUS))
        paste0(paste(FOCUS, collapse = "-"), "_", NAME)
    else NAME
}

# ---------------------------------------------------------------- fitting ----
if (mode == "fit") {
    fds <- loadFraserDataSet(WD, NAME)
    message("=== FIT: ", ncol(fds), " samples, ", nrow(fds), " junctions ===")

    # The fit does no counting, but it does choose the testing threshold, and
    # --min-count is pinned. A fit run at a different value than the cohort was
    # built with silently changes which junctions are tested and therefore every
    # p value and FDR -- it landed only in runs.tsv, where nobody compares it.
    # Scoped to min_count on purpose: every other pinned row is a counting
    # parameter this run never uses (see check_run_params).
    check_run_params(only = "min_count", critical = "min_count", record = FALSE)

    # A fit must start from the COUNTS object, never from another fit's output.
    # `passedExpression` is written by filterExpressionAndVariability and is
    # absent from anything build/add produced, so it identifies that state
    # exactly. Fitting an already-filtered object is not merely redundant: the
    # private-novel table is derived from the junctions the filter REJECTS, so
    # on a pre-filtered input it is silently empty -- which is how job 30219333
    # replaced 1,395 private-novel junctions with a bare header.
    if ("passedExpression" %in% colnames(mcols(fds, type = "j")) && !ALLOWFILT) {
        stop("savedObjects/", NAME, " already carries filter results, so it is ",
             "the OUTPUT of an earlier fit, not a counts object.\n",
             "  This workdir predates the counts/fit split. Re-run `add` with ",
             "the same samplesheet to rebuild the full junction set -- it needs ",
             "no BAM passes, the non-split cache is reused -- then fit again.\n",
             "  --allow-filtered-input fits it anyway; the private-novel table ",
             "will be empty.")
    }
    if (!is.null(FOCUS)) {
        miss <- setdiff(FOCUS, samples(fds))
        if (length(miss))
            stop("--results-samples names sample(s) not in this cohort: ",
                 paste(miss, collapse = ", "))
        message("reporting on ", length(FOCUS), " of ", ncol(fds), " sample(s): ",
                paste(FOCUS, collapse = ", "))
    }
    # The private-novel table is written well before write_results_tables, so
    # the destination has to exist by now.
    dir.create(res_dir(NULL), recursive = TRUE, showWarnings = FALSE)
    if (!is.null(OUTDIR)) {
        message("results -> ", normalizePath(OUTDIR))
        # With a scope the cohort tables are prefixed with it (cohort_prefix)
        # and concurrent runs cannot collide. Without one there is nothing to
        # disambiguate them with, and --out-dir is documented as a directory
        # many samples share.
        if (is.null(FOCUS))
            warning("--out-dir without --results-samples: the cohort-wide ",
                    "tables are named ", NAME, "_* with nothing to tell one ",
                    "run's from another's. If more than one sample writes to ",
                    "this directory they WILL overwrite each other. Pass ",
                    "--results-samples, or give each run its own --out-dir.",
                    immediate. = TRUE)
    }
    N_IN <- nrow(fds)

    # Everything the fit writes from here on goes to savedObjects/<name>_fit.
    # This is deliberately BEFORE calculatePSIValues rather than at the closing
    # save: FRASER() persists assays to disk as it computes them, and those
    # writes -- not the final save -- are what truncated the counts object's
    # rawCountsJ.h5 to the filtered junctions.
    name(fds) <- FITNAME
    message("counts object savedObjects/", NAME, " is read-only for this run; ",
            "the fit writes savedObjects/", FITNAME)

    fds <- calculatePSIValues(fds)

    # Annotate the filter WITHOUT subsetting, so the junctions it rejects are
    # still reachable for the private-novel table, then subset by hand using
    # FRASER's own idiom (applyExpressionFilters_jaccard / applyVariabilityFilters
    # both do fds[<flag>, by="psi5"]). Equivalent to filter=TRUE: passedExpression
    # and passedVariability are per-row statistics over a junction's own counts,
    # so computing variability on the full set gives the same values as computing
    # it on the expression-filtered subset.
    fds    <- filterExpressionAndVariability(fds, minExpressionInOneSample = MINCNT,
                                             quantile = FQUANT,
                                             quantileMinExpression = QMINEXPR,
                                             filter = FALSE)
    passed <- mcols(fds, type = "j")$passed
    message("junctions passing filter (rowMaxs(K) >= ", MINCNT,
            if (QMINEXPR > 0)
                paste0(" & rowQuantiles(N, ", FQUANT, ") >= ", QMINEXPR)
            else "  [quantile term DISABLED, --quantile-min-expression 0]",
            "): ", sum(passed), " of ", length(passed))
    if (QMINEXPR <= 0)
        warning("--quantile-min-expression 0: the cohort-coverage half of the ",
                "expression filter is OFF. Single-carrier junctions whose splice ",
                "sites no other sample covers now enter the fit and are tested ",
                "against a column of zeros -- including STRUCTURAL zeros (sites ",
                "never counted, col4=0 junctions whose strand was never ",
                "resolved), which the model cannot distinguish from real ",
                "absence. Expect many more tests, a heavier BY correction, and ",
                "a class of 'outliers' that are artefacts of the frozen ",
                "universe. Results are NOT comparable to a run at the default 10.",
                immediate. = TRUE)

    priv <- NULL
    if (PRIVTAB) {
        # Resolved here, independently of the fds's own (later) annotation --
        # see annotate_gene_symbol()'s doc comment for why: this table is
        # written before fds gets annotated, and the strand-recovered private
        # rows never enter fds at all. build_txdb_and_map() is cheap when
        # reusing an existing cache ("Reusing txdb: ..."), so resolving it a
        # second time here just reloads it, it does not rebuild anything.
        PRIV_TXDB_USE <- TXDB
        if (is.null(PRIV_TXDB_USE) && is.null(GTF)) {
            auto <- find_txdb(WD)
            if (!is.na(auto)) PRIV_TXDB_USE <- auto
        }
        PRIV_ANNO <- NULL
        if (!is.null(GTF) || !is.null(PRIV_TXDB_USE))
            PRIV_ANNO <- build_txdb_and_map(gtf = GTF, dir = WD, rebuild = REBUILD,
                                            txdb_file = PRIV_TXDB_USE, map_file = GMAP)
        priv <- write_private_novel_table(fds, passed, UNIVERSE,
                                  file.path(res_dir(NULL),
                                            paste0(cohort_prefix(), "_private_novel.tsv")),
                                  min_count = PRIVMIN, max_other = PRIVOTH,
                                  anno = PRIV_ANNO,
                                  extra_file = file.path(WD, "strand_recovered_private.tsv"))
    }

    # A private-novel junction is by definition carried by exactly one sample,
    # so the per-sample table is a FILTER of the cohort one, not a second
    # computation. The cohort-wide table stays on disk regardless: it is the
    # context you need to judge whether a junction is really private, and it is
    # already paid for.
    if (PRIVTAB && !is.null(FOCUS)) {
        for (sid in FOCUS) {
            d <- res_dir(sid)
            dir.create(d, recursive = TRUE, showWarnings = FALSE)
            f   <- file.path(d, paste0(sid, "_private_novel.tsv"))
            sub <- if (is.null(priv)) private_novel_empty() else priv[sampleID == sid]
            fwrite(sub, f, sep = "\t")
            message("private novel for ", sid, ": ", nrow(sub), "  ->  ", f)
        }
    }
    rm(priv); invisible(gc())

    fds <- fds[passed, by = "psi5"]

    # Gene annotation MUST run BEFORE FRASER(), for two independent reasons:
    #  (1) FRASER() calls calculatePadjValues() internally, which derives the
    #      GENE-level p values (pValueGene/padjustGene) by aggregating the
    #      junction p values per hgnc_symbol. Without that mcol it only WARNS
    #      ("Gene-level pvalues could not be calculated as column hgnc_symbol
    #      was not found") and silently skips them -- annotating after the fit
    #      leaves the gene stats permanently absent unless calculatePadjValues
    #      is re-run by hand.
    #  (2) results() reads hgncSymbol and annotatedJunction off the fds at
    #      results() time, so annotating any later would not put them in the
    #      tables either.
    # Placed AFTER the filter subset so we annotate the ~200k passing junctions
    # rather than the full ~1.7M. Nothing here depends on the fit:
    # annotateRangesWithTxDb reads rowRanges only, and
    # annotateIntronReferenceOverlap reads rowRanges + K(fds, type="j").
    # Conversely the fit never clobbers the annotation: calculatePSIValues only
    # READS rowData, and fit/calculatePvalues write assays, so the mcols survive.
    # A rebuild re-reads the GTF, so it cannot run off a prebuilt txdb alone.
    # Caught here rather than inside build_txdb_and_map(), which would fail
    # obscurely on basename(NULL) after the annotation had already started.
    if (isTRUE(REBUILD) && is.null(GTF))
        stop("--rebuild-txdb needs --gtf (rebuilding parses the GTF)")

    # Explicit flags win; otherwise fall back to the shared cache. Auto-discovery
    # only fires when NEITHER flag is given, so --gtf keeps its exact old meaning.
    TXDB_USE <- TXDB
    if (is.null(TXDB_USE) && is.null(GTF)) {
        auto <- find_txdb(WD)
        if (!is.na(auto)) {
            message("no --txdb/--gtf given; using default annotation: ", auto)
            TXDB_USE <- auto
        }
    }

    HAVE_ANNO <- !is.null(GTF) || !is.null(TXDB_USE)
    if (HAVE_ANNO) {
        message("\n=== gene annotation ===")
        anno <- build_txdb_and_map(gtf = GTF, dir = WD, rebuild = REBUILD,
                                   txdb_file = TXDB_USE, map_file = GMAP)
        fds  <- annotate_fds_genes(fds, anno$txdb, anno$map)
    } else {
        message("\nNOTE: no --gtf/--txdb given and no prebuilt txdb found at\n    ",
                path.expand(DEFAULT_TXDB), "\n  -> results will have no hgncSymbol, no ",
                "annotatedJunction and NO gene-level p values.\n  ",
                "parse_fraser.R needs hgncSymbol; annotatePotentialImpact() needs ",
                "annotatedJunction.\n  Build the shared cache once with ",
                "build_fraser_txdb.R --gtf <annotation.gtf>,\n  or pass --txdb / ",
                "--gtf, or set $FRASER_TXDB, or retrofit with annotate_fraser_genes.R.")
    }

    # q: --q overrides everything. Else useOHT analytic q (skips the slow grid
    # search) as a starting point, with a fallback since OHT can degenerate on
    # small/low-rank matrices. If --q-search is also given (below), its winner
    # then REPLACES this OHT-based q as the one the fit actually uses.
    # NB (FRASER 2.6.x API): estimateBestQ() RETURNS the fds (q stored in its
    # hyperParams slot), NOT an integer -- extract with bestQ(). And it defaults
    # plot=TRUE, which renders a ggplot via plotEncDimSearch and CRASHES in a
    # headless batch job ("no method for coercing this S4 class to a vector") ->
    # must pass plot=FALSE.
    if (!is.null(QOPT)) {
        q <- as.integer(QOPT)
        message("using q override: ", q)
    } else {
        q <- tryCatch({
            fds <- estimateBestQ(fds, type = "jaccard", useOHT = TRUE, plot = FALSE)
            bestQ(fds, type = "jaccard")
        }, error = function(e) {
            qf <- max(2L, min(10L, ncol(fds) - 1L))
            message("estimateBestQ(useOHT) failed (", conditionMessage(e),
                    "); falling back to q=", qf)
            qf
        })
        message("estimated q (OHT): ", q)
    }

    # --q-search: runs the injected-outlier AROC grid, i.e. estimateBestQ()
    # with useOHT=FALSE. useOHT=TRUE (used above by default) is the fast
    # analytic Gavish-Donoho estimate; it exists specifically to skip the
    # slower injected-outlier grid search. That grid search is STILL in
    # FRASER 2.6.x, under the same estimateBestQ() function with
    # useOHT=FALSE -- verified by reading the installed package source,
    # since this is what classic FRASER's optimHyperParams() became. When
    # --q was NOT given (no explicit override), this block's winner REPLACES
    # the OHT-based `q` above as the one the fit actually uses -- that's the
    # whole point of asking for the slower, more thorough search. An explicit
    # --q still wins over both. Operates on a throwaway copy (estimateBestQ
    # (useOHT=FALSE) internally calls injectOutliers() on a local fds_copy,
    # but only writes hyperParams onto the object it RETURNS -- capturing
    # that into fds_q, not reassigning fds, keeps the real counts/assays fds
    # untouched).
    if (QSEARCH) {
        message("\n=== q search: injected-outlier AROC grid (useOHT=FALSE) ===")
        save_enc_dim_plot <- function(path, fds_, plotType) {
            png(path, width = 1200, height = 800, res = 150)
            on.exit(dev.off(), add = TRUE)
            print(plotEncDimSearch(fds_, type = "jaccard", plotType = plotType))
        }
        tryCatch({
            fds_q <- estimateBestQ(fds, type = "jaccard", useOHT = FALSE,
                                   plot = FALSE, BPPARAM = BPP)
            # hyperParams()'s GETTER already reduces to the single winning row
            # by default (`ans[aroc == max(aroc)][1]` when an aroc column is
            # present, `ans[oht == TRUE][1]` for the OHT path) -- verified in
            # source AND empirically on a toy fds: with all=FALSE this returns
            # exactly one row regardless of how many q's were actually tried.
            # bestQ() (`hyperParams(fds,type)[1,q]`) is therefore correct as
            # written for BOTH branches; it is filtering on an already-reduced
            # table, not picking an arbitrary row 1 out of the full grid. The
            # full sweep -- what this diagnostic exists to write out -- needs
            # `all=TRUE` explicitly, confirmed to return every candidate q on
            # the same toy fds where all=FALSE gave only the winner.
            grid   <- hyperParams(fds_q, type = "jaccard", all = TRUE)
            best_q <- bestQ(fds_q, type = "jaccard")
            best   <- grid[q == best_q][1]
            tbl <- file.path(res_dir(NULL), paste0(cohort_prefix(), "_q_search.tsv"))
            fwrite(grid[order(q)], tbl, sep = "\t")
            plotpath <- file.path(res_dir(NULL), paste0(cohort_prefix(), "_q_search.png"))
            tryCatch(save_enc_dim_plot(plotpath, fds_q, "auc"),
                     error = function(e) message("q search plot failed: ",
                                                  conditionMessage(e)))
            if (is.null(QOPT)) {
                message("q search: best q by AROC = ", best_q, " (aroc=",
                        signif(best$aroc, 4), ") -> USING for the fit ",
                        "(was q=", q, " from OHT)  ->  ", tbl, ", ", plotpath)
                q <- best_q
            } else {
                message("q search: best q by AROC = ", best_q, " (aroc=",
                        signif(best$aroc, 4), "); fit still uses --q ",
                        "override q=", q, "  ->  ", tbl, ", ", plotpath)
            }
        }, error = function(e) {
            message("q search failed (", conditionMessage(e),
                    ") -- continuing the fit with q=", q)
        })
        # Free here too, same reasoning as OHT: the OHT spectrum was already
        # computed as a side effect of the block above whenever QOPT was not
        # given, i.e. exactly the case --q-search is meant to be paired with.
        if (is.null(QOPT)) {
            svpath <- file.path(res_dir(NULL), paste0(cohort_prefix(), "_q_oht_spectrum.png"))
            tryCatch(save_enc_dim_plot(svpath, fds, "sv"),
                     error = function(e) message("OHT spectrum plot failed: ",
                                                  conditionMessage(e)))
            message("OHT spectrum plot -> ", svpath)
        }
    }

    fds <- FRASER(fds, q = c(jaccard = q), implementation = "PCA", iterations = 15,
                  BPPARAM = BPP)

    fds <- saveFraserDataSet(fds)

    # Now that a real p value exists for anything that turned out testable,
    # patch it into the private-novel files written earlier (pre-fit, for
    # crash-safety -- see write_private_novel_table's doc comment). Cheap: it
    # only rereads/rewrites the small private-novel tables, never the fds.
    if (PRIVTAB) {
        enrich_private_novel_pvalues(fds, file.path(res_dir(NULL),
                                     paste0(cohort_prefix(), "_private_novel.tsv")))
        if (!is.null(FOCUS))
            for (sid in FOCUS)
                enrich_private_novel_pvalues(fds, file.path(res_dir(sid),
                                             paste0(sid, "_private_novel.tsv")))
    }

    #' Write the results tables for one reporting scope.
    #'
    #' scope = NULL  -> the whole cohort, into WD/<name>_*        (old behaviour)
    #' scope = "SID" -> that sample alone, into WD/samples/SID/SID_*
    #'
    #' Per-sample is the point of --results-samples: after an `add` only the new
    #' sample needs reporting, and results(sampleIDs=) makes every pass smaller
    #' by a factor of nSamples without changing a single number (see the note on
    #' FOCUS above). It also keeps each sample's deliverables in their own
    #' directory, so a later run for a different sample cannot overwrite them.
    write_results_tables <- function(fds, scope = NULL) {
        dir <- res_dir(scope)
        if (is.null(scope)) {
            pre <- cohort_prefix(); ids <- samples(fds)
        } else {
            pre <- scope; ids <- scope
        }
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        pth <- function(suffix) file.path(dir, paste0(pre, suffix))

        # Two junction tables, matching DROP: the unfiltered one is what
        # parse_fraser.R takes (its nominal-p gate needs the sub-threshold rows
        # that all=FALSE discards).
        allOut  <- pth("_results_all.tsv.gz")
        res_all <- as.data.table(results(fds, sampleIDs = ids,
                                         aggregate = FALSE, all = TRUE))
        fwrite(res_all, allOut, sep = "\t")
        message(sprintf("all junctions: %d  ->  %s", nrow(res_all), allOut))
        # Free it before anything else runs: at cohort scope this is ~5M rows x
        # 20 cols, and leaving it resident is half of why job 30219333 went past
        # 64G once the gene passes started.
        rm(res_all); invisible(gc())

        out <- pth("_results.tsv")
        res <- results(fds, sampleIDs = ids, aggregate = FALSE, all = FALSE,
                       deltaPsiCutoff = DPSI)
        fwrite(as.data.table(res), out, sep = "\t")
        n_out <- length(res); rm(res); invisible(gc())

        # Gene-level tables, mirroring DROP's resultTableGene_full/_aberrant
        # (08_extract_results_FraseR.R:70). aggregate=TRUE is what actually EMITS
        # the pValueGene/padjustGene columns computed during the fit -- with
        # aggregate=FALSE they are never written to the table
        # (AllGenerics.R:823). Each row is still a real junction: FRASER picks
        # the gene's most aberrant one (lowest pValue, ties by largest
        # |deltaPsi|), so the junction columns stay interpretable.
        # collapse=FALSE keeps one row per (gene, sample). Only possible when
        # annotation ran -- without hgnc_symbol there are no gene-level p values.
        geneOut <- NA_character_; n_gene <- NA_integer_
        if (HAVE_ANNO && GENETAB) {
            # Aberrant table FIRST: it is small and it is what gets triaged, so
            # writing it before the heavy all=TRUE pass means an OOM there
            # cannot cost it. (In job 30219333 the order of the two junction
            # tables saved them for exactly this reason.)
            geneOut  <- pth("_results_gene.tsv")
            res_gene <- results(fds, sampleIDs = ids, aggregate = TRUE,
                                collapse = FALSE, all = FALSE,
                                deltaPsiCutoff = DPSI)
            fwrite(as.data.table(res_gene), geneOut, sep = "\t")
            n_gene <- length(res_gene)
            rm(res_gene); invisible(gc())

            # Opt-in (--gene-all-table): aggregate=TRUE + all=TRUE materialises
            # dense gene x sample p-value matrices AND re-sorts the exploded
            # gene-junction table inside every forked worker. This is the call
            # that OOM-killed a worker at 64G AT COHORT SCOPE; restricted to one
            # sample it is a fraction of that, but give it ~96G, or fewer
            # --threads, before enabling it for the whole cohort.
            if (GENEALL) {
                geneAllOut   <- pth("_results_gene_all.tsv.gz")
                res_gene_all <- as.data.table(results(fds, sampleIDs = ids,
                                    aggregate = TRUE, collapse = FALSE, all = TRUE))
                fwrite(res_gene_all, geneAllOut, sep = "\t")
                message(sprintf("all genes: %d  ->  %s", nrow(res_gene_all), geneAllOut))
                rm(res_gene_all); invisible(gc())
            }
        }
        list(scope = if (is.null(scope)) NAME else scope, dir = dir,
             focused = !is.null(scope),
             results = out, all = allOut, gene = geneOut,
             n_outliers = n_out, n_gene = n_gene)
    }

    scopes <- if (is.null(FOCUS)) list(NULL) else as.list(FOCUS)
    summ   <- lapply(scopes, function(sc) write_results_tables(fds, sc))

    rec <- base_record(mode, ncol(fds),
                       junctions_in = N_IN, junctions_passed = sum(passed),
                       q = q, filter_quantile = FQUANT,
                       quantile_min_expression = QMINEXPR,
                       focus = if (is.null(FOCUS)) "<all>" else FOCUS)
    write_run_record(rec)
    # Same record beside each sample's tables: those directories get handed on
    # for triage, and a p value is only interpretable against the cohort and
    # parameters that produced it. Every add refits, so this moves.
    # Keyed on the SCOPE, not on the directory: under --out-dir the cohort and
    # per-sample tables share one directory, so comparing paths to WD would
    # either skip the per-sample record or emit a cohort-scoped one.
    for (x in summ)
        if (isTRUE(x$focused))
            fwrite(rec, file.path(x$dir, paste0(x$scope, "_fit_info.tsv")), sep = "\t")

    message("\n=== done ===")
    for (x in summ) {
        message("[", x$scope, "] outliers: ", x$n_outliers, "  ->  ", x$results)
        if (HAVE_ANNO) {
            message("[", x$scope, "] gene-level outliers: ", x$n_gene, "  ->  ", x$gene)
            message("[", x$scope, "] feed ", basename(x$all), " to parse_fraser.R with -r")
        }
    }
}
