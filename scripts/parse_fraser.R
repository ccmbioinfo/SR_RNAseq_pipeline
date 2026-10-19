#!/usr/bin/env Rscript
# Parse FRASER2 results for a single individual. Parallel to parse_outrider.R.
#
# Steps:
#   1. Read the FRASER2 results table (.tsv, or an .rds holding the same
#      data.frame) and extract the rows for the input sample.
#   2. Clean the messy hgncSymbol column: FRASER/DROP joins each junction to
#      overlapping genes as a ';'-separated list with make-unique '_N' suffixes
#      and 'NA' fillers (e.g. 'GUSBP11_2;GUSBP11', 'DEPDC5;NA_34608'). Strip the
#      suffix, drop NA/blank, uniquify, then EXPLODE to one row per
#      (junction x clean gene). Junctions with no real gene stay as one row with
#      gene = NA.
#   3. Add gene-level annotations (matched by gene SYMBOL) from the hg38 gene
#      annotation table. NA-gene rows get no annotation.
#   4. Read the individual's HPO terms and, per gene, report how many overlap the
#      gene's associated HPO terms (names), a yes/no flag, and a Jaccard score.
#   5. (optional) If a model_variants.tsv (-V) is supplied, add per-junction
#      COORDINATE-based rare-variant columns (no gene/ENSG match needed):
#        cand_var_n / cand_variants  variants within --cand-flank bp of either
#                                    splice site (start or end) of the junction
#        near_var_n / near_variants  variants within --near-flank bp of either
#                                    splice site, excluding the cand set
#        cand_models                 inheritance model(s) among cand variants
#        cand_max_impact             most severe VEP IMPACT among cand variants
#      Each variant cell carries a dSS=<bp> distance-to-nearest-splice-site token.
#   5b. (optional) If --fds-file/--fds-dir is given, deltaPsi_raw's cohort
#       row-median psi is read directly off the fds's own jaccard assay
#       (every sample the fit ever saw) instead of from -r. Without it, -r
#       must itself be a multi-sample export (e.g. a cohort-wide
#       cohort_results_all.tsv.gz) or deltaPsi_raw is 0 for every row -- a
#       --results-samples-scoped fit output (one row per junction, one
#       sample) has no cohort to take a median over on its own.
#   6. Write the annotated, sorted table, as two files off one -o path:
#        <base>.tsv.gz        every row, gzipped (archive / audit copy)
#        <base>.filtered.tsv  the shortlist (see "filtered subset" below)
#      Pass -o <name>.xlsx to ALSO get the shortlist as a workbook (frozen header
#      + autofilter), produced by tsv_to_xlsx.py so no Excel R package is needed.
#      Default output is .tsv -- no workbook unless you ask for one.
#
# Usage:
#   Rscript parse_fraser.R -r <fraser_results.tsv> -s <sampleID> \
#       -p <hpo_terms.txt> -a <gene_annotation.tsv> [-o <output.xlsx>] \
#       [-V <model_variants.tsv> --cand-flank 10 --near-flank 250 --topk 20]
#   Rscript parse_fraser.R --help

suppressMessages({
  library(tidyverse)
  library(optparse)
})

## ---- arguments -------------------------------------------------------------
option_list <- list(
  # -r/-s are REQUIRED and deliberately have no default: they used to default to
  # the chr22 test fixtures, so a batch loop with an empty variable would quietly
  # emit a plausible-looking file for the wrong sample instead of failing.
  make_option(c("-r", "--results"), type = "character", default = NULL,
              help = "REQUIRED. FRASER2 results: a .tsv/.tsv.gz, or an .rds holding the same data.frame"),
  make_option(c("-s", "--sample"), type = "character", default = NULL,
              help = "REQUIRED. Sample ID to extract (matches the sampleID column)"),
  make_option(c("-p", "--hpo"), type = "character", default = NULL,
              help = "Individual's HPO terms file (one HP:xxxxxxx id per line); omit to leave the HPO columns NA"),
  make_option(c("-a", "--annotation"), type = "character",
              default = "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_annotation/gtf115_annotation_table_Nov3.tsv",
              help = "Gene-level annotation table (.tsv, joined on gene_symbol) [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output path; .tsv (default) or .xlsx, which additionally converts the filtered table via tsv_to_xlsx.py [default: <sample>_FRASER_results_annotated.tsv]"),
  make_option(c("-V", "--variants"), type = "character", default = NULL,
              help = "model_variants.tsv (VEP-expanded, model-tagged variant table); if given, per-junction coordinate-based variant columns are added [default: none]"),
  make_option(c("-N", "--private-novel"), type = "character", default = NULL,
              help = "fraser_cohort.R's <name>_private_novel.tsv (or the per-sample <ID>_private_novel.tsv); if given, the SAME gene annotation (-a), blacklist, artifact classification, HPO overlap (-p) and variant proximity (-V) applied to the filtered junctions are attached to it and written to --private-novel-out. Independent of -r/-s -- this table has its own schema (one clean gene per row already, mostly no pValue) [default: none]"),
  make_option("--private-novel-out", type = "character", default = NULL,
              help = "Output path for the annotated private-novel table [default: <basename of -N minus .tsv>_annotated.tsv]"),
  make_option(c("--cand-flank"), type = "double", default = 10,
              help = "bp flank around EACH splice site for 'candidate' variants [default: %default]"),
  make_option(c("--near-flank"), type = "double", default = 250,
              help = "bp flank around EACH splice site for 'nearby' variants [default: %default]"),
  make_option(c("--topk"), type = "integer", default = 20,
              help = "max variants listed per cell [default: %default]"),
  make_option("--tier-p", type = "double", default = 1e-3,
              help = "pass_tier outlier gate: max pValue [default: %default]"),
  make_option("--tier-dpsi", type = "double", default = 0.1,
              help = "pass_tier outlier gate: min |deltaPsi| [default: %default]"),
  make_option("--tier-sai", type = "double", default = 0.2,
              help = "pass_tier scored gate: min SpliceAI [default: %default]"),
  make_option("--tier-sai-strong", type = "double", default = 0.5,
              help = "pass_tier concordant gate: SpliceAI that alone counts as splice-disrupting [default: %default]"),
  make_option("--tier-cadd", type = "double", default = 20,
              help = "pass_tier scored gate: min CADD [default: %default]"),
  make_option("--tier-revel", type = "double", default = 0.5,
              help = "pass_tier scored gate: min REVEL [default: %default]"),
  make_option("--artifact-biotypes", type = "character",
              default = paste0("^(IGH|IGK|IGL|TRA|TRB|TRG|TRD)[VDJ]([0-9]|OR|[IVX]+-)|",
                               "^(IGHA[12]|IGHG[1-4]|IGHE|IGHM|IGHD|IGKC|IGLC[1-7]|",
                               "TRAC|TRBC[12]|TRGC[12]|TRDC)$"),
              help = "regex on gene SYMBOL (checked against every gene in the junction's un-exploded hgncSymbol_all, not just the row's own exploded hgncSymbol) flagged pass_tier=potential_artifact and dropped from the filtered output. Matches V/D/J segment genes (locus prefix + V/D/J + digit/OR/roman-numeral-family, e.g. IGKV3-20, TRBV20-1, IGLVI-70, IGLVVI-22-1) and named constant genes (e.g. IGHG1, TRAC) -- i.e. immunoglobulin/TCR genes, somatically rearranged in lymphoid tissue. Deliberately narrower than a bare locus-prefix match, which would also catch unrelated real genes like TRAF1/TRAK1/TRAPPC1/TRDN/TRDMT1/IGHMBP2/IGLL1 [default: %default]"),
  make_option("--artifact-genes", type = "character", default = "^(HLA-|MIC[A-G]$)",
              help = "regex on gene symbol flagged pass_tier=potential_artifact: classical MHC class I/II genes, which span several biotypes so they need a symbol match rather than a biotype match [default: %default]"),
  make_option("--dark-genes", type = "character", default = NULL,
              help = "Alliance dark-gene panel BED (chrom/start/end/comma-separated gene names); default: annotation_files/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed, or $DARK_GENES. Use 'none' to skip."),
  make_option("--no-mhc-by-position", action = "store_true", default = FALSE,
              help = "do NOT use the dark panel's HLA intervals to flag MHC mapping artifacts positionally; fall back to the --artifact-genes symbol match alone"),
  make_option("--keep-nogene", action = "store_true", default = FALSE,
              help = "keep junctions with no assigned gene in the filtered output [default: drop them]"),
  make_option("--lenient-p", type = "double", default = 0.01,
              help = "rescue gate: max pValue for HPO/candidate-variant rows [default: %default]"),
  make_option("--lenient-dpsi", type = "double", default = 0.05,
              help = "rescue gate: min |deltaPsi| for HPO/candidate-variant rows [default: %default]"),
  make_option("--blacklist", type = "character", default = NULL,
              help = "ENCODE blacklist BED (.bed/.bed.gz); default: <assembly>-blacklist.v2.bed.gz in blacklist_regions/ beside this script, or $FRASER_BLACKLIST. Use 'none' to skip."),
  make_option("--assembly", type = "character", default = "hg38",
              help = "genome assembly for the default blacklist BED (hg19/hg38) [default: %default]"),
  make_option("--blacklist-minoverlap", type = "integer", default = 5,
              help = "bp around each splice site that must overlap a blacklist region to flag it; matches FRASER's flagBlacklistRegions [default: %default]"),
  # --- optional fds-backed features (need FRASER; --fds-* alone or with --txdb) --
  make_option("--txdb", type = "character", default = NULL,
              help = "txdb.db; with an --fds-* option, adds potentialImpact/causesFrameshift/UTR_overlap to the FILTERED output. Omit to skip (keeps deps to tidyverse+optparse)."),
  make_option("--fds-file", type = "character", default = NULL,
              help = "path to fds-object.RDS, e.g. <WD>/savedObjects/cohort_fit/fds-object.RDS. Used alone, this makes deltaPsi_raw's cohort median come from the fds's own jaccard assay (every sample the fit ever saw) instead of from -r -- the only way to get a real deltaPsi_raw when -r is a single-sample (--results-samples-scoped) export. Add --txdb too for splicing-consequence annotation as well."),
  make_option("--fds-dir", type = "character", default = NULL,
              help = "DROP 'datasets' dir (use with --dataset/--annotation instead of --fds-file)"),
  make_option("--dataset", type = "character", default = NULL,
              help = "DROP dataset name, e.g. LCL"),
  make_option("--annotation-name", type = "character", default = NULL,
              help = "DROP annotation name, e.g. v115"),
  make_option("--impact-blacklist", action = "store_true", default = FALSE,
              help = "also add blacklist_fraser via FRASER's own flagBlacklistRegions, to cross-check the native 'blacklist' column")
)
opt <- parse_args(OptionParser(option_list = option_list))

# opt[["x"]] not opt$x: optparse omits NULL-default options, and $ partial-matches.
results_file   <- opt[["results"]]
sample         <- opt[["sample"]]
hpo_file       <- opt[["hpo"]]
gene_anno_file <- opt[["annotation"]]

missing_req <- c(if (is.null(results_file)) "-r/--results",
                 if (is.null(sample))       "-s/--sample")
if (length(missing_req))
  stop("missing required argument(s): ", paste(missing_req, collapse = ", "),
       "\n  usage: parse_fraser.R -r <fraser_results.tsv[.gz]> -s <sampleID> ",
       "[-p <hpo.txt>] [-V <model_variants.tsv>] [-o <out.tsv|out.xlsx>]",
       "\n  run with --help for the full option list.")
if (!file.exists(path.expand(results_file)))
  stop("results file not found: ", results_file)
# the Snakefile passes the literal "NULL" when no VCF is configured -- that means
# "not given", not a path, so drop it before the existence check below
if (!is.null(opt$variants) && toupper(opt$variants) %in% c("NULL", "NA", "")) opt$variants <- NULL
# -V given but unreadable used to be skipped silently, leaving every variant column NA
if (!is.null(opt$variants) && nzchar(opt$variants) && !file.exists(opt$variants))
  stop("-V/--variants file not found: ", opt$variants, " (relative paths resolve against ", getwd(), ")")
if (!nzchar(sample))
  stop("-s/--sample is empty -- check the variable in your batch loop")
outfile        <- if (!is.null(opt$output)) opt$output else paste0(sample, "_FRASER_results_annotated.tsv")

## ---- read FRASER results ---------------------------------------------------
# Either a plain results .tsv/.tsv.gz (the usual FRASER export) or an .rds
# holding the same data.frame. FRASER's results table is already one row per
# junction x sample, so no results()-style expansion is needed (unlike OUTRIDER).
if (grepl("\\.rds$", results_file, ignore.case = TRUE)) {
  res_all <- as.data.frame(readRDS(results_file))
} else {
  con <- if (grepl("\\.gz$", results_file, ignore.case = TRUE))
           gzfile(results_file) else results_file
  res_all <- read.delim(con, sep = "\t", header = TRUE, quote = "",
                        check.names = FALSE, stringsAsFactors = FALSE,
                        na.strings = c("NA", ""))
}

if (!sample %in% res_all$sampleID) {
  stop(sprintf("Sample '%s' not found in %s. Available samples:\n%s",
               sample, results_file, paste(unique(res_all$sampleID), collapse = ", ")))
}

## ---- optional early fds load (deltaPsi_raw's cohort median + later, -----
## ---- if --txdb is also given, splicing-consequence annotation) --------
# Resolved here rather than only at the splicing-consequence block further
# down, because deltaPsi_raw's cohort median can use the SAME fds independently
# of --txdb: the fds already holds every sample's psi, computed once at fit
# time, regardless of how many samples -r itself was scoped to
# (--results-samples writes a single-sample results_all.tsv.gz, but the fds
# behind it was fit on the whole cohort -- see fraser-private-strand-recovery
# memory). Sourced lazily so the default (no --fds-*) path still pulls in
# only tidyverse+optparse.
fds_file <- opt[["fds-file"]]
fds_dir  <- opt[["fds-dir"]]
have_fds <- !is.null(fds_file) || !is.null(fds_dir)
fds_early <- NULL
if (have_fds) {
  fia_helper <- file.path(dirname(normalizePath(sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE)[1]),
                   mustWork = FALSE)), "fraser_impact_annotation.R")
  if (!file.exists(fia_helper))
    stop("--fds-file/--fds-dir given but fraser_impact_annotation.R was not ",
         "found beside this script (looked at ", fia_helper, ")")
  source(fia_helper)
  fds_early <- load_fds_for_impact(fds_file, fds_dir, opt[["dataset"]],
                                   opt[["annotation-name"]])
}

## ---- cohort-wide row-median psi, for the simple (non-model) deltaPsi -------
# FRASER itself computes this pre-fit, per junction, as `psi - rowMedian(psi)`
# across ALL cohort samples (calculatePSIValue.R::calculateDeltaPsiValue) --
# unlike the model-based `deltaPsi` column (observed - autoencoder-fitted psi),
# this needs no fit, just the cohort's raw psiValue for each junction.
#
# TWO SOURCES, preferring the fds when given:
#  - fds (--fds-file/--fds-dir): assay(fds,"jaccard") is the psi/jaccard value
#    matrix for EVERY junction the fit counted x EVERY sample it was fit on --
#    always the true cohort, however -r itself was scoped.
#  - -r (the old default, unchanged): results(fds, all=TRUE) run with no
#    sampleIDs restriction carries every cohort sample's row per junction, so
#    the median can be taken from res_all directly, before subsetting to -s.
#    This collapses to a median-of-one (delta always 0) if -r itself is a
#    single-sample export -- e.g. a --results-samples-scoped fit output.
n_cohort <- length(unique(res_all$sampleID))
if (!is.null(fds_early)) {
  jr_all    <- rowRanges(fds_early, type = "j")
  jkey_fds  <- paste(seqnames(jr_all), start(jr_all), end(jr_all),
                     strand(jr_all), sep = "\r")
  psi_mat   <- as.matrix(SummarizedExperiment::assay(fds_early, "jaccard"))
  n_fds     <- ncol(psi_mat)
  psi_median <- stats::setNames(matrixStats::rowMedians(psi_mat, na.rm = TRUE),
                                jkey_fds)
  message(sprintf(
    "Row-median psi computed from the fds directly (--fds-file/--fds-dir), across %d cohort sample(s) -- independent of the %d sample(s) in -r.",
    n_fds, n_cohort))
} else {
  jkey_all <- paste(res_all$seqnames, res_all$start, res_all$end, res_all$strand,
                    sep = "\r")
  psi_median <- tapply(res_all$psiValue, jkey_all, median, na.rm = TRUE)
  if (n_cohort == 1) {
    warning("Only 1 sample present in -r and no --fds-file/--fds-dir given; ",
            "deltaPsi_raw will be 0 for every row (no cohort to take a median ",
            "over). Pass the full results_all.tsv, or --fds-file pointing at ",
            "the fitted fds, if you want this column to mean anything.")
  }
  message(sprintf("Row-median psi computed across %d cohort sample(s) (-r).",
                  n_cohort))
}

res <- res_all[res_all$sampleID == sample, , drop = FALSE]
message(sprintf("Extracted %d junction rows for sample '%s'.", nrow(res), sample))

## ---- clean hgncSymbol and explode to one row per (junction x gene) ---------
# Split on ';', strip the make-unique '_N' suffix, drop 'NA'/blank, uniquify.
clean_symbols <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character(0))
  parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  parts <- sub("_[0-9]+$", "", parts)            # make-unique suffix from DROP
  parts <- parts[parts != "" & parts != "NA"]
  unique(parts)
}
# Vectorized: does ANY gene in a junction's raw (un-exploded) hgncSymbol string
# match `pattern`? Used for the IG/TR artifact check -- testing the whole raw
# string (every gene the junction overlaps) rather than one exploded row's own
# single gene means the classification is a property of the JUNCTION, not of
# which gene a particular row happens to represent, and it doesn't need a
# successful gene_biotype join to work.
any_symbol_matches <- function(raw, pattern) {
  vapply(raw, function(x) any(grepl(pattern, clean_symbols(x))), logical(1))
}
# Keep a length-1 NA for junctions with no real gene so the row survives unnest.
res$gene <- lapply(res$hgncSymbol, function(x) {
  g <- clean_symbols(x); if (length(g) == 0) NA_character_ else g
})
res <- tidyr::unnest(res, gene)
res$gene <- as.character(res$gene)
message(sprintf("Exploded to %d (junction x gene) rows (%d with a named gene).",
                nrow(res), sum(!is.na(res$gene))))

# deltaPsi_raw: psiValue minus the cohort's row-median psi for that junction
# (computed above, before the sample subset). Exploding to one row per gene
# only duplicates rows -- the junction key (seqnames/start/end/strand) is
# unchanged, so the lookup is still correct per row.
jkey <- paste(res$seqnames, res$start, res$end, res$strand, sep = "\r")
res$deltaPsi_raw <- res$psiValue - psi_median[jkey]

## ---- gene-level annotation (join on gene_symbol) ---------------------------
#' Read + dedupe the shared gene-annotation table and left_join it onto `df`
#' by `df[[gene_col]]`. Factored out (2026-09-22) so the private-novel table
#' (--private-novel below) can attach the SAME columns via the SAME join,
#' instead of a second implementation.
join_gene_annotation <- function(df, gene_col, anno_file) {
  gene_anno <- read.delim(anno_file, sep = "\t", header = TRUE,
                          quote = "", stringsAsFactors = FALSE,
                          check.names = FALSE, na.strings = c("NA", ""))
  # Collapse to one annotation row per gene_symbol so the join stays 1:1.
  gene_anno <- gene_anno %>%
    filter(!is.na(gene_symbol) & gene_symbol != "") %>%
    group_by(gene_symbol) %>%
    summarise(
      across(c(RefSeq_nuc, Ensembl_nuc, MANE_status),
             ~ paste(unique(.x), collapse = ";")),
      across(-c(RefSeq_nuc, Ensembl_nuc, MANE_status),
             ~ .x[1]),
      .groups = "drop"
    )
  by <- setNames("gene_symbol", gene_col)
  left_join(df, gene_anno, by = by)
}

res <- join_gene_annotation(res, "gene", gene_anno_file)

## ---- individual HPO terms --------------------------------------------------
ind_hpo <- character(0)
if (!is.null(hpo_file) && nzchar(hpo_file) && file.exists(hpo_file)) {
  ind_hpo <- readLines(hpo_file, warn = FALSE)
  ind_hpo <- trimws(gsub("_x000D_|\r", "", ind_hpo))
  ind_hpo <- unique(ind_hpo[grepl("^HP:", ind_hpo)])
}
have_hpo <- length(ind_hpo) > 0

## ---- HPO overlap between individual and each gene --------------------------
#' Per-row HPO overlap against `ind_hpo` (the individual's own HP: id vector),
#' given the row's OWN `hpo_id`/`hpo_name` (already joined from the annotation
#' table, `|`-separated, order-aligned). Factored out (2026-09-22) so the
#' private-novel table reuses this identically. Returns a data.frame of the
#' four output columns, same length as hpo_id/hpo_name.
compute_hpo_overlap <- function(hpo_id, hpo_name, ind_hpo, label = "individual") {
  n <- length(hpo_id)
  if (!length(ind_hpo)) {
    message(sprintf("No HPO terms provided for '%s'; HPO overlap columns set to NA.", label))
    return(data.frame(overlapping_hpo_terms = NA_character_,
                      n_hpo_overlap = NA_integer_,
                      hpo_overlap = NA_character_,
                      hpo_jaccard = NA_real_)[rep(1, n), ])
  }
  message(sprintf("'%s' has %d HPO terms.", label, length(ind_hpo)))

  gene_hpo_ids   <- strsplit(ifelse(is.na(hpo_id),   "", hpo_id),   "\\|")
  gene_hpo_names <- strsplit(ifelse(is.na(hpo_name), "", hpo_name), "\\|")

  n_gene_hpo <- vapply(gene_hpo_ids,
                       function(g) length(unique(trimws(g[trimws(g) != ""]))),
                       integer(1))

  overlap_terms <- character(n)
  n_overlap     <- integer(n)
  for (i in seq_len(n)) {
    ids <- trimws(gene_hpo_ids[[i]])
    nms <- trimws(gene_hpo_names[[i]])
    k   <- min(length(ids), length(nms))
    if (k == 0) next
    ids <- ids[seq_len(k)]; nms <- nms[seq_len(k)]
    keep <- ids != "" & ids %in% ind_hpo & !duplicated(ids)
    if (any(keep)) {
      overlap_terms[i] <- paste(nms[keep], collapse = ";")
      n_overlap[i]     <- sum(keep)
    }
  }

  hpo_overlap <- ifelse(n_overlap > 0, "yes", "no")
  union_size  <- length(ind_hpo) + n_gene_hpo - n_overlap
  hpo_jaccard <- ifelse(union_size > 0, n_overlap / union_size, 0)
  data.frame(overlapping_hpo_terms = overlap_terms, n_hpo_overlap = n_overlap,
            hpo_overlap = hpo_overlap, hpo_jaccard = hpo_jaccard)
}

hpo_cols <- compute_hpo_overlap(res$hpo_id, res$hpo_name, ind_hpo, sample)
res$overlapping_hpo_terms <- hpo_cols$overlapping_hpo_terms
res$n_hpo_overlap         <- hpo_cols$n_hpo_overlap
res$hpo_overlap            <- hpo_cols$hpo_overlap
res$hpo_jaccard             <- hpo_cols$hpo_jaccard

## ---- rare-variant annotation (optional, coordinate-based) ------------------
# Both tiers are purely coordinate-based off the junction (no gene/ENSG match):
#   cand : within --cand-flank bp of EITHER splice site (start or end)
#   near : within --near-flank bp of EITHER splice site, minus the cand set
# Both windows are measured from the splice sites themselves, so the interior of
# a long intron is NOT swept in (that would be deep-intronic, not regulatory).
#
# Factored into compute_variant_columns() (2026-09-22) so the private-novel
# table (--private-novel below) gets the SAME variant proximity/scoring
# instead of a second implementation. Takes the DISTINCT junctions (one row
# per seqnames/start/end/strand) and returns them augmented with the 11
# variant/helper columns; the caller left_joins the result back onto its own
# (possibly gene-exploded) table by junction key.
VARIANT_COLS <- c("cand_var_n", "cand_variants", "near_var_n", "near_variants",
                  "cand_models", "cand_max_impact", "helper_sai", "helper_cadd",
                  "helper_revel", "helper_clinpath", "helper_splice")
empty_variant_cols <- function(n) data.frame(
  cand_var_n = rep(NA_integer_, n), cand_variants = rep(NA_character_, n),
  near_var_n = rep(NA_integer_, n), near_variants = rep(NA_character_, n),
  cand_models = rep(NA_character_, n), cand_max_impact = rep(NA_character_, n),
  helper_sai = rep(NA_real_, n), helper_cadd = rep(NA_real_, n),
  helper_revel = rep(NA_real_, n), helper_clinpath = rep(FALSE, n),
  helper_splice = rep(FALSE, n))

compute_variant_columns <- function(juncs, vt_file, cand_bp, near_bp, top_k) {
  win_bp <- max(near_bp, cand_bp)   # scan window; guards --cand-flank > --near-flank

  # VEP consequence severity, most -> least severe (rank = index).
  VEP_ORDER <- c(
    "transcript_ablation","splice_acceptor_variant","splice_donor_variant",
    "stop_gained","frameshift_variant","stop_lost","start_lost",
    "transcript_amplification","feature_elongation","feature_truncation",
    "inframe_insertion","inframe_deletion","missense_variant",
    "protein_altering_variant","splice_donor_5th_base_variant",
    "splice_region_variant","splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant","incomplete_terminal_codon_variant",
    "start_retained_variant","stop_retained_variant","synonymous_variant",
    "coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant",
    "3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant",
    "NMD_transcript_variant","non_coding_transcript_variant",
    "coding_transcript_variant","upstream_gene_variant","downstream_gene_variant",
    "TFBS_ablation","TFBS_amplification","TF_binding_site_variant",
    "regulatory_region_ablation","regulatory_region_amplification",
    "regulatory_region_variant","intergenic_variant","sequence_variant")
  MISS_RANK <- length(VEP_ORDER) + 1L
  IMPACT_RANK <- c(HIGH = 0L, MODERATE = 1L, LOW = 2L, MODIFIER = 3L)

  csq_rank <- function(conseq) {
    if (is.na(conseq) || conseq == "") return(MISS_RANK)
    r <- match(strsplit(conseq, "&", fixed = TRUE)[[1]], VEP_ORDER)
    r[is.na(r)] <- MISS_RANK
    min(r)
  }
  norm_chrom <- function(x) {
    x <- ifelse(grepl("^chr", x), x, paste0("chr", x))
    ifelse(x == "chrMT", "chrM", x)
  }
  max_spliceai <- function(s) {
    if (is.na(s) || s == "" || s == ".") return(NA_real_)
    best <- NA_real_
    for (rec in strsplit(s, ",", fixed = TRUE)[[1]]) {
      f <- strsplit(rec, "|", fixed = TRUE)[[1]]
      if (length(f) >= 6) {
        v <- suppressWarnings(as.numeric(f[3:6]))
        v <- v[!is.na(v)]
        if (length(v)) best <- max(best, v, na.rm = TRUE)
      }
    }
    best
  }
  blank <- function(x) is.na(x) | x == "." | x == ""
  fmt_hgvs <- function(hgvsc, hgvsp) {
    c1 <- ifelse(blank(hgvsc), "", hgvsc)
    pp <- ifelse(blank(hgvsp), "", sub("^[^:]*:", "", hgvsp))
    pp <- gsub("%3D", "=", pp, fixed = TRUE)
    ifelse(pp != "", ifelse(c1 != "", paste(c1, pp), pp), c1)
  }

  ## read the variant table and coerce fields
  vt <- read.delim(vt_file, sep = "\t", header = TRUE, quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))
  for (c in c("hgvsc", "hgvsp"))
    if (is.null(vt[[c]])) vt[[c]] <- "."
  vt$chrom     <- norm_chrom(vt$chrom)
  vt$vid       <- paste(vt$chrom, vt$pos, vt$ref, vt$alt, sep = "_")
  vt$popmax    <- suppressWarnings(as.numeric(vt$popmax))
  vt$cadd      <- suppressWarnings(as.numeric(vt$cadd))
  vt$revel     <- suppressWarnings(as.numeric(vt$revel))
  vt$spliceaiN <- vapply(as.character(vt$spliceai), max_spliceai, numeric(1))
  vt$clinvar   <- ifelse(is.na(vt$clinvar) | vt$clinvar == ".", "", vt$clinvar)
  vt$gt        <- ifelse(is.na(vt$gt) | vt$gt == "", ".", vt$gt)
  vt$rank      <- vapply(as.character(vt$consequence), csq_rank, integer(1))

  message(sprintf("      %d variant rows, %d distinct variants",
                  nrow(vt), length(unique(vt$vid))))

  ## variant master: one row per variant, keeping the worst-consequence
  ## transcript as the representative (its consequence/impact/HGVS is displayed).
  vt <- vt[order(vt$vid, vt$rank), ]
  master <- vt[!duplicated(vt$vid),
               c("vid","chrom","pos","ref","alt","gt","popmax","cadd",
                 "spliceaiN","clinvar","revel","hgvsc","hgvsp",
                 "consequence","impact")]
  master$models <- tapply(vt$model, vt$vid,
                          function(m) paste(sort(unique(m)), collapse = ","))[master$vid]

  ## precompute the base display string per variant (without dSS, which is
  ## junction-specific and appended per junction). Space-separated: phased GT
  ## '0|1' carries '|', so '|' can't be a delimiter.
  trunc8 <- function(x) ifelse(nchar(x) <= 8, x, paste0(substr(x, 1, 8), "…"))
  add <- function(suf, cond, txt) ifelse(cond & !is.na(cond),
                                         ifelse(suf == "", txt, paste(suf, txt)), suf)
  m_head <- sprintf("%s:%d:%s>%s GT=%s %s", master$chrom, master$pos,
                    trunc8(master$ref), trunc8(master$alt), master$gt, master$models)
  m_hgvs <- fmt_hgvs(master$hgvsc, master$hgvsp)
  m_suf <- rep("", nrow(master))
  m_suf <- add(m_suf, !is.na(master$popmax),    sprintf("pmax=%.2g", master$popmax))
  m_suf <- add(m_suf, !is.na(master$cadd),      sprintf("CADD=%.1f", master$cadd))
  m_suf <- add(m_suf, !is.na(master$spliceaiN), sprintf("SAI=%.2f",  master$spliceaiN))
  m_suf <- add(m_suf, !is.na(master$revel),     sprintf("REVEL=%.2f",master$revel))
  m_suf <- add(m_suf, master$clinvar != "",     paste0("CLIN=", master$clinvar))
  m_base_cell <- gsub(" +", " ",
                      trimws(paste(m_head, master$consequence, master$impact,
                                   m_hgvs, m_suf)))
  m_sai0 <- ifelse(is.na(master$spliceaiN), -1, master$spliceaiN)
  midx   <- setNames(seq_len(nrow(master)), master$vid)

  ## position index (per chrom, ascending) for the window scans
  mo <- order(master$chrom, master$pos)
  pos_by_chrom <- split(master$pos[mo], master$chrom[mo])
  vid_by_chrom <- split(master$vid[mo], master$chrom[mo])

  juncs <- juncs %>%
    mutate(chrom = norm_chrom(as.character(seqnames)),
           start = as.integer(start), end = as.integer(end))

  # chromosome naming must agree after normalisation (FRASER seqnames are
  # Ensembl-style "1", the variant table UCSC-style "chr1"); a silent mismatch
  # leaves every junction with 0 variants.
  j_chroms <- unique(juncs$chrom); v_chroms <- names(pos_by_chrom)
  shared <- intersect(j_chroms, v_chroms)
  message(sprintf("      chromosome check: %d/%d junction chromosomes found in the variant table (junction e.g. %s; variant e.g. %s)",
                  length(shared), length(j_chroms),
                  paste(head(j_chroms, 3), collapse = ","), paste(head(v_chroms, 3), collapse = ",")))
  if (length(shared) == 0)
    stop("no junction chromosome matches the variant table after normalisation -- check chromosome naming in -r vs -V")
  missing_main <- setdiff(intersect(j_chroms, paste0("chr", c(1:22, "X"))), v_chroms)
  if (length(missing_main))
    warning("junction chromosomes absent from the variant table: ", paste(missing_main, collapse = ","))
  nj <- nrow(juncs)
  j_cand_n <- integer(nj); j_cand <- character(nj)
  j_near_n <- integer(nj); j_near <- character(nj)
  j_models <- character(nj); j_impact <- character(nj)
  j_sai <- rep(NA_real_, nj); j_cadd <- rep(NA_real_, nj); j_revel <- rep(NA_real_, nj)
  j_clinpath <- logical(nj); j_splice <- logical(nj)

  # order variants within a cell and build the topk string with dSS tokens.
  build_cell <- function(vids, dss) {
    if (length(vids) == 0) return(list(n = 0L, str = ""))
    ii  <- midx[vids]
    ord <- order(dss, -m_sai0[ii])              # nearest splice site, then SpliceAI
    ii  <- ii[ord]; dss <- dss[ord]
    cells <- paste0(m_base_cell[ii], sprintf(" dSS=%d", dss))
    list(n = length(vids), str = paste(head(cells, top_k), collapse = ";"))
  }
  worst_impact <- function(vids) {
    if (length(vids) == 0) return("")
    imp <- master$impact[midx[vids]]
    ir  <- IMPACT_RANK[imp]; ir[is.na(ir)] <- 9L
    imp[which.min(ir)]
  }

  for (j in seq_len(nj)) {
    ch <- juncs$chrom[j]; s <- juncs$start[j]; e <- juncs$end[j]
    pv <- pos_by_chrom[[ch]]; vv <- vid_by_chrom[[ch]]
    if (is.null(pv)) next
    sel <- which(pmin(abs(pv - s), abs(pv - e)) <= win_bp)  # near EITHER splice site
    if (length(sel) == 0) next
    p    <- pv[sel]; ids <- vv[sel]
    dss  <- pmin(abs(p - s), abs(p - e))                  # dist to nearest SS
    is_c <- dss <= cand_bp
    is_n <- !is_c & dss <= near_bp
    cv <- build_cell(ids[is_c], dss[is_c])
    nv <- build_cell(ids[is_n], dss[is_n])
    j_cand_n[j] <- cv$n; j_cand[j] <- cv$str
    j_near_n[j] <- nv$n; j_near[j] <- nv$str
    if (cv$n > 0) {
      ci <- midx[ids[is_c]]
      j_models[j] <- paste(sort(unique(unlist(
                        strsplit(master$models[ci], ",")))), collapse = ",")
      j_impact[j] <- worst_impact(ids[is_c])
      j_sai[j]   <- suppressWarnings(max(master$spliceaiN[ci], na.rm = TRUE))
      j_cadd[j]  <- suppressWarnings(max(master$cadd[ci],      na.rm = TRUE))
      j_revel[j] <- suppressWarnings(max(master$revel[ci],     na.rm = TRUE))
      cvv <- master$clinvar[ci]
      j_clinpath[j] <- any(grepl("athogenic", cvv, ignore.case = TRUE) &
                           !grepl("onflicting", cvv, ignore.case = TRUE))
      j_splice[j] <- any(grepl("splice", master$consequence[ci], ignore.case = TRUE))
    }
  }

  juncs$cand_var_n <- j_cand_n; juncs$cand_variants <- j_cand
  juncs$near_var_n <- j_near_n; juncs$near_variants <- j_near
  juncs$cand_models <- j_models; juncs$cand_max_impact <- j_impact
  j_sai[is.infinite(j_sai)]     <- NA
  j_cadd[is.infinite(j_cadd)]   <- NA
  j_revel[is.infinite(j_revel)] <- NA
  juncs$helper_sai <- j_sai; juncs$helper_cadd <- j_cadd; juncs$helper_revel <- j_revel
  juncs$helper_clinpath <- j_clinpath; juncs$helper_splice <- j_splice
  juncs$chrom <- NULL
  juncs
}

# Columns created up-front (NA) so they always exist for the select() below,
# regardless of whether -V was given.
res <- dplyr::bind_cols(res, empty_variant_cols(nrow(res)))

if (!is.null(opt$variants) && nzchar(opt$variants)) {
  message(sprintf("Adding coordinate-based variant columns from %s ...", opt$variants))
  juncs <- res %>% distinct(seqnames, start, end, strand)
  juncs <- compute_variant_columns(juncs, opt$variants, as.integer(opt[["cand-flank"]]),
                                   as.integer(opt[["near-flank"]]), opt$topk)

  # drop the placeholder cols, then join the computed ones back by junction key.
  res <- res %>%
    dplyr::select(-all_of(VARIANT_COLS)) %>%
    mutate(start = as.integer(start), end = as.integer(end)) %>%
    left_join(juncs, by = c("seqnames", "start", "end", "strand"))

  message(sprintf("      %d junctions with a candidate (splice-site) variant",
                  sum(juncs$cand_var_n > 0, na.rm = TRUE)))
}

## ---- dark-gene panel (Alliance / PacBio long-read panel) -------------------
# 374 capture-target intervals, hg38, over loci short reads align badly in
# (segmental duplications, paralogy, GC extremes, repeats) -- the PacBio "dark
# genes" panel. Two uses, both here so the BED is read once:
#   (a) the `dark_gene` column        -- informational, see below
#   (b) the positional MHC artifact rule -- see the pass_tier block
#
# CRITERION: a junction is IN a dark region when EITHER SPLICE SITE falls inside
# one. Measured alternatives on a real 247,642-row table (2026-09-02):
#   either splice site      10,502 (4.24%)   <- this one
#   both splice sites        6,806 (2.75%)
#   junction interval overlaps 15,161 (6.12%)
#   fully contained          6,694 (2.70%)
#   either site +/-500bp    10,772 (4.35%)
# Why "either site": a split read is anchored on both sides of the junction and
# aligned to NEITHER the intron nor anything between -- so what makes the call
# unreliable is hard sequence AT AN ANCHOR. Interval-overlap is the tempting
# choice and the wrong one: it nearly doubles the shortlist (22 vs 12 rows) by
# counting junctions whose intron merely SPANS a dark region, where no read
# aligns at all. Requiring both sites misses the likeliest mismapping case, a
# cryptic junction called off one paralogous anchor. No +/- window is used
# (minoverlap 0): regions are gene-scale (median 31 kb) so a 500 bp pad moves
# 4.24% -> 4.35% and changes nothing on the shortlist. That is the difference
# from the ENCODE blacklist below, whose regions are small enough to need one.
#
# The gene SYMBOL is deliberately NOT the primary key: the panel names `GBA`
# where Ensembl 115 says `GBA1` (and `FAM189B`, also renamed), so a symbol join
# silently loses one of the highest-value dark loci -- along with 20 of 704
# panel symbols (2.8%) that no longer reconcile. Coordinates do not rot.
# Symbol match is kept only as a weaker secondary signal, tagged "gene-only".
DARK_DEFAULT <- paste0("/hpf/largeprojects/tgnode/sandbox/yliang_analysis/",
                       "pacbio_longread/annotation_files/",
                       "Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed")

find_dark_genes <- function() {
  cands <- character(0)
  if (nzchar(Sys.getenv("DARK_GENES"))) cands <- c(cands, Sys.getenv("DARK_GENES"))
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f))
    cands <- c(cands, file.path(dirname(normalizePath(f, mustWork = FALSE)),
                                "dark_genes", basename(DARK_DEFAULT)))
  cands <- path.expand(c(cands, DARK_DEFAULT))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

# -> data.frame(chrom, start, end, names) in 1-based inclusive coordinates, with
# Ensembl-style contigs, so it joins to the results table without a seqlevels
# dance. BED is 0-based half-open; the results table is 1-based inclusive.
read_dark_bed <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  d <- utils::read.table(con, sep = "\t", header = FALSE, quote = "",
                         comment.char = "#", stringsAsFactors = FALSE)
  if (ncol(d) < 4)
    stop("dark-gene BED needs 4 columns (chrom, start, end, gene names): ", path)
  d <- d[, 1:4]; names(d) <- c("chrom", "start0", "end", "names")
  d$chrom <- sub("^chr", "", d$chrom)
  d$chrom[d$chrom == "M"] <- "MT"
  d$start <- as.integer(d$start0) + 1L
  d$end   <- as.integer(d$end)
  d[, c("chrom", "start", "end", "names")]
}

#' Region names hit by either of two query points, per row.
#' Vectorised per region (374 of them, tens per contig) rather than per row, so
#' this stays a handful of vector ops -- no GenomicRanges, keeping the dep set
#' at tidyverse+optparse.
#' @return character vector, "" where nothing is hit, ";"-joined when several.
dark_hits <- function(regions, qchrom, qstart, qend) {
  out <- character(length(qchrom))
  for (ch in intersect(unique(qchrom), unique(regions$chrom))) {
    idx <- which(qchrom == ch)
    if (!length(idx)) next
    rr <- regions[regions$chrom == ch, , drop = FALSE]
    s <- qstart[idx]; e <- qend[idx]
    for (k in seq_len(nrow(rr))) {
      hit <- (!is.na(s) & s >= rr$start[k] & s <= rr$end[k]) |
             (!is.na(e) & e >= rr$start[k] & e <= rr$end[k])
      if (!any(hit)) next
      j <- idx[hit]
      out[j] <- ifelse(nzchar(out[j]), paste(out[j], rr$names[k], sep = ";"),
                       rr$names[k])
    }
  }
  out
}

DARK <- NULL; DARK_SYMS <- character(0); dark_file <- NA_character_
{
  dk <- opt[["dark-genes"]]
  if (!is.null(dk) && identical(tolower(dk), "none")) {
    message("Dark-gene annotation skipped (--dark-genes none).")
  } else {
    dark_file <- if (!is.null(dk)) dk else find_dark_genes()
    if (is.na(dark_file) || !file.exists(path.expand(dark_file))) {
      warning("Dark-gene BED not found; 'dark_gene' column left as NA. ",
              "Pass --dark-genes <file> or set $DARK_GENES.")
      dark_file <- NA_character_
    } else {
      dark_file <- path.expand(dark_file)
      DARK <- read_dark_bed(dark_file)
      DARK_SYMS <- unique(trimws(unlist(strsplit(DARK$names, ","))))
      DARK_SYMS <- DARK_SYMS[nzchar(DARK_SYMS) & DARK_SYMS != "NA"]
      message(sprintf("Dark-gene panel: %d regions, %d gene names (%s).",
                      nrow(DARK), length(DARK_SYMS), basename(dark_file)))
    }
  }
}

## ---- pass_tier: outlier x variant-quality triage tier ----------------------
# Ordinal; the best tier a row qualifies for. Orthogonal to hpo_overlap (combine
# the two for the phenotype funnel). Non-outliers get "" so the autofilter can
# hide them. Concordance for splicing = the candidate is splice-disrupting (a
# splice_* consequence, or SpliceAI >= --tier-sai-strong); deltaPsi direction is
# NOT gated (exon-skip = loss, cryptic = gain, both valid).
#   1_concordant  outlier + a splice-disrupting candidate variant
#   2_scored      outlier + a candidate variant passing a score gate
#   3_outlier     significant splicing outlier, no qualifying variant
#   potential_artifact
#                 any gene in the junction's un-exploded hgncSymbol matches
#                 --artifact-biotypes (immunoglobulin / T-cell-receptor V,D,J,C
#                 segments, somatically rearranged in lymphoid tissue, so
#                 "aberrant" junctions there are V(D)J recombination, not
#                 splicing) OR gene symbol matches
#                 --artifact-genes (classical MHC class I/II: hyper-polymorphic,
#                 so reads mismap against the single reference haplotype, and
#                 the MIC paralogs multi-map onto each other -- either way the
#                 junction is a mapping artifact). MHC genes span protein_coding,
#                 pseudogene and lncRNA biotypes, hence a symbol match rather
#                 than a biotype match. Overrides any tier above.
#' IG/TR-symbol + MHC-symbol + MHC-by-position artifact classification.
#' Factored out (2026-09-22) so the private-novel table (--private-novel
#' below) reuses this EXACT logic instead of a second copy. IG/TR segments
#' are identified by a gene-SYMBOL regex tested against the junction's whole
#' un-exploded `hgnc_all` (every gene it overlaps), so a multi-gene junction
#' is classified consistently across all its exploded rows and the check
#' doesn't depend on a successful gene_biotype join. HLA genes span
#' protein_coding, pseudogene and lncRNA biotypes, so MHC needs its own
#' row-local symbol regex on `gene`, PLUS a positional check against the
#' dark-gene panel's HLA-named intervals: the symbol regex alone misses MHC
#' junctions labelled with a co-captured pseudogene/snRNA neighbour (on a
#' real shortlist, 6 of 12 dark-region rows were HLA-B/C, HLA-DRB1 or
#' HLA-DQB1 junctions carrying labels like USP8P1/RPL3P2/RNU1-61P/MTCO3P1/
#' MIR3135B, and every one escaped as 3_outlier without it). The dark panel's
#' own HLA-named regions are the coordinate list -- no second one to maintain:
#'     chr6:31,266,572-31,359,304  HLA-B,MIR6891,HLA-C,LINC02571   92.7 kb
#'     chr6:32,575,766-32,592,730  HLA-DRB1                        17.0 kb
#'     chr6:32,656,355-32,670,655  HLA-DQB1-AS1,HLA-DQB1           14.3 kb
#' 124 kb total = 2.5% of the chr6:28.5-33.4 Mb MHC block, considered and
#' rejected as too broad (it would take NOTCH4, C4A/C4B and TNF with it).
classify_artifact <- function(gene, hgnc_all, seqnames, start, end,
                              bt_re, gn_re, dark, no_mhc_by_position = FALSE,
                              label = "") {
  pfx <- if (nzchar(label)) paste0(label, ": ") else ""
  n <- length(gene)
  is_bt <- rep(FALSE, n); is_gn <- rep(FALSE, n)
  if (nzchar(bt_re)) is_bt <- any_symbol_matches(hgnc_all, bt_re)
  if (nzchar(gn_re)) is_gn <- !is.na(gene) & grepl(gn_re, gene)

  is_mhc <- rep(FALSE, n)
  if (!is.null(dark) && !isTRUE(no_mhc_by_position)) {
    hla <- dark[grepl("(^|,)\\s*HLA-", dark$names), , drop = FALSE]
    if (nrow(hla)) {
      is_mhc <- nzchar(dark_hits(hla, sub("^chr", "", as.character(seqnames)),
                                 as.integer(start), as.integer(end)))
      message(sprintf("%sMHC by position: %d region(s) from the dark panel, %d row(s) hit.",
                      pfx, nrow(hla), sum(is_mhc)))
    }
  }

  is_artifact <- is_bt | is_gn | is_mhc
  if (any(is_artifact))
    message(sprintf("%sFlagged %d rows as potential_artifact (%d IG/TR symbol ~ %s, %d gene ~ %s, %d MHC by position).",
                    pfx, sum(is_artifact), sum(is_bt), bt_re, sum(is_gn), gn_re, sum(is_mhc)))
  list(is_bt = is_bt, is_gn = is_gn, is_mhc = is_mhc, is_artifact = is_artifact)
}

{
  tp <- opt[["tier-p"]]; tdpsi <- opt[["tier-dpsi"]]
  tsai <- opt[["tier-sai"]]; tsais <- opt[["tier-sai-strong"]]
  tcadd <- opt[["tier-cadd"]]; trevel <- opt[["tier-revel"]]
  is_out <- !is.na(res$pValue) & !is.na(res$deltaPsi) & res$pValue < tp & abs(res$deltaPsi) >= tdpsi
  conc <- is_out & ((res$helper_splice %in% TRUE) |
                    (!is.na(res$helper_sai) & res$helper_sai >= tsais))
  scored <- is_out & (
    (!is.na(res$helper_sai)   & res$helper_sai   >= tsai)   |
    (!is.na(res$helper_cadd)  & res$helper_cadd  >= tcadd)  |
    (!is.na(res$helper_revel) & res$helper_revel >= trevel) |
    (res$helper_clinpath %in% TRUE) |
    (!is.na(res$cand_max_impact) & res$cand_max_impact %in% c("HIGH", "MODERATE")))
  res$pass_tier <- ifelse(!is_out, "",
                   ifelse(conc, "1_concordant",
                   ifelse(scored, "2_scored", "3_outlier")))

  bt_re <- opt[["artifact-biotypes"]]
  gn_re <- opt[["artifact-genes"]]
  art <- classify_artifact(res$gene, res$hgncSymbol, res$seqnames, res$start, res$end,
                           bt_re, gn_re, DARK, opt[["no-mhc-by-position"]])
  is_bt <- art$is_bt; is_gn <- art$is_gn; is_mhc <- art$is_mhc
  is_artifact <- art$is_artifact
  res$pass_tier[is_artifact] <- "potential_artifact"
}

## ---- blacklist-region flag -------------------------------------------------
# Reimplementation of FRASER::flagBlacklistRegions with identical semantics: a
# junction is flagged when EITHER splice site falls within --blacklist-minoverlap
# bp of an ENCODE blacklist region (Amemiya, Kundaje & Boyle 2019). FRASER builds
# a window around the donor ([start-mo, start+mo-1]) and the acceptor
# ([end-mo+1, end+mo]) and flags on a hit in either. Done here on plain intervals
# so the script keeps its tidyverse+optparse-only dependency set -- no
# GenomicRanges/rtracklayer needed. Factored into functions (2026-09-22) so the
# private-novel table (--private-novel below) reuses this EXACT implementation
# instead of a second copy -- there was briefly a duplicate in fraser_cohort.R,
# reverted in favour of routing private-novel annotation through this script.
find_blacklist <- function(assembly) {
  cands <- character(0)
  if (nzchar(Sys.getenv("FRASER_BLACKLIST"))) cands <- c(cands, Sys.getenv("FRASER_BLACKLIST"))
  cands <- path.expand(c(cands,
    "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/hg38-blacklist.v2.bed"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

#' Flag rows whose start or end splice site falls within `minoverlap` bp of a
#' blacklist region. Pure function of its coordinate vectors -- no dependency
#' on `res`/`opt` -- so any table with seqnames/start/end can reuse it.
flag_blacklist_rows <- function(seqnames, start, end, bl_file, minoverlap) {
  con <- if (grepl("\\.gz$", bl_file)) gzfile(bl_file) else bl_file
  bl <- utils::read.table(con, sep = "\t", header = FALSE, quote = "",
                          comment.char = "#", stringsAsFactors = FALSE)[, 1:3]
  names(bl) <- c("chrom", "start0", "end")
  # BED is 0-based half-open; the results table is 1-based inclusive. Strip the
  # "chr" prefix so Ensembl-style seqnames (1, 22, MT) join to UCSC-style BED.
  bl$chrom <- sub("^chr", "", bl$chrom)
  bl$chrom[bl$chrom == "M"] <- "MT"
  bl$start <- bl$start0 + 1L
  mo <- as.integer(minoverlap)
  qchrom <- sub("^chr", "", as.character(seqnames))
  hit <- logical(length(seqnames))
  for (ch in unique(qchrom)) {
    b   <- bl[bl$chrom == ch, , drop = FALSE]
    idx <- which(qchrom == ch)
    if (nrow(b) == 0L || length(idx) == 0L) next
    o <- order(b$start); s <- b$start[o]; cmax <- cummax(b$end[o])
    # sorted by start with a running max end: a window [qs,qe] overlaps some
    # region iff the last region starting at/before qe ends at/after qs.
    any_ov <- function(qs, qe) {
      j  <- findInterval(qe, s)
      ok <- j > 0L
      out <- logical(length(qs))
      out[ok] <- cmax[j[ok]] >= qs[ok]
      out
    }
    st <- as.integer(start[idx]); en <- as.integer(end[idx])
    hit[idx] <- any_ov(st - mo, st + mo - 1L) | any_ov(en - mo + 1L, en + mo)
  }
  hit
}

#' Resolve --blacklist/--assembly to a real BED path (or NA) and apply
#' flag_blacklist_rows(). Warns (never errors) and returns all-NA on a
#' missing/skipped file -- same degrade-gracefully convention as find_txdb()
#' etc. `label` prefixes the message/warning so the private-novel call site
#' (below) is distinguishable from the main table's in the log.
annotate_blacklist <- function(seqnames, start, end, bl_opt, assembly, minoverlap,
                               label = "") {
  pfx <- if (nzchar(label)) paste0(label, ": ") else ""
  bl_file <- if (!is.null(bl_opt)) bl_opt else find_blacklist(assembly)
  if (!is.null(bl_opt) && identical(tolower(bl_opt), "none")) {
    message(pfx, "Blacklist annotation skipped (--blacklist none).")
    return(rep(NA, length(seqnames)))
  }
  if (is.na(bl_file) || !file.exists(path.expand(bl_file))) {
    warning(pfx, "Blacklist BED not found; 'blacklist' column left as NA. ",
            "Pass --blacklist <file> or set $FRASER_BLACKLIST.")
    return(rep(NA, length(seqnames)))
  }
  bl_file <- path.expand(bl_file)
  hit <- flag_blacklist_rows(seqnames, start, end, bl_file, minoverlap)
  message(sprintf("%sBlacklist: %d of %d rows overlap a region (%s, minoverlap=%d).",
                  pfx, sum(hit), length(hit), basename(bl_file), as.integer(minoverlap)))
  hit
}

# NB: opt[["blacklist"]], not opt$blacklist -- optparse omits NULL defaults from
# the list, and R's $ would then partial-match "blacklist-minoverlap".
res$blacklist <- annotate_blacklist(res$seqnames, res$start, res$end,
                                    opt[["blacklist"]], opt$assembly,
                                    opt[["blacklist-minoverlap"]])

## ---- dark_gene flag --------------------------------------------------------
# One column, not four. Values:
#   ""                     no overlap
#   "GBA,FAM189B"          a splice site is inside a panel region; the value is
#                          the region's OWN gene list, so a row labelled with a
#                          co-captured neighbour still reads correctly (e.g.
#                          gene USP8P1 -> "HLA-B,MIR6891,HLA-C,LINC02571")
#   "PNKP, gene-only"      the gene is a panel target but neither splice site is
#                          inside a captured interval -- weaker, and flagged as
#                          such rather than silently equated with the above
# Informational ONLY: it does NOT feed pass_tier and does NOT drop rows. Dark
# genes are the inverse of the IG/TR/MHC artifact rule -- those are loci where a
# real finding is near-impossible, these are loci curated BECAUSE they carry real
# findings that short reads miss (NPRL3, ANKRD11, PMS2, SMN1/2, STRC, GBA1...).
# Folding them into potential_artifact would delete them from the shortlist,
# which is the opposite of why the panel exists. If you later want an artifact
# rule out of this, the defensible one is "dark AND recurrent across the cohort",
# and cohort recurrence is not computed here yet.
res$dark_gene <- NA_character_
if (!is.null(DARK)) {
  hits <- dark_hits(DARK, sub("^chr", "", as.character(res$seqnames)),
                    as.integer(res$start), as.integer(res$end))
  sym_only <- !nzchar(hits) & !is.na(res$gene) & res$gene %in% DARK_SYMS
  hits[sym_only] <- paste0(res$gene[sym_only], ", gene-only")
  res$dark_gene <- hits
  message(sprintf("Dark genes: %d of %d rows in a panel region, %d more by gene name only.",
                  sum(nzchar(hits)) - sum(sym_only), length(hits), sum(sym_only)))
}

## ---- exon annotation -------------------------------------------------------
# Which exons of the gene each junction touches.
#
# FRASER's start/end are the INTRON BODY: the flanking exons end at start-1 and
# begin at end+1. Verified exactly on 100% of fully annotated junctions, so the
# anchors are matched by exact arithmetic -- no maxgap fudge (unlike
# clinical_pipeline/scripts/junction_outlier.R, whose maxgap=5 can mis-assign in
# genes with closely spaced exons).
#
# Numbers come from the pre-extracted exon table built by
# scripts/exon_annotation/extract_exons.py: Ensembl's own strand-aware
# exon_number on the MANE Select transcript where the gene has one, else a
# per-gene union model of all its exons. Both refs are present for every gene,
# which is what lets a boundary that is novel w.r.t. MANE but used by another
# isoform ("alt-") be told apart from one that is novel everywhere ("novel-").
#
# STRAND: on the minus strand the exon at end+1 is the DONOR, not the acceptor.
find_exon_table <- function() {
  cands <- character(0)
  if (nzchar(Sys.getenv("FRASER_EXONS"))) cands <- Sys.getenv("FRASER_EXONS")
  cands <- path.expand(c(cands, file.path(
    "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files",
    "Homo_sapiens.GRCh38.115.exons.tsv.gz")))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

res$exons_affected <- NA_character_
{
  ex_file <- find_exon_table()
  qgene <- as.character(res$gene)
  qok   <- !is.na(qgene) & trimws(qgene) != "" & qgene != "NA"
  if (is.na(ex_file)) {
    warning("Exon table not found; 'exons_affected' left as NA. Build it with ",
            "scripts/exon_annotation/extract_exons.py, or set $FRASER_EXONS.")
  } else if (!any(qok)) {
    message("No named genes; 'exons_affected' left as NA.")
  } else {
    ex <- read.delim(gzfile(ex_file), sep = "\t", header = TRUE, quote = "",
                     stringsAsFactors = FALSE)
    ex <- ex[ex$gene_symbol %in% unique(qgene[qok]), , drop = FALSE]

    qchr <- sub("^chr", "", as.character(res$seqnames))
    qs   <- as.integer(res$start)
    qe   <- as.integer(res$end)
    qst  <- as.character(res$strand)
    # one work unit per unique (junction x gene); gene-exploded rows repeat a
    # junction, and the same junction/gene pair can recur across the table.
    key  <- paste(qgene, qchr, qs, qe, qst, sep = "\r")
    ui   <- which(qok & !duplicated(key))

    # MANE where available, else the union model. Decided per (gene, contig,
    # strand), not per gene: repeat families (Metazoa_SRP, snoRNA classes) reuse
    # one symbol at many loci, and a MANE transcript at one locus says nothing
    # about a junction at another.
    locus <- paste(ex$gene_symbol, ex$chrom, ex$strand, sep = "\r")
    mane_loci <- unique(locus[ex$ref == "MANE"])
    uref <- ifelse(paste(qgene[ui], qchr[ui], qst[ui], sep = "\r") %in% mane_loci,
                   "MANE", "union")

    # Exact point lookups. Within one transcript (and within a merged union
    # model) exons do not overlap, so a position identifies at most one exon.
    mk <- function(g, r, c, st, p) paste(g, r, c, st, p, sep = "\r")
    end_key   <- mk(ex$gene_symbol, ex$ref, ex$chrom, ex$strand, ex$end)
    start_key <- mk(ex$gene_symbol, ex$ref, ex$chrom, ex$strand, ex$start)
    ug <- qgene[ui]; uc <- qchr[ui]; us <- qs[ui]; ue <- qe[ui]; ust <- qst[ui]

    lef <- ex$exon_number[match(mk(ug, uref,    uc, ust, us - 1L), end_key)]
    rig <- ex$exon_number[match(mk(ug, uref,    uc, ust, ue + 1L), start_key)]
    lea <- ex$exon_number[match(mk(ug, "union", uc, ust, us - 1L), end_key)]
    ria <- ex$exon_number[match(mk(ug, "union", uc, ust, ue + 1L), start_key)]
    # n_exons is per (gene, ref, contig, strand) for the same reason.
    tot <- ex$n_exons[match(paste(ug, uref, uc, ust, sep = "\r"),
                            paste(ex$gene_symbol, ex$ref, ex$chrom, ex$strand,
                                  sep = "\r"))]

    # Skipped exons: those falling entirely inside the intron. Exons of a single
    # transcript / union model are non-overlapping, so sorting by start also
    # sorts by end -- the qualifying exons are then a CONTIGUOUS range, found
    # with two findInterval() calls instead of an all-pairs comparison.
    skipped <- rep(NA_character_, length(ui))
    gkey <- paste(ex$gene_symbol, ex$ref, ex$chrom, ex$strand, sep = "\r")
    egrp <- split(seq_len(nrow(ex)), gkey)
    qkey <- paste(ug, uref, uc, ust, sep = "\r")
    for (g in unique(qkey)) {
      ei <- egrp[[g]]
      if (is.null(ei)) next
      o  <- order(ex$start[ei]); ei <- ei[o]
      es <- ex$start[ei]; ee <- ex$end[ei]; en <- ex$exon_number[ei]
      qi <- which(qkey == g)
      if (is.unsorted(ee)) {
        # Should not happen -- extract_exons.py merges the union model so a
        # group's exons never overlap. Fall back rather than trust the shortcut.
        skipped[qi] <- vapply(qi, function(k) {
          h <- en[es >= us[k] & ee <= ue[k]]
          if (length(h)) paste(sort(h), collapse = ",") else NA_character_
        }, character(1))
        next
      }
      a  <- findInterval(us[qi] - 1L, es) + 1L   # first exon starting at/after us
      b  <- findInterval(ue[qi], ee)             # last exon ending at/before ue
      hit <- which(a <= b)
      if (length(hit))
        skipped[qi[hit]] <- vapply(hit, function(k)
          paste(sort(en[a[k]:b[k]]), collapse = ","), character(1))
    }

    # Label each side by strand: '+' donor is the exon at start-1, '-' donor is
    # the exon at end+1. A side with no match is "alt-" when some other isoform
    # has that boundary, else "novel-".
    lab <- function(n, alt, what)
      ifelse(!is.na(n), paste0("exon", n),
             ifelse(!is.na(alt), paste0("alt-", what), paste0("novel-", what)))
    plus <- ust == "+"
    don <- ifelse(plus, lab(lef, lea, "donor"),    lab(rig, ria, "donor"))
    acc <- ifelse(plus, lab(rig, ria, "acceptor"), lab(lef, lea, "acceptor"))

    lbl <- paste0(don, " -> ", acc,
                  ifelse(is.na(skipped), "", paste0(", skips exon", skipped)),
                  " (of ", tot, ifelse(uref == "union", ", union", ""), ")")
    lbl[is.na(tot)] <- NA_character_   # gene absent from the exon table
    res$exons_affected[qok] <- lbl[match(key[qok], key[ui])]

    nres <- sum(!is.na(res$exons_affected))
    message(sprintf(paste("Exon annotation: %d of %d rows assigned (%d unique",
                          "junction x gene; %d used the union model, %d skip an exon)."),
                    nres, nrow(res), length(ui), sum(uref == "union"),
                    sum(!is.na(skipped))))
  }
}

## ---- ranking helper flags --------------------------------------------------
nonblank <- function(x) !is.na(x) & trimws(as.character(x)) != ""
res$has_disease_anno <- ifelse(
  nonblank(res$omim_phenotype) | nonblank(res$Orphanet) |
  nonblank(res$clinGen_disease_label) | nonblank(res$hpo_id),
  "yes", "no")

## ---- pvalue_rank: full-table significance rank, one value per JUNCTION -----
# Ranks every JUNCTION in the full table purely by pValue -- no hpo_overlap
# tie-break (unlike `rank` below) -- so it answers "how significant is this
# call in the whole sample, independent of phenotype match". Computed once per
# UNIQUE junction (dedup on seqnames/start/end/strand) so a junction mapped to
# several genes gets ONE consistent value across all its exploded rows, not a
# different one per gene row. Ties share the same (minimum) rank --
# dplyr::min_rank(), i.e. competition ranking (1,1,3,4,...). Present in BOTH
# outputs (full + filtered), since it's computed here on the full table.
jkey_rows   <- paste(res$seqnames, res$start, res$end, res$strand, sep = "\r")
junc_pval   <- tapply(res$pValue, jkey_rows, function(x) x[1])
junc_rank   <- stats::setNames(dplyr::min_rank(junc_pval), names(junc_pval))
res$pvalue_rank <- junc_rank[jkey_rows]

## ---- arrange and write -----------------------------------------------------
# Sort priority:
#   1. any HPO overlap (yes first)
#   2. pValue (smallest first), then |deltaPsi| (largest first)
# `rank` (the filtered report's own row-order column) is assigned LATER, on
# res_filt only, once the filtered subset is known -- see "filtered subset"
# below.
res <- res %>%
  dplyr::select(pvalue_rank,
                seqnames, start, end, width, strand, gene, hgncSymbol,
                exons_affected, pass_tier,
                hpo_overlap, n_hpo_overlap, hpo_jaccard,
                overlapping_hpo_terms, has_disease_anno,
                omim_phenotype, omim_inheritance, Orphanet, clinGen_disease_label,
                type, pValue, padjust, deltaPsi, psiValue, deltaPsi_raw,
                counts, totalCounts, meanCounts, meanTotalCounts,
                nonsplitCounts, nonsplitProportion, nonsplitProportion_99quantile,
                cand_var_n, cand_variants, near_var_n, near_variants,
                cand_models, cand_max_impact,
                annotatedJunction, blacklist, dark_gene,
                everything()) %>%
  dplyr::select(-dplyr::any_of(c("sampleID", "gene_id", "hgnc_symbols",
                                 "helper_sai", "helper_cadd", "helper_revel",
                                 "helper_clinpath", "helper_splice",
                                 "GO.term.accession", "GO.domain"))) %>%
  arrange(desc(hpo_overlap == "yes"),
          pValue,
          desc(abs(deltaPsi)))

## ---- filtered subset -------------------------------------------------------
# The full table is always written (gzipped); this is the shortlist meant to be
# opened in Excel. A row is kept when it has a gene (unless --keep-nogene) and:
#   - it is NOT an artifact locus, and either carries a pass_tier (the outlier
#     gate: --tier-p / --tier-dpsi) or clears the softer --lenient-p /
#     --lenient-dpsi gate AND is phenotype- or variant-supported (hpo_overlap
#     == yes, or a candidate splice-site variant); or
#   - it IS an artifact locus on the MHC side (gene symbol ~ --artifact-genes,
#     or the positional MHC-by-dark-panel hit) and clears the strict outlier
#     gate -- kept and still labelled potential_artifact, since a significant
#     call at a locus known to mismap is still worth a look, just flagged; or
#   - it IS an artifact locus on the IG/TR side (any gene in the junction's
#     un-exploded hgncSymbol ~ --artifact-biotypes, i.e. V(D)J recombination)
#     and clears the strict outlier gate AND has hpo_overlap == yes -- IG/TR
#     "outliers" are near-universal background otherwise (measured on 8 real
#     LCL samples: 0-1 of several hundred significant IG/TR rows per sample
#     carry any HPO overlap, vs 4-7 on the MHC side), so the HPO match is
#     required, not just a rescue.
# NB: is_igtr/is_mhc_side are recomputed here (not reusing is_bt/is_gn/is_mhc
# from the pass_tier block above) because those were computed BEFORE the
# arrange() a few lines up reordered res -- reusing them here would silently
# misalign flags with rows. any_symbol_matches(res$hgncSymbol, ...) is a fresh
# COLUMN read, so it's safe to recompute post-arrange exactly like is_artifact
# (which reads the res$pass_tier column, also safe for the same reason).
{
  is_artifact <- res$pass_tier %in% "potential_artifact"
  bt_re       <- opt[["artifact-biotypes"]]
  is_igtr     <- is_artifact & nzchar(bt_re) &
                 any_symbol_matches(res$hgncSymbol, bt_re)
  is_mhc_side <- is_artifact & !is_igtr

  has_gene <- nonblank(res$gene) & !res$gene %in% "NA"
  gene_ok  <- opt[["keep-nogene"]] | has_gene

  tp <- opt[["tier-p"]]; tdpsi <- opt[["tier-dpsi"]]
  is_out    <- !is.na(res$pValue) & !is.na(res$deltaPsi) &
               res$pValue < tp & abs(res$deltaPsi) >= tdpsi
  lenient   <- !is.na(res$pValue) & !is.na(res$deltaPsi) &
               res$pValue < opt[["lenient-p"]] &
               abs(res$deltaPsi) >= opt[["lenient-dpsi"]]
  supported <- (res$hpo_overlap %in% "yes") |
               (!is.na(res$cand_var_n) & res$cand_var_n > 0)

  keep_nonartifact <- !is_artifact & gene_ok & (nonblank(res$pass_tier) | (lenient & supported))
  keep_mhc         <- is_mhc_side & gene_ok & is_out
  keep_igtr        <- is_igtr & gene_ok & is_out & (res$hpo_overlap %in% "yes")

  keep <- keep_nonartifact | keep_mhc | keep_igtr
  res_filt <- res[keep, , drop = FALSE]
  if (any(is_artifact))
    message(sprintf("Artifact rows admitted to filtered output: %d MHC-side (of %d), %d IG/TR-side (of %d), still labelled potential_artifact.",
                    sum(keep_mhc), sum(is_mhc_side), sum(keep_igtr), sum(is_igtr)))

  # `rank`: the filtered report's OWN row order, 1-based, JUNCTION-based (a
  # junction spanning several kept gene-rows gets ONE rank, not a different
  # number per row) -- res_filt already inherits its display order from res's
  # arrange() above (a boolean subset never reorders), so this just numbers
  # the DISTINCT junctions in that already-correct order: match() against
  # unique() gives every row the 1-based position of its junction's FIRST
  # occurrence, i.e. how many distinct junctions have appeared so far.
  jkey_filt <- paste(res_filt$seqnames, res_filt$start, res_filt$end,
                     res_filt$strand, sep = "\r")
  res_filt$rank <- match(jkey_filt, unique(jkey_filt))
  res_filt <- dplyr::relocate(res_filt, rank, .before = 1)
}

## ---- splicing consequences on the filtered subset (optional) ---------------
# Only when --txdb + an --fds-* option are given -- --fds-* ALONE is no longer
# an error (2026-09-21): it's independently legitimate for deltaPsi_raw's
# cohort median above, which needed no txdb at all. Reuses fds_early (already
# loaded above, if given) rather than opening the fds a second time. Applied
# to the FILTERED rows only -- DROP likewise runs annotatePotentialImpact on
# the aberrant subset, never on the all=TRUE table, because it is expensive.
{
  txdb_path <- opt[["txdb"]]
  if (!is.null(txdb_path) && have_fds && nrow(res_filt) > 0) {
    message("Adding splicing consequences to the filtered rows ...")
    txdb <- AnnotationDbi::loadDb(path.expand(txdb_path))
    res_filt <- add_splice_impact(res_filt, fds_early, txdb,
                                  add_blacklist = isTRUE(opt[["impact-blacklist"]]),
                                  assembly = opt$assembly)
    # add_splice_impact() appends its columns; move them next to the other
    # variant/impact interpretation columns (right after cand_max_impact, i.e.
    # after the nonsplit-PSI block) rather than leaving them stranded at the
    # far right of the sheet.
    res_filt <- res_filt %>%
      dplyr::relocate(dplyr::any_of(c("potentialImpact", "causesFrameshift",
                                      "UTR_overlap")),
                      .after = dplyr::any_of("cand_max_impact"))
  } else if (!is.null(txdb_path) && !have_fds) {
    stop("--txdb needs an fds too: pass --fds-file, or --fds-dir with ",
         "--dataset and --annotation-name")
  }
  # --fds-* with no --txdb: fine now -- it was already used above for
  # deltaPsi_raw's cohort median; splicing-consequence annotation is simply
  # skipped, same as giving neither.
}

## ---- drop the constant 'type' column ---------------------------------------
# FRASER2 has exactly one psi type (jaccard), so this column is always the
# same literal value -- dropped from the output here, not earlier, because
# add_splice_impact() (just above) requires it as an input column (it reads
# df$type to call FRASER's own addPotentialImpactLabels with the right psi
# type) whenever --txdb + an --fds-* option are given.
res      <- dplyr::select(res,      -dplyr::any_of("type"))
res_filt <- dplyr::select(res_filt, -dplyr::any_of("type"))

## ---- rename gene -> hgncSymbol, hgncSymbol -> hgncSymbol_all ---------------
# `gene` was the cleaned, EXPLODED one-gene-per-row column; the original
# `hgncSymbol` was FRASER/DROP's raw ';'-joined, un-exploded, multi-gene-per-
# junction string. Renamed here, at the very last step, so every internal
# computation above (gene annotation join, pass_tier, the IG/TR artifact
# check, has_gene, ...) keeps referring to them by their original names --
# only the OUTPUT header changes. dplyr::rename() resolves both renames
# against the ORIGINAL column names simultaneously (verified: no collision
# from the shared "hgncSymbol" name on both sides of the swap).
res      <- dplyr::rename(res,      hgncSymbol = gene, hgncSymbol_all = hgncSymbol)
res_filt <- dplyr::rename(res_filt, hgncSymbol = gene, hgncSymbol_all = hgncSymbol)

# Write the table. A .tsv output is written directly; a .xlsx output is written
# as a sibling TSV and converted by tsv_to_xlsx.py (shipped next to this script),
# so no Excel-writing R package is needed. The .xlsx keeps a frozen header row +
# autofilter. Requires python3 with pandas + xlsxwriter/openpyxl.
find_converter <- function() {
  cands <- character(0)
  a <- commandArgs(trailingOnly = FALSE)                 # Rscript exposes --file=
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f))
    cands <- c(cands, file.path(dirname(normalizePath(f, mustWork = FALSE)),
                                "tsv_to_xlsx.py"))
  if (nzchar(Sys.getenv("TSV_TO_XLSX"))) cands <- c(Sys.getenv("TSV_TO_XLSX"), cands)
  cands <- path.expand(c(cands, "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/SR_RNAseq_pipeline/scripts/tsv_to_xlsx.py"))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}
write_tsv <- function(df, path) {
  write.table(df, path, col.names = TRUE, row.names = FALSE,
              sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}
write_tsv_gz <- function(df, path) {
  con <- gzfile(path, "w")
  on.exit(close(con), add = TRUE)
  write.table(df, con, col.names = TRUE, row.names = FALSE,
              sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

# Three outputs off one -o path, e.g. -o SAMPLE.xlsx gives
#   SAMPLE.tsv.gz        every row, gzipped (the archive/audit copy)
#   SAMPLE.filtered.tsv  the shortlist
#   SAMPLE.xlsx          the same shortlist, converted (only if -o ends in .xlsx)
base     <- tools::file_path_sans_ext(outfile)
full_gz  <- paste0(base, ".tsv.gz")
filt_tsv <- paste0(base, ".filtered.tsv")
is_xlsx  <- grepl("\\.xlsx$", outfile, ignore.case = TRUE)

## ---- number formatting for the human-facing (filtered) output -------------
# Applied to res_filt only -- res (the .tsv.gz) is explicitly the archive/audit
# copy, kept at full precision on purpose (see the comment above write_tsv_gz).
#
# Most numeric columns (deltaPsi, psiValue, meanCounts, meanTotalCounts,
# nonsplitCounts, nonsplitProportion, nonsplitProportion_99quantile) are
# already written pre-rounded to 2dp by FRASER itself (AllGenerics.R:691-704),
# so the blanket round(.,3) below is a no-op for them. FOUR columns are NOT
# pre-rounded and genuinely carry many decimal digits:
#   pValue, padjust   raw beta-binomial p-values, e.g. 1.2352e-09
#   hpo_jaccard       n_hpo_overlap / union_size, a raw ratio of integers,
#                     never rounded (e.g. 1/7 = 0.142857142857143)
#   deltaPsi_raw      psiValue (2dp) minus the cohort row-median psi; when
#                     that median comes from --fds-file's raw jaccard assay
#                     (the normal case for a Snakemake-pipeline parse -- see
#                     the "cohort-wide row-median psi" block above) it is full
#                     double precision, not 2dp
# For these four: round to 3 decimals; a naive round(x, 3) would zero out any
# real value below 5e-4 (e.g. a real 1.2e-09 becomes "0.000", indistinguishable
# from a true zero or NA), so anything with |x| < 0.001 (and x != 0) is shown
# in scientific notation instead. Turning them to character first means the
# blanket round() pass below skips them (mutate() sees each step's own
# updates sequentially).
fmt_small <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x != 0 & abs(x) < 0.001,
                formatC(x, format = "e", digits = 2),
                as.character(round(x, 3))))
}
res_filt <- res_filt %>%
  dplyr::mutate(
    dplyr::across(dplyr::any_of(c("pValue", "padjust", "hpo_jaccard",
                                  "deltaPsi_raw")),
                  fmt_small),
    dplyr::across(dplyr::where(is.numeric), ~ round(.x, 3))
  )

write_tsv_gz(res, full_gz)
write_tsv(res_filt, filt_tsv)
if (is_xlsx) {
  conv <- find_converter()
  if (is.na(conv))
    stop("output is .xlsx but tsv_to_xlsx.py was not found (set $TSV_TO_XLSX or ",
         "keep it beside this script). Filtered TSV written to ", filt_tsv)
  py <- if (nzchar(Sys.which("python3"))) "python3" else "python"
  st <- system2(py, c(shQuote(conv), shQuote(filt_tsv), shQuote(outfile)))
  if (st != 0)
    stop("tsv_to_xlsx.py failed (exit ", st, "); filtered TSV kept at ", filt_tsv)
}

hpo_note <- if (have_hpo) {
  sprintf("%d rows with HPO overlap", sum(res$hpo_overlap == "yes", na.rm = TRUE))
} else {
  "HPO columns set to NA (no terms provided)"
}
message(sprintf("Wrote %d rows (full, %s) to %s.", nrow(res), hpo_note, full_gz))
message(sprintf("Wrote %d rows (filtered, %.2f%% of full) to %s%s.",
                nrow(res_filt), 100 * nrow(res_filt) / max(nrow(res), 1), filt_tsv,
                if (is_xlsx) paste0(" and ", outfile) else ""))

## ---- optional: annotate fraser_cohort.R's private-novel table -------------
# Independent of -r/-s's own table -- the private-novel table already has ONE
# clean gene per row (hgnc_symbol, no explode needed) and mostly no pValue
# (that is WHY it exists: those junctions were never tested). Reuses the SAME
# join_gene_annotation()/compute_hpo_overlap()/annotate_blacklist()/
# classify_artifact()/compute_variant_columns() functions used above for the
# filtered junctions, so the two tables can never drift apart -- this is what
# fraser_cohort.R's private_novel_empty() doc comment always intended
# ("Expect mapping artifacts to dominate until parse_fraser.R's artifact
# filters are applied"), rather than a second implementation living in
# fraser_cohort.R itself (tried, then reverted, same day).
if (!is.null(opt[["private-novel"]])) {
  n_file <- opt[["private-novel"]]
  if (!file.exists(n_file)) stop("--private-novel file not found: ", n_file)
  message(sprintf("Annotating private-novel table %s ...", n_file))

  pn <- read.delim(n_file, sep = "\t", header = TRUE, quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))
  if (!"hgnc_symbol" %in% names(pn))
    stop("--private-novel table has no hgnc_symbol column -- wrong file?")
  if ("sampleID" %in% names(pn) && length(unique(pn$sampleID)) > 1) {
    n0 <- nrow(pn)
    pn <- pn[pn$sampleID == sample, , drop = FALSE]
    message(sprintf("  private-novel: %d/%d rows for sample '%s'.", nrow(pn), n0, sample))
  }

  # Gene annotation -- same table, same join as the filtered junctions.
  pn <- join_gene_annotation(pn, "hgnc_symbol", gene_anno_file)

  # HPO overlap -- same individual HPO list as -p above.
  hpo_cols_pn <- compute_hpo_overlap(pn$hpo_id, pn$hpo_name, ind_hpo, sample)
  pn$overlapping_hpo_terms <- hpo_cols_pn$overlapping_hpo_terms
  pn$n_hpo_overlap         <- hpo_cols_pn$n_hpo_overlap
  pn$hpo_overlap           <- hpo_cols_pn$hpo_overlap
  pn$hpo_jaccard           <- hpo_cols_pn$hpo_jaccard

  # Blacklist -- same BED/resolution as the filtered junctions.
  pn$blacklist <- annotate_blacklist(pn$seqnames, pn$start, pn$end,
                                     opt[["blacklist"]], opt$assembly,
                                     opt[["blacklist-minoverlap"]],
                                     label = "private-novel")

  # Artifact classification. This table has only ONE gene per row already
  # (no explode), so hgnc_symbol doubles as both the row-local symbol (MHC
  # check) and the "every gene the junction overlaps" string (IG/TR check) --
  # there is no separate raw/un-exploded column to read here.
  art_pn <- classify_artifact(pn$hgnc_symbol, pn$hgnc_symbol, pn$seqnames,
                              pn$start, pn$end, opt[["artifact-biotypes"]],
                              opt[["artifact-genes"]], DARK,
                              opt[["no-mhc-by-position"]], label = "private-novel")
  pn$potential_artifact <- art_pn$is_artifact

  # Candidate-variant proximity -- optional, same -V file as the filtered
  # junctions, same columns created up-front so they always exist.
  pn <- dplyr::bind_cols(pn, empty_variant_cols(nrow(pn)))
  if (!is.null(opt$variants) && nzchar(opt$variants)) {
    juncs_pn <- pn %>% distinct(seqnames, start, end, strand)
    juncs_pn <- compute_variant_columns(juncs_pn, opt$variants,
                                        as.integer(opt[["cand-flank"]]),
                                        as.integer(opt[["near-flank"]]), opt$topk)
    pn <- pn %>%
      dplyr::select(-all_of(VARIANT_COLS)) %>%
      mutate(start = as.integer(start), end = as.integer(end)) %>%
      left_join(juncs_pn, by = c("seqnames", "start", "end", "strand"))
  }
  pn <- pn %>% dplyr::select(-dplyr::any_of(c("helper_sai", "helper_cadd", "helper_revel",
                                       "helper_clinpath", "helper_splice")))

  # Same top-level sort priority as the filtered junctions: HPO-matched rows
  # first (the single most powerful triage signal, and orthogonal to
  # everything else). potential_artifact rows STAY in the table (this is a
  # raw evidence list, not a filtered shortlist -- see fraser_cohort.R's own
  # doc comment: "rank it on read support, gene context and variant
  # proximity"), just sorted after non-artifact rows within each HPO group;
  # ties broken by read count descending, matching fraser_cohort.R's own
  # setorder(dt, -count).
  pn <- pn %>% arrange(desc(hpo_overlap == "yes"), potential_artifact, desc(count))

  n_out <- if (!is.null(opt[["private-novel-out"]])) opt[["private-novel-out"]]
           else paste0(tools::file_path_sans_ext(n_file), "_annotated.tsv")
  write_tsv(pn, n_out)
  message(sprintf(
    "Wrote %d annotated private-novel row(s) (%d potential_artifact, %d blacklisted, %d HPO-matched) to %s.",
    nrow(pn), sum(pn$potential_artifact, na.rm = TRUE),
    sum(pn$blacklist, na.rm = TRUE),
    sum(pn$hpo_overlap == "yes", na.rm = TRUE), n_out))
}
