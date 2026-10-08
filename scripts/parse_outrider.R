#!/usr/bin/env Rscript
# Parse OUTRIDER results for a single individual.
#
# Steps:
#   1. Read the OUTRIDER result .rds and extract the rows for the input sample.
#      The rds may be a results data.frame/data.table already, OR a fitted
#      OutriderDataSet object (then OUTRIDER::results(ods, all=TRUE) is run).
#   2. Add gene-level annotations (matched by gene id) from the hg38 gene
#      annotation table. If a -g coords file is supplied, also add a gene_coord
#      column (chr:start-end from the GTF/coords file).
#   3. Read the individual's HPO terms and, for each gene, report how many of
#      the individual's HPO terms overlap the gene's associated HPO terms (as
#      term NAMES), plus a yes/no flag for whether there is any overlap.
#   4. (optional) If a model_variants.tsv from trio_rare_filter.py is supplied
#      (-V) with a gene-coords file (-g), add per-gene rare-variant columns
#      (cand_var_n/cand_variants/near_var_n/near_variants/cand_models/
#      cand_max_impact) after l2fc. This folds in outrider_add_variants.py, so no
#      bcftools is needed here: trio_rare_filter.py has already done the VEP work.
#   5. Write the annotated, sorted table. Default output is .xlsx (frozen header
#      + autofilter), produced by writing a sibling TSV and converting it with
#      tsv_to_xlsx.py (so no Excel R package is needed); a .tsv output path skips
#      the conversion and writes the TSV directly.
#
# Usage:
#   Rscript parse_outrider.R -r <outrider_result.rds> -s <sampleID> \
#       -p <hpo_terms.txt> -a <gene_annotation.tsv> [-o <output.tsv>] \
#       [-V <model_variants.tsv> -g <gene_coords.tsv> --flank 10 --topk 20]
#   Rscript parse_outrider.R --help

suppressMessages({
  library(tidyverse)
  library(optparse)
})

## ---- arguments -------------------------------------------------------------
option_list <- list(
  make_option(c("-r", "--rds"), type = "character",
              default = "LCL_TG_OUTRIDER_Mar16.rds",
              help = "OUTRIDER result .rds: a results data.frame/data.table, or a fitted OutriderDataSet (results(ods,all=TRUE) is run) [default: %default]"),
  make_option(c("-s", "--sample"), type = "character",
              default = "DSK_023_03",
              help = "Sample ID to extract from the results [default: %default]"),
  make_option(c("-p", "--hpo"), type = "character",
              default = "DSK_023_03_geneCEP152_RL_LR_hpo.txt",
              help = "Individual's HPO terms file (one HP:xxxxxxx id per line) [default: %default]"),
  make_option(c("-a", "--annotation"), type = "character",
              default = "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_annotation/gtf115_annotation_table_Nov3.tsv",
              help = "Gene-level annotation table (.tsv, gene_id matches OUTRIDER geneID) [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output path; .xlsx (via tsv_to_xlsx.py) or .tsv by extension [default: <sample>_OUTRIDER_results_annotated.xlsx]"),
  make_option(c("-V", "--variants"), type = "character", default = NULL,
              help = "model_variants.tsv from trio_rare_filter.py (VEP-expanded, model-tagged variant table); if given, per-gene rare-variant columns are added [default: none]"),
  make_option(c("-g", "--coords"), type = "character", default = "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_coords.GRCh38.115.tsv",
              help = "gene coordinates for the 'nearby variant' window and the dark-gene overlap: a compact 6-col TSV (ensg,chrom,start,end,strand,name) or an Ensembl GTF. Default: annotation_files/gene_coords.GRCh38.115.tsv, or $GENE_COORDS. Use 'none' to skip."),
  make_option("--dark-genes", type = "character", default = "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed",
              help = "Alliance dark-gene panel BED (chrom/start/end/comma-separated gene names); default: /hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed, or $DARK_GENES. Use 'none' to skip."),
  make_option(c("--flank"), type = "double", default = 10,
              help = "kb flank around the gene body for 'nearby' variants [default: %default]"),
  make_option(c("--topk"), type = "integer", default = 20,
              help = "max variants listed per cell [default: %default]"),
  make_option("--tier-p", type = "double", default = 0.05,
              help = "pass_tier outlier gate: max pValue [default: %default]"),
  make_option("--tier-z", type = "double", default = 2,
              help = "pass_tier outlier gate: min |zScore| [default: %default]"),
  make_option("--tier-sai", type = "double", default = 0.2,
              help = "pass_tier scored gate: min SpliceAI [default: %default]"),
  make_option("--tier-cadd", type = "double", default = 20,
              help = "pass_tier scored gate: min CADD [default: %default]"),
  make_option("--tier-revel", type = "double", default = 0.5,
              help = "pass_tier scored gate: min REVEL [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

rds_file       <- opt$rds
sample         <- opt$sample
hpo_file       <- opt$hpo
gene_anno_file <- opt$annotation
outfile        <- if (!is.null(opt$output)) opt$output else paste0(sample, "_OUTRIDER_results_annotated.xlsx")

## ---- read OUTRIDER results -------------------------------------------------
# The rds may be either (a) a results data.frame/data.table already, or (b) a
# fitted OutriderDataSet object. In case (b), run OUTRIDER::results(ods, all=TRUE)
# to get the full per-gene results table (all=TRUE keeps every gene, not just the
# aberrant calls, so HPO/variant annotation and ranking see all rows).
raw <- readRDS(rds_file)
if (methods::is(raw, "OutriderDataSet")) {
  if (!requireNamespace("OUTRIDER", quietly = TRUE))
    stop("rds is an OutriderDataSet but the OUTRIDER package is not installed")
  message("rds is an OutriderDataSet; running OUTRIDER::results(ods, all=TRUE) ...")
  suppressMessages(library(OUTRIDER))
  res_all <- as.data.frame(OUTRIDER::results(raw, all = TRUE))
} else {
  res_all <- as.data.frame(raw)
}

if (!sample %in% res_all$sampleID) {
  stop(sprintf("Sample '%s' not found in %s. Available samples:\n%s",
               sample, rds_file, paste(unique(res_all$sampleID), collapse = ", ")))
}

# Extract the input sample's results.
res <- res_all[res_all$sampleID == sample, , drop = FALSE]
message(sprintf("Extracted %d gene rows for sample '%s'.", nrow(res), sample))

## ---- gene-level annotation -------------------------------------------------
gene_anno <- read.delim(gene_anno_file, sep = "\t", header = TRUE,
                        quote = "", stringsAsFactors = FALSE,
                        check.names = FALSE, na.strings = c("NA", ""))

# Collapse to one annotation row per gene_id so the join stays 1:1 (the table
# has a handful of duplicated gene_ids that would otherwise inflate the result).
# Transcript-level columns are concatenated; everything else takes the first value.
gene_anno <- gene_anno %>%
  group_by(gene_id) %>%
  summarise(
    across(c(RefSeq_nuc, Ensembl_nuc, MANE_status),
           ~ paste(unique(.x), collapse = ";")),
    across(-c(RefSeq_nuc, Ensembl_nuc, MANE_status),
           ~ .x[1]),
    .groups = "drop"
  )

# Strip any Ensembl version suffix so the ids match (e.g. ENSG000...1 -> ENSG000...).
res$geneID_nov <- sub("\\..*$", "", res$geneID)
gene_anno$gene_id_nov <- sub("\\..*$", "", gene_anno$gene_id)

res <- left_join(res, gene_anno, by = c("geneID_nov" = "gene_id_nov"))

## ---- gene coordinates (chr:start-end) --------------------------------------
# The gene-annotation table carries no coordinates, so they come from the -g
# coords file (a compact ensg,chrom,start,end TSV or an Ensembl GTF). Added
# whenever -g is supplied, independent of -V; the parsed table is reused below
# for the 'nearby variant' window. norm_chrom/read_gene_coords are defined here
# (hoisted above the variant block) so both sections share them.
norm_chrom <- function(x) {
  x <- ifelse(grepl("^chr", x), x, paste0("chr", x))
  ifelse(x == "chrMT", "chrM", x)
}
# gene coordinates: compact 6-col TSV (ensg,chrom,start,end,strand,name) or GTF.
read_gene_coords <- function(path) {
  con <- file(path, "r"); first <- NULL
  repeat {
    l <- readLines(con, n = 1, warn = FALSE)
    if (length(l) == 0) break
    if (!startsWith(l, "#")) { first <- l; break }
  }
  close(con)
  if (is.null(first)) stop("empty coords file: ", path)
  f1 <- strsplit(first, "\t")[[1]]
  compact <- length(f1) >= 6 && grepl("^ENSG", f1[1]) && grepl("^[0-9]+$", f1[3])
  if (compact) {
    d <- read.delim(path, header = FALSE, comment.char = "#", sep = "\t",
                    quote = "", stringsAsFactors = FALSE)
    d <- d[grepl("^ENSG", d[[1]]), , drop = FALSE]     # tolerate a header row
    data.frame(ensg  = sub("\\..*$", "", d[[1]]), chrom = norm_chrom(d[[2]]),
               start = as.integer(d[[3]]), end = as.integer(d[[4]]),
               stringsAsFactors = FALSE)
  } else {
    # GTF: pull 'gene' feature lines via awk (files can be multi-GB).
    cmd <- sprintf("awk -F'\\t' '$3==\"gene\"{print $1\"\\t\"$4\"\\t\"$5\"\\t\"$9}' %s",
                   shQuote(path))
    g <- read.delim(pipe(cmd), header = FALSE, sep = "\t", quote = "",
                    stringsAsFactors = FALSE)
    data.frame(ensg  = sub('.*gene_id "([^".]+).*', "\\1", g[[4]]),
               chrom = norm_chrom(g[[1]]),
               start = as.integer(g[[2]]), end = as.integer(g[[3]]),
               stringsAsFactors = FALSE)
  }
}

# The coords file now has a default, resolved the same way as the blacklist in
# parse_fraser.R: $GENE_COORDS, then the shared annotation_files copy. It used to
# be supplied by hand or not at all, which left gene_coord empty on most runs --
# and the dark-gene overlap below NEEDS coordinates to be correct (see there).
GENE_COORDS_DEFAULT <- paste0("/hpf/largeprojects/tgnode/sandbox/yliang_analysis/",
                              "pacbio_longread/annotation_files/",
                              "gene_coords.GRCh38.115.tsv")
find_gene_coords <- function() {
  cands <- character(0)
  if (nzchar(Sys.getenv("GENE_COORDS"))) cands <- c(cands, Sys.getenv("GENE_COORDS"))
  cands <- path.expand(c(cands, GENE_COORDS_DEFAULT))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

coords <- NULL
res$gene_coord <- NA_character_
{
  cf <- opt[["coords"]]
  coords_file <- if (!is.null(cf) && identical(tolower(cf), "none")) NA_character_
                 else if (!is.null(cf) && nzchar(cf)) cf
                 else find_gene_coords()
  if (!is.na(coords_file) && file.exists(path.expand(coords_file))) {
    coords_file <- path.expand(coords_file)
    coords <- read_gene_coords(coords_file)
    ci0 <- match(sub("\\..*$", "", res$geneID), coords$ensg)
    res$gene_coord <- ifelse(is.na(ci0), NA_character_,
                             sprintf("%s:%d-%d", coords$chrom[ci0],
                                     coords$start[ci0], coords$end[ci0]))
    message(sprintf("Added gene_coord for %d/%d genes from %s.",
                    sum(!is.na(res$gene_coord)), nrow(res), basename(coords_file)))
  } else if (!is.null(cf) && identical(tolower(cf), "none")) {
    message("Gene coordinates skipped (--coords none).")
  } else {
    warning("Gene coordinates not found; gene_coord and the positional half of ",
            "dark_gene will be NA. Pass --coords <file> or set $GENE_COORDS.")
  }
}

## ---- dark-gene panel (Alliance / PacBio long-read panel) -------------------
# Capture-target intervals over loci short reads align badly in (segmental
# duplications, paralogy, GC extremes, repeats). Kept in step with the same
# block in parse_fraser.R -- same file, same column name, same vocabulary.
#
# CRITERION, and it differs from FRASER's on purpose. A FRASER row is a junction:
# a point pair, in or out. An OUTRIDER row is a GENE, and its statistic is a
# count summed over the whole gene body, so the question is not "in or out" but
# HOW MUCH of the counted region is dark. That matters here: only 66% of the
# panel's named genes lie entirely inside their capture interval, and the 10th
# percentile of coverage is 0.093 -- a 9%-covered gene should not read the same
# as a 100%-covered one. So: flag on any overlap, and carry the covered fraction
# in the cell.
#
# The gene SYMBOL is not the primary key, for the same reason as in
# parse_fraser.R: the panel names `GBA` where Ensembl 115 says `GBA1`, so a
# symbol join silently drops it (20 of 704 panel symbols, 2.8%, no longer
# reconcile). Symbol match is the fallback only, tagged "gene-only".
DARK_DEFAULT <- paste0("/hpf/largeprojects/tgnode/sandbox/yliang_analysis/",
                       "pacbio_longread/annotation_files/",
                       "Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed")
find_dark_genes <- function() {
  cands <- character(0)
  if (nzchar(Sys.getenv("DARK_GENES"))) cands <- c(cands, Sys.getenv("DARK_GENES"))
  cands <- path.expand(c(cands, DARK_DEFAULT))
  hit <- cands[file.exists(cands)]
  if (length(hit)) hit[1] else NA_character_
}

res$dark_gene <- NA_character_
{
  dk <- opt[["dark-genes"]]
  if (!is.null(dk) && identical(tolower(dk), "none")) {
    message("Dark-gene annotation skipped (--dark-genes none).")
  } else {
    dark_file <- if (!is.null(dk)) dk else find_dark_genes()
    if (is.na(dark_file) || !file.exists(path.expand(dark_file))) {
      warning("Dark-gene BED not found; 'dark_gene' column left as NA. ",
              "Pass --dark-genes <file> or set $DARK_GENES.")
    } else {
      dark_file <- path.expand(dark_file)
      d <- read.table(dark_file, sep = "\t", header = FALSE, quote = "",
                      comment.char = "#", stringsAsFactors = FALSE)
      if (ncol(d) < 4)
        stop("dark-gene BED needs 4 columns (chrom, start, end, gene names): ",
             dark_file)
      d <- d[, 1:4]; names(d) <- c("chrom", "start0", "end", "names")
      d$chrom <- norm_chrom(sub("^chr", "", d$chrom))   # match coords' style
      d$start <- as.integer(d$start0) + 1L              # BED 0-based half-open
      d$end   <- as.integer(d$end)
      dsyms <- unique(trimws(unlist(strsplit(d$names, ","))))
      dsyms <- dsyms[nzchar(dsyms) & dsyms != "NA"]

      val <- character(nrow(res)); covf <- rep(NA_real_, nrow(res))
      if (!is.null(coords)) {
        gi <- match(sub("\\..*$", "", res$geneID), coords$ensg)
        gch <- coords$chrom[gi]; gs <- coords$start[gi]; ge <- coords$end[gi]
        ovbp <- rep(0, nrow(res))
        for (ch in intersect(unique(gch[!is.na(gch)]), unique(d$chrom))) {
          idx <- which(!is.na(gch) & gch == ch)
          if (!length(idx)) next
          rr <- d[d$chrom == ch, , drop = FALSE]
          for (k in seq_len(nrow(rr))) {
            # intersection length of the gene body with this capture interval
            lo <- pmax(gs[idx], rr$start[k]); hi <- pmin(ge[idx], rr$end[k])
            w  <- pmax(0, hi - lo + 1)
            hit <- w > 0
            if (!any(hit)) next
            j <- idx[hit]
            ovbp[j] <- ovbp[j] + w[hit]
            val[j] <- ifelse(nzchar(val[j]),
                             paste(val[j], rr$names[k], sep = ";"), rr$names[k])
          }
        }
        glen <- ge - gs + 1
        covf <- ifelse(nzchar(val) & !is.na(glen) & glen > 0,
                       pmin(1, ovbp / glen), NA_real_)
        val <- ifelse(nzchar(val),
                      sprintf("%s (%d%%)", val, as.integer(round(covf * 100))),
                      val)
      }
      # gene named on the panel but no positional overlap (or no coords at all)
      sym_only <- !nzchar(val) & !is.na(res$gene_symbol) &
                  res$gene_symbol %in% dsyms
      val[sym_only] <- paste0(res$gene_symbol[sym_only], ", gene-only")
      res$dark_gene <- val
      message(sprintf("Dark genes: %d regions; %d of %d rows overlap a region, %d more by gene name only.",
                      nrow(d), sum(nzchar(val)) - sum(sym_only), nrow(res),
                      sum(sym_only)))
    }
  }
}

## ---- individual HPO terms --------------------------------------------------
# One HPO id per line; strip whitespace and stray carriage-return artifacts.
ind_hpo <- character(0)
if (!is.null(hpo_file) && nzchar(hpo_file) && file.exists(hpo_file)) {
  ind_hpo <- readLines(hpo_file, warn = FALSE)
  ind_hpo <- trimws(gsub("_x000D_|\r", "", ind_hpo))
  ind_hpo <- unique(ind_hpo[grepl("^HP:", ind_hpo)])
}
have_hpo <- length(ind_hpo) > 0

## ---- HPO overlap between individual and each gene --------------------------
if (have_hpo) {
  message(sprintf("Individual '%s' has %d HPO terms.", sample, length(ind_hpo)))

  # Gene-associated HPO ids and their human-readable names are stored as aligned
  # pipe-separated lists in the annotation's hpo_id / hpo_name columns (same
  # element order). We match on ids, but report the names.
  gene_hpo_ids   <- strsplit(ifelse(is.na(res$hpo_id),   "", res$hpo_id),   "\\|")
  gene_hpo_names <- strsplit(ifelse(is.na(res$hpo_name), "", res$hpo_name), "\\|")

  # Distinct gene-associated HPO terms (set size, for the Jaccard denominator).
  n_gene_hpo <- vapply(gene_hpo_ids,
                       function(g) length(unique(trimws(g[trimws(g) != ""]))),
                       integer(1))

  # For each gene, keep the NAMES whose id is one of the individual's HPO terms.
  # Count overlaps directly (on distinct ids), not by splitting the name string,
  # so a name that itself contains ';' can't corrupt n_hpo_overlap.
  overlap_terms <- character(nrow(res))
  n_overlap     <- integer(nrow(res))
  for (i in seq_len(nrow(res))) {
    ids <- trimws(gene_hpo_ids[[i]])
    nms <- trimws(gene_hpo_names[[i]])
    k   <- min(length(ids), length(nms))          # guard against ragged rows
    if (k == 0) next
    ids <- ids[seq_len(k)]; nms <- nms[seq_len(k)]
    keep <- ids != "" & ids %in% ind_hpo & !duplicated(ids)
    if (any(keep)) {
      overlap_terms[i] <- paste(nms[keep], collapse = ";")
      n_overlap[i]     <- sum(keep)
    }
  }

  res$overlapping_hpo_terms <- overlap_terms
  res$n_hpo_overlap <- n_overlap
  res$hpo_overlap   <- ifelse(res$n_hpo_overlap > 0, "yes", "no")

  # Jaccard similarity between the individual's HPO set and the gene's HPO set:
  #   |intersection| / |union| = overlap / (n_individual + n_gene - overlap)
  union_size <- length(ind_hpo) + n_gene_hpo - res$n_hpo_overlap
  res$hpo_jaccard <- ifelse(union_size > 0, res$n_hpo_overlap / union_size, 0)
} else {
  # No HPO terms provided for this individual: leave the HPO columns as NA.
  message(sprintf("No HPO terms provided for '%s'; HPO overlap columns set to NA.", sample))
  res$overlapping_hpo_terms <- NA_character_
  res$n_hpo_overlap       <- NA_integer_
  res$hpo_overlap         <- NA_character_
  res$hpo_jaccard         <- NA_real_
}

## ---- rare-variant annotation (optional) ------------------------------------
# If a model_variants.tsv from trio_rare_filter.py is supplied, add per-gene
# rare-variant columns (ported from outrider_add_variants.py):
#   cand_var_n / cand_variants  variants VEP-annotated to this gene (ENSG match)
#   near_var_n / near_variants  variants within +/-flank kb of the gene body that
#                               are NOT annotated to it (regulatory/intergenic)
#   cand_models                 inheritance model(s) hit by the gene's variants
#   cand_max_impact             most severe VEP IMPACT among the gene's variants
# Columns are created up-front (NA) so they always exist for the select() below.
res$cand_var_n      <- NA_integer_
res$cand_variants   <- NA_character_
res$near_var_n      <- NA_integer_
res$near_variants   <- NA_character_
res$cand_models     <- NA_character_
res$cand_max_impact <- NA_character_
# hidden helpers for pass_tier (max candidate scores); dropped before output.
res$helper_sai      <- NA_real_
res$helper_cadd     <- NA_real_
res$helper_revel    <- NA_real_
res$helper_clinpath <- FALSE

if (!is.null(opt$variants) && nzchar(opt$variants) && file.exists(opt$variants)) {
  if (is.null(opt$coords) || !file.exists(opt$coords))
    stop("--variants given but --coords (gene coordinates) is missing/not found")
  message(sprintf("Adding rare-variant columns from %s ...", opt$variants))

  flank <- as.integer(opt$flank * 1000)
  top_k <- opt$topk

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

  # Rank of the most severe term in an '&'-joined consequence string.
  csq_rank <- function(conseq) {
    if (is.na(conseq) || conseq == "") return(MISS_RANK)
    r <- match(strsplit(conseq, "&", fixed = TRUE)[[1]], VEP_ORDER)
    r[is.na(r)] <- MISS_RANK
    min(r)
  }
  # spliceai_score: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|... -> max DS.
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
  # Main nomenclature token (vectorized): HGVSc (transcript:c./n. change) plus the
  # p. change when present. '.'/''/NA -> "" so it drops out of the cell cleanly.
  blank <- function(x) is.na(x) | x == "." | x == ""
  fmt_hgvs <- function(hgvsc, hgvsp) {
    c1 <- ifelse(blank(hgvsc), "", hgvsc)
    pp <- ifelse(blank(hgvsp), "", sub("^[^:]*:", "", hgvsp))   # drop ENSP..: prefix
    pp <- gsub("%3D", "=", pp, fixed = TRUE)                    # VEP encodes p.=
    ifelse(pp != "", ifelse(c1 != "", paste(c1, pp), pp), c1)
  }
  ## read the variant table and coerce fields
  vt <- read.delim(opt$variants, sep = "\t", header = TRUE, quote = "",
                   check.names = FALSE, stringsAsFactors = FALSE,
                   na.strings = c("NA", ""))
  for (c in c("hgvsc", "hgvsp"))          # tolerate an older TSV without HGVS
    if (is.null(vt[[c]])) vt[[c]] <- "."
  vt$chrom     <- norm_chrom(vt$chrom)
  vt$vid       <- paste(vt$chrom, vt$pos, vt$ref, vt$alt, sep = "_")
  vt$gene      <- sub("\\..*$", "", ifelse(is.na(vt$gene), "", vt$gene))
  vt$popmax    <- suppressWarnings(as.numeric(vt$popmax))
  vt$cadd      <- suppressWarnings(as.numeric(vt$cadd))
  vt$revel     <- suppressWarnings(as.numeric(vt$revel))
  vt$spliceaiN <- vapply(as.character(vt$spliceai), max_spliceai, numeric(1))
  vt$clinvar   <- ifelse(is.na(vt$clinvar) | vt$clinvar == ".", "", vt$clinvar)
  vt$gt        <- ifelse(is.na(vt$gt) | vt$gt == "", ".", vt$gt)
  vt$rank      <- vapply(as.character(vt$consequence), csq_rank, integer(1))

  message(sprintf("      %d variant rows, %d distinct variants",
                  nrow(vt), length(unique(vt$vid))))

  ## variant master: one row per variant. Coordinate/score fields are identical
  ## across a variant's transcripts, but HGVS is not, so pick the worst-consequence
  ## transcript as the representative (its HGVS is shown for 'nearby' variants).
  vt <- vt[order(vt$vid, vt$rank), ]
  master <- vt[!duplicated(vt$vid),
               c("vid","chrom","pos","ref","alt","gt","popmax","cadd",
                 "spliceaiN","clinvar","revel","hgvsc","hgvsp")]
  master$models <- tapply(vt$model, vt$vid,
                          function(m) paste(sort(unique(m)), collapse = ","))[master$vid]

  # precompute the display strings (space-separated: phased GT '0|1' carries '|')
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
  # nearby cell: head + representative HGVS + scores (collapse gaps from empties)
  m_near_cell <- gsub(" +", " ", trimws(paste(m_head, m_hgvs, m_suf)))
  m_cadd0 <- ifelse(is.na(master$cadd), 0, master$cadd)
  midx <- setNames(seq_len(nrow(master)), master$vid)

  ## gene association: worst consequence per (gene, variant), with that
  ## transcript's HGVS nomenclature
  ga <- vt[vt$gene != "", c("gene","vid","rank","consequence","impact","hgvsc","hgvsp")]
  ga <- ga[order(ga$gene, ga$vid, ga$rank), ]
  ga <- ga[!duplicated(ga[, c("gene","vid")]), ]
  ga_by_gene <- split(ga, ga$gene)

  ## position index (per chrom, ascending) for the nearby-window scan
  mo <- order(master$chrom, master$pos)
  pos_by_chrom <- split(master$pos[mo], master$chrom[mo])
  vid_by_chrom <- split(master$vid[mo], master$chrom[mo])

  # coords was parsed once above (for the gene_coord column); reuse it, but
  # re-read defensively in case gene_coord parsing was skipped.
  if (is.null(coords)) coords <- read_gene_coords(opt$coords)
  cidx <- match(sub("\\..*$", "", res$geneID), coords$ensg)

  ens <- sub("\\..*$", "", res$geneID)
  nr  <- nrow(res)
  v_cand_n <- integer(nr); v_cand <- character(nr)
  v_near_n <- integer(nr); v_near <- character(nr)
  v_models <- character(nr); v_impact <- character(nr)
  v_sai <- rep(NA_real_, nr); v_cadd <- rep(NA_real_, nr)
  v_revel <- rep(NA_real_, nr); v_clinpath <- logical(nr)

  for (i in seq_len(nr)) {
    a <- ga_by_gene[[ens[i]]]
    assoc_vids <- character(0)
    if (!is.null(a)) {
      ii <- midx[a$vid]
      a  <- a[order(a$rank, -m_cadd0[ii]), ]   # severity, then CADD desc
      ii <- midx[a$vid]
      assoc_vids <- a$vid
      # gene-specific nomenclature after the consequence/impact (collapse the gaps
      # left when a token is empty)
      cells <- gsub(" +", " ", trimws(paste(m_head[ii], a$consequence, a$impact,
                                             fmt_hgvs(a$hgvsc, a$hgvsp), m_suf[ii])))
      v_cand_n[i] <- length(cells)
      v_cand[i]   <- paste(head(cells, top_k), collapse = ";")
      v_models[i] <- paste(sort(unique(unlist(strsplit(master$models[ii], ",")))),
                           collapse = ",")
      ir <- IMPACT_RANK[a$impact]; ir[is.na(ir)] <- 9L
      v_impact[i] <- a$impact[which.min(ir)]
      v_sai[i]   <- suppressWarnings(max(master$spliceaiN[ii], na.rm = TRUE))
      v_cadd[i]  <- suppressWarnings(max(master$cadd[ii],      na.rm = TRUE))
      v_revel[i] <- suppressWarnings(max(master$revel[ii],     na.rm = TRUE))
      cvv <- master$clinvar[ii]
      v_clinpath[i] <- any(grepl("athogenic", cvv, ignore.case = TRUE) &
                           !grepl("onflicting", cvv, ignore.case = TRUE))
    }
    ci <- cidx[i]
    if (!is.na(ci)) {
      ch <- coords$chrom[ci]; pv <- pos_by_chrom[[ch]]
      if (!is.null(pv)) {
        sel <- which(pv >= coords$start[ci] - flank & pv <= coords$end[ci] + flank)
        nv  <- setdiff(vid_by_chrom[[ch]][sel], assoc_vids)
        if (length(nv)) {
          jj <- midx[nv]
          nv <- nv[order(-m_cadd0[jj])]
          jj <- midx[nv]
          v_near_n[i] <- length(nv)
          v_near[i]   <- paste(head(m_near_cell[jj], top_k), collapse = ";")
        }
      }
    }
  }

  # genes with no candidate variant get 0 / "" (not NA) once annotation has run;
  # NA is reserved for "no variant table supplied" (the up-front placeholders).
  res$cand_var_n      <- v_cand_n
  res$cand_variants   <- v_cand
  res$near_var_n      <- v_near_n
  res$near_variants   <- v_near
  res$cand_models     <- v_models
  res$cand_max_impact <- v_impact
  v_sai[is.infinite(v_sai)]     <- NA
  v_cadd[is.infinite(v_cadd)]   <- NA
  v_revel[is.infinite(v_revel)] <- NA
  res$helper_sai      <- v_sai
  res$helper_cadd     <- v_cadd
  res$helper_revel    <- v_revel
  res$helper_clinpath <- v_clinpath
  message(sprintf("      %d genes with an associated variant",
                  sum(res$cand_var_n > 0, na.rm = TRUE)))
}

## ---- pass_tier: outlier x variant-quality triage tier ----------------------
# Ordinal; the best tier a row qualifies for. Orthogonal to hpo_overlap (combine
# the two for the phenotype funnel). Non-outliers get "" so the autofilter can
# hide them. Thresholds are all --tier-* flags.
#   1_concordant  outlier + HIGH-impact variant AND under-expression (LoF -> NMD)
#   2_scored      outlier + a candidate variant passing a score gate
#   3_outlier     significant outlier, no qualifying variant
{
  tp <- opt[["tier-p"]]; tz <- opt[["tier-z"]]
  tsai <- opt[["tier-sai"]]; tcadd <- opt[["tier-cadd"]]; trevel <- opt[["tier-revel"]]
  is_out <- !is.na(res$pValue) & !is.na(res$zScore) & res$pValue < tp & abs(res$zScore) > tz
  conc <- is_out & !is.na(res$cand_max_impact) & res$cand_max_impact == "HIGH" & res$zScore < 0
  scored <- is_out & (
    (!is.na(res$helper_sai)   & res$helper_sai   >= tsai)   |
    (!is.na(res$helper_cadd)  & res$helper_cadd  >= tcadd)  |
    (!is.na(res$helper_revel) & res$helper_revel >= trevel) |
    (res$helper_clinpath %in% TRUE) |
    (!is.na(res$cand_max_impact) & res$cand_max_impact %in% c("HIGH", "MODERATE")))
  res$pass_tier <- ifelse(!is_out, "",
                   ifelse(conc, "1_concordant",
                   ifelse(scored, "2_scored", "3_outlier")))
}

## ---- ranking helper flags --------------------------------------------------
# Does the gene have ANY disease annotation (OMIM/Orphanet/ClinGen/HPO)?
nonblank <- function(x) !is.na(x) & trimws(as.character(x)) != ""
res$has_disease_anno <- ifelse(
  nonblank(res$omim_phenotype) | nonblank(res$Orphanet) |
  nonblank(res$clinGen_disease_label) | nonblank(res$hpo_id),
  "yes", "no")

## ---- arrange and write -----------------------------------------------------
# Sort priority:
#   1. any HPO overlap (yes first)
#   2. pValue (smallest first), then |zScore| (largest first)
# The final row order is then frozen into a 1-based `rank` column (col 1), so the
# ranking survives re-sorting the sheet in Excel.
res <- res %>%
  dplyr::select(gene_symbol, geneID, gene_coord, dark_gene, pass_tier,
                hpo_overlap, n_hpo_overlap, hpo_jaccard,
                overlapping_hpo_terms, has_disease_anno,
                omim_phenotype, omim_inheritance, Orphanet, clinGen_disease_label,
                pValue, padjust, zScore, l2fc,
                cand_var_n, cand_variants, near_var_n, near_variants,
                cand_models, cand_max_impact,
                rawcounts, normcounts, meanCorrected,
                aberrant, AberrantBySample, AberrantByGene, padj_rank,
                everything()) %>%
  # drop working/unwanted columns; any_of() tolerates any that are absent
  dplyr::select(-dplyr::any_of(c("sampleID", "geneID_nov", "gene_id",
                                 "theta", "hgnc_symbols",
                                 "helper_sai", "helper_cadd", "helper_revel", "helper_clinpath",
                                 "GO.term.accession", "GO.domain"))) %>%
  arrange(desc(hpo_overlap == "yes"),
          pValue,
          desc(abs(zScore))) %>%
  dplyr::mutate(rank = dplyr::row_number(), .before = 1)

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
write_output <- function(res, outfile) {
  is_xlsx  <- grepl("\\.xlsx$", outfile, ignore.case = TRUE)
  tsv_path <- if (is_xlsx) paste0(tools::file_path_sans_ext(outfile), ".tsv") else outfile
  write.table(res, tsv_path, col.names = TRUE, row.names = FALSE,
              sep = "\t", quote = FALSE, na = "NA")
  if (!is_xlsx) return(invisible(tsv_path))
  conv <- find_converter()
  if (is.na(conv))
    stop("output is .xlsx but tsv_to_xlsx.py was not found (set $TSV_TO_XLSX or ",
         "keep it beside this script). TSV written to ", tsv_path)
  py <- if (nzchar(Sys.which("python3"))) "python3" else "python"
  st <- system2(py, c(shQuote(conv), shQuote(tsv_path), shQuote(outfile)))
  if (st != 0)
    stop("tsv_to_xlsx.py failed (exit ", st, "); intermediate TSV kept at ", tsv_path)
  invisible(outfile)
}
write_output(res, outfile)
hpo_note <- if (have_hpo) {
  sprintf("%d genes with HPO overlap", sum(res$hpo_overlap == "yes", na.rm = TRUE))
} else {
  "HPO columns set to NA (no terms provided)"
}
message(sprintf("Wrote %d rows to %s (%s).", nrow(res), outfile, hpo_note))

