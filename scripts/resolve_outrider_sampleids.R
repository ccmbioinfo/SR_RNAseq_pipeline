#!/usr/bin/env Rscript
# Map an OUTRIDER results table's sampleID values to a "plain" sample ID (the
# basename after the last '/'), because some control-cohort OUTRIDER tables
# carry the FULL source directory as sampleID rather than the bare sample name
# -- e.g. control_cohort/Sept2026/Blood/blood_bestQres.rds: 134 of 175 samples
# have sampleID "/hpf/largeprojects/tgnode/sandbox/rnaseq_repo/ComplexCare/21-8881"
# instead of "21-8881" (Fibroblasts/LCL's *_bestQres.rds do NOT have this -- it
# is inconsistent across tissues/cohort builds, so re-check, don't assume).
# parse_outrider.R's `-s` does an EXACT match against this column and hard-errors
# if it isn't found, so submit_parse_control_cohort.sh looks up the literal
# string to pass here instead of assuming the highQual-list ID is what's in the
# table.
#
# Usage: Rscript resolve_outrider_sampleids.R <rds_file> <out_tsv>
# Writes a 2-col TSV: plain_id<TAB>literal_sampleid

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: resolve_outrider_sampleids.R <rds_file> <out_tsv>")
}
rds_file <- args[1]
out_tsv  <- args[2]

raw <- readRDS(rds_file)
if (methods::is(raw, "OutriderDataSet")) {
  if (!requireNamespace("OUTRIDER", quietly = TRUE)) {
    stop("rds is an OutriderDataSet but the OUTRIDER package is not installed")
  }
  suppressMessages(library(OUTRIDER))
  ids <- unique(as.data.frame(OUTRIDER::results(raw, all = TRUE))$sampleID)
} else {
  ids <- unique(as.data.frame(raw)$sampleID)
}

# basename() passes a slash-free ID through unchanged, and strips a directory
# path down to its last component -- exactly the two cases seen in practice.
plain <- basename(ids)

dup <- unique(plain[duplicated(plain)])
if (length(dup)) {
  warning("Ambiguous plain sample ID(s) resolve to >1 literal sampleID in ",
          rds_file, "; keeping the first occurrence: ", paste(dup, collapse = ", "))
}

map <- data.frame(plain_id = plain, literal_sampleid = ids, stringsAsFactors = FALSE)
map <- map[!duplicated(map$plain_id), , drop = FALSE]

write.table(map, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("Wrote %d sample ID mapping(s) to %s", nrow(map), out_tsv))
