# fraser_impact_annotation.R ---------------------------------------------------
#
# Splicing-consequence annotation for a FRASER results table, wrapping the same
# calls DROP makes in modules/aberrant-splicing-pipeline/FRASER/08_extract_results_FraseR.R:
#     annotatePotentialImpact(result=, txdb=, fds=)
#     flagBlacklistRegions(result=, assemblyVersion=)
#
# Sourced by BOTH parse_fraser.R (--txdb + --fds-*, annotating the filtered
# shortlist inline) and annotate_splice_impact.R (annotating an existing
# .filtered.tsv after the fact). One copy so the two cannot drift.
#
# This file pulls in FRASER/Bioconductor, so it is sourced LAZILY -- only when
# the caller actually asks for impact annotation. parse_fraser.R keeps its
# tidyverse+optparse-only dependency set on the default path.
#
# annotatePotentialImpact() needs the FraserDataSet itself, not just a table: it
# pulls every junction overlapping the query region out of the fds plus their
# counts, because classifying exon skipping vs truncation vs elongation requires
# knowing which neighbouring introns are actually expressed.

.fia_libs <- function() {
  suppressPackageStartupMessages({
    library(FRASER); library(AnnotationDbi); library(data.table)
  })
}

#' Load the fds either directly by file or DROP-style by dir + names.
#'
#' Our own cohort pipeline saves as <WD>/savedObjects/<name>/, DROP as
#' <dir>/savedObjects/<dataset>--<annotation>/, hence both forms.
#' NB since the 2026-08-28 counts/fit split (see fraser-cohort-per-sample-outputs),
#' <name> for a FITTED object is "<name>_fit" (e.g. savedObjects/cohort_fit/),
#' NOT savedObjects/cohort/ -- that one is counts-only (pre-filter, pre-fit,
#' no hgnc_symbol/annotatedJunction/deltaPsi) and a per-sample Snakemake run may
#' have already deleted it by the time this runs. Point --fds-file at the _fit
#' one. Verified 2026-09-03: loadFraserDataSet(file=".../cohort_fit/fds-object.RDS")
#' loads with hgnc_symbol and annotatedJunction present.
load_fds_for_impact <- function(fds_file = NULL, fds_dir = NULL,
                                dataset = NULL, annotation = NULL) {
  .fia_libs()
  if (!is.null(fds_file)) {
    ff <- path.expand(fds_file)
    if (!file.exists(ff)) stop("fds not found: ", ff)
    message("  loading fds from file: ", ff)
    return(loadFraserDataSet(file = ff))
  }
  if (is.null(fds_dir) || is.null(dataset) || is.null(annotation))
    stop("provide fds_file, or all of fds_dir, dataset and annotation")
  nm <- paste(dataset, annotation, sep = "--")
  message("  loading fds ", nm, " from ", fds_dir)
  loadFraserDataSet(dir = path.expand(fds_dir), name = nm)
}

#' Add potentialImpact / causesFrameshift / UTR_overlap to a results table.
#'
#' Junctions are deduplicated before annotation and joined back, so the cost
#' scales with unique junctions rather than gene-exploded rows.
#'
#' @param df   data.frame/data.table with seqnames,start,end,strand,type,annotatedJunction
#' @param add_blacklist  also run FRASER's own flagBlacklistRegions, stored as
#'   blacklist_fraser -- a cross-check against the native `blacklist` column
#'   parse_fraser.R computes, not a replacement for it.
#' @return df with the new columns appended (input row order preserved)
add_splice_impact <- function(df, fds, txdb, add_blacklist = FALSE,
                              assembly = "hg38") {
  .fia_libs()
  dt <- as.data.table(df)
  if (nrow(dt) == 0L) return(df)

  # deltaPsi is as hard a requirement as annotatedJunction: addPotentialImpactLabels()
  # splits annotatedJunction=="both" into increased/reducedUsage on its sign.
  need <- c("seqnames", "start", "end", "strand", "type", "annotatedJunction",
            "deltaPsi")
  absent <- setdiff(need, colnames(dt))
  if (length(absent))
    stop("cannot annotate splice impact: missing column(s) ",
         paste(absent, collapse = ", "),
         "\n  annotatedJunction comes from FRASER's annotateIntronReferenceOverlap;",
         "\n  run fraser_cohort.R fit with --gtf (or annotate_fraser_genes.R) first.")

  # Optional intron-retention refinement. FRASER gates it on only the first two
  # being present but then also reads nonsplitCounts, so carry the trio or none
  # of it -- passing 2 of 3 reproduces exactly the "Object not found" error this
  # subsetting caused in the first place.
  retention <- c("nonsplitProportion", "nonsplitProportion_99quantile",
                 "nonsplitCounts")
  if (!all(retention %in% colnames(dt))) {
    message("  note: ", paste(setdiff(retention, colnames(dt)), collapse = ", "),
            " absent -- skipping the (partial)intronRetention refinement")
    retention <- character(0)
  }

  key   <- c("seqnames", "start", "end", "strand")
  carry <- c(key, "type", "annotatedJunction", "deltaPsi", retention)
  # unique() by the junction key alone: the caller's rows are gene-exploded, so
  # one junction can appear many times with identical stats, and deduping on the
  # full carry set would let a stray difference inflate the join.
  juncs <- unique(dt[, carry, with = FALSE], by = key)
  message(sprintf("  %d unique junctions to annotate", nrow(juncs)))

  before <- colnames(juncs)
  juncs  <- as.data.table(annotatePotentialImpact(result = juncs, txdb = txdb,
                                                  fds = fds))
  if (isTRUE(add_blacklist)) {
    juncs <- as.data.table(flagBlacklistRegions(result = juncs,
                                                assemblyVersion = assembly))
    if ("blacklist" %in% colnames(juncs))
      setnames(juncs, "blacklist", "blacklist_fraser")
  }
  added <- setdiff(colnames(juncs), before)
  if (!length(added)) stop("annotatePotentialImpact() added no new columns")
  message("  added: ", paste(added, collapse = ", "))

  # Drop any pre-existing copy so the join cannot produce .x/.y duplicates.
  dupe <- intersect(added, colnames(dt))
  if (length(dupe)) {
    message("  overwriting pre-existing column(s): ", paste(dupe, collapse = ", "))
    dt[, (dupe) := NULL]
  }
  # merge() reorders BOTH ways, and both have to be undone.
  #   rows:    sort=FALSE is not enough, hence the ..fia_row index.
  #   columns: data.table::merge puts the `by` columns FIRST. That silently
  #            rewrote the caller's layout -- seqnames/start/end/strand jumped to
  #            the front and `rank`, which parse_fraser.R deliberately places at
  #            column 1, slid to 5. Only --txdb runs were affected, so the
  #            documented positions in outlier_schema.html were correct for a
  #            default run and wrong for every annotated one.
  col_order <- setdiff(colnames(dt), "..fia_row")
  dt[, ..fia_row := .I]
  out <- merge(dt, juncs[, c(key, added), with = FALSE], by = key,
               all.x = TRUE, sort = FALSE)
  setorder(out, ..fia_row)
  out[, ..fia_row := NULL]
  want <- c(col_order, added)
  if (!setequal(want, colnames(out)))
    stop("column bookkeeping error in add_splice_impact(): expected ",
         length(want), " columns, merge produced ", ncol(out))
  setcolorder(out, want)

  for (cl in added) {
    v <- out[[cl]]
    n <- if (is.logical(v)) sum(v, na.rm = TRUE) else sum(!is.na(v))
    message(sprintf("  %-22s %d / %d rows", cl, n, nrow(out)))
  }
  as.data.frame(out)
}
