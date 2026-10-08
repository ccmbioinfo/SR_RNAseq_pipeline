# fraser_gene_annotation.R -----------------------------------------------------
#
# Gene annotation for a FraserDataSet, replaying what DROP does in
#   template/Scripts/Pipeline/preprocessGeneAnnotation.R   (GTF -> txdb + mapping)
#   modules/aberrant-splicing-pipeline/FRASER/06_annotate_genes.R  (annotate fds)
#
# Sourced by BOTH fraser_cohort.R (fit mode, --gtf) and annotate_fraser_genes.R
# (retrofitting a cohort that was already fit). Kept in one place so the two
# cannot drift -- the seqlevelsStyle handling in particular is a silent-failure
# trap and must not exist in two copies.
#
# Adds two things the cohort pipeline otherwise lacks:
#   hgnc_symbol       -> surfaces as hgncSymbol in results(); parse_fraser.R needs it
#   annotatedJunction -> required by FRASER::annotatePotentialImpact(), which
#                        hard-errors without it

.fga_libs <- function() {
  suppressPackageStartupMessages({
    library(GenomicFeatures); library(txdbmaker); library(rtracklayer)
    library(AnnotationDbi);   library(GenomeInfoDb); library(data.table)
  })
}

#' Resolve where a GTF's txdb + gene-name mapping are cached.
#'
#' DEFAULT IS BESIDE THE GTF, not in the workdir. makeTxDbFromGFF on a full
#' Ensembl annotation is slow and memory-hungry, and NOTHING about it depends on
#' the cohort -- so it is built once next to the annotation and every workdir
#' just loads it. Files are named after the GTF, so several annotations can
#' share one directory. Falls back to `dir` (the workdir) if the GTF's directory
#' is not writable.
#'
#' @return list(txdb_file=, map_file=, dir=, stem=)
txdb_paths <- function(gtf, dir = NULL) {
  gtf  <- path.expand(gtf)
  stem <- sub("\\.gz$", "", basename(gtf), ignore.case = TRUE)
  stem <- sub("\\.(gtf|gff3|gff)$", "", stem, ignore.case = TRUE)
  gdir <- dirname(gtf)
  base <- gdir
  if (file.access(gdir, 2) != 0) {
    if (is.null(dir))
      stop("cannot write next to the GTF (", gdir, ") and no fallback dir given")
    message("NOTE: ", gdir, " is not writable -> caching txdb in ", dir)
    base <- dir
  }
  list(txdb_file = file.path(base, paste0(stem, ".txdb.db")),
       map_file  = file.path(base, paste0(stem, ".gene_name_mapping.tsv")),
       dir = base, stem = stem)
}

#' Build (or load) a txdb + gene-name mapping from a GTF.
#'
#' Build it ONCE with build_fraser_txdb.R; every fit afterwards just loads it.
#' Both artefacts are needed downstream -- the TxDb for FRASER's
#' annotateRangesWithTxDb()/annotateIntronReferenceOverlap(), the mapping table
#' for the gene_id -> gene_name lookup -- so they are always resolved as a pair.
#'
#' @param gtf       GTF path. Optional when the cache already exists.
#' @param dir       fallback cache dir (the workdir) if the GTF dir is read-only.
#' @param rebuild   rebuild even if the cache exists.
#' @param txdb_file explicit path to a prebuilt txdb.db (skips path resolution).
#' @param map_file  explicit path to the gene-name mapping; defaults to the
#'                  sibling of `txdb_file` (<stem>.gene_name_mapping.tsv).
#' @return list(txdb=<TxDb>, map=<data.table>, txdb_file=, map_file=)
build_txdb_and_map <- function(gtf = NULL, dir = NULL, rebuild = FALSE,
                               txdb_file = NULL, map_file = NULL) {
  .fga_libs()

  if (is.null(txdb_file)) {
    if (is.null(gtf))
      stop("build_txdb_and_map(): need either `gtf` or `txdb_file`")
    p <- txdb_paths(gtf, dir)
    txdb_file <- p$txdb_file
    if (is.null(map_file)) map_file <- p$map_file
  } else {
    txdb_file <- path.expand(txdb_file)
    if (is.null(map_file)) {
      map_file <- sub("\\.txdb\\.db$", ".gene_name_mapping.tsv", txdb_file)
      # a txdb.db that does not follow our naming (e.g. DROP's plain "txdb.db")
      if (identical(map_file, txdb_file))
        map_file <- file.path(dirname(txdb_file), "gene_name_mapping.tsv")
    } else {
      map_file <- path.expand(map_file)
    }
  }

  cached <- file.exists(txdb_file) && file.exists(map_file)
  if (!cached && is.null(gtf))
    stop("no cached annotation found at\n  ", txdb_file, "\n  ", map_file,
         "\nBuild it once with:\n",
         "  Rscript --vanilla build_fraser_txdb.R --gtf <annotation.gtf>\n",
         "or pass --gtf to build it inline.")
  if (!is.null(gtf)) {
    gtf <- path.expand(gtf)
    if (!file.exists(gtf)) stop("GTF not found: ", gtf)
  }
  dir.create(dirname(txdb_file), recursive = TRUE, showWarnings = FALSE)

  if (file.exists(txdb_file) && !isTRUE(rebuild)) {
    message("Reusing txdb: ", txdb_file)
    txdb <- loadDb(txdb_file)
  } else {
    message("Building txdb from ", basename(gtf), " (slow, minutes) ...")
    txdb <- makeTxDbFromGFF(gtf)
    txdb <- keepStandardChromosomes(txdb)
    tmp <- tempfile(); saveDb(txdb, tmp)
    file.copy(tmp, txdb_file, overwrite = TRUE); unlink(tmp)
    txdb <- loadDb(txdb_file)
    message("  wrote ", txdb_file)
  }

  if (file.exists(map_file) && !isTRUE(rebuild)) {
    message("Reusing gene-name mapping: ", map_file)
    map <- fread(map_file)
  } else {
    message("Extracting gene-name mapping from the GTF ...")
    gtf_dt <- as.data.table(rtracklayer::import(gtf))
    if (!"gene_name" %in% colnames(gtf_dt)) gtf_dt[, gene_name := gene_id]
    if ("gene_biotype" %in% colnames(gtf_dt)) setnames(gtf_dt, "gene_biotype", "gene_type")
    gtf_dt <- gtf_dt[type == "gene"]
    cols <- intersect(c("seqnames", "start", "end", "strand",
                        "gene_id", "gene_name", "gene_type"), colnames(gtf_dt))
    gtf_dt <- gtf_dt[, cols, with = FALSE]
    # Make gene_names unique exactly as DROP does. THIS is where the "_2" /
    # "_34608" suffixes in hgncSymbol come from, which parse_fraser.R strips.
    gtf_dt[, N := seq_len(.N), by = gene_name]
    gtf_dt[, gene_name_orig := gene_name]
    gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = "_")]
    gtf_dt[, N := NULL]
    # rtracklayer::import yields seqnames/strand as FACTORS; coerce BEFORE the
    # write so the file is plain unquoted text (fwrite quotes factor columns),
    # byte-identical to what DROP's fread-based path produces.
    for (cc in intersect(c("seqnames", "strand"), colnames(gtf_dt)))
      gtf_dt[[cc]] <- as.character(gtf_dt[[cc]])
    fwrite(gtf_dt, map_file, sep = "\t", na = "NA")
    map <- gtf_dt
    message(sprintf("  %d genes -> %s", nrow(map), map_file))
  }

  # Defensive, for the CACHE-REUSE path. fwrite(na="NA") happens to quote the
  # character columns, so fread currently returns seqnames as character even for
  # Ensembl "1".."22" -- but that is a side effect of the quoting, not a
  # guarantee, and an unquoted mapping table (hand-edited, or written by another
  # tool) would come back INTEGER, on which seqlevelsStyle<- misbehaves silently.
  # Idempotent for the fresh-build path. Cheap insurance now that the cache is
  # shared across workdirs and the reuse path is the normal one.
  for (cc in intersect(c("seqnames", "strand"), colnames(map)))
    map[[cc]] <- as.character(map[[cc]])

  list(txdb = txdb, map = map, txdb_file = txdb_file, map_file = map_file)
}

#' Annotate a FraserDataSet with gene symbols + reference-overlap.
#'
#' Must run AFTER the fit but BEFORE results() -- annotatedJunction is only
#' carried into a results table if it is on the fds when results() is called.
#'
#' @return the annotated fds (NOT saved; caller decides where/whether to save)
annotate_fds_genes <- function(fds, txdb, map) {
  .fga_libs()
  # Match chromosome naming across fds / txdb / mapping. Our cohort is Ensembl
  # ("22"); a UCSC-style GTF is "chr22" and EVERY overlap would silently return
  # nothing -- no error, just an empty annotation. DROP does the same coercion.
  style <- seqlevelsStyle(fds)[1]
  message("Matching seqlevelsStyle to the fds: ", style)
  suppressWarnings({
    seqlevelsStyle(map$seqnames) <- style
    seqlevelsStyle(txdb)         <- style
  })

  message("annotateRangesWithTxDb (gene symbols) ...")
  fds <- annotateRangesWithTxDb(fds, txdb = txdb, orgDb = map,
                                feature = "gene_name",
                                featureName = "hgnc_symbol",
                                keytype = "gene_id")

  message("annotateIntronReferenceOverlap (annotatedJunction) ...")
  fds <- annotateIntronReferenceOverlap(fds, txdb)

  # Report the fill rate: a near-zero rate is the signature of a seqlevelsStyle
  # or gene_id-format mismatch, which is otherwise completely silent.
  sym <- tryCatch(as.character(mcols(rowRanges(fds, type = "j"))$hgnc_symbol),
                  error = function(e) character(0))
  if (length(sym)) {
    n_ok <- sum(!is.na(sym) & nzchar(sym))
    message(sprintf("  hgnc_symbol set on %d / %d junctions (%.1f%%)",
                    n_ok, length(sym), 100 * n_ok / length(sym)))
    if (n_ok == 0)
      warning("hgnc_symbol is empty for EVERY junction -- check that the GTF's ",
              "chromosome naming and gene_id format match the fds.")
  }
  fds
}
