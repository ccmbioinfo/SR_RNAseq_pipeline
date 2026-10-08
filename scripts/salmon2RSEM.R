#!/usr/bin/env Rscript

# ============================================================
# salmon_to_rsem.R
#
# Convert a single Salmon quant.sf file to an RSEM-like
# gene-level results file using tximport.
#
# Output columns:
#   gene_id
#   transcript_id(s)
#   length
#   effective_length
#   expected_count
#   TPM
#
# Usage:
#
# Rscript salmon_to_rsem.R \
#     --tx2gene /path/to/tx2gene.txt \
#     --salmon_dir /path/to/salmon_sample \
#     --output_dir /path/to/output
#
# Output:
#   <sampleID>_salmon_rsem.results
#
# ============================================================

suppressPackageStartupMessages({
    library(optparse)
    library(tximport)
    library(readr)
    library(dplyr)
})

# ------------------------------------------------------------
# Command-line options
# ------------------------------------------------------------

option_list <- list(

    make_option(
        c("-t", "--tx2gene"),
        type = "character",
        help = "Path to tx2gene.txt",
        metavar = "FILE"
    ),

    make_option(
        c("-s", "--salmon_dir"),
        type = "character",
        help = "Directory containing quant.sf",
        metavar = "DIR"
    ),

    make_option(
        c("-o", "--output_dir"),
        type = "character",
        help = "Output directory",
        metavar = "DIR"
    ),
    make_option(
      c("-c", "--countsFrom"),
      type = "character",
      help = "Argument passed to countsFromAbundance, default set to 'no', can be set to 'lengthScaledTPM' or 'scaledTPM'",
      default = "no"
    ),
    make_option(
      c("--genelist"),
      type = "character",
      help = "List of gene IDs exist in gtf files, only these genes are reported. Default to '/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_id_name_biotype_gtf115_new.tsv'",
      default = "/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_id_name_biotype_gtf115_new.tsv"
    )
)

parser <- OptionParser(
    option_list = option_list,
    description = "Convert Salmon quant.sf to RSEM-like gene-level output"
)

opt <- parse_args(parser)

# ------------------------------------------------------------
# Check arguments
# ------------------------------------------------------------

if (is.null(opt$tx2gene) ||
    is.null(opt$salmon_dir) ||
    is.null(opt$output_dir)) {

    print_help(parser)

    stop(
        "\n--tx2gene, --salmon_dir, and --output_dir are required.",
        call. = FALSE
    )
}

# ------------------------------------------------------------
# Define paths
# ------------------------------------------------------------

tx2gene_file <- opt$tx2gene
salmon_dir <- opt$salmon_dir
output_dir <- opt$output_dir

salmon_file <- file.path(
    salmon_dir,
    "quant.sf"
)

genelist_file <- opt$genelist

# ------------------------------------------------------------
# Determine sample ID
# ------------------------------------------------------------

sample_id <- basename(dirname(
    normalizePath(
        salmon_dir,
        mustWork = TRUE
    )
))

output_file <- file.path(
    output_dir,
    paste0(
        sample_id,
        "_salmon_rsem.results"
    )
)

# ------------------------------------------------------------
# Check inputs
# ------------------------------------------------------------

if (!file.exists(tx2gene_file)) {
    stop(
        "tx2gene file does not exist: ",
        tx2gene_file
    )
}

if (!dir.exists(salmon_dir)) {
    stop(
        "Salmon directory does not exist: ",
        salmon_dir
    )
}

if (!file.exists(salmon_file)) {
    stop(
        "Cannot find quant.sf in: ",
        salmon_dir
    )
}

if (!dir.exists(output_dir)) {
    dir.create(
        output_dir,
        recursive = TRUE
    )
}

# ------------------------------------------------------------
# Read tx2gene
# ------------------------------------------------------------

cat(
    "\nReading tx2gene:\n  ",
    tx2gene_file,
    "\n",
    sep = ""
)

tx2gene <- read.table(
    tx2gene_file,
    col.names = c("TXNAME","GENEID")
)

# Check columns
required_tx2gene <- c(
    "TXNAME",
    "GENEID"
)

missing_tx2gene <- setdiff(
    required_tx2gene,
    colnames(tx2gene)
)

if (length(missing_tx2gene) > 0) {
    stop(
        "tx2gene.txt is missing required columns: ",
        paste(
            missing_tx2gene,
            collapse = ", "
        )
    )
}

# ------------------------------------------------------------
# Remove transcript versions
# ------------------------------------------------------------

tx2gene <- tx2gene %>%
    mutate(
        TXNAME = sub(
            "\\.[0-9]+$",
            "",
            TXNAME
        )
    ) %>%
    select(
        TXNAME,
        GENEID
    )

# Remove duplicate transcript mappings
tx2gene <- tx2gene %>%
    distinct(
        TXNAME,
        .keep_all = TRUE
    )

tx2gene_df <- as.data.frame(tx2gene)

# ------------------------------------------------------------
# Read Salmon quant.sf
#
# We read this separately because we need:
#   - transcript IDs
#   - transcript lengths
#   - effective lengths
#   - transcript -> gene relationships
# ------------------------------------------------------------

cat(
    "\nReading Salmon quant.sf:\n  ",
    salmon_file,
    "\n",
    sep = ""
)

salmon <- read_tsv(
    salmon_file,
    show_col_types = FALSE
)

required_salmon <- c(
    "Name",
    "Length",
    "EffectiveLength",
    "TPM",
    "NumReads"
)

missing_salmon <- setdiff(
    required_salmon,
    colnames(salmon)
)

if (length(missing_salmon) > 0) {
    stop(
        "quant.sf is missing required columns: ",
        paste(
            missing_salmon,
            collapse = ", "
        )
    )
}

# Remove transcript versions
salmon <- salmon %>%
    mutate(
        transcript_id = sub(
            "\\.[0-9]+$",
            "",
            Name
        )
    )

# ------------------------------------------------------------
# Map transcripts to genes
# ------------------------------------------------------------

salmon <- salmon %>%
    left_join(
        tx2gene,
        by = c(
            "transcript_id" = "TXNAME"
        )
    ) %>%
    rename(
        gene_id = GENEID
    )

# ------------------------------------------------------------
# Identify transcripts missing from tx2gene
# ------------------------------------------------------------

missing_transcripts <- salmon %>%
    filter(is.na(gene_id)) %>%
    select(
        transcript_id
    ) %>%
    distinct() %>%
    arrange(transcript_id)

missing_file <- file.path(
    output_dir,
    paste0(
        sample_id,
        "_transcripts_missing_tx2gene.txt"
    )
)

if(nrow(missing_transcripts) >0){
  write_tsv(
    missing_transcripts,
    missing_file,
    col_names = FALSE
  )
}


n_total <- nrow(salmon)

n_missing <- sum(
    is.na(salmon$gene_id)
)

cat(
    "Total transcripts: ",
    n_total,
    "\n",
    sep = ""
)

cat(
    "Transcripts without gene ID: ",
    n_missing,
    "\n",
    sep = ""
)

if (n_missing > 0) {
    warning(
        n_missing,
        " transcripts were not found in tx2gene.txt."
    )
}

# Keep mapped transcripts
salmon_mapped <- salmon %>%
    filter(
        !is.na(gene_id)
    )

# ------------------------------------------------------------
# Import Salmon using tximport
# ------------------------------------------------------------

cat(
    "\nRunning tximport...\n"
)

files <- salmon_file
names(files) <- sample_id

txi <- tximport(
    files,
    type = "salmon",
    tx2gene = tx2gene_df,
    ignoreTxVersion = TRUE,
    countsFromAbundance = opt$countsFrom
)

# ------------------------------------------------------------
# Transcript lists for each gene
# ------------------------------------------------------------

transcript_lists <- salmon_mapped %>%
    group_by(gene_id) %>%
    summarise(
        `transcript_id(s)` = paste(
            transcript_id,
            collapse = ","
        ),
        .groups = "drop"
    )

# ------------------------------------------------------------
# Gene-level length and effective length
#
# RSEM defines both as isoform-fraction-weighted averages of the
# transcript values, so the weight is relative abundance (TPM), not
# read count -- NumReads is proportional to TPM * EffectiveLength and
# would bias both columns toward long isoforms.
#
# effective_length comes straight from tximport, which does exactly
# that weighting and substitutes the unweighted mean for genes with
# zero abundance. Weighting by hand there gives 0/0 = NaN, which is
# most of the annotation on a typical sample.
# ------------------------------------------------------------

eff_lengths <- data.frame(

    gene_id = rownames(txi$length),

    effective_length = as.numeric(
        txi$length[, sample_id]
    ),

    stringsAsFactors = FALSE

)

# `length` gets the same weighting, applied to the un-effective lengths,
# with the same fallback for genes carrying no signal.
gene_lengths <- salmon_mapped %>%
    group_by(gene_id) %>%
    summarise(

        length = if (sum(TPM, na.rm = TRUE) > 0) {
            weighted.mean(
                Length,
                w = TPM,
                na.rm = TRUE
            )
        } else {
            mean(
                Length,
                na.rm = TRUE
            )
        },

        .groups = "drop"
    )

# ------------------------------------------------------------
# Gene-level results from tximport
# ------------------------------------------------------------

gene_results <- data.frame(

    gene_id = rownames(txi$counts),

    expected_count = as.numeric(
        txi$counts[, sample_id]
    ),

    TPM = as.numeric(
        txi$abundance[, sample_id]
    ),

    stringsAsFactors = FALSE

)

# ------------------------------------------------------------
# Add transcript IDs and lengths
# ------------------------------------------------------------

gene_results <- gene_results %>%
    left_join(
        transcript_lists,
        by = "gene_id"
    ) %>%
    left_join(
        gene_lengths,
        by = "gene_id"
    ) %>%
    left_join(
        eff_lengths,
        by = "gene_id"
    )

# ------------------------------------------------------------
# Reorder columns
# ------------------------------------------------------------

gene_results <- gene_results %>%
    select(
        gene_id,
        `transcript_id(s)`,
        length,
        effective_length,
        expected_count,
        TPM
    )

# ------------------------------------------------------------
# Round to RSEM's output precision
#
# RSEM writes the length columns to 2 decimals. expected_count is
# rounded to whole reads so the column is usable directly by
# count-based models, which is the point of the lengthScaledTPM
# default above.
# filter for only genes in the gtf file
# ------------------------------------------------------------

cat(
  "\nFiltering for genes in: ",
  genelist_file,
  "\n",
  sep = ""
)

used_genes <- as.character(read.table(genelist_file, header = F, sep = "\t")[,1])

n_missing <- length(gene_results$gene_id[!gene_results$gene_id %in% used_genes])

cat(
  "\n",
  n_missing,
  " genes removed, not in provided gene list\n",
  sep = ""
)

gene_results <- gene_results %>%
    mutate(
        length           = round(length, 2),
        effective_length = round(effective_length, 2),
        expected_count   = round(expected_count)
    ) %>% 
  filter(gene_id %in% used_genes)



# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------

cat(
    "\nWriting output:\n  ",
    output_file,
    "\n",
    sep = ""
)

write_tsv(
    gene_results,
    output_file
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

cat(
    "\n========================================\n",
    "Finished successfully!\n",
    "========================================\n",
    "Sample ID: ",
    sample_id,
    "\n",
    "Genes written: ",
    nrow(gene_results),
    "\n",
    "Output: ",
    output_file,
    "\n",
    sep = ""
)
