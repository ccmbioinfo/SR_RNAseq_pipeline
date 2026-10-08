suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(OUTRIDER))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(optparse))

option_list = list(
  make_option(c("-s", "--sample_salmon"), type="character", help="sample salmon output"),
  make_option(c("-c", "--cohort_salmon"), type="character", help="cohort salmon output"),
  make_option(c("--txdb"), type="character", help="path to txdb"),
  make_option(c("--tx2gene"), type="character", help="path to tx2gene file"),
  make_option(c("--q_value"), type="integer", help="q for OUTRIDER"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output path;[default: <sample>_OUTRIDER_results.rds]")
)

opt_parser<-OptionParser(option_list=option_list)
opt<-parse_args(opt_parser)

cohort_quant_files_path <- opt$cohort_salmon
sample_quant_files_path <- opt$sample_salmon

quant_files <- scan(
    cohort_quant_files_path,
    what = "character"
)

quant_files <- c(quant_files, sample_quant_files_path)

# tx2gene file
tx2gene_file <- opt$tx2gene

tx2gene <- read.delim(
    tx2gene_file,
    header = TRUE,
    stringsAsFactors = FALSE
)

txi <- tximport(
    files = quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = T
)

all_count <- txi$counts

#remove all zero rows
all_count <- all_count[
    rowSums(all_count) > 0,
]

# Remove version from gene IDs
rownames(all_count) <- sub("\\..*$", "", rownames(all_count))

#remove globin
keep_gene <- read.csv('/hpf/largeprojects/tgnode/sandbox/yliang_analysis/pacbio_longread/annotation_files/gene_id_name_biotype_gtf115_new.tsv',
                      sep='\t',header =F)
keep_gene <- keep_gene[,1]

genes_to_remove <- c('ENSG00000206172','ENSG00000188536','ENSG00000244734') #remove globin
all_count <- all_count[
  !rownames(all_count) %in% genes_to_remove,
]

all_count <- all_count[!grepl("ERCC-|SIRV", row.names(all_count)),,drop = FALSE]

all_count <- all_count[
  rownames(all_count) %in% keep_gene,
]

all_count <- na.omit(all_count)

dim(all_count)

colnames(all_count) <- basename(dirname(dirname(quant_files)))
#colnames(all_count) <- dirname(gsub('/hpf/largeprojects/tgnode/sandbox/rnaseq_repo/TG_RNAseq/','',dirname(quant_files)))

ods <- OutriderDataSet(countData = round(all_count))
txdb <- loadDb(opt$txdb)
implementation <- 'autoencoder'

fpkmCutoff <- 1
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE,
                        fpkmCutoff=fpkmCutoff, addExpressedGenes=TRUE)

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

## subset filtered
ods <- ods[mcols(ods)$passedFilter,]

register(MulticoreParam(10))
ods <- estimateSizeFactors(ods)
ods <- controlForConfounders(ods, q=opt$q_value, implementation=implementation)
ods <- fit(ods)
ods <- computePvalues(ods)
ods <- computeZscores(ods)
all_res <- results(ods, all = TRUE)
saveRDS(all_res,opt$output)
