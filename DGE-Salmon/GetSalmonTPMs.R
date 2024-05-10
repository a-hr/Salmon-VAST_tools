require("tidyverse")
require("tximport")
require("rtracklayer")
require("BiocParallel")

paramMulti <- MulticoreParam(workers = 10)
paramSerial <- SerialParam()
register(paramSerial)


args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Provide a metric to extract TPM counts", call.=FALSE)
}

project_dir <- args[2]
salmon_out <- args[3]
design_tab <- args[4]
gtf <- args[5]
outname <- args[6]

setwd(project_dir)

# Prepare annotation ------------------------------------------------------

# Get transcripts
human_gtf <- rtracklayer::import(gtf)
transcripts <- human_gtf %>% as.data.frame(.) %>%
  dplyr::select(c(gene_id, transcript_id)) %>%
  dplyr::rename(
    GENEID = gene_id, 
    TXNAME = transcript_id) %>%
  dplyr::mutate(GENEID = str_remove_all(
    GENEID, ".[0-9]*$")) %>%
  dplyr::distinct(TXNAME, GENEID) %>%
  dplyr::filter(!is.na(TXNAME))

# Get gene names
gene_names <-  human_gtf %>% as.data.frame(.) %>%
  dplyr::select(c(gene_id, gene_name)) %>% 
  dplyr::rename(
    GENEID = gene_id,
    GENENAME = gene_name) %>%
  dplyr::mutate(GENEID = str_remove_all(
    GENEID, ".[0-9]*$")) %>%
  dplyr::distinct(GENEID, GENENAME)

# Final annotation
# human_annotation <- transcripts %>% 
#   dplyr::left_join(., gene_names, by="GENEID")

rm(human_gtf)

# Prepare data ------------------------------------------------------------

design <- read_delim(
    design_tab, 
    delim = "\t",
    escape_double = FALSE, 
    trim_ws = TRUE
)

files <- file.path(project_dir, salmon_out, design$SampleName , "quant.sf")
names(files) <- paste0(design$SampleName)

# Prepare counts ----------------------------------------------------------

# Raw TPM counts
txi.inf <- tximport(files, type = "salmon", tx2gene = transcripts, dropInfReps=TRUE, countsFromAbundance=args[1])
counts_df <- txi.inf$counts %>%
  as.data.frame() %>%
  rownames_to_column("GENEID") %>%
  dplyr::left_join(., gene_names, by="GENEID")

write.table(
  counts_df, 
  str_interp("${outname}_TPM_COUNTS_AGGREGATED-${args[1]}.tab"),
  sep = "\t", quote = FALSE, 
  row.names = FALSE
)
rm(txi.inf)
