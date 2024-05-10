require("tidyverse")
require("itertools")
require("DESeq2")
require("sleuth")
require("tximport")
require("biomaRt")
require("hash")
require("plotly")
require("ggcorrplot")
require("FactoMineR")
require("factoextra")
require("pheatmap")
require("reticulate")
require("wasabi")
require("GenomicFeatures")
require("apeglm")
require("vsn")
require("gprofiler2")

# ---- Inputs ----
project_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Biodonostia/GB/CGGA_TPMs"
gtf_file <- 'gencode.v41.primary_assembly.annotation.gtf'
design_tab <- "design/design_CGGA.tab"
tpms_path <- "TPMs/CGGA_TPM_COUNTS_AGGREGATED-lengthScaledTPM.tab.gz"

# Select groups for comparison
group_a <- "Astrocytoma_G2"
group_b <- "Glioblastoma"

# Specify minimum number of reads across counts data
nreads <- 10
# Specify minimum log2FC for the regulated genes
fold_change <- 2

# ---- Prepare gene annotation ----
setwd(project_dir)

human_gtf <- rtracklayer::import(
  gtf_file)
gene_names <-  human_gtf %>% as.data.frame(.) %>%
  dplyr::select(c(gene_id, gene_name)) %>% 
  dplyr::rename(
    GENEID = gene_id,
    GENENAME = gene_name) %>%
  dplyr::mutate(GENEID = str_remove_all(
    GENEID, ".[0-9]*$")) %>%
  dplyr::distinct(GENEID, GENENAME)
rm(human_gtf)


# ---- Set up DGE analysis ----

# Prepare design
design <- read_delim(
    design_tab,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE
)
unique(design$GroupName)

selected_samples <- design %>%
  dplyr::filter(grepl(str_interp("${group_a}|${group_b}"), GroupName)) %>%
  dplyr::filter(GroupName != "EXCLUDE")
samples <- selected_samples$SampleName

selected_design <- design %>%
  dplyr::filter(grepl(str_interp("${group_a}|${group_b}"), GroupName))
selected_design$GroupName <- factor(selected_design$GroupName, levels = c(group_a, group_b))
selected_design$GroupName  <- relevel(selected_design$GroupName, ref=group_a)


# Load counts data
counts <- read_delim(
    tpms_path, 
    delim = "\t", 
    escape_double = FALSE,
    trim_ws = TRUE
    ) %>%
    dplyr::select(-GENENAME) %>%
    column_to_rownames("GENEID")

# Select samples for comparison
counts_matrix <- counts[, samples] %>% as.matrix() %>%
  round(.)

# Create DESeq object
model <- as.formula("~GroupName")
deseq_object <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = selected_design,
  design = model
)

# Filter counts by minimum number of reads
keep <- rowSums(counts(deseq_object)) >= nreads
deseq_object <- deseq_object[keep,]
deseq_object <- DESeq(deseq_object)

# Select contrast
contrast <- resultsNames(deseq_object)[2]
print(contrast)

contrast_dir <- contrast %>% str_replace("GroupName_", "lengthScaledTPM_")
print(contrast_dir)
if (!dir.exists(contrast_dir)){
  dir.create(contrast_dir)
} else {
  print("Directory already exists!")
}
setwd(contrast_dir)

# Perform GDE between two groups
objectLFC <- lfcShrink(deseq_object, 
                       coef=contrast, 
                       type="apeglm")

vsd <- vst(deseq_object, blind=FALSE)
ntd <- normTransform(deseq_object)

# Extract normalised counts
normalised_counts <- round(counts(deseq_object, normalized=TRUE)) %>%
  as.data.frame(.) %>%
  rownames_to_column("GENEID")

# Visualize Dist matrix for samples
counts_sampleDists <- dist(t(assay(ntd)))
counts_sampleDistMatrix <- as.matrix(counts_sampleDists)
dev.off()
pdf(file.path(
  str_interp("SAMPLES_NORM_COUNTS-${contrast_dir}.pdf"))
)
pheatmap(counts_sampleDistMatrix,
         clustering_distance_rows=counts_sampleDists,
         clustering_distance_cols=counts_sampleDists,
)
dev.off()

# Visualize variance stabilizing transformation for samples
vst_sampleDists <- dist(t(assay(vsd)))
vst_sampleDistMatrix <- as.matrix(vst_sampleDists)
pdf(file.path(
  str_interp("SAMPLES_VST-${contrast_dir}.pdf"))
)
pheatmap(vst_sampleDistMatrix,
         clustering_distance_rows=vst_sampleDists,
         clustering_distance_cols=vst_sampleDists,
)
dev.off()

# Visualize the Cooks distance for sample
pdf(file.path(
  str_interp("COOKS_DIST-${contrast_dir}.pdf")
), paper = "a4r"
)
boxplot(log10(assays(deseq_object)[["cooks"]]), range=0, las=2)
dev.off()

# dispersion estimates
pdf(file.path(
  str_interp("DISP_ESTS-${contrast_dir}.pdf"))
)
plotDispEsts(deseq_object)
dev.off()

# Get GDE data
gde_results <- objectLFC %>%
  as.data.frame(.) %>% 
  rownames_to_column("GENEID") %>%
  dplyr::left_join(., gene_names, by="GENEID") %>%
  dplyr::mutate(
    Significant = ifelse(
      padj <= 0.05, "Yes", "No"),
    Regulated = ifelse(
      log2FoldChange >= fold_change | 
        log2FoldChange <= -fold_change, "Yes", "No"),
    DiffRegulated = ifelse(
      Significant == "Yes" & Regulated == "Yes", "Yes", "No")
  )

# Get GDE data for significantly regulated genes
sig_gde_results <- gde_results %>% dplyr::filter(Significant == "Yes")

# Get GDE data for significantly regulated genes passing log2FC threshold
gde_regulated <- gde_results %>% dplyr::filter(DiffRegulated == "Yes")

# Volcano plot
de <- gde_results
min_fc <- round(min(gde_results$log2FoldChange) - .75)
max_fc <- round(max(gde_results$log2FoldChange) + .75)
max_pvalue <- max(-log10(gde_results$pvalue)[is.finite(-log10(gde_results$pvalue))])

de$diffexpressed <- "NO"
de$diffexpressed[de$log2FoldChange > fold_change & de$padj < 0.05] <- "UP"
de$diffexpressed[de$log2FoldChange < -fold_change & de$padj < 0.05] <- "DOWN"

volcano_plot <- ggplot(data=de, aes(x=log2FoldChange, 
                                    y=-log10(pvalue), col=diffexpressed)) +
  geom_point(alpha=0.5, size=3) + 
  theme_minimal() +
  scale_color_manual(values=c("red", "darkgrey", "blue")) +
  geom_vline(xintercept=c(-fold_change, fold_change), col="lightgrey") +
  geom_hline(yintercept=-log10(0.05), col="lightgrey") +
  ylim(c(0, max_pvalue)) +
  xlim(c(min_fc, max_fc)) +
  theme(legend.position="top") +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major = element_line(colour = "grey", linewidth = 0.1),
    panel.ontop = TRUE,
    plot.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    axis.ticks.length = unit(0.3, "cm"),
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text.x = element_text(size = 15, face = "bold"),
    panel.spacing = unit(1.25, "lines"),
    legend.text = element_text(size=14),
    legend.title = element_text(size=14),
  )

volcano_plot

ggsave(filename = file.path(str_interp("VOLCANO_PLOT_FC-${contrast_dir}-${fold_change}.pdf")),
       device = "pdf",
       plot = volcano_plot, width = 250, height = 250,  units = "mm")

# Save data
write_delim(
  normalised_counts, 
  file.path(str_interp("NORM_COUNTS-${contrast_dir}.tsv")),
  delim = "\t"
)

write_delim(
  gde_results, 
  file.path(
    str_interp("GDE_RESULTS-${contrast_dir}_NREADS-${nreads}.tsv")),
  delim = "\t"
)

write_delim(
  sig_gde_results, 
  file.path(
    str_interp("GDE_SIG_RESULTS-${contrast_dir}_NREADS-${nreads}.tsv")),
  delim = "\t"
)

write_delim(
  gde_regulated, 
  file.path(
    str_interp("GDE_DIFF_RESULTS-${contrast_dir}_NREADS-${nreads}_FC-${fold_change}.tsv")),
  delim = "\t"
)

write_delim(
  selected_design, 
  file.path(
    str_interp("DESIGN-${contrast_dir}.tsv")),
  delim = "\t"
)


# ---- Gene Ontology ----
background <- gde_results$GENEID

genes_up <- gde_regulated %>% dplyr::filter(log2FoldChange > 0)
genes_down <- gde_regulated %>% dplyr::filter(log2FoldChange < 0)
ids_up <- genes_up$GENEID
ids_down <- genes_down$GENEID

gostres_up <- gost(query = ids_up, 
                organism = "hsapiens",
                multi_query = TRUE, 
                significant = TRUE, 
                user_threshold = 0.05, 
                correction_method = "fdr", 
                domain_scope = "custom", 
                custom_bg = background, 
                numeric_ns = "", 
                highlight = TRUE,
                sources = c("GO:BP", "GO:MF", "GO:CC")
                )

gostres_down <- gost(query = ids_down, 
                   organism = "hsapiens",
                   multi_query = TRUE, 
                   significant = TRUE, 
                   user_threshold = 0.05, 
                   correction_method = "fdr", 
                   domain_scope = "custom", 
                   custom_bg = background, 
                   numeric_ns = "", 
                   highlight = TRUE,
                   sources = c("GO:BP", "GO:MF", "GO:CC")
                   )

results_up <- gostres_up$result %>% as.data.frame() %>%
  dplyr::select(-parents)
results_up <- apply(results_up,2,as.character)

results_down <- gostres_down$result %>% as.data.frame() %>%
  dplyr::select(-parents)
results_down <- apply(results_down,2,as.character)

write_delim(
  results_up %>% as.data.frame(),
  file.path(
    str_interp("GOTERMS-${contrast_dir}-Up.tsv")),
  delim = "\t"
)

write_delim(
  results_down %>% as.data.frame(), 
  file.path(
    str_interp("GOTERMS-${contrast_dir}-Down.tsv")),
  delim = "\t"
)

