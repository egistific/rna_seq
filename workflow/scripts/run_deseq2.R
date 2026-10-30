suppressPackageStartupMessages({
    library(optparse)
    library(tximport)
    library(readr)
    library(DESeq2)
    library(ggplot2)
})

option_list <- list(
    make_option("--quantdir", type="character"),
    make_option("--samples", type="character"),
    make_option("--tx2gene", type="character"),
    make_option("--formula", type="character"),
    make_option("--condition-col", type="character", dest="condition_col"),
    make_option("--reference-level", type="character", dest="reference_level"),
    make_option("--fdr", type="double", default=0.05),
    make_option("--min-count", type="integer", dest="min_count", default=10),
    make_option("--lfc-threshold", type="double", dest="lfc_threshold", default=1),
    make_option("--known-de-genes", type="character", dest="known_de_genes", default=""),
    make_option("--results", type="character"),
    make_option("--ma-plot", type="character", dest="ma_plot"),
    make_option("--pca-plot", type="character", dest="pca_plot"),
    make_option("--significant-genes", type="character", dest="significant_genes"),
    make_option("--validation", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

# read samples file to give dataframe
samples_df <- as.data.frame(read_tsv(opt$samples, show_col_type = FALSE))
print(opt$samples)
print(samples_df)
# assign rownames to rows of dataframe
rownames(samples_df) <-samples_df$sample_id
# build one filepath per sample and save in a vector
quant_filepaths <- file.path(opt$quantdir, samples_df$sample_id, 'quant.sf')
# create a named character vector of filepaths
names(quant_filepaths) <-samples_df$sample_id
# check where there are any missing salmon output files
# missing is a named character vector containing filepaths of the missing files
missing <- quant_filepaths[!file.exists(quant_filepaths)]
if (length(missing) > 0) {
    stop("Missing salmon output files for: ", paste(names(missing), collapse=", "))
}

# create transcript and gene mapping file as a lookup table in dataframe format
tx2gene_df <- as.data.frame(read_tsv(opt$tx2gene, col_names=c("tx_name", "gene_id"), show_col_type=FALSE))

# read all salmon files and collapse the transcript level values to gene level
txi <- tximport(quant_filepaths, type="salmon", tx2gene=tx2gene_df, dropInfReps=TRUE)
# set donor and condition as factors
for (v in all.vars(as.formula(opt$formula))) {
    samples_df[[v]] <- factor(samples_df[[v]])
    }
# assign a certain condition to be the baseline
samples_df[[opt$condition_col]] <- relevel(samples_df[[opt$condition_col]], ref=opt$reference_level)

dds <- DESeqDataSetFromTximport(txi, colData=samples_df, design=as.formula(opt$formula))
# pre-filtering: remove low count genes
dds <- dds[rowSums(counts(dds)) > opt$min_count, ]
dds <- DESeq(dds)

res <- results(dds, alpha=opt$fdr)
summary(res)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
# reorder columns so "gene" comes first
res_df <- res_df[, c("gene", setdiff(colnames(res_df), "gene"))]
res_df <- res_df[order(res_df$padj), ]
write.table(res_df, opt$results, sep="\t", quote=FALSE, row.names=FALSE)

png(opt$ma_plot, width=900, height=700, res=120)
# x-axis: average gene across all samples and conditions; y-axis: log2 fold change
plotMA(res, alpha=opt$fdr, main="MA plot")
dev.off()

# let data have constant variance and normalize for library size
# vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE)
# data labelled by condition
pca_data <- plotPCA(rld, intgroup=opt$condition_col, returnData=TRUE)
# get the fraction of total variance explained by each component
percent_var <- round(100 * attr(pca_data, "percentVar"))
p <- ggplot(pca_data, aes(x=PC1, y=PC2, color=.data[[opt$condition_col]])) + geom_point(size=4) + 
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  # ggtitle("PCA of variance-stabilized counts") +
  ggtitle("PCA of rlog-transformed normalized counts") +
  theme_bw()

ggsave(opt$pca_plot, p, width=7, height=5, dpi=120)

lines <- c(
    sprintf("Total genes tested (non-NA padj): %d", sum(!is.na(res_df$padj))),
    sprintf("Significant genes (padj < %s): %d", opt$fdr, sum(res_df$padj < opt$fdr, na.rm=TRUE)), 
    ""
)

sig_genes_df <- res_df[!is.na(res_df$padj) & res_df$padj < opt$fdr & abs(res_df$log2FoldChange) > opt$lfc_threshold, ]
write.table(sig_genes_df, opt$significant_genes, sep="\t", quote=FALSE, row.names=FALSE)

known_genes <- trimws(strsplit(opt$known_de_genes, ",")[[1]])
known_genes <- known_genes[known_genes != ""]

for (g in known_genes) {
    if (g %in% res_df$gene) {
        row <- res_df[res_df$gene==g, ]
        status <- if (!is.na(row$padj) && row$padj < opt$fdr) "PASS" else "FAIL"
        lines <- c(lines, sprintf("%s: padj=%.3e log2F=%.3f -- %s", g, row$padj, row$log2FoldChange, status))
    } else {
        lines <- c(lines, sprintf("%s: NOT FOUND in results; check gene ID / tx2gene mapping", g))
    }
}

writeLines(lines, opt$validation)
writeLines(lines)
print("deseq2 done!")