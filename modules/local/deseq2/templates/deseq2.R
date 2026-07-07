#!/usr/bin/env Rscript

# ---- Nextflow-injected variables ----
# When Nextflow runs a template, it substitutes \$variables with real values.
prefix        <- "${meta.id}"
counts_file   <- "${counts}"
samplesheet   <- "${samplesheet}"

suppressPackageStartupMessages({
    library(DESeq2)
})

# ---- 1. Read the count matrix ----
# Columns: gene_id, gene_name, then one column per sample.
raw <- read.delim(counts_file, header = TRUE, sep = "\t", check.names = FALSE)
rownames(raw) <- raw\$gene_id
# Drop the two annotation columns; keep only numeric sample counts.
counts <- raw[, !(colnames(raw) %in% c("gene_id", "gene_name")), drop = FALSE]
# DESeq2 requires INTEGER counts; tximport length-scaled values are decimals.
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# ---- 2. Read sample metadata (which sample is which condition) ----
ss <- read.csv(samplesheet, header = TRUE, stringsAsFactors = FALSE)
# Keep only samples present in the count matrix, in the same column order.
ss <- ss[match(colnames(counts), ss\$sample), ]
coldata <- data.frame(
    row.names = ss\$sample,
    condition = factor(ss\$condition)
)

# ---- 3. Build the DESeq2 dataset with design = ~ condition ----
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = coldata,
    design    = ~ condition
)

# ---- 4. Run the differential expression analysis ----
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res\$padj), ]  # sort by adjusted p-value (most significant first)

# ---- 5. Write the results table ----
res_df <- as.data.frame(res)
res_df\$gene_id <- rownames(res_df)
write.table(res_df, paste0(prefix, ".de_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 6. PCA plot (do samples cluster by condition?) ----
# Variance-stabilising transform makes counts comparable for PCA.
vsd <- tryCatch(vst(dds, blind = TRUE),
                error = function(e) varianceStabilizingTransformation(dds, blind = TRUE))
png(paste0(prefix, ".pca.png"), width = 800, height = 600)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv <- round(100 * attr(pcaData, "percentVar"))
plot(pcaData\$PC1, pcaData\$PC2,
     col = as.integer(pcaData\$condition) + 1, pch = 19, cex = 2,
     xlab = paste0("PC1: ", pv[1], "% variance"),
     ylab = paste0("PC2: ", pv[2], "% variance"),
     main = "PCA of samples")
text(pcaData\$PC1, pcaData\$PC2, labels = rownames(pcaData), pos = 3, cex = 0.8)
legend("topright", legend = levels(pcaData\$condition),
       col = seq_along(levels(pcaData\$condition)) + 1, pch = 19)
dev.off()

# ---- 7. Volcano plot (fold-change vs significance) ----
png(paste0(prefix, ".volcano.png"), width = 800, height = 600)
with(res_df, {
    plot(log2FoldChange, -log10(pvalue),
         pch = 20, cex = 0.6, col = "grey50",
         xlab = "log2 fold change", ylab = "-log10(p-value)",
         main = "Volcano plot")
    sig <- !is.na(padj) & padj < 0.05
    points(log2FoldChange[sig], -log10(pvalue[sig]), pch = 20, cex = 0.7, col = "red")
})
abline(h = -log10(0.05), lty = 2, col = "blue")
dev.off()

# ---- 8. Assemble a self-contained HTML report ----
top <- head(res_df[!is.na(res_df\$padj), ], 20)
top_rows <- paste0(
    "<tr><td>", top\$gene_id, "</td><td>",
    round(top\$log2FoldChange, 3), "</td><td>",
    signif(top\$padj, 3), "</td></tr>", collapse = "\n"
)
n_sig <- sum(!is.na(res_df\$padj) & res_df\$padj < 0.05)

# Link the PNGs by filename (report is a folder: HTML + PNGs). No extra deps needed.
pca_img     <- paste0(prefix, ".pca.png")
volcano_img <- paste0(prefix, ".volcano.png")

html <- paste0(
"<!DOCTYPE html><html><head><meta charset='utf-8'><title>DE report: ", prefix, "</title>",
"<style>body{font-family:sans-serif;max-width:900px;margin:2em auto;padding:0 1em}",
"table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccc;padding:6px 10px;text-align:left}",
"th{background:#f4f4f4}img{max-width:100%;border:1px solid #eee;margin:1em 0}</style></head><body>",
"<h1>Differential expression report</h1>",
"<p><b>Dataset:</b> ", prefix, " &nbsp; <b>Genes tested:</b> ", nrow(res_df),
" &nbsp; <b>Significant (padj &lt; 0.05):</b> ", n_sig, "</p>",
"<h2>PCA</h2><img src='", pca_img, "'>",
"<h2>Volcano</h2><img src='", volcano_img, "'>",
"<h2>Top 20 genes by adjusted p-value</h2>",
"<table><tr><th>gene_id</th><th>log2 fold change</th><th>padj</th></tr>",
top_rows, "</table></body></html>"
)
writeLines(html, paste0(prefix, ".de_report.html"))

# ---- 9. Version reporting (for the versions topic) ----
writeLines(c(
    '"DESEQ2":',
    paste0("    bioconductor-deseq2: ", as.character(packageVersion("DESeq2")))
), "versions.yml")
