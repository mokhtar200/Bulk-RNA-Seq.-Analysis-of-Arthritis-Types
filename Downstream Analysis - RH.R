# Load necessary libraries
library(DESeq2)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)
library(AnnotationDbi)
library(biomaRt)
library(ReactomePA)
library(EnhancedVolcano)
library(pheatmap)

setwd("D:/rheumatoid/")

# Load count matrix and metadata
rh_count = read.table("D:/rheumatoid/GSE89408_GEO_count_matrix_rename.txt.gz", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
metadata = read.csv("D:/rheumatoid/metadata_1.csv", row.names = 1)

# Match colnames and metadata
colnames(rh_count) <- rownames(metadata)
metadata$disease <- gsub("[^a-zA-Z0-9._]", "_", metadata$disease)

# Filter low-expressed genes
keep <- rowSums(rh_count >= 10) >= 2
rh_count <- rh_count[keep, ]
rh_count <- round(rh_count)
metadata$disease <- factor(metadata$disease)

# DESeq2 pipeline
dds <- DESeqDataSetFromMatrix(countData = rh_count, colData = metadata, design = ~ disease)
vsd <- vst(dds, blind = TRUE)

# PCA plot
pcaData <- plotPCA(vsd, intgroup = "disease", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = disease)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of RNA-Seq Samples (VST Transformed)")

# DESeq2 analysis
dds <- DESeq(dds)
res <- results(dds)
res <- lfcShrink(dds, coef=2, type="apeglm")
res_sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >= 1), ]
res_sig <- na.omit(res_sig)
write.csv(as.data.frame(res), file = "DEGs_deseq2.csv")

# MA plot
plotMA(res, ylim = c(-5, 5))

# Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                title = 'Disease_Compare')

# Heatmap of top 50 DEGs
top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), 50)
write.csv(top_genes, file = "top50_genes.csv", row.names = TRUE)
mat <- assay(vsd)[top_genes, ]
mat <- t(scale(t(mat)))
annotation_col <- as.data.frame(metadata["disease"])
pheatmap(mat,
         annotation_col = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         fontsize_row = 6,
         main = "Top 50 Differentially Expressed Genes")

# Functional enrichment (GO + KEGG)
# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(res_sig)
entrez_ids <- mapIds(org.Hs.eg.db,
                     keys = gene_symbols,
                     column = "ENTREZID",
                     keytype = "SYMBOL",
                     multiVals = "first")
entrez_ids <- na.omit(entrez_ids)

# GO enrichment
ego <- enrichGO(gene         = entrez_ids,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.2,
                readable     = TRUE)

# KEGG enrichment
ekegg <- enrichKEGG(gene         = entrez_ids,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05)

# Visualize enrichment
barplot(ego, showCategory = 10, title = "GO Biological Process")
dotplot(ekegg, showCategory = 10, title = "KEGG Pathways")
emapplot(pairwise_termsim(ego), showCategory = 10)
cnetplot(ego, showCategory = 5, foldChange = res_sig$log2FoldChange[names(entrez_ids)])

# Save enrichment results
write.csv(as.data.frame(ego), file = "GO_Enrichment.csv")
write.csv(as.data.frame(ekegg), file = "KEGG_Enrichment.csv")
