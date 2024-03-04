library(DESeq2)
library(ggplot2)
library(pheatmap)
library(Rsubread)
library(dplyr)
library(apeglm)

fcdata <- read.table("/home/majnusova/all/projects/jotnarlogs/results/aurantiochytrium/aurantiochytrium_combined_featurecounts.txt", sep="\t", header=T, skip = 1, row.names="Geneid", comment.char = "")
head(fcdata)

# extracing only those columns that contain the number of reads mapped
fcdata <- fcdata[ ,6:ncol(fcdata)]
head(fcdata)

# ediding column names (absolute paths -> names of the bam files)
colnames(fcdata) <- gsub("^.+\\.(\\w+)\\.bam$", "\\1", colnames(fcdata))
head(fcdata)

fcdata <- as.matrix(fcdata)
head(fcdata)

group <- factor(c(rep("motile", 3), rep("nonmotile", 3)))
condition <- factor(c(rep("poor_medium", 3), rep("rich_medium", 3)))

coldata <- data.frame(row.names=colnames(fcdata), group, condition) 
coldata

dds <-  DESeqDataSetFromMatrix(countData = fcdata,
                              colData = coldata,
                              design = ~ group)

dds$group <- relevel(dds$group, ref="nonmotile")

dds <- DESeq(dds)
dds

# list the coefficients
resultsNames(dds)

res <- results(dds, name="group_motile_vs_nonmotile")

res <- lfcShrink(dds, coef="group_motile_vs_nonmotile", type="apeglm")

rld = rlog(dds, blind=T)

resOrdered <- res[order(res$padj),]
resOrderedDF <- as.data.frame(resOrdered) 
write.csv(resOrderedDF, file ="aurantio_resordered.csv") 

counts_ = counts(dds, normalized=T) 
write.csv(counts_, file="aurantio_norm_counts_.csv")



resSig <- subset(res, padj < 0.05) 
resSig <- subset(resSig, abs(log2FoldChange) > 0.6) 
select <- rownames(resSig[ order(resSig$padj), ], ) 

mat = assay(rld)[select,] 
mat = t(scale(t(mat))) # z-score
nrow(mat)
head(mat)

write.csv(mat, file="aurantio_significant_genes_heatmap.csv")

# Heatmap: genes of our interest
anno <- as.data.frame(colData(rld)[c("group", "condition")])

p <- pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE, 
              cluster_cols = TRUE, annotation_col = anno, show_colnames = TRUE,
              labels_row = ifelse(rownames(mat) == "fgenesh1_pg.22_#_86", "RAW",
                         ifelse(rownames(mat) == "e_gw1.9.395.1", "RAQ",
                         ifelse(rownames(mat) == "estExt_Genemark1.C_2_t10283", "GOR3P", "")))) 
pdf("aurantio_significant_genes_heatmap.pdf", height = 300)
p
dev.off()

write.csv(mat, file="aurantio_allgenes_padj<1_nologfoldchangefilter_heatmap.csv")

# all gene names
anno <- as.data.frame(colData(rld)[c("group", "condition")])

p <- pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE, 
              cluster_cols = TRUE, annotation_col = anno, 
              show_colnames = TRUE, labels_row = rownames(mat))
pdf("aurantio_allgenes_padj<1_nologfoldchangefilter_heatmap.pdf", height = 1500)
p
dev.off()

# MA plot
resLFC <- lfcShrink(dds, coef="group_motile_vs_nonmotile", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-4,4), alpha=0.05)

genes_to_highlight <- c("fgenesh1_pg.22_#_86", "e_gw1.9.395.1", "estExt_Genemark1.C_2_t10283")
labels_for_genes <- c("RAW", "RAQ", "GOR3P")

indices <- match(genes_to_highlight, rownames(resLFC))

# Plot points for the highlighted genes
points(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], col="red", pch=20, cex=2)

# Add labels for the highlighted genes
text(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], labels=labels_for_genes, pos=4, cex=0.8)
library(dplyr)
library(ggrepel)install.packages("ggrepel")
