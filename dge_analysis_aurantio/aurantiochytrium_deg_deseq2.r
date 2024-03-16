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
group
condition <- factor(c(rep("poor_medium", 3), rep("rich_medium", 3)))
condition

# to use dds, one needs to have count matrix cts and a table of sample information called coldata
# the design indicates how to model the samples (to measure the effect of the condition), controlling for batch differences. 
# The two factor variables batch and condition should be columns of coldata.
coldata <- data.frame(row.names=colnames(fcdata), group, condition) #row.names=colnames(fcdata) assigns the column names of the fcdata matrix as row names for colData, linking the sample metadata to the corresponding count data for each sample
coldata

dds <-  DESeqDataSetFromMatrix(countData = fcdata,
                              colData = coldata,
                              design = ~ group)
dds

# jde zde pouzit?

# By default, R will choose a reference level for factors based on alphabetical order. 
# If you never tell the DESeq2 functions which level you want to compare against (e.g. which level represents the control group), the comparisons will be based on the alphabetical order of the level
dds$group <- relevel(dds$group, ref="nonmotile")
dds

dds <- DESeq(dds)
dds

# list the coefficients
resultsNames(dds)

res <- results(dds, name="group_motile_vs_nonmotile") # results of the differential expression test for the comparison between "motile" and "nonmotile" 
res

# shrinkage je dulezity pro snizeni logfoldu u genu s nizkymi counts: vyradi to tedy radu false positive genu
# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
# shrink log fold changes association with condition:
res <- lfcShrink(dds, coef="group_motile_vs_nonmotile", type="apeglm")
res

rld = rlog(dds, blind=T)
rld

resOrdered <- res[order(res$padj),] # ordering results by padj in ascending order = most statistically significant genes (those with the smallest adjusted p-values) will come first.
resOrderedDF <- as.data.frame(resOrdered) 
write.csv(resOrderedDF, file ="aurantio_resordered.csv")

counts_ = counts(dds, normalized=T)
write.csv(counts_, file="aurantio_norm_counts_.csv")

# resSig subsetuje pouze geny, kde je vyrazna zmena; v tomhle pripade je absolute log2Foldchange vetsi nez 1 = ta zmena musi byt alespon dvakrat vice nebo mene, log2(2) je totiz 1

resSig <- subset(res, padj < 0.1) # subset of res containing only the results that are statistically significant at an adjusted p-value threshold of 0.01
#resSig <- subset(resSig, abs(log2FoldChange) > 1) # the subset is further refined to include only those genes with a log2 fold change greater than 1 or less than -1, which typically signifies genes that are at least 2-fold upregulated or downregulated
select <- rownames(resSig[ order(resSig$padj), ], ) 

mat = assay(rld)[select,] 
mat = t(scale(t(mat)))
nrow(mat)
head(mat)

write.csv(mat, file="aurantio_heatmap_allgenes.csv")

# Create an annotation data frame
anno <- as.data.frame(colData(rld)[c("group", "condition")])

p <- pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE, 
              cluster_cols = TRUE, annotation_col = anno, show_colnames = TRUE,
              labels_row = ifelse(rownames(mat) == "fgenesh1_pg.22_#_86", "RAW",
                         ifelse(rownames(mat) == "gw1.9.395.1", "RAQ",
                         ifelse(rownames(mat) == "e_gw1.2.260.1", "GOR3P", "")))) 

pdf("aurantio_significant_genes_heatmap.pdf", height = 300)
p
dev.off()


write.csv(mat, file="aurantio_heatmap_allgenes.csv")

# Create an annotation data frame
anno <- as.data.frame(colData(rld)[c("group", "condition")])

# Generate the heatmap with labels for all genes
p <- pheatmap(mat, cluster_rows = TRUE, show_rownames = TRUE, 
              cluster_cols = TRUE, annotation_col = anno, 
              show_colnames = TRUE, labels_row = rownames(mat))

# Save the heatmap to a PDF file
pdf("aurantio_significant_genes_heatmap_allnames.pdf", height = 1200)
p
dev.off()


 plotMA(res, ylim=c(-4,4))

# Highlight specific genes of interest with their correct labels
genes_to_highlight <- c("fgenesh1_pg.22_#_86", "gw1.9.395.1", "e_gw1.2.260.1")
labels_for_genes <- c("RAW", "RAQ", "GOR3P")

# Find the indices of the genes to highlight in the results
indices <- match(genes_to_highlight, rownames(resLFC))

# Plot points for the highlighted genes
points(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], col="red", pch=20, cex=2)

# Add labels for the highlighted genes
text(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], labels=labels_for_genes, pos=4, cex=0.8)


##Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes
# shrink log fold changes association with condition:
resLFC <- lfcShrink(dds, coef="group_motile_vs_nonmotile", type="apeglm")
resLFC
plotMA(resLFC, ylim=c(-4,4))

genes_to_highlight <- c("fgenesh1_pg.22_#_86", "gw1.9.395.1", "e_gw1.2.260.1")
labels_for_genes <- c("RAW", "RAQ", "GOR3P")

indices <- match(genes_to_highlight, rownames(resLFC))

# Plot points for the highlighted genes
points(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], col="red", pch=20, cex=2)

# Add labels for the highlighted genes
text(resLFC$baseMean[indices], resLFC$log2FoldChange[indices], labels=labels_for_genes, pos=4, cex=0.8)


