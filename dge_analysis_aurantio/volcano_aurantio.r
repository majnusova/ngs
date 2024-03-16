library(ggrepel)

library(tidyverse)

deseq_res <- read_csv("/home/majnusova/all/projects/jotnarlogs/scripts/aurantio_resordered.csv") #do prvniho sloupce dopsat gene
head(deseq_res)

# Parse gene that are differently up or down regulated 
deseq_res$diffexpressed <- "NO"

deseq_res$diffexpressed[deseq_res$log2FoldChange > 0.6 & deseq_res$padj < 0.05 ] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.06, set as "DOWN"
deseq_res$diffexpressed[deseq_res$log2FoldChange < -0.6 & deseq_res$padj < 0.05 ] <- "DOWN"


deseq_res$delabel <- NA

diff_ord <- deseq_res %>%
        filter(diffexpressed != 'NO' & abs(log2FoldChange)>1) %>%
        arrange(padj)
diff_genes_na <- rbind(head(diff_ord %>% filter(log2FoldChange > 0), 3), head(diff_ord %>% filter(log2FoldChange < 0), 3))
deseq_res$delabel[deseq_res$gene %in% c(diff_genes_na$gene)] <- deseq_res$gene[deseq_res$gene %in% c(diff_genes_na$gene)]

gene_nickname_map <- setNames(
  c("RAW", "RAQ", "GOR3P", "Rab23", "Rab28", "Rabl2", "RJL", "IFT27", "Arl3", "Arl6", "IFT22", "D2LIC"),
  c("fgenesh1_pg.22_#_86", "e_gw1.9.395.1", "estExt_Genemark1.C_2_t10283", "fgenesh1_pg.9_#_334", 
    "fgenesh1_pm.7_#_316", "e_gw1.16.392.1", "e_gw1.10.693.1","fgenesh1_pg.11_#_69", "e_gw1.5.671.1", 
    "fgenesh1_pg.4_#_88", "gw1.1.1510.1", "fgenesh1_pg.29_#_26" )
)

# Přiřazení přezdívek genům
deseq_res$delabel <- ifelse(deseq_res$gene %in% names(gene_nickname_map), gene_nickname_map[deseq_res$gene], NA)


# prvních 10 hodnot sloupce delabel
#print(head(deseq_res$delabel, 10))

# deseq_res %>% write.csv(paste0( "unpaired_batch/volcano_lable.csv"))
#plot annotated volcano plot

pdf(paste0("volcano_correctmodels.pdf"))
#cols <- c("DOWN" = "#0072ba", "NO" = "#bfbfbf", "UP" = "#aa27a9")
cols <- c("DOWN" = "#23619b", "NO" = "#bfbfbf", "UP" = "#ae242e")
#p <- ggplot(data=deseq_res, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=delabel )) +
p <- ggplot(data=deseq_res, aes(x=log2FoldChange, y=-log10(padj), colour =factor(diffexpressed), fill = factor(diffexpressed), label=delabel )) +
        labs(fill="Differentially Expressed", title = "Motile vs. Non-motile Stage")+
        geom_point(shape = 21, alpha = 0.4, size = 3) + #alpha je pruhlednost bodu
        geom_label_repel(colour = "black", fill = "#e5fbe5", box.padding = 0.4, max.overlaps = 200, force=20, segment.size=0.3) +
        theme_bw() +
    theme(plot.title=element_text(face="bold", hjust=0.5)) +

        scale_color_manual(
                values = cols,
                aesthetics = c("colour", "fill"),
#                values=c("blue", "grey", "red"),
                ) +
        xlab("Log2FoldChange") +
        ylab("-log10(padj)") + 
        geom_vline(xintercept = -0.6, colour="grey50", linetype = "longdash") +
        geom_vline(xintercept = 0.6, colour="grey50", linetype = "longdash") +
        geom_hline(yintercept = -log10(0.05), colour="grey50", linetype = "longdash") + 
        theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.ticks = element_line(colour = "grey70", size = 0.2),
                legend.position = "none",
            
        )
p
dev.off()

Rab23:
(base) majnusova@Swift:~/all/projects/jotnarlogs/scripts$ grep "fgenesh1_pm.7_#_316" aurantio_resordered.csv 
"fgenesh1_pm.7_#_316",213.392346254487,4.63967312470034,0.32215930281018,6.84814508020044e-48,5.027472167284e-47
IFT22:
grep "fgenesh1_pg.11_#_69" aurantio_resordered.csv 
"fgenesh1_pg.11_#_69",720.332351300818,4.2668639286455,0.253475215376877,1.18075560883926e-64,1.4554184669023e-63
RJL:
grep "e_gw1.10.693.1" aurantio_resordered.csv 
"e_gw1.10.693.1",206.428212338803,4.8151729962228,0.366451143649792,2.29320950427372e-40,1.30963475796302e-39
