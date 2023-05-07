
library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)
library(cowplot)
library(patchwork)
library(doParallel)
library(ggVennDiagram)
library(tidytext)
library(viridis)
library(MyEllipsefit)
library(princurve)

source('./util_funcs.R')

getPcaMetaData <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## ggplot 
# L.atac.sds <-  readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_sds_data_YR.rds')
# 
# p <- ggplot(L.atac.sds, aes(x=PC_1,y=PC_2)) +
#   geom_point(aes(fill = phase, color =  phase), shape=21, size = 1)+ 
#   scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
#   scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
#   geom_point(aes(x=sc1[order(pt.shifted.scaled)],y=sc2[order(pt.shifted.scaled)]), col = 'black',fill = "black", 
#              shape=21, size = 0.01, stroke = 1) +
#   # geom_path(aes(x=sc1[order(cell.ord)],y=sc2[order(cell.ord)]), col = 'black', size = 1.5, 
#   #           arrow = arrow(type = "open", angle = 30, ends = "first",length = unit(0.25, "inches")))+
#   # geom_point(aes(x=sc1[order(pt.shifted.scaled)][1],y=sc2[order(pt.shifted.scaled)][1]), col = 'black',fill = "black", 
#   #            shape=21, size = 7, stroke = 1.5)+
#   theme_bw(base_size = 14) +
#   theme(axis.text.y = element_text( size = 16, face="bold", color = "black"),
#         axis.text.x = element_text(size = 16, face="bold", color = "black"),
#         plot.title = element_text(size=18, face = "bold.italic", color = 'black', hjust = 0.5),
#         axis.title.x = element_text(size=18, face="bold"),
#         axis.title.y = element_text(size=18, face="bold"),
#         legend.title  = element_text(size = 18, face = "bold"),
#         legend.text = element_text(size =  "16", face = "bold", color = "black"),
#         strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
#   ylab('PC2') + xlab('PC1') +
#   ggtitle("scATAC") 
# 
# plot(p) + guides(color = guide_legend(override.aes = list(size = 5)))
# 
# 
# ggsave("../Output/toxo_cdc/ME49_59/figures_paper/pca_scATAC.pdf", 
#        width = 6, height = 6, dpi = 300)

## Pseudo time
sc.rna.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_genes_expr_pt.rds')
S.O.rna  <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
rna.sds.data <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_sds_data.rds')


L <- wiskerPlot(S.O.rna)
L$pc$PC_2 <- L$pc$PC_2 * -1
L$fit$s[,2] <- L$fit$s[,2]  * -1 

pdf(file="../Output/toxo_cdc/ME49_59/figures_paper/whiskers_intra_rna_no_curve.pdf",
    width=6, height=6)


par(mar = c(5, 5, 4, 4) + 0.1)
plot(x = -40:12, y = -26:26, type = 'n',xlab = '', ylab = '',
     lwd = 2, cex.lab = 2, cex.main = 1.5, cex.axis = 1.5)
grid(13,13, lwd = 1, col = "lightgray", lty = 1)
box(which = "plot", lty = "solid")
whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
color = rep(NA, length=length(rna.sds.data$phase))
color[which(rna.sds.data$phase=="G1.a")] = "#b6232a"
color[which(rna.sds.data$phase=="G1.b")] = "#ed7202"
color[which(rna.sds.data$phase=="S")] = "#caae05"
color[which(rna.sds.data$phase=="M")] = "#6f883a"
color[which(rna.sds.data$phase=="C")] = "#b138ee"
points(L$pc$PC_1, L$pc$PC_2, cex = 0.5, col = color, pch = 20)
points(rna.sds.data$sc1[rna.sds.data$cell.ord],-rna.sds.data$sc2[rna.sds.data$cell.ord], cex = 0.2, col = 'black')
          
# grid(nx = NULL, ny = NULL,
#      lty = 1, col = "gray", lwd = 1)

          
dev.off()


## ggplot 

S.O.rna@meta.data$spp2 <- S.O.rna@meta.data$spp
intra.pca <- getPcaMetaData(S.O.rna)

## intra - rna  - pca
p1  <- ggplot(intra.pca, aes(x= PC_1,y=-PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("intra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p1)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_rna_intra.pdf", 
       plot = p1, width = 6, height = 6, dpi = 300)


##
p1  <- ggplot(intra.pca, aes(x= -UMAP_1,y=-UMAP_2)) +
  geom_point(aes(fill = phase, color = phase), 
             #color = 'blue', 
             shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("intra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))

plot(p1)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/Umap_rna_intra.pdf", 
       width = 6, height = 6, dpi = 300)


## density 
p  <- ggplot(intra.pca, aes(x= PC_1,y= -PC_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("intra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


p



ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/pca_intra_density.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## density 
p  <- ggplot(intra.pca, aes(x= -UMAP_1,y= -UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("intra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


p



ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/umap_intra_density.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## atac 

## Pseudo time

sc.atac.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_genes_expr_pt.rds')
S.O.atac <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')
atac.sds.data <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_sds_data.rds')

pca = S.O.atac@reductions$pca
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
varExplained[1:3]

L <- wiskerPlot(S.O.atac)
L$pc$PC_2 <- L$pc$PC_2 * -1
L$fit$s[,2] <- L$fit$s[,2] * -1

pdf(file="../Output/toxo_cdc/ME49_59/figures_paper/whiskers_atac.pdf",
    width=6, height=6)


par(mar = c(5, 5, 4, 4) + 0.1)
plot(x = -35:10, y = -25:20, type = 'n',  xlab = '', ylab = '',  #xaxt = "n", yaxt = "n", axes=FALSE,  
     lwd = 2, cex.lab = 1.5, cex.main = 2, cex.axis = 1.5)
grid(10,10, lwd = 1, col = "lightgray", lty = 1)
box(which = "plot", lty = "solid")
whiskers(as.matrix(L$pc[,c(1,2)]), L$fit$s, col = "gray")
color = rep(NA, length=length(atac.sds.data$phase))
color[which(atac.sds.data$phase=="G1.a")] = "#b6232a"
color[which(atac.sds.data$phase=="G1.b")] = "#ed7202"
color[which(atac.sds.data$phase=="S")] = "#caae05"
color[which(atac.sds.data$phase=="M")] = "#6f883a"
color[which(atac.sds.data$phase=="C")] = "#b138ee"
points(atac.sds.data$PC_1, -atac.sds.data$PC_2, cex = 0.5, col = color, pch = 20)
points(atac.sds.data$sc1[atac.sds.data$cell.ord],-atac.sds.data$sc2[atac.sds.data$cell.ord], cex = 0.2, col = 'black')
          
dev.off()
          

## ggplot

S.O.atac@meta.data$spp2 <- "atac"
atac.pca <- getPcaMetaData(S.O.atac)


p1  <- ggplot(atac.pca, aes(x= PC_1,y=-PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("intra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p1)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_atac_intra.pdf", 
       plot = p1, width = 6, height = 6, dpi = 300)


p1  <- ggplot(atac.pca, aes(x= UMAP_1,y=UMAP_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("atac") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p1)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/umap_atac_intra.pdf", 
       plot = p1, width = 6, height = 6, dpi = 300)

## Density 

p  <- ggplot(atac.pca, aes(x= UMAP_1,y= UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("atac") + 
  guides(color = guide_legend(override.aes = list(size = 7)))

p



ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/umap_atac_density.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## extra
S.O.intra.extra <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_extra_lables_not_anchored_list.rds')
names(S.O.intra.extra) <- c("intra", "extra")
S.O.extra <- S.O.intra.extra$extra
S.O.extra@meta.data$spp2 <- S.O.extra@meta.data$spp
extra.pca <- getPcaMetaData(S.O.extra)
extra.pca$phase <- factor(extra.pca$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))

## extra - umap
p2  <- ggplot(extra.pca, aes(x= PC_1,y=-PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("extra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p2)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/pca_rna_extra_not_anchored.pdf", 
       plot = p2, width = 6, height = 6, dpi = 300)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_rna_extra_not_anchored.pdf", 
       width = 6, height = 6, dpi = 300)

# extra density 
p  <- ggplot(extra.pca, aes(x= UMAP_1,y= -UMAP_2)) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + ggtitle("extra") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


p



ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/umap_extra_density.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## box plot 
rna.sds.data <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_sds_data.rds')

p <- ggplot(rna.sds.data, aes(x=phase, y=pt.shifted.scaled, color  = phase)) + 
  geom_boxplot(size = 0.7) +
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  theme_bw() +
  xlab('phase') + ylab("pseudo-time") +
  theme(axis.text.y = element_text( size = 18, face="bold", color = "black"),
        axis.text.x = element_text(size = 18, face="bold", color = "black"),
        plot.title = element_text(size=20, face = "bold.italic", color = 'black', hjust = 0.5),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.title  = element_text(size = 20, face = "bold"),
        legend.text = element_text(size =  "18", face = "bold", color = "black"),
        strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
  ggtitle("scRNA")

p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/boxPlot_scRNA.pdf", 
       plot = p,  width = 5, height = 5, dpi = 300)


## atac box plot

atac.sds.data <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_sds_data.rds')
p <- ggplot(atac.sds.data, aes(x=phase, y=pt.shifted.scaled, color  = phase)) + 
  geom_boxplot(size = 0.7) +
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  theme_bw() +
  xlab('phase') + ylab("pseudo-time") +
  theme(axis.text.y = element_text( size = 18, face="bold", color = "black"),
        axis.text.x = element_text(size = 18, face="bold", color = "black"),
        plot.title = element_text(size=20, face = "bold.italic", color = 'black', hjust = 0.5),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        legend.title  = element_text(size = 20, face = "bold"),
        legend.text = element_text(size =  "18", face = "bold", color = "black"),
        strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
  ggtitle("scATAC")

p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/boxPlot_scATAC.pdf", 
       plot = p,  width = 5, height = 5, dpi = 300)


## cross correlation 
cc.dat <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_sc_atac_cross_cor_lag.rds')
p <- ggplot(cc.dat, aes(x = ccs)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1,
               linetype = 1,
               colour = 2)+
  theme_bw() + 
  theme(plot.title = element_text(face = "bold.italic", size = 18),
        axis.title = element_text(face = "bold", size = 14)) +
  theme(axis.text.y = element_text( size = 16, face="bold", color = "black"),
        axis.text.x = element_text(size = 16, face="bold", color = "black"),
        plot.title = element_text(size=18, face = "bold.italic", color = 'black', hjust = 0.5),
        axis.title.x = element_text(size=18, face="bold"),
        axis.title.y = element_text(size=18, face="bold"),
        legend.title  = element_text(size = 18, face = "bold"),
        legend.text = element_text(size =  "16", face = "bold", color = "black"),
        strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid"))+
  ggtitle("scRNA & scATAC cross-correlation")

p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/sc_rna_sc_atac_cross_corr_lag.png", 
       plot = p, height = 6, width = 6, dpi = 300)

#### 
Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')
colnames(Intra.markers.sig) <- gsub("cluster", "phase", colnames(Intra.markers.sig)) 
ss <- Intra.markers.sig %>% group_by(phase) %>% summarise(num.DEG = n())
ss$phase <- factor(ss$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))

p <- ggplot(data = ss, aes(x=phase, y=num.DEG, fill = phase, color = phase)) +
  geom_bar(stat="identity",size = 2, width = 0.75, aes(fill = phase))+
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  geom_text(aes(label=num.DEG), vjust=1.5, color="white", size=8, fontface = 'bold')+
  # facet_grid(. ~ spp, scales = "free", space='free',labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) +
  ggtitle("cell cycle markers")
#theme(legend.position = "none")
p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/bar_plot_cell_cycle_markers.pdf", 
       plot = p, height = 6, width = 6, dpi = 300, )


## transition 

L.trans.atac <- readRDS('../Input/toxo_cdc/rds/atac_based_transition_points_v2.rds')
atac_peaks.dat <- L.trans.atac$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
atac_peaks.dat$value[atac_peaks.dat$value < 0] <- 0
atac_peaks.dat$value[atac_peaks.dat$value > 6] <- 6


atac_peaks.dat$drivs <- factor(atac_peaks.dat$drivs, levels = c('y', 'yp'))
p1  <- ggplot(atac_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1.2)+ 
  scale_color_manual(values = c("y" = "blue3","yp" ='chartreuse4'))+
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.atac$transition.points$x, linetype=2, color = 'brown', size = 1) + 
  facet_grid(drivs~., scales = 'free')+
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", colour = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white", size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 18, face="bold", angle = 0))  +
  ggtitle("atac - transition") +
  theme(plot.title = element_text(size=18, face = "bold.italic", color = 'brown'),
    axis.title.x = element_text(size=20, face="bold", hjust = 1),
    axis.title.y = element_text(size=20, face="bold")) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=16, face="bold"),
    legend.text = element_text(colour="black", size=16, face="bold"))


plot(p1)


ggsave("../Output/toxo_cdc/ME49_59/figures_paper/atac_trnas_points.pdf", 
       plot = p1,
       height = 5, width = 7, dpi = 300)


L.trans.rna <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_based_transition_points_v2.rds')

rna_peaks.dat <- L.trans.rna$spline.fit.peaks.smooth %>% 
  transmute(g = x, y = s0, yp = s1) %>% pivot_longer(-g, names_to = 'drivs', values_to = 'value')
rna_peaks.dat$value[rna_peaks.dat$value < 0] <- 0
rna_peaks.dat$value[rna_peaks.dat$value > 6] <- 6


rna_peaks.dat$drivs <- factor(rna_peaks.dat$drivs, levels = c('y', 'yp'))
p2  <- ggplot(rna_peaks.dat, aes(x= g,y=value)) +
  geom_path(aes(color = drivs),alpha = 0.8, size = 1.2)+ 
  scale_color_manual(values = c("y" = "blue3","yp" ='chartreuse4'))+
  theme_bw(base_size = 14) +
  geom_vline(xintercept=L.trans.rna$transition.points$x, linetype=2, color = 'brown', size = 1) + 
  facet_grid(drivs~., scales = 'free')+
  ylab('peak time') + xlab('genes') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", colour = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", colour = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white", size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 18, face="bold", angle = 0))  +
  ggtitle("rna - transition") +
  theme(plot.title = element_text(size=18, face = "bold.italic", color = 'brown'),
        axis.title.x = element_text(size=20, face="bold", hjust = 1),
        axis.title.y = element_text(size=20, face="bold")) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=16, face="bold"),
    legend.text = element_text(colour="black", size=16, face="bold"))


plot(p2)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/rna_trnas_points.pdf", 
       plot = p2,
       height = 5, width = 7, dpi = 300)


## map to PCA 
rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_rna_atac_trnasition_v2.rds')
rna_sub@reductions[["pca"]]@cell.embeddings[,2] <- -1 * rna_sub@reductions[["pca"]]@cell.embeddings[,2]

getPcaMetaData.trans <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$phase,
                          transition.rna = S.O@meta.data$transition.rna, 
                          transition.atac = S.O@meta.data$transition.atac)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

#fb9062
#ff93ac
rna_sub@meta.data$spp2 <- "rna"
rna_sub_pca <- getPcaMetaData.trans(rna_sub)

## rna peak transition
p1  <- ggplot(rna_sub_pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = transition.rna,
    color = transition.rna
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + 
  #ggtitle("Transition") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p1)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_rna_transition.pdf", 
       plot = p1, width = 7, height = 6, dpi = 300)

plot_ly(umap.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
        color =  umap.data$inferred.phase, 
        size = 0.2)

pca.data <- FetchData(object = rna_sub, vars = c("PC_1", "PC_2", "PC_3", "phase", "transition.rna"))
pca.data$phase <- factor(pca.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
pca.data$transition.rna <- factor(pca.data$transition.rna, levels = c("T1", "T2", "T3", "T4"))
plot_ly(pca.data, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.data$phase, 
        colors = c("#b6232a",'#ed7202', '#caae05', '#6f883a', '#b138ee'),
        size = 0.4)


p2  <- ggplot(rna_sub_pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = transition.atac,
    color = transition.atac
  ), #color = 'blue',
  shape=21, size = 1)+
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +

  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) +
  #ggtitle("atac-transition") +
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p2)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_atac_transition.pdf",
       plot = p2, width = 7, height = 6, dpi = 300)


p3  <- ggplot(rna_sub_pca, aes(x= PC_1,y=PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + 
  #ggtitle("inferred-phase") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p3)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_inferred_phase.pdf", 
       plot = p3, width = 7, height = 6, dpi = 300)

## 

getPcaMetaData.atac.trans <- function(S.O){
  
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, 
                          phase = S.O@meta.data$phase,
                          #transition.rna = S.O@meta.data$transition.rna, 
                          transition.atac = S.O@meta.data$transition.atac)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data) 
}

atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_atac_atac_trnasition_v2.rds')
atac_sub@meta.data$spp2 <- "atac"
atac_sub_pca <- getPcaMetaData.atac.trans(atac_sub)

p4  <- ggplot(atac_sub_pca, aes(x= PC_1,y=-PC_2)) +
  geom_point(aes(#fill = lable.prob,
    fill = transition.atac,
    color = transition.atac
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  
  theme_bw(base_size = 14) +
  ylab('PC_2') + xlab('PC_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + 
  #ggtitle("atac-transition") + 
  guides(color = guide_legend(override.aes = list(size = 7)))

plot(p4)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/atac_PCA_atac_transition.pdf.pdf", 
       plot = p4, width = 7, height = 6, dpi = 300)


pca.data <- FetchData(object = atac_sub, vars = c("PC_1", "PC_2", "PC_3", "phase", "transition.atac"))
pca.data$phase <- factor(pca.data$phase, levels = c("G1.a", "G1.b", "S", "M", "C"))
pca.data$transition.atac <- factor(pca.data$transition.atac, levels = c("T1", "T2", "T3", "T4"))
plot_ly(pca.data, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.data$transition.atac, 
        size = 0.4)

plot_ly(pca.data, 
        x = ~PC_1, y = ~PC_2, z = ~PC_3,
        color =  pca.data$phase, 
        colors = c("#b6232a",'#ed7202', '#caae05', '#6f883a', '#b138ee'),
        size = 0.4)



## markers within each transition group 
prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1_ME49.xlsx') 
GO.MJ <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation.xlsx")

prod.desc <- left_join(prod.desc, GO.MJ, by = "TGME49")

rna.sig.markers.rna.trans <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
rna.sig.markers.rna.trans <- left_join(rna.sig.markers.rna.trans, prod.desc, by = c("GeneID" = "TGME49"))
rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>%
  mutate(category = ifelse(str_detect(ProductDescription, "hypothetical protein"), "hypo.", "others"))

ss.rna <- rna.sig.markers.rna.trans %>% group_by(cluster, category) %>% summarise(num.DEG = n())
sum(ss.rna$num.DEG)


ss.rna$cluster <- factor(ss.rna$cluster, levels = c('T1', 'T2', 'T3', 'T4'))


p <- ggplot(data = ss.rna, 
            aes(x=cluster, y=num.DEG, fill = category, color = cluster)) +
  geom_bar(stat="identity",position="stack",size = 2, width = 0.75, 
           aes(fill = cluster, color = category))+
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  geom_text(aes(label=num.DEG), vjust=1.5, color="white", size=10, fontface = 'bold')+
  # facet_grid(. ~ spp, scales = "free", space='free',labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 18))+
  theme(panel.spacing = unit(1.5, "lines")) +
  ggtitle("scRNA markers - rna transition")
#theme(legend.position = "none")
p

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_numbers_rna_transition.pdf",
       plot=p,
       width = 6, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


transition.markers.sig <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_sig_atac_trans_v2.rds')
ss.atac <- transition.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss.atac$cluster <- factor(ss.atac$cluster, levels = c('T1', 'T2', 'T3', 'T4'))
ss.atac

p <- ggplot(data = ss.atac, 
            aes(x=cluster, y=num.DEG, fill = cluster, color = cluster)) +
  geom_bar(stat="identity",size = 2, width = 0.75, aes(fill = cluster))+
  scale_color_manual(values = c("T1" = "#C2401F", 'T2' = '#caae05', 'T3' = '#6f883a', 'T4' = '#b138ee')) +
  scale_fill_manual(values = c("T1" = "#C2401F", 'T2' = '#caae05', 'T3' = '#6f883a', 'T4' = '#b138ee')) +
  geom_text(aes(label=num.DEG), vjust=1.5, color="white", size=10, fontface = 'bold')+
  # facet_grid(. ~ spp, scales = "free", space='free',labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 18))+
  theme(panel.spacing = unit(1.5, "lines")) +
  ggtitle("scRNA markers - atac transition")
#theme(legend.position = "none")
p

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_numbers_atac_transition.pdf",
       plot=p,
       width = 6, height = 5,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## heatmaps
 ## run the heatmaps_trans.r

## expr and access profile of matched genes withing transition of rna and atac
trans.list <- readRDS( "../Input/toxo_cdc/rds_ME49_59/sc_rna_atac_mus_trans_matched_genes_df_plot.rds")
pp <- grid.arrange(grobs = trans.list, ncol = 2)
ggsave(plot = pp , "../Output/toxo_cdc/ME49_59/figures_paper/dtw_clust_matched_rna_atac_transition.pdf",
       width = 12, height = 14, dpi = 300)


## AP2XII-8 expression profile 

# intra
S.O.rna  <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
DefaultAssay(S.O.rna) <- "RNA"
Idents(S.O.rna) <- 'phase'
levels(S.O.rna) <- c("G1.a", "G1.b", "S", "M", "C")

VlnPlot(S.O.rna, "TGME49-250800",assay = "RNA", log = F, slot = "data") 


p1 <- VlnPlot(S.O.rna, "TGME49-250800") +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  theme(axis.text.x  = element_text(face = "bold", size = 18, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(face = "bold", size = 16),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(face = "bold", size = 18),
        plot.title = element_text(face = "bold", size = 18)) + 
  theme(legend.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16))+
  theme(legend.position = "none")+
  theme(plot.title = element_text(size = 18, face = "bold", colour = "black"))+
  ggtitle("i.c.")+
  ylab("expr")
  
p1
# FeaturePlot(S.O.rna, "TGME49-250800", split.by = 'orig.ident',
#                 cols = c("grey", "blue"), reduction = "pca") 
# VlnPlot(S.O.rna, "TGME49-250800")


# extra
S.O.intra.extra <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_extra_lables_not_anchored_list.rds')
names(S.O.intra.extra) <- c("intra", "extra")
S.O.extra <- S.O.intra.extra$extra
DefaultAssay(S.O.extra) <- "RNA"
Idents(S.O.extra) <- 'phase'
levels(S.O.extra) <- c("G1.a", "G1.b", "S", "M", "C")

p2 <- VlnPlot(S.O.extra, "TGGT1-250800") +
  scale_fill_manual(values = c("G1.a" = "#b6232a","G1.b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  theme(axis.text.x  = element_text(face = "bold", size = 16, angle = 0, hjust = 0.5),
        axis.text.y  = element_text(face = "bold", size = 18),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold", size = 18)) + 
  ylab("") +
  theme(legend.text = element_text(face = "bold", size = 16),
        legend.title = element_text(face = "bold", size = 16),
        legend.position = "none")+
  theme(plot.title = element_text(size = 18, face = "bold", colour = "black"))+
  ggtitle("e.c.")

p2

p <- wrap_plots(list(p1, p2))
p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/AP2XII8_expresseion_intra_extra_violin.pdf", 
       plot = p, height = 4, width = 7, dpi = 300)
#prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
#TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## Map to ME49
#counts = S.O@assays$RNA@counts
#rownames(S.O.extra@assays$RNA@counts) <- gsub('_', '-', TGGT1_ME49$TGME49[match(gsub('-', '_', rownames(S.O.extra@assays$RNA@counts)), TGGT1_ME49$TGGT1)])


Idents(S.O.extra) <- 'phase'
levels(S.O.extra) <- c("G1.a", "G1.b", "S", "M", "C")
FeaturePlot(S.O.extra, "TGGT1-250800", split.by = 'orig.ident',
                  cols = c("grey", "blue"), reduction = "umap") 

VlnPlot(S.O.extra,"TGGT1-250800")


## TSS 
tss.df <- readRDS("../Input/toxo_cdc/rds_ME49_59/motif_dist_from_TSS.rds")
motif.name <- c(paste("T1.motif", 1:3, sep = "."), 
               paste("T2.motif", 4:7, sep = "."), 
               paste("T3.motif", 8:10, sep = "."), 
               paste("T4.motif", 11, sep = "."))

tss.df <- lapply(1:length(motif.name), function(i){
  
  tss <- tss.df[[i]] %>% mutate(motif = motif.name[i])
  
})
names(tss.df) <-  motif.name

tss.df.all <- Reduce("rbind", tss.df) 

tss.df.all <- tss.df.all %>% filter(!motif =="T3.motif.10" ) # not a significant motif

p <- ggplot(tss.df.all, aes(x = tss.dist)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 0.7,linetype = 1,colour = 2) +
  theme_bw() + 
  facet_wrap(~ motif, scales='free', ncol = 3) +
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"),
        axis.text.y =  element_text(size = 10, face = "bold", colour = "black"),
        axis.text.x =  element_text(size = 10, face = "bold", colour = "black", angle = 35))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines"))
p  

ggsave(plot = p , "../Output/toxo_cdc/ME49_59/figures_paper/motifs_dist_tss.pdf",
       width = 9, height = 8, dpi = 300)

#genes.with.motif <- readRDS("../Input/toxo_cdc/rds_ME49_59/motif_genes_dtw_clust_df_plot.rds")

genes.with.motif <- readRDS("../Input/toxo_cdc/rds_ME49_59/motif_genes_dtw_clust_df_plot_new_strandness.rds")

genes.with.motif.list <- lapply(genes.with.motif, "[[", 1) 
genes.with.motif.list <- genes.with.motif.list[-10]
motif.name <- c(paste("T1.motif", 1:3, sep = "."), 
                paste("T2.motif", 4:7, sep = "."), 
                paste("T3.motif", 8:9, sep = "."), 
                paste("T4.motif", 10, sep = "."))

genes.with.motif.list <- lapply(1:length(genes.with.motif.list), function(i){
  tmp <- genes.with.motif.list[[i]] %>% mutate(motif = motif.name[i])
})

plot_rna_trends_dtw <- function(df){
  p  <- ggplot(df, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 16, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ motif, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=16, face="bold", hjust = 1),
      axis.title.y = element_text(size=12, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  
}
#plot_rna_trends_dtw(genes.with.motif.list[[1]])

pp.list <- lapply(1:length(genes.with.motif.list), function(i) {
  p <- plot_rna_trends_dtw(genes.with.motif.list[[i]]) 
  p
  
})
pp <- grid.arrange(grobs = pp.list, ncol = 4)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/genes_with_motif_cluster_rna_new_strandness.pdf", 
       plot = pp, height = 15, width = 14, dpi = 300)

## cut and run clusters

plot_rna_trends <- function(df){
  p  <- ggplot(df, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 16, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=16, face="bold", hjust = 1),
      axis.title.y = element_text(size=12, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  
}

sc.rna.sc.atac.joint.long <- readRDS('../Input/toxo_cdc/rds_ME49_59/cut_run_sc_rna_sc_atac_joint_dtw_clust_6k_orig.rds')
p1 <- plot_rna_trends(sc.rna.sc.atac.joint.long) + ggtitle("cut & run targets")

plot(p1)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/cut_run_targets_clusters_6k.pdf",
       plot=p1,
       width = 6, height = 8,
       units = "in"
)
sc.rna.sc.atac.joint.long <- readRDS('../Input/toxo_cdc/rds_ME49_59/cut_run_sc_rna_sc_atac_joint_dtw_clust_4k_orig.rds')
p1 <- plot_rna_trends(sc.rna.sc.atac.joint.long) + ggtitle("cut & run targets")

plot(p1)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/cut_run_targets_clusters_4k.pdf",
       plot=p1,
       width = 6, height = 8,
       units = "in"
)

## 

## bar plot - DEGs KDvsWT
KD.vs.WT.phase.marker.sig.desc <- readRDS("../Input/toxo_cdc/rds_ME49_59/KD_vs_WT_markers_sig_NO_phase.rds")
KD.vs.WT.phase.marker.stat <- KD.vs.WT.phase.marker.sig.desc %>% group_by(cluster, dir, comp.group) %>% 
  summarise(genes = list(GeneID), num.DEG = n()) 
KD.vs.WT.phase.marker.stat$dir <- factor(KD.vs.WT.phase.marker.stat$dir, levels = c("activated", "repressed"))


KD.vs.WT.phase.marker.sig.desc %>% group_by(dir) %>% summarise(n())



p <- ggplot(data = KD.vs.WT.phase.marker.stat, aes(x=dir, y=num.DEG, fill = dir, color = dir)) +
  geom_bar(stat="identity",size = 2, width = 0.75, aes(fill = dir))+
  scale_color_manual(values = c("activated" = "#6565bf","repressed" ='#ee5d6c')) +
  scale_fill_manual(values = c("activated" = "#6565bf","repressed" ='#ee5d6c')) +
  geom_text(aes(label=num.DEG), vjust=1.5, color="white", size=10, fontface = 'bold')+
  theme_bw()+ xlab('')+
  theme(
    axis.text.x = element_text(size = 20, face = "bold", colour = "black"),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside", 
        panel.border = element_blank()) +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  #theme(legend.position = "none")
p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/KD_vs_WT_DEGs_bar_plot.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)


## cut& run peaks assigned to genes
peak.genes.bed.merged.bed <- readRDS("../Input/toxo_cdc/rds_ME49_59/Cut_run_peak_gene_assigned_final.rds")

KD.vs.WT.phase.marker.sig.desc <- readRDS("../Input/toxo_cdc/rds_ME49_59/KD_vs_WT_markers_sig_NO_phase.rds")
KD.vs.WT.phase.marker.stat <- KD.vs.WT.phase.marker.sig.desc %>% group_by(cluster, dir, comp.group) %>% 
  summarise(genes = list(GeneID), num.DEG = n()) 
KD.vs.WT.phase.marker.stat$dir <- factor(KD.vs.WT.phase.marker.stat$dir, levels = c("activated", "repressed"))

venn.list <- list(cutRun = peak.genes.bed.merged.bed$gene_name,
                  activated = unlist(KD.vs.WT.phase.marker.stat$genes[1]),
                  repressed = unlist(KD.vs.WT.phase.marker.stat$genes[2]))


pdf("../Output/toxo_cdc/ME49_59/figures_paper/KD_DEGs_cut_run_overlap.pdf",
    height = 6, width = 6)
ggvenn::ggvenn(venn.list, fill_color = c("orange", "#8080FA", "#fc6c85"),
               set_name_color = c("orange","#8080FA", "#fc6c85" ),
               set_name_size = 8,
               fill_alpha = 0.9,
               text_color = "black",
               text_size = 8,
               show_percentage = FALSE)


dev.off()
venn <- Venn(venn.list)
data <- process_data(venn)

p <- ggplot() +
  geom_sf(aes(fill = count), data = venn_region(data)) +
  geom_sf(aes(color = id), data = venn_setedge(data), show.legend = FALSE) +
  geom_sf_label(aes(label = count), data = venn_region(data),size =4,  alpha = 0.5) +
  geom_sf_text(aes(label = name), color=c("cutRun" = "orange","activated" ="#8080FA",
                                          'repressed' = '#fc6c85'),
               fontface = "bold", size= 3,nudge_y = 0,nudge_x=0,
               data = venn_setlabel(data))+
  geom_sf(color=c("cutRun" = "orange","activated" ="#8080FA",
                  'repressed' = '#fc6c85'), 
          size = 1, data = venn_setedge(data), show.legend = FALSE) + 
  theme_void()+ 
  theme(legend.position = "none")+
  theme(plot.margin=unit(c(0,0,0,0),"cm"))+
  #ggtitle(clusters[i])+
  theme(plot.title = element_text(size = "22", face = "bold.italic"))
p

listforVenn <- readRDS("./rds/cell_cycle_markers_shared_across_speciies_Venn_diag.rds")
ps <- lapply(1:length(clusters), function(i){
  
  
  
  
})
grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 4)

## plot multiomics
getPcaMetaData.multiome <- function(S.O){
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2, UMAP_3 = UMAP_3)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), 
                          spp = S.O@meta.data$spp2, phase = S.O@meta.data$inferred.phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  return(meta.data)  
}

## old multiomics 
S.O.integrated <- readRDS("../Input/toxo_cdc/rds_ME49_59/S.O.integrated_rna_WT_KD_atac_WT_KD_reference_rna_WT_inferred_phases.rds")

S.O.integrated@reductions$pca@cell.embeddings[,2] <- -1 * S.O.integrated@reductions$pca@cell.embeddings[,2]
S.O.integrated$orig.ident <- factor(S.O.integrated$orig.ident, levels = c("scRNA", "scRNA_KD", "scATAC", "scATAC_KD"))
S.O.integrated$inferred.phase <- factor(S.O.integrated$inferred.phase, levels = c("G1a", "G1b", "S", "M", "C"))

DefaultAssay(S.O.integrated)


S.O.rna.KD <- subset(S.O.integrated, ident = "scRNA.KD")
S.O.rna.KD@meta.data$spp2 <- "rna.KD"
pca.KD <- getPcaMetaData.multiome(S.O.rna.KD)

p1  <- ggplot(pca.KD, aes(x=UMAP_2 ,y=UMAP_3)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1a" = "#b6232a","G1b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1a" = "#b6232a","G1b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + 
  ggtitle("scRNA.KD") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p1)


## WT 

S.O.rna.WT<- subset(S.O.integrated, ident = "scRNA")
S.O.rna.WT@meta.data$spp2 <- "rna.WT"
pca.WT <- getPcaMetaData.multiome(S.O.rna.WT)

p2  <- ggplot(pca.WT, aes(x=UMAP_2 ,y=UMAP_3)) +
  geom_point(aes(#fill = lable.prob,
    fill = phase,
    color = phase
  ), #color = 'blue', 
  shape=21, size = 1)+ 
  scale_color_manual(values = c("G1a" = "#b6232a","G1b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee')) +
  scale_fill_manual(values = c("G1a" = "#b6232a","G1b" ='#ed7202', 'S' = '#caae05', 'M' = '#6f883a', 'C' = '#b138ee'))+
  
  theme_bw(base_size = 14) +
  ylab('UMAP_2') + xlab('UMAP_1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) + 
  ggtitle("scRNA.WT") + 
  guides(color = guide_legend(override.aes = list(size = 7)))


plot(p2)

pp <- p2 / p1
pp
# ggsave("../Output/toxo_cdc/ME49_59/figures_paper/PCA_sRNA_WT_KD.pdf", 
#        height = 12, width = 6, dpi = 300, plot = pp)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/UMAP_sRNA_WT_KD_multiomics_old.pdf", 
       height = 8, width = 14, dpi = 300, plot = pp)

