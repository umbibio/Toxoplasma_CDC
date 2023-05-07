library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)


source('./util_funcs.R')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## from the review paper
Cyclic.AP2s <- read.xlsx("../Input/Toxo_genomics/genes/microarray_cyclic_timing2.xlsx")
Cyclic.AP2s <- Cyclic.AP2s %>% dplyr::select(GeneID, ProductDescription, cyclic) %>% distinct() %>%
  filter(cyclic == 1) %>% filter(grepl("^AP2",ProductDescription ))

Cyclic.AP2s$Name <- unlist(lapply(str_split(Cyclic.AP2s$ProductDescription,pattern = " "), "[[", 5))
write.xlsx(Cyclic.AP2s, "../Input/Toxo_genomics/genes/Cyclic_AP2s_review_paper.xlsx")

## 

rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

## Clustering AP2s

sc.rna.AP2 <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% Cyclic.AP2s$GeneID ]
sc.atac.AP2 <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% Cyclic.AP2s$GeneID ]

sc.rna.AP2.markers.hc_dtw <- dtwClustCurves(sc.rna.AP2, nclust = 4L)
sc.atac.AP2.markers.hc_dtw <- dtwClustCurves(sc.atac.AP2, nclust = 4L)

plot(sc.rna.AP2.markers.hc_dtw, type = 'sc')
plot(sc.atac.AP2.markers.hc_dtw, type = 'sc')


## GGplot cluster graphs

sc.rna.long.AP2 <- inner_join(sc.rna.mu.scale, Cyclic.AP2s, by = c('GeneID')) 
sc.rna.long.AP2 <- sc.rna.long.AP2 %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()

sc.atac.long.AP2 <- inner_join(sc.atac.mu.scale, Cyclic.AP2s, by = c('GeneID')) 
sc.atac.long.AP2 <- sc.atac.long.AP2 %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()


sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.AP2), cluster = cutree(sc.rna.AP2.markers.hc_dtw, k = 4))
sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.AP2), cluster = cutree(sc.atac.AP2.markers.hc_dtw, k = 4))

sc.rna.long.clust <- inner_join(sc.rna.long.AP2, sc.rna.clust.info, by = 'GeneID')
sc.atac.long.clust <- inner_join(sc.atac.long.AP2, sc.atac.clust.info, by = 'GeneID')

sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust, 
                                   by = c("time", "GeneID", "Name"))
colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
  pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"), 
               names_to = 'data', values_to = 'normExpr') 


sc.rna.sc.atac.joint.long$label <- NA
sc.rna.sc.atac.joint.long$label[which(sc.rna.sc.atac.joint.long$time == 3)] <-
  sc.rna.sc.atac.joint.long$Name[which(sc.rna.sc.atac.joint.long$time == 3)]

sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)

saveRDS(sc.rna.sc.atac.joint.long, '../Input/toxo_cdc/rds_ME49_59/Cyclic_33_AP2s_sc_rna_sc_atac_dtw_4_clust.rds')
#saveRDS(sc.rna.sc.atac.joint.long, '../Input/toxo_cdc/rds_ME49_59/Cyclic_33_AP2s_sc_rna_sc_atac_dtw_5_clust.rds')

plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = Name,),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    geom_text_repel(aes(label = label), size = 3.4, fontface = "bold",
                    box.padding = unit(0.6, "lines"),
                    max.overlaps = 300,
                    #segment.angle = 180,
                    nudge_x = 0.25, 
                    nudge_y = 0.25,
                    hjust=0.25,
                    #nudge_x=0.25, 
                    segment.size = 0.1,
                    na.rm = TRUE)+ 
    facet_grid(cluster.RNA~data, scales = 'free', space = 'free') +
    
    
    #ggtitle(titles[i]) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

p1 <- plot_rna_atac_trends(sc.rna.sc.atac.joint.long)

plot(p1)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/AP2_33_Cyclic_4_Clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/AP2_33_Cyclic_5_Clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)

