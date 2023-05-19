
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

##  DEGs
DEG.sig <- readRDS( '../Input/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')

## New

BCs <- read.xlsx("../Input/Toxo_genomics/genes/BC gene list 021723.xlsx")
BCs <- BCs[-1,] %>% dplyr::select(TGGT1, Name )
BCs <- left_join(BCs, TGGT1_ME49, by = "TGGT1")
BCs$TGME49 <- gsub("_", "-", BCs$TGME49)
colnames(BCs) <- gsub("Name", "ProductDescription", colnames(BCs))


tab <- BCs
sig.tab <- tab[tab$TGME49 %in% c(DEG.sig$gene),]

# scRNA and scATAC
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

na.ind <- which(apply(sc.rna.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.rna.dtw.wide <- sc.rna.dtw.wide[,-na.ind]
}


sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.atac.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.atac.dtw.wide <- sc.atac.dtw.wide[,-na.ind]
}

sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


## Clustering BCs

k <- 4
sc.rna.tab <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% sig.tab$TGME49 ]
sc.atac.tab <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% sig.tab$TGME49 ]

sc.rna.tab.markers.hc_dtw <- dtwClustCurves(sc.rna.tab, nclust = k)
sc.atac.tab.markers.hc_dtw <- dtwClustCurves(sc.atac.tab, nclust = k)

plot(sc.rna.tab.markers.hc_dtw, type = 'sc')
plot(sc.atac.tab.markers.hc_dtw, type = 'sc')


## GGplot cluster graphs

sc.rna.long.tab <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'TGME49')) 
sc.rna.long.tab <- sc.rna.long.tab %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()

sc.atac.long.tab <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'TGME49')) 
sc.atac.long.tab <- sc.atac.long.tab %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()


sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.tab), cluster = cutree(sc.rna.tab.markers.hc_dtw, k = k))
sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.tab), cluster = cutree(sc.atac.tab.markers.hc_dtw, k = k))

sc.rna.long.clust <- inner_join(sc.rna.long.tab, sc.rna.clust.info, by = 'GeneID')
sc.atac.long.clust <- inner_join(sc.atac.long.tab, sc.atac.clust.info, by = 'GeneID')

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


saveRDS(sc.rna.sc.atac.joint.long, '../Input/toxo_cdc/rds_ME49_59/BC_sc_rna_sc_atac_joint_dtw_clust_new_update.rds')

sc.rna.sc.atac.joint.long <- readRDS('../Input/toxo_cdc/rds_ME49_59/BC_sc_rna_sc_atac_joint_dtw_clust_new_update.rds')
sc.rna.sc.atac.joint.long$data <- factor(sc.rna.sc.atac.joint.long$data, 
                                         levels = c("scRNA", "scATAC"))
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
    geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
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

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/BCs_dtw_4_clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)


