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
prod.desc  <- read.xlsx('../Input_YR/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input_YR/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

##  DEGs
DEG.sig <- readRDS( '../Input_YR/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')

## New
gene.fam <- read.xlsx("../Input/toxo_cdc/gene_families/gene_fam_KZ.xlsx")
gene.fam <- left_join(gene.fam, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

# scRNA and scATAC
rna_sub <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_KZ///toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

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



stats <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/all_genes_cyclic_timing.rds')

##
stats.cyclic.rna <- stats %>% dplyr::filter(rna.cyclic == 1)


## Filter out familes with ony one member
fam.stat <- gene.fam %>% group_by(Family) %>% summarise(num.mem = n())
gene.fam <- gene.fam %>% filter(!(Family %in% fam.stat$Family[fam.stat$num.mem <= 4]))
gene.fam <- left_join(gene.fam, stats, by = c('GeneID' = 'GeneID.y'))
gene.fam <- gene.fam %>% na.omit()
gene.fam$Family <- factor(gene.fam$Family)

fam.ord <- gene.fam %>% group_by(Family) %>% summarise(m = median(rna.peak.time))
gene.fam$Family <- factor(gene.fam$Family, levels = fam.ord$Family[sort(fam.ord$m, index.return = T)$ix])

p1 <- ggplot(data = gene.fam, aes(x = Family, y = rna.peak.time, color = Family)) + 
  geom_boxplot() + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('peak expression time') + xlab('Family') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
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


p1

ggsave(filename='../Output_KZ/figures/peak_expr_dist_gene_fam.pdf',
       plot=p1,
       width = 6, height = 6,
       units = "in")
    




p2 <- ggplot(data = gene.fam, aes(x = Family, y = atac.peak.time, color = Family)) + 
  geom_boxplot() + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('peak access time') + xlab('Family') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
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



p2
ggsave(filename='../Output_KZ/figures/peak_access_dist_gene_fam.pdf',
       plot=p2,
       width = 6, height = 6,
       units = "in")






## Both
gene.fam.long <- gene.fam %>% dplyr::select(Family, GeneID, rna.peak.time, atac.peak.time) %>%
  pivot_longer(-c('Family', 'GeneID'), names_to = 'mode', values_to = 'peak_time')
p1 <- ggplot(data = gene.fam.long, aes(x = Family, y = peak_time, color = mode)) + 
  geom_boxplot(width = 0.5, position = position_dodge(width=1)) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + 
  geom_jitter(shape=16,  position = position_jitterdodge(0.1)) + 
  #geom_jitter(shape=16,  position=position_jitter(0.2)+position_dodge(width=1)) + 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('peak expression time') + xlab('Family') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(legend.position = c(0.75, 0.2),
    #legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


p1

ggsave(filename='../Output_KZ/figures/peak_expr_access_dist_gene_fam.pdf',
       plot=p1,
       width = 7, height = 6,
       units = "in")





fams <- unique(gene.fam$Family) %>% sort()

ks <- c(2, 2, 2, 4, 2, 2, 2, 2, 2, 2, 2,2)
sc.rna.sc.atac.clusts <- lapply(1:length(fams), function(i){
  print(i)
  print(fams[i])
  tmp <- gene.fam %>% dplyr::filter(Family == fams[i])
  tmp$TGME49 <- gsub('_', '-', tmp$TGME49)
  sc.rna.tab <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tmp$TGME49 ]
  sc.atac.tab <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tmp$TGME49 ]
  
  sc.rna.tab.markers.hc_dtw <- dtwClustCurves(sc.rna.tab, nclust = ks[i])
  sc.atac.tab.markers.hc_dtw <- dtwClustCurves(sc.atac.tab, nclust = ks[i])
  
  #plot(sc.rna.tab.markers.hc_dtw, type = 'sc')
  #plot(sc.atac.tab.markers.hc_dtw, type = 'sc')

  return(list(sc.rna.tab.markers.hc_dtw, sc.atac.tab.markers.hc_dtw))
})


names(sc.rna.sc.atac.clusts) <- fams

ks <- c(1, 4, 3, 4, 4, 1, 3, 3, 2, 2, 2,2)
sc.rna.sc.atac.joint.long.list <- lapply(1:length(sc.rna.sc.atac.clusts), function(i){
  print(i)
  print(fams[i])
  tmp <- gene.fam %>% dplyr::filter(Family == fams[i])
  tmp$TGME49 <- gsub('_', '-', tmp$TGME49)
  sc.rna.tab <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tmp$TGME49 ]
  sc.atac.tab <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tmp$TGME49 ]
  
  
  sc.rna.tab.markers.hc_dtw <- sc.rna.sc.atac.clusts[[i]][[1]]
  sc.atac.tab.markers.hc_dtw <- sc.rna.sc.atac.clusts[[i]][[2]]
  
  sc.rna.long.tab <- inner_join(sc.rna.mu.scale, tmp, by = c('GeneID' = 'TGME49')) 
  sc.rna.long.tab <- sc.rna.long.tab %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()
  
  sc.atac.long.tab <- inner_join(sc.atac.mu.scale, tmp, by = c('GeneID' = 'TGME49')) 
  sc.atac.long.tab <- sc.atac.long.tab %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()
  
  
  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.tab), cluster = cutree(sc.rna.tab.markers.hc_dtw, k = ks[i]))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.tab), cluster = cutree(sc.atac.tab.markers.hc_dtw, k = ks[i]))
  
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
  
  return(sc.rna.sc.atac.joint.long)
})
names(sc.rna.sc.atac.joint.long.list) <- fams

## GGplot cluster graphs
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = Name,),alpha = 0.8, linewidth = 0.8)+ 
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

#out.dir <- '../Input/toxo_cdc/rds/'

out.pic.dir <- '../Output_KZ/figures/gene_fam_scaled//'
lapply(1:length(sc.rna.sc.atac.joint.long.list), function(i){
  #out.file <- paste0(out.dir, names(sc.rna.tab.markers.hc_dtw)[i], 'sc_rna_sc_atac_joint_dtw_clust_new_update.rds')
  #saveRDS(c.rna.tab.markers.hc_dtw[[i]], out.file)
  sc.rna.sc.atac.joint.long.list[[i]]$data <- factor(sc.rna.sc.atac.joint.long.list[[i]]$data, 
                                           levels = c("scRNA", "scATAC"))
  
  p <- plot_rna_atac_trends(sc.rna.sc.atac.joint.long.list[[i]])
  num.clust <- length(unique(sc.rna.sc.atac.joint.long.list[[i]]$cluster.RNA))
  out.pic  <- paste0(out.pic.dir, fams[i],'_', ks[i], '_clusters.pdf')
  ggsave(filename=out.pic,
         plot=p,
         width = 10, height = 8,
         units = "in"
  )
  
})
  


