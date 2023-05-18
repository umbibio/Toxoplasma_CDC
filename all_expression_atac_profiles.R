
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

# expression and atac profile of one gene at a time
plot_trends <- function(my.GeneID, sc.rna.spline.fits,sc.atac.spline.fits ){
  
  
  ## Turn the data into wide format (time by gene) and center & scale each gene
  sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = F, scale = F)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  my.rna <- sc.rna.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  my.atac <- sc.atac.mu.scale %>% dplyr::filter(GeneID == my.GeneID)
  
  p1  <- ggplot(my.rna , aes(x= x,y=expr)) +
    geom_line(color = 'blue',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('rna', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  
  p2  <- ggplot(my.atac , aes(x= x,y=expr)) +
    geom_line(color = 'red',alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    ggtitle(paste('atac', my.GeneID)) + 
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'black'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) 
  
  p <- grid.arrange(p1, p2)
  
  return(p)
}


## this function plots rna and atac profile of genes in a table of interest

plot_rna_atac <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, size = 20, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 20, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(. ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=20, face="bold", hjust = 1),
      axis.title.y = element_text(size=20, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  return(p)
  
}


## this function extracts the expression and accessibility profile of genes of interest
## need to give rna and atac splines and genes of interest as input
## set the scale T/F
get_rna_atac_profile <- function(rna.splines, atac.splines, genes.tab, scale = T) {


  sc.rna.dtw.wide <- rna.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  sc.atac.dtw.wide <- atac.splines %>% 
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
    mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = scale)) %>%
    as.data.frame()
  
  
  sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
    pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')
  
  
  ## table of genes to plot their expression 
  tab.genes <- data.frame(TGME49 = gsub("_", "-", genes.tab$gene_name), 
                          Name = HC.peaks$ProductDescription)
  
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.rna.long <- sc.rna.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name)  %>% distinct()
  
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab.genes, by = c('GeneID' = 'TGME49')) 
  sc.atac.long <- sc.atac.long %>% 
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Name) %>% distinct()
  
  
  sc.rna.sc.atac.joint <- inner_join(sc.rna.long, sc.atac.long, 
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "scATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
    pivot_longer(-c('time', "GeneID", "Name"), 
                 names_to = 'data', values_to = 'normExpr') 
  
  return(sc.rna.sc.atac.joint.long)
  
}


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## table of genes of interest
## cut and run intersection and global KD vs WT
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

## global 
tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

## phase based 
tab <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"


HC.peaks <- tab %>% filter(Category == "ribosomal")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))


expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                      atac.splines = sc.atac.spline.fits, 
                                      genes.tab = HC.peaks, scale = F)

## if you wan the atac profiles as well, then do not filter data == "scRNA
expr.tab <- expr.atac.tab %>% filter(data == "scRNA")
p1 <- plot_rna_atac(expr.tab)
#p1 <- p1 + ggtitle("ribosomals")
p1 
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_Global_KD_vs_WT_ribosomal_No_clustering.pdf", 
       plot = p1, height = 3, width = 6, dpi = 300)

plot_trends("TGME49-227600",sc.rna.spline.fits, sc.atac.spline.fits)

## table of genes of interest
## cut and run intersection and phase based KD vs WT 
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"

HC.peaks <- HC.peaks %>% filter(Category == "ribosomal")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                      atac.splines = sc.atac.spline.fits, 
                                      genes.tab = HC.peaks, scale = T)

expr.tab <- expr.atac.tab %>% filter(data == "scRNA")
p1 <- plot_rna_atac(expr.tab)
p1 <- p1 + ggtitle("ribosomals")
p1 
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_phase_based_KD_vs_WT_ribosomal_No_clustering.pdf", 
       plot = p1, height = 4, width = 6, dpi = 300)


