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

## scDATA

rna_sub <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes_1.2.rds')

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




plot_trends <- function(my.GeneID){
  
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




gene <- "TGME49-283900"
gene <- 'TGME49-232970'
#gene <- gsub("_", "-",gene)
p <- plot_trends(gene)
p



stats.cyclic.AP2 <- stats %>% dplyr::filter(rna.cyclic == 1 & grepl("AP2 domain transcription factor", ProductDescription))
AP2s.c <- stats.cyclic.AP2$GeneID

for(i in 1:length(unique(sc.rna.sc.atac.joint.long$GeneID))){
  gene <- unique(sc.rna.sc.atac.joint.long$GeneID)[i]
  p <- plot_trends(gene)
  Sys.sleep(0.6)
  
}




regs <- read.xlsx('../Input/toxo_cdc/gene_families/regulatory_toxo.xlsx')
regs <- left_join(regs, stats, by = c('GeneID' = 'GeneID.y'))
cyclic.regs <- regs %>% dplyr::filter(rna.cyclic == 1)

DefaultAssay(rna_sub) <- 'RNA'
x <- FetchData(rna_sub, vars = cyclic.regs$GeneID.y)
dim(x)


library(ggfortify)
df <- t(x)

vars <- apply(df, 2, sd)
rm.ind <- which(vars == 0)
df <- df[,-rm.ind]
pca_res <- prcomp(df, scale. = TRUE)

df <- as.data.frame(df)
df$GeneID <- gsub('-1', '', gsub('\\.', '-', rownames(df)))
df$class <- regs$Family[match(df$GeneID, regs$GeneID.y)]

autoplot(pca_res, data = df, colour = 'class', label = F, label.size = 3)
