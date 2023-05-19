
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
library(bigmemory)


source('./util_funcs.R')
source('./util_funcs_YR.R')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


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



trans.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")
trans.all <- do.call("rbind", trans.list)
trans.all <- trans.all %>% dplyr::select(GeneID, cluster.RNA.ordered, group,Name) %>% distinct()
trans.all$trans.cluster.rna <- paste(trans.all$group, trans.all$cluster.RNA.ordered, sep = "_")
colnames(trans.all) <- gsub("GeneID", "gene_name", colnames(trans.all))
colnames(trans.all) <- gsub("Name", "ProductDescription", colnames(trans.all))
trans.all.sum <- trans.all %>% 
  group_by(trans.cluster.rna) %>% summarise(genes = list(gsub("-", "_", gene_name)), total = n())


trans.clust.rna.list  <- split(trans.all, f = trans.all$trans.cluster.rna)
atac.list <- c()

out.dir <- "../OutPut/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/"
i <- 16
k <- 2
atac.list <- lapply(1:length(trans.clust.rna.list), function(i) {
  my.df <- trans.clust.rna.list[[i]]
  df <- clust.atac.df(my.df, num.clust = k)
  df$cluster.ATAC <- gsub(" ","", df$cluster.ATAC)
  df$group <- names(trans.clust.rna.list)[i]
  df$trans.rna.atac.clust <- paste(df$group, df$cluster.ATAC, sep = "_")
  df <- left_join(df, prod.desc, by = "GeneID")
  p1 <- plot_atac_trand(df) +
    ggtitle(names(trans.clust.rna.list[i]))
  p1
  
  # ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/T1C3.pdf", 
  #        plot = p1, height = 3, width = 5, dpi = 300)
  df.Go.tab <- df %>% dplyr::select(GeneID, trans.cluster.rna, trans.rna.atac.clust) %>% distinct()
  write.xlsx(df.Go.tab, paste(out.dir, paste(unique(df.Go.tab$trans.cluster.rna), ".xlsx", sep = ""), sep = ""))

  my.df.GO <- my.df %>% select(gene_name, trans.cluster.rna) %>% distinct()
  write.xlsx(my.df.GO, paste(out.dir, paste(unique(my.df$trans.cluster.rna), ".xlsx", sep = ""), sep = ""))
  
})

## T1, C1


expr.atac.tab <- get_rna_atac_profile(rna.splines = sc.rna.spline.fits, 
                                      atac.splines = sc.atac.spline.fits, 
                                      genes.tab = my.df, scale = T)

## if you wan the atac profiles as well, then do not filter data == "scRNA
expr.tab <- expr.atac.tab %>% filter(data == "scATAC")
p1 <- plot_rna_atac(expr.tab)
p1 
ggsave("../OutPut/toxo_cdc/ME49_59/figures_paper/T1C1C1.pdf", 
       plot = p1, height = 3, width = 5, dpi = 300)


prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )


in.dir <- "../OutPut/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters_manual/"
files <- gsub(".xlsx", "", list.files(in.dir))


for (i in 1: length(files)){
  
  tab <- read.xlsx(paste(in.dir, paste(files[i], ".xlsx", sep = ""), sep = ""))
  colnames(tab)[1] <- "gene_name"
  tab$gene_name <- gsub("-", "_", tab$gene_name)
  tab <- inner_join(tab, prod.desc, by = c("gene_name" = "TGME49"))
  
  write.xlsx(tab, paste(in.dir, paste(files[i], "desc.xlsx", sep = "_") , sep = ""))
  
}
