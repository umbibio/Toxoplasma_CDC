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


## scDATA

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


## 

## This function performs dynamic time warping clustering for both rna and atac profiles
## it gets a table of genes, 
## the table should include the product description (gene ID and gene Name) 
## you can specify the number of clusters you are looking for (cannot be less than 2)

clust.df <- function(tab , num.clust) {
  
  k <- num.clust
  sc.rna <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% tab$gene_name ]
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac<- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% colnames(sc.rna)]
  
  sc.rna.markers.hc_dtw <- dtwClustCurves(sc.rna, nclust = k)
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  
  tab <- tab[tab$gene_name %in% colnames(sc.rna),]
  
  sc.rna.long <- inner_join(sc.rna.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.rna.long <- sc.rna.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription)  %>% distinct()

  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, Name = ProductDescription) %>% distinct()


  sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna), cluster = cutree(sc.rna.markers.hc_dtw, k = k))
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))

  sc.rna.long.clust <- inner_join(sc.rna.long, sc.rna.clust.info, by = 'GeneID')
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')

  sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust,
                                     by = c("time", "GeneID", "Name"))
  colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
  
  sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
    pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"),
                 names_to = 'data', values_to = 'normExpr')

  sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)
  sc.rna.sc.atac.joint.long$cluster.ATAC <- paste('C', sc.rna.sc.atac.joint.long$cluster.ATAC)
  
  return(sc.rna.sc.atac.joint.long)
  
}


## plot the expression and accessibility of genes within each cluster
# facet by rna cluster
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

## plot the expression and accessibility of genes within each cluster
#facet by atac cluster
plot_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

###########################################################################
###########################################################################

## genes with the corresponding motif in each transition cluster
## in v2 strand info is disabled
# input.dir <- "../Input/toxo_cdc/rds_ME49_59/bed_files/rna_transitions_v2/motif_beds_genes/" 
# file.names <- list.files(input.dir, "*_GO.xlsx")

# in v4 strandness option in getfasta has been enabled

input.dir <- "../Input/toxo_cdc/bed_files/rna_transitions_v4/genes_bed_BAMM/"
file.names <- list.files(input.dir, ".xlsx")

genes.list <- lapply(1:length(file.names), function(i){
  
  df <- read.xlsx(paste(input.dir, file.names[i],  sep = ""))
  df$V13 <- gsub("_", "-", df$V13)
  df <- df %>% distinct(V13, .keep_all = T)
  colnames(df) <- gsub("V13", "gene_name", colnames(df))
  return(df)
  
})
#names(genes.list) <- gsub(".xlsx", "", gsub("_GO", "", file.names))
names(genes.list) <- gsub(".xlsx", "", file.names)


pp.list <- lapply(1:length(genes.list), function(i) {
  
  my.df <- genes.list[[i]]
  df <- clust.df(my.df, num.clust = 4)
  df <- df[df$data == "scRNA",]
  p1 <- plot_rna_atac_trends(df) +
    ggtitle(names(genes.list[i]))
  p1
  
})
names(pp.list) <- names(genes.list)


saveRDS(pp.list, "../Input/toxo_cdc/rds_ME49_59/motif_genes_dtw_clust_df_plot_new_strandness.rds")

## T3_motif3 is not significant , so remove it 
pp <- grid.arrange(grobs = pp.list[-10], ncol = 4)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/genes_in_each_rna_transition_with_motif_dtw_clust.pdf", 
        plot = pp,
        height = 14,width = 14, dpi = 300)
ggsave( "../Output/toxo_cdc/ME49_59/figures/T1_motif_3_genes_dtw_clust.pdf", 
        plot = pp.list[[3]],
        height = 6,width = 6, dpi = 300)


########################################################
## coexpressed genes in each transition (rna based) ####
########################################################
## Fig 5 B

rna.trans.marker.genes <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
rna.trans.marker.genes <- rna.trans.marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()
rna.trans.marker.genes %>% group_by(phase) %>% summarise(n())

prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1"))
prod.desc <- prod.desc %>% 
  transmute(GeneID = gsub("_", "-", TGME49), ProductDescription = ProductDescription) %>% na.omit()

rna.trans.marker.genes <- left_join(rna.trans.marker.genes, prod.desc, by = "GeneID")
colnames(rna.trans.marker.genes)[1] <- "gene_name"

rna.trans.marker.genes.list <- split(rna.trans.marker.genes, rna.trans.marker.genes$phase) 

trans.list <- c()
trans.list <- lapply(1:length(rna.trans.marker.genes.list), function(i) {
  
  my.df <- rna.trans.marker.genes.list[[i]]
  df <- clust.df(my.df, num.clust = 4)
  df$group <- names(rna.trans.marker.genes.list)[i]
  #df <- df[df$data == "scRNA",]
  p1 <- plot_rna_atac_trends(df) +
    ggtitle(names(rna.trans.marker.genes.list[i]))
  p1
  #tt <- list(df, p1)
})

names(trans.list) <- names(rna.trans.marker.genes.list)
saveRDS(trans.list, "../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list.rds")


pp <- grid.arrange(grobs = trans.list, ncol = 2)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/rna_markers_rna_transitions_dtw_clust_list.pdf", 
        plot = pp,
        height = 12,width = 14, dpi = 300)


trans.clust.df <- lapply(trans.list, "[[", 1) 
trans.clust.df.all <- do.call("rbind", trans.clust.df)
trans.clust.df.all <- trans.clust.df.all %>% dplyr::select(GeneID, cluster.RNA, group) %>% distinct()
trans.clust.df.all.sum <- trans.clust.df.all %>% 
  group_by(group, cluster.RNA) %>% summarise(genes = list(gsub("-", "_", GeneID)), total = n())

write.xlsx(trans.clust.df.all, "../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_clust_based.xlsx")
write.xlsx(trans.clust.df.all.sum, "../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_clust_based_sum.xlsx") 
#trans.clust.df.all.sum <- split(trans.clust.df.all.sum, f = trans.clust.df.all.sum$group)


##########################################################
# intersection of 4 CUT&RUN + DEG (KD_vs_WT, phase based)
##########################################################
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"


### cluster up_reg genes
HC.peaks <- HC.peaks %>% filter(dir == "up_reg")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
#HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) 
p
ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_up_reg_KD_vs_WT_phase_based_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

## 
HC.peaks.clust.df <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_up_reg_KD_vs_WT_phase_based_3_clust.xlsx")

###  cluster down_reg genes
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"

HC.peaks <- HC.peaks %>% filter(dir == "down_reg")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
#HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) 
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_KD_vs_WT_phase_based_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)
HC.peaks.clust.df <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_KD_vs_WT_phase_based_3_clust.xlsx")


### cluster ribosomal genes 
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

HC.peaks <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"


HC.peaks <- HC.peaks %>% filter(Category == "ribosomal")
colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 3)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) 
p


ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_KD_vs_WT_ribosomal_3_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)


###########################################################
# intersection of 4 CUT&RUN + DEG (KD_vs_WT, Global)
##########################################################
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

direct <- "down_reg"
HC.peaks <- tab %>% filter(dir == direct)

colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_Global_KD_vs_WT_2_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

HC.peaks.clust.df.down <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df.down, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_Global_KD_vs_WT_2_clust.xlsx")


# up-reg
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

direct <- "up_reg"
HC.peaks <- tab %>% filter(dir == direct)

colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]

p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p

ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_up_reg_Global_KD_vs_WT_2_clust.pdf', 
       plot = p, width = 6, height = 6, dpi = 300)

HC.peaks.clust.df.up <- HC.peaks.clust %>% dplyr::select(GeneID, Name, cluster.RNA) %>% distinct()
write.xlsx(HC.peaks.clust.df, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_up_reg_Global_KD_vs_WT_2_clust.xlsx")


up.down.global <- rbind(HC.peaks.clust.df.up, HC.peaks.clust.df.down)


## ribosomals (which are only present in down-reg list)
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"

tab <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT, ProductDescription.x , Category ) %>% 
  distinct()
names(tab)[9] <- "dir"

HC.peaks <- tab %>% filter(Category == "ribosomal")


colnames(HC.peaks) <- gsub("ProductDescription.x", "ProductDescription", colnames(HC.peaks))

HC.peaks$gene_name <- gsub("_", "-", HC.peaks$gene_name) 
HC.peaks.clust <- clust.df(HC.peaks, num.clust = 2)
HC.peaks.clust <- HC.peaks.clust[HC.peaks.clust$data == "scRNA",]
#p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle(direct)
p <- plot_rna_atac_trends(HC.peaks.clust) + ggtitle("ribosomals")
p



# ribosomals
write.xlsx(HC.peaks.clust, "../Output/toxo_cdc/ME49_59/tables/High_conf_peaks_down_reg_Global_KD_vs_WT_ribosomal_2_clust.xlsx")
ggsave('../Output/toxo_cdc/ME49_59/figures_paper/High_conf_peaks_down_reg_Global_KD_vs_WT_ribosomal_2_clust.pdf',
       plot = p, width = 6, height = 6, dpi = 300)




########################################################
## coaccessible genes in each transition (atac based) ##
########################################################

## we dont use this for paper

da_peaks_sig_genes_trans <- readRDS("../Input/toxo_cdc/rds_ME49_59/da_peaks_sig_genes.rds")
da_peaks_sig_genes_trans <- da_peaks_sig_genes_trans %>% dplyr::select(V11, cluster) %>% distinct()
names(da_peaks_sig_genes_trans) <- c("GeneID", "phase")
da_peaks_sig_genes_trans$GeneID <- gsub("_", "-", da_peaks_sig_genes_trans$GeneID)

prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1"))
prod.desc <- prod.desc %>% 
  transmute(GeneID = gsub("_", "-", TGME49), ProductDescription = ProductDescription) %>% na.omit()

da_peaks_sig_genes_trans <- left_join(da_peaks_sig_genes_trans, prod.desc, by = "GeneID")
colnames(da_peaks_sig_genes_trans)[1] <- "gene_name"

da_peaks_sig_genes_trans.list <- split(da_peaks_sig_genes_trans, da_peaks_sig_genes_trans$phase) 

da_peaks_sig_genes_trans.list <- da_peaks_sig_genes_trans.list[-1]
trans.list <- c()
trans.list <- lapply(1:length(da_peaks_sig_genes_trans.list), function(i) {
  
  my.df <- da_peaks_sig_genes_trans.list[[i]]
  df <- clust.df(my.df, num.clust = 3)
  df$group <- names(da_peaks_sig_genes_trans.list)[i]
  #df <- df[df$data == "scRNA",]
  p1 <- plot_atac_trends(df) +
    ggtitle(names(da_peaks_sig_genes_trans.list[i]))
  p1
  #tt <- list(df, p1)
})

names(trans.list) <- names(da_peaks_sig_genes_trans.list)
saveRDS(trans.list, "../Input/toxo_cdc/rds_ME49_59/atac_markers_atac_transitions_dtw_clust_list_3L.rds")

pp <- grid.arrange(grobs = trans.list, ncol = 2)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/atac_markers_atac_transitions_dtw_clust_list3L.pdf", 
        plot = pp,
        height = 12,width = 14, dpi = 300)


trans.clust.df <- lapply(trans.list, "[[", 1) 
trans.clust.df.all <- do.call("rbind", trans.clust.df)
trans.clust.df.all <- trans.clust.df.all %>% dplyr::select(GeneID, cluster.ATAC, group) %>% distinct()
trans.clust.df.all.sum <- trans.clust.df.all %>% 
  group_by(group, cluster.ATAC) %>% summarise(genes = list(gsub("-", "_", GeneID)), total = n())

write.xlsx(trans.clust.df.all.sum, "../Output/toxo_cdc/ME49_59/tables/atac_peaks_sig_markers_atac_trans_clust_based_sum_3L.xlsx") 



######################################################## 
## genes in the overlap of transitioning in rna and atac 
########################################################

## we dont use this for paper

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(TGGT1_ME49, prod.desc, by = c("TGGT1" = "GeneID"))

DD <- read.xlsx("../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp_lng.xlsx")
DD <- left_join(DD, prod.desc , by = c("GeneID" = "TGME49"))
colnames(DD) <- gsub("gene", "gene_name", colnames(DD))
DD.list <- split(DD, f = DD$cluster.x)

trans.list <- c()
trans.list <- lapply(1:length(DD.list), function(i) {
  
  my.df <- DD.list[[i]]
  df <- clust.df(my.df, num.clust = 4)
  df$group <- names(DD.list)[i]
  #df <- df[df$data == "scRNA",]
  p1 <- plot_rna_atac_trends(df) +
    ggtitle(names(DD.list[i]))
  p1
  #tt <- list(df, p1)
})

names(trans.list) <- names(DD.list)
#saveRDS(trans.list, "../Input/toxo_cdc/rds_ME49_59/sc_rna_atac_mus_trans_overlap_df_plot.rds")
saveRDS(trans.list, "../Input/toxo_cdc/rds_ME49_59/sc_rna_atac_mus_trans_matched_genes_df_plot.rds")

grid.arrange(grobs = trans.list, ncol = 2)

pdf("../Output/toxo_cdc/ME49_59/figures/sc_rna_atac_mus_trans_overlap_clust_all.pdf")
trans.list
dev.off()
 
pp <- grid.arrange(grobs = trans.list, ncol = 2)
ggsave(plot = pp , "../Output/toxo_cdc/ME49_59/figures_paper/dtw_clust_matched_rna_atac_transition.pdf",
    width = 12, height = 14, dpi = 300)

############################################################################################
## rna transition markers with ordered clusters based on the peak (visually inspected)
############################################################################################

plot_rna_atac_trends.ord <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw() +
    #theme_bw(base_size = 16) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 22, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 18, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 22, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA.ordered ~ data, scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=20, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=22, face="bold", hjust = 1),
      axis.title.y = element_text(size=22, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}


rna.trans.marker.genes.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list.rds")
rna.trans.data.list <- lapply(rna.trans.marker.genes.list, "[[", 1)
rna.trans.data <- do.call("rbind", rna.trans.data.list)

rna.trans.data <- rna.trans.data %>% 
  mutate(new.ord.clust = 
           case_when(group == "T1" & cluster.RNA == "C 1" ~ "D4",
                     group == "T1" & cluster.RNA == "C 2" ~ "D2", 
                     group == "T1" & cluster.RNA == "C 3" ~ "D3", 
                     group == "T1" & cluster.RNA == "C 4" ~ "D1",
                     group == "T2" & cluster.RNA == "C 1" ~ "D4",
                     group == "T2" & cluster.RNA == "C 2" ~ "D2", 
                     group == "T2" & cluster.RNA == "C 3" ~ "D1", 
                     group == "T2" & cluster.RNA == "C 4" ~ "D3", 
                     group == "T3" & cluster.RNA == "C 1" ~ "D4",
                     group == "T3" & cluster.RNA == "C 2" ~ "D3",
                     group == "T3" & cluster.RNA == "C 3" ~ "D2", 
                     group == "T3" & cluster.RNA == "C 4" ~ "D1", 
                     group == "T4" & cluster.RNA == "C 1" ~ "D3", 
                     group == "T4" & cluster.RNA == "C 2" ~ "D4", 
                     group == "T4" & cluster.RNA == "C 3" ~ "D1", 
                     group == "T4" & cluster.RNA == "C 4" ~ "D2", 
                     TRUE ~ "NA"))
rna.trans.data$cluster.RNA.ordered <- gsub("\\D", "C", rna.trans.data$new.ord.clust) 
rna.trans.data$data <- factor(rna.trans.data$data, levels = c("scRNA", "scATAC"))
rna.trans.data.list <- split(rna.trans.data, f= rna.trans.data$group)
#saveRDS(rna.trans.data.list, "../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")
rna.trans.data.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list_ordered.rds")
trans.plt <- c()
trans.plt <- lapply(1:length(rna.trans.data.list), function(i) {
  
  my.df <- rna.trans.data.list[[i]]
  
  p1 <- plot_rna_atac_trends.ord(my.df) 
  # +
  #   ggtitle(names(rna.trans.marker.genes.list[i]))
  p1
  #tt <- list(df, p1)
})

pp <- grid.arrange(grobs = trans.plt, ncol = 2)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/rna_markers_rna_transitions_dtw_clust_list_ordered.pdf", 
        plot = pp,
        height = 16,width = 16, dpi = 300)


p1 <-trans.plt[[1]]
p1
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T1_rna_clust_ordered.pdf",
        plot = p1,height = 8,width = 10, dpi = 300)
p2 <-trans.plt[[2]]
p2
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T2_rna_clust_ordered.pdf",
        plot = p2,height = 8,width = 10, dpi = 300)

p3 <-trans.plt[[3]]
p3
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T3_rna_clust_ordered.pdf",
        plot = p3,height = 8,width = 10, dpi = 300)

p4 <-trans.plt[[4]]
p4
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T4_rna_clust_ordered.pdf",
        plot = p4,height = 8,width = 10, dpi = 300)

#################################################################################
