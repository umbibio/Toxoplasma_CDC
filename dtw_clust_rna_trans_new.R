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

clust.atac.df <- function(tab, num.clust = num.clust){
  
  tab <- tab
  k <- num.clust
  
  sc.atac <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% tab$gene_name]
  sc.atac.markers.hc_dtw <- dtwClustCurves(sc.atac, nclust = k)
  tab <- tab[tab$gene_name %in% colnames(sc.atac),]
  sc.atac.long <- inner_join(sc.atac.mu.scale, tab, by = c('GeneID' = 'gene_name'))
  sc.atac.long <- sc.atac.long %>%
    transmute(time = x, GeneID = GeneID, normExpr = expr, group = group, trans.cluster.rna) %>% 
    distinct()
  sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac), cluster = cutree(sc.atac.markers.hc_dtw, k = k))
  sc.atac.long.clust <- inner_join(sc.atac.long, sc.atac.clust.info, by = 'GeneID')
  colnames(sc.atac.long.clust) <- c("time", "GeneID", "normExpr", "group", "trans.cluster.rna", "cluster.ATAC")
  sc.atac.long.clust$cluster.ATAC <- paste('C', sc.atac.long.clust$cluster.ATAC)
  
  return(sc.atac.long.clust)
}

# facet by rna cluster
plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1,face="bold", size = 15, colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, face="bold",size = 15, colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.RNA ~ data , scales = 'free', space = 'free') +
    #facet_grid(data ~ cluster.RNA, scales = 'free', space = 'free') +
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


plot_atac_trand <- function(sc.atac.long.clust){
  
  p  <- ggplot(sc.atac.long.clust, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, face="bold", colour = "black")) +
    theme(strip.background = element_rect(colour="black", fill="white",size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 18, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    facet_grid(cluster.ATAC ~ ., scales = 'free', space = 'free') +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=18, face="bold", hjust = 1),
      axis.title.y = element_text(size=18, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold")) +
    theme(axis.ticks =  element_blank())
  
  return(p)
}


source('./util_funcs.R')

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

## cluster markers of each rna-transition (T1, T2, T3, T4)
trans.list <- c()
k <- 4

trans.list <- lapply(1:length(rna.trans.marker.genes.list), function(i) {
  my.df <- rna.trans.marker.genes.list[[i]]
  df <- clust.df(my.df, num.clust = k)
  df$group <- names(rna.trans.marker.genes.list)[i]
  df$data <- factor(df$data, levels = c("scRNA", "scATAC"))
  p1 <- plot_rna_atac_trends(df) +
    ggtitle(names(rna.trans.marker.genes.list[i]))
  p1
})

names(trans.list) <- names(rna.trans.marker.genes.list)
saveRDS(trans.list,"../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_3_clust_list.rds" )

trans.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list.rds")
names(trans.list) <- c("T1", "T2", "T3", "T4")

pp <- grid.arrange(grobs = trans.list, ncol = 2)
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/rna_markers_rna_transitions_dtw_4_clust_list.pdf",
        plot = pp,height = 18,width = 16, dpi = 300)




trans.clust.df <- lapply(trans.list, "[[", 1) 
trans.clust.df.all <- do.call("rbind", trans.clust.df)
trans.clust.df.all <- trans.clust.df.all %>% dplyr::select(GeneID, cluster.RNA, group) %>% distinct()
trans.clust.df.all$trans.cluster.rna <- paste(trans.clust.df.all$group, trans.clust.df.all$cluster.RNA, sep = "_")
trans.clust.df.all$trans.cluster.rna <- gsub(" ", "", trans.clust.df.all$trans.cluster.rna)
colnames(trans.clust.df.all) <- gsub("GeneID", "gene_name", colnames(trans.clust.df.all))
trans.clust.df.all.sum <- trans.clust.df.all %>% 
  group_by(trans.cluster.rna) %>% summarise(genes = list(gsub("-", "_", gene_name)), total = n())


# 
# plot_rna_atac_trends(trans.list$T1$data)
# T1.C1 <- trans.clust.df.all %>% filter(trans.cluster.rna == "T1_C1")
# colnames(T1.C1) <- gsub("GeneID", "gene_name", colnames(T1.C1))
# 


trans.clust.rna.list  <- split(trans.clust.df.all, f = trans.clust.df.all$trans.cluster.rna)
atac.list <- c()
#i <- 2
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
})

names(atac.list) <- names(trans.clust.rna.list)

saveRDS(atac.list, "../Input/toxo_cdc/rds_ME49_59/atac_dtw_clust.rds")
atac.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/atac_dtw_clust.rds")

pp2 <- grid.arrange(grobs = atac.list, ncol = 4)
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/atac_clust_V2.pdf",
        plot = pp2,height = 13,width = 24, dpi = 300)


atac <- lapply(atac.list, "[[", 1)
atac.df <- do.call("rbind", atac)
atac.df$gene_name <- gsub("-", "_", atac.df$GeneID)
atac.df <- atac.df %>% select(gene_name, trans.rna.atac.clust) %>% 
  group_by(trans.rna.atac.clust) %>% distinct() %>% mutate(total = n())
atac.df.list <- split(atac.df, f = atac.df$trans.rna.atac.clust)

out.dir <- "../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters/"
lapply(1:length(atac.df.list), function(i){
  
  tab <- atac.df.list[i]
 
  name.file<- paste(names(atac.df.list)[i], ".xlsx", sep = "")
  
  write.xlsx(tab, paste(out.dir, name.file))
  
})


# atac.summ <- atac.df %>% select(gene_name,trans.rna.atac.clust ) %>% distinct() %>%
#   group_by(trans.rna.atac.clust) %>% summarise(genes = list(unique(gene_name)), total = n())

# write.xlsx(atac.summ, "../Output/toxo_cdc/ME49_59/tables/atac_clusters_within_rna_tran_dtw_clusters.xlsx")
# tmp <- df %>% dplyr::select(GeneID, cluster.ATAC, ProductDescription) %>% distinct()
# tmp$GeneID <- gsub("-", "_", tmp$GeneID)
# # write.xlsx(tmp, "../Output/toxo_cdc/tabels/atac_T1C1.xlsx")


trans.list <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_transitions_dtw_clust_list.rds")
names(trans.list) <- c("T1", "T2", "T3", "T4")

# pp <- grid.arrange(grobs = trans.list[1], ncol = 1)
# ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/rna_markers_rna_transitions_dtw_4_clust_list.pdf",
#         plot = pp,height = 18,width = 16, dpi = 300)


i <- 1
df <- trans.list[[i]][["data"]]
df$group <- names(trans.list)[i]
df$data <- factor(df$data, levels = c("scRNA", "scATAC"))
df$cluster.RNA <- gsub(" ", "", df$cluster.RNA)
#df$cluster.RNA <- factor(df$cluster.RNA, levels = c("C4", "C2", "C3", "C1")) # T1
#df$cluster.RNA <- factor(df$cluster.RNA, levels = c("C3", "C2", "C4", "C1")) # T2
#df$cluster.RNA <- factor(df$cluster.RNA, levels = c("C4", "C3", "C2", "C1")) # T3
df$cluster.RNA <- factor(df$cluster.RNA, levels = c("C3", "C4", "C1", "C2")) # T3
p1 <- plot_rna_atac_trends(df) +
  ggtitle(names(rna.trans.marker.genes.list[i]))
p1
  


ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/T4_rna_clust_ordered.pdf",
        plot = p1,height = 8,width = 10, dpi = 300)
