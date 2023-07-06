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
DEG.sig <- readRDS( '../Input_KZ/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')

## New
gene.fam <- read.xlsx("../Input/toxo_cdc/gene_families/gene_fam_KZ.xlsx")
gene.fam <- left_join(gene.fam, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

# scRNA and scATAC
rna_sub <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input_KZ/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes_1.1.rds')
sc.atac.spline.fits <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes_1.1.rds')

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


## Direct targets of AP2XII8
targs <- read.xlsx('../Input_KZ/tab/AP2XII8_cut_run_KD.xlsx')
#direct.targs <- targs
direct.targs <- targs %>% dplyr::select(TGME49, intersection_CutRun_dataSets, Global_KD_vs_WT, motif,
                                        KD_vs_WT_phase_based, ProductDescription, new.name, Category)
direct.targs <- direct.targs %>% dplyr::filter(intersection_CutRun_dataSets == 'yes' & !is.na(KD_vs_WT_phase_based))
hyperlopit <- read.xlsx('../Input_KZ/tab/hyperLopit_KZ.xlsx')
direct.targs <- left_join(direct.targs, hyperlopit, by = c('TGME49' = 'GeneID'))
direct.targs <- direct.targs %>% distinct()

CRISPR <- read.xlsx('../Input/toxo_genomics/gene_function/CRIESPAR_SCORE.xlsx', sheet = 4)
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
CRISPR <- left_join(CRISPR, TGGT1_ME49, by = c('ID' = 'TGGT1'))
direct.targs <- left_join(direct.targs, CRISPR, by = 'TGME49')
direct.targs <- direct.targs %>% distinct()

## Clustering direct targets
direct.targs$TGME49 <- gsub('_', '-', direct.targs$TGME49)


## Merge repeated rows (motif is repeated)
direct.targs <- direct.targs %>% group_by(TGME49) %>% mutate(motifs = paste(motif, collapse = ','))
direct.targs <- direct.targs %>% dplyr::select(-motif) %>% distinct()

direct.targs <- direct.targs %>% dplyr::filter(Category == 'ribosomal')
sc.rna.tab <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% unique(direct.targs$TGME49) ]
sc.atac.tab <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% unique(direct.targs$TGME49) ]

sc.rna.tab.markers.hc_dtw <- dtwClustCurves(sc.rna.tab, nclust = 4)
sc.atac.tab.markers.hc_dtw <- dtwClustCurves(sc.atac.tab, nclust = 4)

plot(sc.rna.tab.markers.hc_dtw, type = 'sc')
plot(sc.atac.tab.markers.hc_dtw, type = 'sc')

sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.tab), cluster = cutree(sc.rna.tab.markers.hc_dtw, k = 2))
sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.tab), cluster = cutree(sc.atac.tab.markers.hc_dtw, k = 2))

direct.targs$RNA.clust <- paste('C', sc.rna.clust.info$cluster[match(direct.targs$TGME49, sc.rna.clust.info$GeneID)])
direct.targs$ATAC.clust <- paste('C', sc.atac.clust.info$cluster[match(direct.targs$TGME49, sc.atac.clust.info$GeneID)])

#write.xlsx(direct.targs, '../Output_KZ/tables/AP2XII8_direct_targets_KZ_4clust.xlsx')
write.xlsx(direct.targs, '../Output_KZ/tables/AP2XII8_ribosomal_targets_KZ_2clust.xlsx')

## Plotting curves
sc.rna.long.tab <- inner_join(sc.rna.mu.scale, direct.targs, by = c('GeneID' = 'TGME49')) 
sc.rna.long.tab <- sc.rna.long.tab %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = new.name, 
            RNA.clust = RNA.clust, CRISPR = mean.phenotype, direction = KD_vs_WT_phase_based)  %>% distinct()

sc.atac.long.tab <- inner_join(sc.atac.mu.scale, direct.targs, by = c('GeneID' = 'TGME49')) 
sc.atac.long.tab <- sc.atac.long.tab %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = new.name,
            ATAC.clust = ATAC.clust, CRISPR = mean.phenotype, direction = KD_vs_WT_phase_based)  %>% distinct()




sc.rna.long.tab$Name[is.na(sc.rna.long.tab$Name)] <- sc.rna.long.tab$GeneID[is.na(sc.rna.long.tab$Name)]
sc.atac.long.tab$Name[is.na(sc.atac.long.tab$Name)] <- sc.atac.long.tab$GeneID[is.na(sc.atac.long.tab$Name)]


sc.rna.long.tab$label <- NA
sc.rna.long.tab$label[which(sc.rna.long.tab$time == 3)] <-
  sc.rna.long.tab$Name[which(sc.rna.long.tab$time == 3)]


sc.atac.long.tab$label <- NA
sc.atac.long.tab$label[which(sc.atac.long.tab$time == 3)] <-
  sc.rna.long.tab$Name[which(sc.atac.long.tab$time == 3)]





sc.rna.sc.atac.joint <- inner_join(sc.rna.long.tab, sc.atac.long.tab, 
                                   by = c("time", "GeneID", "Name", "CRISPR", 'label', 'direction'))
colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA",
                                    "CRISPR", "direction", "label", "scATAC", "cluster.ATAC")
sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>%
  pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC", "label", "CRISPR", "direction"),
               names_to = 'data', values_to = 'normExpr')


#ribosomals <- direct.targs$TGME49[direct.targs$Category == 'ribosomal']
#sc.rna.sc.atac.joint.long.ribo <- sc.rna.sc.atac.joint.long %>% dplyr::filter(GeneID %in% ribosomals)

p  <- ggplot(sc.rna.sc.atac.joint.long, aes(x= time,y=normExpr)) +
  geom_path(aes(color = GeneID),alpha = 0.8, linewidth = 0.8)+ 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('normExpr') + xlab('Time') +
  #scale_color_gradient2(low = "blue", midpoint = -0.8, mid = "gray50", high = "red") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  
  coord_cartesian(xlim = c(0,6.5)) +
  # geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
  #                 box.padding = unit(0.6, "lines"),
  #                 max.overlaps = 300,
  #                 #segment.angle = 180,
  #                 nudge_x = 0.25,
  #                 nudge_y = 0.25,
  #                 hjust=0.25,
  #                 #nudge_x=0.25,
  #                 segment.size = 0.1,
  #                 na.rm = TRUE)+
  facet_grid(cluster.ATAC~data, scales = 'free') +
  
  
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'None',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


p

#'../Output_KZ/figures/AP2XII_8_direct_targets_profile_RNA_6clust.pdf'
ggsave(filename='../Output_KZ/figures/ribosomal_2_ATAC_clust.pdf',
       plot=p,
       width = 10, height = 8,
       units = "in")
