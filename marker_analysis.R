library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(openxlsx)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))


atac_sub <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')
rna_sub <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')

### Differential gene expression
Idents(rna_sub) <- 'phase'
DefaultAssay(rna_sub) <- 'RNA'
Intra.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)

Intra.markers$GeneID <- gsub('-', '_', Intra.markers$gene)
Intra.markers.top <- Intra.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = rna_sub, 
            features = Intra.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")


Intra.markers.sig <- Intra.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())

saveRDS(Intra.markers.sig, '../Input_KZ//toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')


## SI figure 1?
phase.cols <- c('#b7222a', '#ee7202', '#caae05', '#70883a', '#b138ee')
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG, fill = cluster)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=8, fontface="bold")+
  theme_minimal() + 
  scale_fill_manual(values = phase.cols) + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold", color= "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(legend.position = 'None',
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) 

p <- plot(p) + ggtitle("scRNA-seq markers - phase-based")

p
ggsave(filename="../Output_KZ/figures/intra_DEGs_phases.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

sum(ss$num.DEG)
###################################################
## marker analysis using inferred transition points 
###################################################


rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_rna_atac_trnasition_v2.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.intra_atac_atac_trnasition_v2.rds')


## Fig 5A
## rna 
## Differential gene expression usin rna transition
## rna markers using rna transition points

Idents(rna_sub) <- 'transition.rna'
rna.markers.rna.trans <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0.0)
rna.markers.rna.trans$GeneID <- gsub('-', '_', rna.markers.rna.trans$gene)

rna.sig.markers.rna.trans <- rna.markers.rna.trans %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(rna.sig.markers.rna.trans)

saveRDS(rna.sig.markers.rna.trans, '../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
write.xlsx(rna.sig.markers.rna.trans, '../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')

ss.rna <- rna.sig.markers.rna.trans %>% group_by(cluster) %>% summarise(num.DEG = n())
sum(ss.rna$num.DEG)


ss.rna$cluster <- factor(ss.rna$cluster, levels = c('T1', 'T2', 'T3', 'T4'))
p <- ggplot(data=ss.rna, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=8, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold", color= "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) 

p <- plot(p) + ggtitle("scRNA-seq markers - rna transition")

ggsave(filename="../Output/toxo_cdc/ME49_59/figures/intra_DEGs_rna_transition.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## table for performing GO term 
rna.sig.markers.rna.trans <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')

rna.sig.markers.rna.trans.sum <- rna.sig.markers.rna.trans %>%
  group_by(cluster) %>% summarise(genes = list(GeneID), total = n())
write.xlsx(rna.sig.markers.rna.trans.sum, "../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_sum_v2.xlsx")



## Differential gene expression usin atac transition
## rna markers using atac transition points

Idents(rna_sub) <- 'transition.atac'
transition.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)
transition.markers$GeneID <- gsub('-', '_', transition.markers$gene)
transition.markers.sig <- transition.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(transition.markers.sig)

saveRDS(transition.markers.sig, '../Input/toxo_cdc/rds_ME49_59/rna_markers_sig_atac_trans_v2.rds')
write.xlsx(transition.markers.sig, '../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')

ss.atac <- transition.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss.atac$cluster <- factor(ss.atac$cluster, levels = c('T1', 'T2', 'T3', 'T4'))
ss.atac

p <- ggplot(data=ss.atac, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=8, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold",color = "black" )) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold", color = "black")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold", hjust = 1),
    axis.title.y = element_text(size=18, face="bold")
  ) 

p <- plot(p) + ggtitle("scRNA-seq markers - atac transitions ")
p
ggsave(filename="../Output/toxo_cdc/ME49_59/figures/tg_Intra_deg_numbers_atac_transition.pdf",
       plot=p,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## ATAC
## atac markers using atac transitions  (this is based on gene activity)
## we do not use this analysis 
Idents(atac_sub) <- 'transition.atac'
DefaultAssay(atac_sub) <- "RNA"

atac.markers.atac.trans <- FindAllMarkers(object = atac_sub, only.pos = T, min.pct = 0.0)
atac.markers.atac.trans$GeneID <- gsub('-', '_', atac.markers.atac.trans$gene)

atac.sig.markers.atac.trans <- atac.markers.atac.trans %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)
dim(atac.sig.markers.atac.trans)

# saveRDS(atac.sig.markers.atac.trans, '../Input/toxo_cdc/rds_ME49_59/atac_markers_atac_trns_sig_v2.rds')
# write.xlsx(atac.sig.markers.atac.trans, '../Output/toxo_cdc/ME49_59/tables/atac_markers_sig_atac_trns_v2.xlsx')
# 
# ss.atac <- atac.sig.markers.atac.trans %>% group_by(cluster) %>% summarise(num.DEG = n())
# sum(ss.atac$num.DEG)


### intersect the rna marker genes (rna and atac transitions were used to perform DEG)
## we do not use this analysis 
rna.sig.markers.rna.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')
rna.sig.markers.atac.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')

rna.markers.rna.trans.sum <- rna.sig.markers.rna.trans %>% 
  group_by(cluster) %>% summarise( genes = list(GeneID), total = n())

rna.markers.atac.trans.sum <- rna.sig.markers.atac.trans %>%
  group_by(cluster) %>% summarise(genes = list(GeneID), total = n())

XX <- left_join(rna.markers.rna.trans.sum, rna.markers.atac.trans.sum, by = "cluster")
XX <- XX %>% rowwise() %>%
  mutate(overlap = list(intersect(genes.x, genes.y)), 
         overlap.num = length(intersect(genes.x, genes.y)),
         diff = list(setdiff(genes.x, genes.y)), 
         diff.num = length(setdiff(genes.x, genes.y)))

#write.xlsx(XX, "../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp.xlsx")


## contingency 

rna.sig.markers.rna.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_trns_v2.xlsx')
rna.sig.markers.rna.trans$data <- "rna.transition"
#rna.sig.markers.rna.trans <- rna.sig.markers.rna.trans %>% dplyr::select(GeneID, cluster, data) 
rna.sig.markers.rna.trans$cluster.trans <- paste(rna.sig.markers.rna.trans$cluster, "rna", sep = ".")

rna.sig.markers.atac.trans <- read.xlsx('../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_atac_trns_v2.xlsx')
rna.sig.markers.atac.trans$data <- "atac.transition"
#rna.sig.markers.atac.trans <- rna.sig.markers.atac.trans %>% dplyr::select(GeneID, cluster , data)
rna.sig.markers.atac.trans$cluster.trans <- paste(rna.sig.markers.atac.trans$cluster, "atac", sep = ".")

DD <- left_join(rna.sig.markers.rna.trans, rna.sig.markers.atac.trans, by = "GeneID" ,  multiple = "all") %>% na.omit()
DD.tab <- table(DD$cluster.trans.x, DD$cluster.trans.y)

DD.overlap <- DD %>% dplyr::filter(cluster.x == cluster.y)
write.xlsx(DD.overlap, "../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp_lng.xlsx")
write.xlsx(DD.tab,"../Output/toxo_cdc/ME49_59/tables/rna_markers_sig_rna_atac_trans_ovlp_contingency.xlsx" )

matched.genes <- DD.tab %>% data.frame()
p <- ggplot(matched.genes, aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile(color = "white",
            lwd = 1,
            linetype = 1) + 
  geom_text(aes(label = Freq), color = "white", size = 8, fontface = "bold") +
  theme_bw() +
  scale_fill_continuous(high = "#132B43", low = "#56B1F7") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold.italic", color = "black")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold.italic", color = "black")) +
  theme(legend.text = element_text(size = 14, face = "bold"),
        legend.title = element_text(size = 18, face = "bold")) +
  theme(
    plot.title = element_text(size=20, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=20, face="bold", hjust = 0.5),
    axis.title.y = element_text(size=20, face="bold", hjust = 0.5)
  ) +
  ylab("") + xlab("")+
  theme(legend.position = "none")+
  guides(color = guide_legend(override.aes = list(size = 7)))

p
ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/matched_DEGs_rna_atac_transitions.pdf", 
        plot = p, height = 8, width = 8, dpi = 300)

