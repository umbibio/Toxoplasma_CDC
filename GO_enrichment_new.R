
## GO term 
## rna transition markers GO term - Slide 5 and 8 of toxo_cdc_figures_V7

in.dir <- '../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_sum_GO_toxodb_OutPut/'
in.dir <- '../Output/toxo_cdc/ME49_59/tables/rna_sig_markers_rna_trans_clust_based_sum/'

## Cut&run and DEGs upon AP2XII8 KD
in.dir <- '../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks/GO/' ## phase based
in.dir <- '../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks_globalDEG/GO/'

all.files <- list.files(in.dir)

all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  cluster <- strsplit(nn, split = '_')[[1]][2]
  GF <- strsplit(nn, split = '_')[[1]][1]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$cluster <- cluster
  all.clust.items <- c(all.clust.items, list(tmp))
}


all.clust.items <- do.call(rbind, all.clust.items)


write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_GO_term.xlsx')
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_clust_based_GO_term.xlsx')

## 
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_phase_based_DEG.xlsx')
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG.xlsx')



filtered.Go <- all.clust.items
colnames(all.clust.items) <- gsub("P-value", "pval", colnames(all.clust.items))


filtered.Go <- all.clust.items %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 30) %>% arrange(cluster, Benjamini)



filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))
filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=rev(unique(filtered.Go$Name)))
colnames(filtered.Go) <- gsub(" ", "_", colnames(filtered.Go))



write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_GO_term_filt.xlsx')
write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/rna_sig_markers_rna_trans_clust_based_GO_term_filt.xlsx')

write.xlsx(filtered.Go, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_phase_based_DEG_filt.xlsx')
write.xlsx(all.clust.items, '../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG_filt.xlsx')

# filtered.Go.long <- filtered.Go %>%
#   mutate(gene_name = strsplit(Result_gene_list, ",")) %>%
#   unnest(gene_name)
########################
## transition color ####
########################

filtered.Go$cluster <- factor(filtered.Go$cluster, levels = sort(unique(filtered.Go$cluster)))

filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=rev(unique(filtered.Go$Name)))
colnames(filtered.Go) <- gsub(" ", "_", colnames(filtered.Go))

p <- ggplot(filtered.Go, aes(x = cluster, y = Name,fill = cluster, color = cluster)) + 
  geom_point(aes( size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  scale_fill_manual(values = c("T1" = "#ff9a00", 'T2' = '#9ca820', 'T3' = '#615FB1', 'T4' = '#8f139f')) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold", colour = "black")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10)) + 
  xlab("cluster") + ylab("GO term") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"))

plot(p)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_sig_markers_rna_trans_GO_term_top.pdf",
       plot=p ,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)
##############################
## transition cluster based ##
##############################
filtered.Go <- all.clust.items 
# filtered.Go$cluster <- factor(filtered.Go$cluster, 
#                               levels = c("T1C4", "T1C2", "T1C3","T1C1", 
#                                          "T2C3", "T2C2", "T2C4", "T2C1",
#                                          "T3C4", "T3C3", "T3C2", "T3C1", 
#                                          "T4C3", "T4C4", "T4C1", "T4C2"))

filtered.Go <- filtered.Go %>% arrange(cluster, Benjamini) %>% distinct() %>%
  group_by(cluster) %>% mutate(rank = row_number()) %>%
  #dplyr::filter(pval < 0.05 )
  dplyr::filter(Benjamini < 0.1 & rank < 5) %>% arrange(cluster, Benjamini)

filtered.Go$ID <- factor(filtered.Go$ID, level=unique(filtered.Go$ID))
filtered.Go$Name <- factor(filtered.Go$Name, level=rev(unique(filtered.Go$Name)))
colnames(filtered.Go) <- gsub(" ", "_", colnames(filtered.Go))

filtered.Go <- filtered.Go %>%  arrange(cluster, Benjamini)
# transition color clust based
p <- ggplot(filtered.Go, aes(x = cluster, y = Name,fill = cluster, color = cluster)) + 
  geom_point(aes( size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("T1C1" = "#ff9a00", 'T1C2' = '#ff9a00', 'T1C3' = '#ff9a00', 'T1C4' = '#ff9a00',
                                'T2C1' = '#9ca820', 'T2C2' = '#9ca820', 'T2C3' = '#9ca820', 'T2C4' = '#9ca820', 
                                'T3C1' = '#615FB1', 'T3C2' = '#615FB1', 'T3C3' = '#615FB1', 'T3C4' = '#615FB1',
                                'T4C1' = '#8f139f', 'T4C2' = '#8f139f', 'T4C3' = '#8f139f', 'T4C4' = '#8f139f')) +
  scale_fill_manual(values = c("T1C1" = "#ff9a00", 'T1C2' = '#ff9a00', 'T1C3' = '#ff9a00', 'T1C4' = '#ff9a00',
                               'T2C1' = '#9ca820', 'T2C2' = '#9ca820', 'T2C3' = '#9ca820', 'T2C4' = '#9ca820', 
                               'T3C1' = '#615FB1', 'T3C2' = '#615FB1', 'T3C3' = '#615FB1', 'T3C4' = '#615FB1',
                               'T4C1' = '#8f139f', 'T4C2' = '#8f139f', 'T4C3' = '#8f139f', 'T4C4' = '#8f139f')) +
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold", colour = "black")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 9, face="bold", colour = "black")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  theme(plot.title = element_text(size = 10)) + 
  xlab("cluster") + ylab("GO term") +
  theme(axis.title = element_text(size = 14, colour = "black", face = "bold"))

plot(p)
ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/rna_sig_markers_rna_trans_clust_based_GO_term.pdf",
       plot=p ,
       width = 10, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




#################################################
### cut&run and global/phase based DEGs #########
#################################################

p <- ggplot(filtered.Go, aes(x = cluster, y = Name)) + 
  geom_point(aes(colour = cluster, size =  -log(Benjamini))) +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c')) +
  scale_fill_manual(values = c("downReg" = "#6565bf","upReg" ='#ee5d6c'))+
  #scale_colour_gradient(limits=c(0, 0.01), low="red") +
  ylab(NULL) + xlab(NULL) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face="bold")) + 
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 10, face="bold")) +
  #theme(legend.position="none") +
  theme(strip.background = element_rect(colour="black", fill="white", 
                                        size=0.5, linetype="solid")) + guides(color = FALSE)+
  #ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
  ggtitle("High Conf CutRun peaks Global DEGs in KD_vs_WT")+
  theme(plot.title = element_text(size = 9))

plot(p)

ggsave(filename="../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_phase_based_DEG_GO_term.pdf",
       plot=p ,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


ggsave(filename="../Output/toxo_cdc/ME49_59/tables/all_GO/HighConfPeaks_globalDEG_GO_term.pdf",
       plot=p ,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


