

## 1. peaks in the intersection Overlap with up and down regulated genes in WT vs KD 
tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
names(tab)[1] <- "chr"


## DEGs KD vs WT phase based 

# phase based
KD.vs.WT.phase <- tab %>% dplyr::select(gene_name, KD_vs_WT_phase_based) %>% 
  filter(KD_vs_WT_phase_based != "NA") %>% distinct()

# global 
KD.vs.WT.glob <- tab %>% dplyr::select(gene_name, Global_KD_vs_WT) %>% 
  filter(Global_KD_vs_WT != "NA") %>% distinct()


# cutRun genes in intersection of 4 data
intrsct.peaks  <- tab %>% dplyr::select(gene_name, intersection) %>% 
  distinct() %>% filter(intersection == "yes")


DEGs <- KD.vs.WT.glob
DEGs <- KD.vs.WT.phase

venn.list <- list(KD.vs.WT = unique(DEGs$gene_name), 
                  High.Conf.Peaks = intrsct.peaks$gene_name)

p <- ggVennDiagram(venn.list, set_size = 6, label_size = 8) 
p


ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_Conf_CutRun_peaks_phase_based_DEGs_KD_vs_WT_Venn.pdf",
       plot = p, width = 6, height = 6, dpi = 300)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_Conf_CutRun_peaks_global_DEGs_KD_vs_WT_Venn.pdf",
       plot = p, width = 6, height = 6, dpi = 300)


##################################################

# phase based
HC.peaks <- tab %>% 
  filter(intersection == "yes" & KD_vs_WT_phase_based %in% c('down_reg', 'up_reg', 'modulated') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, KD_vs_WT_phase_based,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"

# global
HC.peaks <- tab %>% 
  filter(intersection == "yes" & Global_KD_vs_WT %in% c('down_reg', 'up_reg') ) %>% 
  dplyr::select(chr, start_peak, end_peak, V4, V5, V11,gene_name,intersection, Global_KD_vs_WT,ProductDescription.x , Category ) %>% 
  distinct()
names(HC.peaks)[9] <- "dir"


## summary and bar plot
HC.peaks.stat <- HC.peaks %>% group_by(dir, Category) %>%
  summarise(total = n())

HC.peaks.stat <- HC.peaks.stat %>% 
  mutate(Color  = ifelse(dir == "down_reg", "#8080FA", "#fc6c85"))

HC.peaks.stat <- HC.peaks.stat %>%
  mutate( ## for this you will need to remove the whitespace 
    Category = stringr::str_trim(Category),
    newcolor = ifelse(grepl("ribo", Category), alpha(Color, .5), Color)
  ) 

#
sum(HC.peaks.stat$total)

p <- ggplot(HC.peaks.stat, aes(dir, total)) +
  geom_col(aes(fill = I(newcolor), color = Category),
           position = position_stack(reverse = FALSE),
           ##remove the border 
           linewidth = 0) +
  geom_text(aes(label = total, group = Category),size = 8, color = "black",
            fontface = 'bold', position = position_stack(vjust = .5, reverse = TRUE)) +
  ## change the legend fill manually 
  guides(color = guide_legend(
    reverse = FALSE,
    override.aes = list(fill = c("grey25", alpha("grey25", .5))))) +
  #theme(legend.position = "none") +
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_text(size = 14, face = "bold", colour = "black"),
    axis.ticks = element_blank())+
  theme(strip.background=element_rect(fill='black', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
        plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 16))+
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  #theme(legend.position = "none")
p

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_Conf_CutRun_peaks_phase_based_DEGs_KD_vs_WT_ribo_stacked.pdf", 
       plot = p, width = 6, height = 4, dpi = 300)

ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_Conf_CutRun_peaks_global_DEGs_KD_vs_WT_ribo_stacked.pdf", 
       plot = p, width = 6, height = 4, dpi = 300)

## 
# p <- ggplot(data=HC.peaks.stat, aes(x=dir, y=total, fill=Category)) +
#   geom_bar(stat="identity", color="black", position=position_dodge())+
#   geom_text(aes(label=total), vjust=1.15, color="black",
#             position = position_dodge(1), size=6) +
#   theme_bw()+
#   theme(strip.background=element_rect(fill='black', color = 'black'),
#         panel.spacing = unit(1.5, "lines"), 
#         strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 14,face = 'bold'),
#         plot.title = element_text(size=16, face = "bold.italic", color = 'black'),
#         axis.title.x = element_text(size=16, face="bold"),
#         axis.title.y = element_text(size=16, face="bold"),
#         axis.text = element_text(size = 14, face = "bold"))  +
#   theme(legend.text = element_text(face = "bold", size = 14),
#         legend.title = element_text(face = "bold", size = 14))
# 
# p
# 
# ggsave("../Output/toxo_cdc/ME49_59/figures_paper/High_Conf_CutRun_peaks_phase_based_DEGs_KD_vs_WT_ribo_stacked.pdf", 
#        plot = p, width = 4, height = 4, dpi = 300)



## write fasta file (cut and run regions) for high conf peaks 
HC.peaks.list <- split(HC.peaks, f = HC.peaks$dir)
#DEG.type <- "KD_vs_WT_phase_based"
DEG.type <- "KD_vs_WT_global"

#out.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks/bedFasta/"
out.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/HighConfPeaks_globalDEG/bedFasta/"
cluster.bed.list <- lapply(1:length(HC.peaks.list), function(i){
  
  bed.name <- paste(names(HC.peaks.list)[i], paste(DEG.type, ".bed", sep = ""), sep = "_")
  fasta.name <- paste(names(HC.peaks.list)[i], paste(DEG.type, ".fasta", sep = ""), sep = "_")
  tmp <- HC.peaks.list[[i]] 
  cluster.bed <- tmp %>% select(chr, start_peak, end_peak ,V5 ,  V11, gene_name,  dir)
  
  write.table(cluster.bed, paste(out.dir, bed.name, sep = ""), sep = "\t", quote = F, row.names = F, col.names = F)
  write.xlsx(cluster.bed, paste(out.dir, gsub('.bed', '.xlsx', bed.name), sep = "") )

  cluster.fasta <- bedtoolsr::bt.getfasta(fi = "../Input/toxo_cdc/Genome/ToxoDB-59_TgondiiME49_Genome.fasta",
                                          bed = tmp,
                                          fo =paste(out.dir, fasta.name, sep = ""))
  return(cluster.bed)
})


## summary of motif 1 and motif 2 

tab <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")


tab.down.reg <- tab %>% 
  dplyr::select(gene_name, has.motif, motif, Global_KD_vs_WT, ProductDescription.x, Category) %>%
  filter(has.motif == "yes" & Global_KD_vs_WT == "down_reg" & motif %in% c("motif_1", "motif_2")) %>% distinct()

tab.down.reg <- tab.down.reg %>% group_by(gene_name) %>% mutate(motif.list = list(motif)) 
tab.down.reg <- tab.down.reg %>% rowwise() %>%
  mutate(which_motif = ifelse(length(unlist(motif.list)) > 1 , "both" , "one")) %>% distinct()

ind.one <- which(tab.down.reg$which_motif == "one")
tab.down.reg$which_motif[ind.one] <- tab.down.reg$motif[ind.one]

tab.sum <- tab.down.reg %>% group_by(which_motif, Category, Global_KD_vs_WT) %>% 
  summarise(genes = list(unique(gene_name)), total = list(length(unique(gene_name))))

tab.down.reg.uniq <- tab.down.reg %>%
  dplyr::select(gene_name, has.motif, Global_KD_vs_WT, ProductDescription.x, Category, which_motif) %>%
  distinct()
write.xlsx(tab.down.reg.uniq, "../Output/toxo_cdc/ME49_59/tables/global_KD_vs_WT_down_reg_motif1_motif2_occurence.xlsx")

m1.v2 <- tab.down.reg %>% filter(motif %in% c("motif_1"))
m2.v2 <- tab.down.reg %>% filter(motif %in% c( "motif_2"))

ven.list <- list(motif_1 = unique(m1.v2$gene_name), motif_2 = unique(m2.v2$gene_name))
p <- ggVennDiagram(ven.list, set_size = 6, label_size = 8)
p <- p + ggtitle("down_reg") +theme(
  plot.title = element_text(size=20, face = "bold.italic", color = 'red'))

p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/global_KD_vs_WT_down_reg_motif1_motif2_occurence_venn.pdf", 
       plot = p, height = 6, width = 6, dpi = 300)



## Proximity of motif 1 and motif 2

CutRunPeaksMotif <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_motif_info_NO_frag_filt.rds")
CutRunPeaksMotif <- do.call("rbind", CutRunPeaksMotif)

# CutRunPeaksMotif <- CutRunPeaksMotif %>% 
#   dplyr::select( V13, ProductDescription, has.motif, motif) %>% distinct() 

CutRunPeaksMotif <- CutRunPeaksMotif %>% 
  dplyr::select(V1, V2, V3, V13, ProductDescription, has.motif, motif) %>% distinct() 

m1 <- CutRunPeaksMotif %>% filter(motif %in% c("motif_1"))
m2 <- CutRunPeaksMotif %>% filter(motif %in% c( "motif_2"))

ven.list <- list(m1 = unique(m1$V13), m2 = unique(m2$V13))
ggVennDiagram(ven.list, set_size = 6, label_size = 8)

df <- full_join(m1, m2, by = "V13", relationship = "many-to-many")
df <- df %>% filter(has.motif.x == "yes" & has.motif.y == "yes")
df <- df %>% mutate(dist = abs(V2.y - V2.x))
df <- df %>% group_by(V13) %>% summarise(min.dist = min(dist))

hist(df$min.dist)

p <- ggplot(df, aes(x = min.dist)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 0.7,linetype = 1,colour = 2) +
  theme_bw() 
p
