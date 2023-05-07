library(tidyverse)
library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidytext)
library(seqinr)

peak.genes <- read.xlsx("../Input/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.xlsx")

sc.rna.sig.markers <- readRDS("../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds")


colnames(sc.rna.sig.markers) <- gsub("cluster", "phase", colnames(sc.rna.sig.markers))
sc.rna.sig.markers <- sc.rna.sig.markers %>% transmute(GeneID, phase) %>% distinct()
sc.rna.sig.markers$GeneID <- gsub("-", "_", sc.rna.sig.markers$GeneID)
sc.rna.sig.markers.stat <- sc.rna.sig.markers %>% group_by(phase) %>% summarise(total = n())
sc.rna.sig.markers.stat

sc.rna.sig.markers.peaks <- 
  right_join(sc.rna.sig.markers, peak.genes, by = c("GeneID" = "gene_name" )) %>% na.omit()


sc.rna.sig.markers.peaks$phase <- gsub("\\.", "_", sc.rna.sig.markers.peaks$phase)
colnames(sc.rna.sig.markers.peaks) <- gsub("V1.x", "chr" , 
                                           gsub("V11", "strand", 
                                                colnames(sc.rna.sig.markers.peaks)))
sc.rna.sig.markers.peaks <- sc.rna.sig.markers.peaks %>% 
  dplyr::select(chr, start_peak, end_peak, strand, phase, GeneID)
sc.rna.sig.markers.peaks.list <- split(sc.rna.sig.markers.peaks, 
                                       f = sc.rna.sig.markers.peaks$phase)



## write fasta file (atac regions) for rna transition markers to perform motif search 
## fasta file have been generated in terminal with strandness option enabled

options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")

out.dir <- "../Input/toxo_cdc/rds_ME49_59/bed_files/rna_transitions_v4/"

cluster.bed.list <- lapply(1:length(sc.rna.sig.markers.peaks.list), function(i){
  
  bed.name <- paste(names(sc.rna.sig.markers.peaks.list)[i], "bed", sep = ".")
  fasta.name <- paste(names(sc.rna.sig.markers.peaks.list)[i], "fasta", sep = ".")
  tmp <- sc.rna.sig.markers.peaks.list[[i]] 
  cluster.bed <- tmp %>% dplyr::transmute(V1 = chr, V2 = start_peak, V3 = end_peak,  V4 = GeneID ,V5 = ".",  V6 = strand)
  
  #write.table(cluster.bed, paste(out.dir, bed.name, sep = ""), sep = "\t", quote = F, row.names = F, col.names = F)
  # cluster.fasta <- bedtoolsr::bt.getfasta(fi = "../Input/toxo_cdc/Genome/ToxoDB-59_TgondiiME49_Genome.fasta",
  #                                         bed = tmp,
  #                                         fo =paste(out.dir, fasta.name, sep = ""))
  return(cluster.bed)
})

##### 
## read in the bed file for each motif identified by BAM
## overlap with the gene.peaks table 
## to see which gene has that specific motif 
## this version is without specifying strand, this gives the bed file corresponding to the occurence of motif 

prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1_ME49.xlsx') 
peak.genes <- read.xlsx("../Input/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.xlsx")

in.dir <- "../Input/toxo_cdc/bed_files/rna_transitions_v2/motif_beds_BAMM/"
out.dir <- "../Input/toxo_cdc/bed_files/rna_transitions_v2/motif_beds_genes/"


bed.files <- list.files(in.dir)
bed.names <- gsub(".bed", "", bed.files)

options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
bed.genes.list <- lapply(1:length(bed.files), function(i){
  
  bed.tmp <- read.table(paste(in.dir, bed.files[i], sep = "/"), header = F, quote = NULL, sep = "\t")
  ## intersect with the peak gene 
  bed.genes <- bedtoolsr::bt.intersect(a = bed.tmp, b = peak.genes, wo = T)
  bed.genes <- bed.genes %>% group_by(V13) %>% mutate(motif_freq = n())
  bed.genes <- left_join(bed.genes, prod.desc, by = c("V13" = "TGME49"))
  #bed.genes <- bed.genes %>% distinct(GeneID, .keep_all = T)
  bed.genes <- bed.genes %>% select(V1, V2, V3, V4, V5, V6, V12,V13, motif_freq, ProductDescription)
  
  #write.xlsx(bed.genes, paste(out.dir, paste(bed.names[i], "genes.xlsx", sep = "_"), sep =""))
  
  # write.table(bed.genes, paste(out.dir, paste(bed.names[i], "genes.txt", sep = "_"), sep =""),
  #             sep = "\t", quote = F, row.names = F, col.names = F)
  # 
  return(bed.genes)
  
})

names(bed.genes.list) <- bed.names

saveRDS(bed.genes.list,"../Input/toxo_cdc/rds_ME49_59/motif_dist_from_TSS.rds")

## distance of motifs relative to TSS
gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*', V1))
gtf.filt$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt$V9)))
gtf.filt$gene_name <- gsub("\"", "", gtf.filt$gene_name)

tss.df <- lapply(1:length(bed.genes.list), function(i){
  
  tmp <- bed.genes.list[[i]]
  tmp.gtf <- left_join(tmp, gtf.filt, by = c("V13" = "gene_name") )
  tmp.gtf <- tmp.gtf %>% 
    filter(V3.y == "transcript") %>% mutate(tss = ifelse((V7 == "+"), V4.y, V5.y))
  tmp.gtf <- tmp.gtf %>% mutate(tss.dist = tss - V2.x)
  tmp.gtf$group <- names(bed.genes.list)[i]
  
  
  return(tmp.gtf)
  
})

names(tss.df) <- bed.names
plt.title <- sub("^([^-]+)_([^-]+)_", "\\1-\\2.", bed.names)

tss.plt <- lapply(1:length(tss.df), function(i){
  
  p <- ggplot(tss.df[[i]], aes(x = tss.dist)) + 
    geom_histogram(aes(y = ..density..),
                   colour = 1, fill = "white") +
    geom_density(lwd = 0.7,linetype = 1,colour = 2) +
    theme_bw() + 
    ggtitle(plt.title[[i]])
  
  return(p)
})

pp <- grid.arrange(grobs = tss.plt, ncol = 3)

ggsave( "../Output/toxo_cdc/ME49_59/figures_paper/motifs_dist_tss.pdf", 
        plot = pp,
        height = 8,width = 8, dpi = 300)
