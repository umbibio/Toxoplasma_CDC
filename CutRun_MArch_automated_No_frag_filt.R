library(tidyverse)
library(bedtoolsr)
library(openxlsx)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(edgeR)
library(dtwclust)
library(geomtextpath)
library(bigmemory)
require(gridExtra)
library(grid)
library(ggVennDiagram)


source('./util_funcs.R')

#######################################################
########### Peak Gene Assignment (CUT&RUN) ############
#######################################################


# all cut and run samples old(AP2XII8) and new batch(multiple concentration of antibody)
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )

get_peak_genes_assign <- function(gtf, peaks, qval = qval){
  
  # sort cut and run peaks 
  CutRun <- peaks %>% filter(V9 > -log10(qval))
  CutRun <- CutRun %>% dplyr::select(V1, V2, V3)
  peaks.all.sort <- CutRun %>% arrange(V1, as.numeric(V2), as.numeric(V3))
  peaks.all.sort$V4 <- paste(paste(peaks.all.sort$V1, peaks.all.sort$V2, sep = ":"),peaks.all.sort$V3 ,sep = "-" )
  
  
  
  # prep gtf file
  gtf.filt <- gtf %>% dplyr::filter(!grepl('KE.*',gtf$V1))
  ## Remove the first Exon from transcripts.
  gtf.exon <- gtf.filt %>% dplyr::filter(V3 == 'exon')
  gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
  parse.str <- strsplit(gtf.exon$V9, split = ' ')
  inds <- unlist(lapply(parse.str , function(x) which(grepl("gene_id", x)) + 1))
  gtf.exon$gene_name <- gsub(";", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))
  gtf.exon$gene_name <- gsub("\"", "", gtf.exon$gene_name)
  gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                   multiple.exon = ifelse(n() > 1, T, F))
  ## Remove the exon1, but keep the Intron 1 , build exon2ton
  gtf.exon.2Ton <- gtf.exon %>% mutate(V10 = ifelse(multiple.exon & V7 == '-', min(V4), 
                                                    ifelse(multiple.exon & V7 == '+', min(V5), V4)),
                                       V11 = ifelse(multiple.exon & V7 == '-', max(V4), 
                                                    ifelse(multiple.exon & V7 == '+', max(V5), V5))) %>%
    mutate(V4 = V10, V5 = V11) %>% 
    dplyr::select(-c(exon.ord,multiple.exon, V10, V11) ) %>% distinct()
  
  
  # peak-gene assignement
  
  ## Overlap with peaks and filter peaks that are entirely within the genes.
  ## Overlapping peaks with exon2Ton and check to see if it is entirely within the gene. 
  ## Then from the sorted peaks we throw out all peaks entirely within the gene & not overlapping with exon1, 
  ## These peaks should not be assigned to any peaks.
  
  options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
  
  peak.genes.ovrlp <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.exon.2Ton, wo = T)
  peak.genes.filt <- peak.genes.ovrlp %>% dplyr::filter(V8  <= V2 & V9 >= V3)
  peak.filt <- peaks.all.sort[!(peaks.all.sort$V4 %in%  peak.genes.filt$V4), ]
  peak.filt.sort <- peak.filt %>% arrange(V1, as.numeric(V2), as.numeric(V3))
  
  
  ## filter gtf for transcripts only to get the coordinates of start and end of gene
  gtf.filt.trn <- gtf.filt %>% filter(V3 == "transcript")
  gtf.filt.trn$gene_name <- gsub("\\;.*", "", gsub("transcript_id ", "", gsub("-t.*", "", gtf.filt.trn$V9)))
  gtf.filt.trn$gene_name <- gsub("\"", "", gtf.filt.trn$gene_name)
  
  
  ## Filter for first exon coordinates (exon1 coordinates)
  tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
  tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
  gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
  gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)
  
  
  ## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)
  
  peaks.genes.dist <- bedtoolsr::bt.closest(a = peak.filt.sort, b = gtf.exon1.sort, D = "b", k = 5)
  parse.str2 <- strsplit(peaks.genes.dist$V13, split = ';')
  peaks.genes.dist$gene_name  <- unlist(lapply(parse.str2, '[[' , 3))
  peaks.genes.dist.trns <- left_join(peaks.genes.dist, gtf.filt.trn, by = "gene_name")
  
  ## V16 is the distance of the peak to the exon 1 
  ## we need to overcome the issue with the  ones with  dist = 0
  ## on pos strand V3.x (end of peak) should not exceed V5.y (end of transcript/exon_n)
  ## on neg strand V2.x (start of peak) is not less than V4.y (beggining of the transcript/exon_1)
  
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "+" & V3.x > V5.y))
  
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% filter(!(V16 == 0 & V11 == "-" & V2.x < V4.y))
  
  ## V16 <= 0 means the peak is at upstream 
  ## Find closest gene among top 5 that is upstreaam (min V16)
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 <= 0)
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% group_by(V4.x) %>% 
    mutate(V17 = V16[which.min(abs(V16))])
  
  
  ## Filter the rest
  peaks.genes.dist.trns <- peaks.genes.dist.trns %>% dplyr::filter(V16 == V17)
  
  ## filter the ones that are too far (2000 bp)
  peaks.genes.dist.trns.filt <- peaks.genes.dist.trns %>% dplyr::filter(abs(V16) < 2000)
  
  
  # merge multiple peaks assigned to a single gene
  # the duplicated peaks are the bidirectioonal peaks 
  
  peak.genes <- peaks.genes.dist.trns.filt
  peak.genes <- peak.genes %>% dplyr::select(V1.x, V2.x, V3.x, V11, gene_name) 
  peak.genes.bed.merged <- peak.genes %>% arrange(V2.x) %>% 
    group_by(gene_name) %>% mutate(start_peak = V2.x[which.min(V2.x)], end_peak = V3.x[which.max(V3.x)])  %>% 
    mutate(V4 = ".", V5 = ".")
  
  peak.genes.bed.merged.bed <- peak.genes.bed.merged %>% dplyr::select(V1.x, start_peak, end_peak, V4, V5, V11, gene_name) %>%
    distinct(gene_name, .keep_all = T)
  
  peak.genes.bed.merged.bed <- left_join(peak.genes.bed.merged.bed, prod.desc, by = c("gene_name" = "TGME49"))
  
  # only peaks iinformation to be loaded into IGV
  peak.merged.bed <- peak.genes.bed.merged.bed %>% 
    ungroup() %>% dplyr::select(V1.x,start_peak, end_peak)
  
  
  return(list(peak.gene.all = peaks.genes.dist.trns.filt, 
              peak.gene.merged = peak.genes.bed.merged, 
              peak.gene.merged.bed = peak.genes.bed.merged.bed, 
              peak.gene.merged.IGV =  peak.merged.bed))
  
}


#############################################################################
### individually processing the CUT&RUN 
#############################################################################

## read peaks called by macs2 and process them individually 

in.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")

# read all caled peaks (new AP2_TY (new) vs 4 controls + MiseqA_2mm (old) vs 4 controls)
all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- gsub("\\.narrowPeak", "", files)


# filter significant peaks
qval <- 0.05
all.peaks.filt <- lapply(1:length(all.peaks), function(i){
  cutRun <- all.peaks[[i]]
  cutRun <- cutRun %>% filter(V9 > -log10(qval))
})
names(all.peaks.filt) <- gsub("\\.narrowPeak", "", files)


# peak gene assignment 
gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05
all.peak.genes <- lapply(1:8, function(i){
  peaks.df <- all.peaks[[i]]
  tmp <- get_peak_genes_assign(gtf = gtf, peaks = peaks.df, qval = qval)
  return(tmp)
})
names(all.peak.genes) <- gsub("\\.narrowPeak", "", files[1:8])

saveRDS(all.peak.genes, "../Input/toxo_cdc/rds_ME49_59/all_new_Old_cut_run_peak_gene_assigned_qval_0.05_NO_frag_filt.rds")

## Summary table 
CutRunPeaks.num <- unlist(lapply(all.peaks, nrow), "[[", 1)
CutRunPeaks.num.filt <- unlist(lapply(all.peaks.filt, nrow), "[[", 1)
CutRunPeaks.num.peak.gene <-  unlist(lapply(lapply(all.peak.genes, "[[", 3), nrow))
venn.table <- data.frame(total_cut_run_peaks = CutRunPeaks.num,
                         total_cut_run_peaks_qval = CutRunPeaks.num.filt, 
                         total_peak_gene_assigned = CutRunPeaks.num.peak.gene)


## intersection of 4 new data
new.peaks.venn.list <- list(AP2XII8_vs_AP2XII8_IgG1 = unique(all.peak.genes$`AP2XII-8_Ty_S4_vs_AP2XII-8_IgG1_peaks`$peak.gene.merged.bed$gene_name),
                            AP2XII8_vs_RH_IgG1 = unique(all.peak.genes$`AP2XII-8_Ty_S4_vs_RH_IgG1_S1_peaks`$peak.gene.merged.bed$gene_name),
                            AP2XII8_vs_RH_Neg = unique(all.peak.genes$`AP2XII-8_Ty_S4_vs_RH_Negative_S2_peaks`$peak.gene.merged.bed$gene_name),
                            AP2XII8_vs_RH_Ty = unique(all.peak.genes$`AP2XII-8_Ty_S4_vs_RH_Ty_S2_peaks`$peak.gene.merged.bed$gene_name))

p <- ggVennDiagram(new.peaks.venn.list)
p
ggsave("../Output/toxo_cdc/ME49_59/figures_paper/all_new_cut_run_peak_gene_assigned_overlap_venn_qval_0.05_NO_filt.pdf", 
       plot = p, height = 10, width = 10, dpi = 300)



##############################################
############ union of the 4 data sets ########
##############################################

## concatenating peaks (union of 4 new data sets peaks)
in.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")

all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- gsub("\\.narrowPeak", "", files)


qval <- 0.05
all.peaks.filt <- lapply(1:length(all.peaks), function(i){
  cutRun <- all.peaks[[i]]
  cutRun <- cutRun %>% filter(V9 > -log10(qval))
})
names(all.peaks.filt) <- gsub("\\.narrowPeak", "", files)

tmp <- do.call("rbind",all.peaks.filt[1:4])

# peak gene union new data sets
gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

# AP2_Ty vs 4 controls
all.peaks.new <- all.peaks[1:4]
all.peaks.new <- do.call("rbind", all.peaks.new)
peak.genes.union.new <- get_peak_genes_assign(gtf,all.peaks.new, qval)

saveRDS(peak.genes.union.new, "../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")

# write the cut and run regions in bed format for motif search
write.table(peak.genes.union.new$peak.gene.merged.bed[1:6], 
            "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/bedFasta/peak_genes_union_0_05_NO_frag_filt.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)

# this file is used for view in IGV browser
peak.genes.bed.Only <- peak.genes.union.new$peak.gene.merged.IGV %>% distinct()

write.table(peak.genes.bed.Only, 
            "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/peak_genes_assigned/peak_genes_union_0_05_NO_frag_filt_IGV_distinct.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)




write.table(peak.genes.union.new$peak.gene.merged.IGV, 
            "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/peak_genes_assigned/peak_genes_union_0_05_NO_frag_filt_IGV.bed",
            sep = "\t", quote = F, row.names = F, col.names = F)





##############################################################
## add info of intersection of 4 new data to the mega table 
##############################################################

# individually processed data
all.peak.genes <- readRDS("../Input/toxo_cdc/rds_ME49_59/all_new_Old_cut_run_peak_gene_assigned_qval_0.05_NO_frag_filt.rds")
all.peak.genes.indv <- lapply(all.peak.genes[1:4], "[[", 3) 
all.peak.genes.indv.df <- do.call("rbind", all.peak.genes.indv)
unique(all.peak.genes.indv.df$gene_name) # 970 genes

intrsct.genes <- data.frame(gene_name = Reduce(intersect, lapply(all.peak.genes.indv, function(x) x$gene_name)),
                            intersection = "yes")

# union of peaks
peak.genes.union.new <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")
peak.genes.union.new.df <- peak.genes.union.new$peak.gene.merged.bed

peak.genes.union.new.df.intrsct <- dplyr::left_join(peak.genes.union.new.df, intrsct.genes, by = "gene_name")
peak.genes.union.new.df.intrsct$cutRun.peaks <- "yes"

# genes in the intersection
tmp2 <- peak.genes.union.new.df.intrsct %>% filter(intersection == "yes")

saveRDS(peak.genes.union.new.df.intrsct, "../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info.rds")

##################################################################
# peak gene union old data set- We are not using this
##################################################################

gtf.file <- "../Genomes/ToxoDB-59_TgondiiME49.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
qval <- 0.05

# Miseq_A vs 4 controls
all.peaks.old <- all.peaks[4:8]
#all.peaks.old <- all.peaks[7] # miseqA_vs_RH_neg
all.peaks.old <- do.call("rbind", all.peaks.old)
peak.genes.union.old <- get_peak_genes_assign(gtf,all.peaks.old, qval)

ven.list <- list(new.union = unique(peak.genes.union.new$peak.gene.merged.bed$gene_name), 
                 MiseqA_vs_RH_neg = unique(peak.genes.union.old$peak.gene.merged.bed$gene_name))
ggVennDiagram(ven.list)


## table 
CutRunPeaks.num <- unlist(lapply(all.peaks[1:4], nrow), "[[", 1)
CutRunPeaks.num.filt <- unlist(lapply(all.peaks.filt[1:4], nrow), "[[", 1)
peak.gene.union.new.num <-  rep(nrow(peak.genes.union.new$peak.gene.merged.bed),4)
venn.table <- data.frame(total_cut_run_peaks = CutRunPeaks.num,
                         total_cut_run_peaks_qval = CutRunPeaks.num.filt, 
                         peak.gene.union.new.num = peak.gene.union.new.num)



############### motif search ################
## to perform motif seach we can use the bed file (regions of peak)
## get fasta file using the following command
## upload fasta files on BAMM tool 
# https://bammmotif.soedinglab.org/job/denovo/

## fasta
#bedtools getfasta -s -fi  ../../../Genome/ToxoDB-59_TgondiiME49_Genome.fasta -bed CutRun_KD_targets.bed -fo CutRun_KD_targets.fasta
#bedtools getfasta -fi  ../../../../Genome/ToxoDB-59_TgondiiME49_Genome.fasta  -bed peak_genes_union_0_05_NO_frag_filt.bed -fo peak_genes_union_0_05_NO_frag_filt.fasta



############### which genes have which motifs?? #############
## take the bed files return from BAMM motif search, 
## each identified motif will have one bed file representing regions with the occurence of motif

## overlap the regions containing motifs with cut and run peak gene assigned,
## to get the genes that have motif

## TGME49_253900
prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1_ME49.xlsx') 
peak.genes.union <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info.rds")
peak.genes.union <- peak.genes.union[1:7]

in.dir <- '../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/BAMM_analysis/allPeaksUnion/bed_BAMM/'
out.dir <- in.dir

bed.files <- list.files(in.dir, pattern = ".bed")
bed.names <- gsub(".bed", "", bed.files)


options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
bed.genes.list <- lapply(1:length(bed.files), function(i){
  
  bed.tmp <- read.table(paste(in.dir, bed.files[i], sep = "/"), header = F, quote = NULL, sep = "\t")
  
  ## intersect with the peak gene 
  bed.genes <- bedtoolsr::bt.intersect(a = bed.tmp, b = peak.genes.union, wo = T)
  bed.genes <- left_join(bed.genes, prod.desc, by = c("V13" = "TGME49"))
  bed.genes <- bed.genes %>% dplyr::select(V1, V2, V3, V4, V5, V6, V12,V13,  ProductDescription)
  bed.genes$has.motif <- "yes"
  bed.genes$motif <- gsub("peak_genes_union_0_05_NO_frag_filt_", "", bed.names[i])
  
  write.xlsx(bed.genes, paste(out.dir, paste(bed.names[i], "genes.xlsx", sep = "_"), sep =""))
  
  return(bed.genes)
  
})

names(bed.genes.list) <- bed.names

saveRDS(bed.genes.list, "../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_motif_info_NO_frag_filt.rds")


CutRunPeaksMotif <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_motif_info_NO_frag_filt.rds")
CutRunPeaksMotif <- do.call("rbind", CutRunPeaksMotif)
CutRunPeaksMotif <- CutRunPeaksMotif %>% 
  dplyr::select(V13, ProductDescription, has.motif, motif) %>% distinct() 


# cut and run peaks (union of new data sets)
peak.genes.union <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info.rds")

## add motif info to the table
CutRun.all.info <- left_join(peak.genes.union, CutRunPeaksMotif, by = c("gene_name" = "V13"))
saveRDS(CutRun.all.info,"../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info_motif_info.rds" )


##################################################
## Add chip info HDAC3, MORC and AP2XII-2 chip seq
##################################################

AP2XII2.Chip <- read.xlsx("../Input/Toxo_genomics/genes/HDAC_AP2XII2_MORC_chip.xlsx", sheet = 1)
AP2XII2.Chip$chip <- "AP2XII-2"
HDAC3.chip <- read.xlsx("../Input/Toxo_genomics/genes/HDAC_AP2XII2_MORC_chip.xlsx", sheet = 2)
HDAC3.chip$chip <- "HDAC3"
MORC.chip <- read.xlsx("../Input/Toxo_genomics/genes/HDAC_AP2XII2_MORC_chip.xlsx", sheet = 3)
MORC.chip$chip <- "Morc"

all.chip <- rbind(AP2XII2.Chip, HDAC3.chip, MORC.chip)
all.chip <- all.chip %>% dplyr::select(Name, chip) %>% distinct()
 

## GCN5 chip 
## modifying gene ids that have been updated on toxodb.org

GCN5 <- read.xlsx('../Input/Toxo_genomics/genes/GCN5_Chip_Chip.xlsx') # all 3 reps of chip-chip
GCN5.targ <- GCN5
GCN5.targ$ID <- gsub("\\_0", "_2", GCN5.targ$ID)
GCN5.targ$ID <- gsub("\\_1", "_3", GCN5.targ$ID)
GCN5.targ <- GCN5.targ %>% dplyr::select(ID, Description)
colnames(GCN5.targ) <- c("TGME49ID", 'Description')
GCN5.targ$chip <- "GCN5"
GCN5.targ$TGME49ID <- gsub("TGME49_319310", "TGME49_319312", 
                           gsub("TGME49_314950", "TGME49_314955", 
                                gsub("TGME49_295900", "TGME49_295850", 
                                     gsub("TGME49_292930", "TGME49_292935", 
                                          gsub("TGME49_285270", "TGME49_285272", 
                                               gsub("TGME49_281410", "TGME49_281400", 
                                                    gsub("TGME49_222950", "TGME49_222948", GCN5.targ$TGME49ID)))))))  

GCN5.targ <- GCN5.targ %>% dplyr::select(-Description)
colnames(GCN5.targ)[1] <- "Name"

all.chip <- rbind(all.chip,GCN5.targ )
saveRDS(all.chip, "../Input/toxo_cdc/rds_ME49_59/all_public_chip_seq.rds")


## add chip info to cut and run table
CutRun <- readRDS( "../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info_motif_info.rds")
CutRun.peaks.motif.chip <- left_join(CutRun, all.chip, by = c("gene_name"= "Name"))



write.xlsx(CutRun.peaks.motif.chip,"../Output/toxo_cdc/ME49_59/tables/union_peaks_cut_run_motif_NO_frag_filt_motif_chip.xlsx")
saveRDS(CutRun.peaks.motif.chip,"../Input/toxo_cdc/rds_ME49_59/union_peaks_cut_run_motif_NO_frag_filt_motif_chip.rds")



########### Overlap of cut and run peaks  with atac- venn diagram  ########
########### after peak - gene assignment ##################################

## this is not the right way the next section is what we need #############
##  atac peaks (after peak gene assignments)

## use the unique peaks that have been assigned to genes 
peak.genes.bed.merged.bed <- read.table("../Input/toxo_scATAC_MJ_ME49_59/peak_gene_assigned_final.bed", sep = "\t")
peak.genes.bed.merged.bed$peakLoc <- paste(peak.genes.bed.merged.bed$V1, paste(peak.genes.bed.merged.bed$V2, peak.genes.bed.merged.bed$V3, sep = "-"), sep = ":")
peak.genes.bed.merged.bed <- peak.genes.bed.merged.bed %>% dplyr::select(-c(V6, V7)) %>% distinct()
nrow(peak.genes.bed.merged.bed)

# cut and run (after peak gene assignment) 
peak.genes.union.new <- readRDS( "../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt.rds")
cut.run.peaks <- peak.genes.union.new$peak.gene.merged.bed[1:5] %>% distinct()
nrow(cut.run.peaks)
  # cut.run.peaks <- peak.genes.union.new$peak.gene.merged.bed %>% ungroup() %>% 
#   dplyr::select(V1.x, start_peak, end_peak) %>% distinct()

# overlap atac and cut RUN
options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
peak.cutRun.atac.ovrlp <- bedtoolsr::bt.intersect(a = cut.run.peaks, b = peak.genes.bed.merged.bed, wb = T)

#peak.cutRun.atac.ovrlp <- peak.cutRun.atac.ovrlp %>% distinct(V7, .keep_all = T)


library(VennDiagram)
pdf(file = "../Output/toxo_cdc/ME49_59/figures_paper/cut_run_union_peaks_overlap_atac_peaks_venn.pdf")
venn.plot <- draw.pairwise.venn(
  area1 = nrow(peak.genes.bed.merged.bed),
  area2 = nrow(cut.run.peaks),
  cross.area = nrow(peak.cutRun.atac.ovrlp),
  #category = c("ATAC", "C&R"),
  fill = c("#F19F39", "#469C2C"),
  lty = "blank",
  cex = 4,
  cat.cex = 3,
  #cat.pos = c(285, 105),
  #cat.dist = 0.09,
  #cat.just = list(c(-1, -1), c(1, 1)),
  #ext.pos = 30,
  #ext.dist = -0.05,
  ext.length = 0.9,
  ext.line.lwd = 2.5,
  #ext.line.lty = "dashed"
)
grid.draw(venn.plot);

dev.off()

########### Overlap of cut and run peaks  with atac- venn diagram  ########
########### befor peak - gene assignment ##################################

in.dir <- "../Input/toxo_cdc/cutNrun/all_macs2_old_new_batch_NO_Filter/"
files <- list.files(in.dir, pattern = ".narrowPeak")

all.peaks <- list()
for (f in files) {
  tmp <- read.table(paste(in.dir, f, sep = ''), header=F, sep="\t", quote = NULL)
  tmp <- tmp %>% filter(!grepl("KE.*", V1)) 
  all.peaks <- c(all.peaks, list(tmp))
}
names(all.peaks) <- gsub("\\.narrowPeak", "", files)


qval <- 0.05
all.peaks.filt <- lapply(1:length(all.peaks), function(i){
  cutRun <- all.peaks[[i]]
  cutRun <- cutRun %>% filter(V9 > -log10(qval))
})
names(all.peaks.filt) <- gsub("\\.narrowPeak", "", files)
all.peaks.filt.df <- do.call("rbind",all.peaks.filt[1:4])
#tmp <- get_peak_genes_assign(all.peaks.filt.df, gtf = gtf, qval = 0.05)


rownames(all.peaks.filt.df) <- NULL
all.peaks.filt.df <- all.peaks.filt.df %>% transmute(V1 = V1, V2 = V2, V3 = V3, V4 = V4) 
cut.run.peaks <- all.peaks.filt.df 

#  atac peaks (before peak gene assignments)

Tg_ATAC <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_ATAC_peak.rds')

peak.regions <- data.frame(ATAC_peak_region = rownames(Tg_ATAC@assays$peaks@data))
peak.regions <- peak.regions %>% filter(!grepl("KE.*", ATAC_peak_region))
peak.regions.bed  <- data.frame(do.call(rbind, strsplit(peak.regions$ATAC_peak_region,"-")))
peak.regions.bed$X1 <- paste(peak.regions.bed$X1, peak.regions.bed$X2, sep = "_")
peaks.all.sort.atac <- peak.regions.bed  %>% dplyr::select(X1, X3, X4) %>%  arrange(X1, as.numeric(X3), as.numeric(X4))
names(peaks.all.sort.atac) <- c("V1", "V2", "V3")
peaks.all.sort.atac$V4 <- paste(paste(peaks.all.sort.atac$V1, peaks.all.sort.atac$V2, sep = ":"),peaks.all.sort.atac$V3 ,sep = "-" )
peak.genes.bed.merged.bed <- peaks.all.sort.atac

# overlap atac and cut RUN
options(bedtools.path = "/Users/kourosh.zarringhalam/miniconda3/bin/")
peak.cutRun.atac.ovrlp <- bedtoolsr::bt.intersect(a = cut.run.peaks, b = peak.genes.bed.merged.bed, wb = T)


library(VennDiagram)
pdf(file = "../Output/toxo_cdc/ME49_59/figures_paper/cut_run_union_peaks_overlap_atac_peaks_venn_v2.pdf",
    width = 12, height = 12)
venn.plot <- draw.pairwise.venn(
  area1 = nrow(cut.run.peaks),
  area2 = nrow(peak.genes.bed.merged.bed),
  cross.area = nrow(peak.cutRun.atac.ovrlp),
  #category = c("ATAC", "C&R"),
  fill = c("#469C2C","#F19F39"),
  lty = rep("solid", 2),
  lwd = 6,
  col = c("darkgreen", "darkorange"),
  cex = 5.5,
  cat.cex = 3,
  #cat.pos = c(285, 105),
  #cat.dist = 0.09,
  #cat.just = list(c(-1, -1), c(1, 1)),
  #ext.pos = 30,
  #ext.dist = -0.05,
  ext.length = 0.9,
  ext.line.lwd = 2.5,
  #ext.line.lty = "dashed"
)
grid.draw(venn.plot);

dev.off()
