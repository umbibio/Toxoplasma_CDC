
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidyverse)
library(RColorBrewer)
library(ggraph)
library(graphlayouts)
library(Seurat)
library(patchwork)
library(hdf5r)
library(GenomeInfoDbData)
library(GenomicRanges)
library(BioCircos)
library(Biostrings)


# cusomized Genome
ME49.gtf  <- read.delim("../Input/Toxo_genomics/genome/ToxoDB-59_TgondiiME49.gtf", header = F)
ME49.gtf <- ME49.gtf %>% filter(str_detect( V1, "chr"))
ME49.trans <- ME49.gtf %>% filter(V3 == "transcript")
ME49.trans$GeneID <- gsub(";", "", unlist(lapply(strsplit(ME49.trans$V9, " "), "[[", 4)))
ME49.trans$GeneID <- gsub("_", "-", ME49.trans$GeneID)
ME49.trans$V1 <- gsub("TGME49_", "", ME49.trans$V1)

ME49.fasta <- readDNAStringSet("../Input/toxo_genomics/genome/ToxoDB-59_TgondiiME49_Genome.fasta")
chrs <- names(ME49.fasta)[grep("TGME49_chr", names(ME49.fasta))]
chr.len <- data.frame(chr = gsub(" ", "", unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 1))),
                      len = as.numeric(gsub('length=', '', unlist(lapply(strsplit(chrs, split = '\\|'), `[[`, 4)))))
chr.len$chr <- unlist(gsub("TGME49_", "", chr.len$chr))

Chr <- chr.len$chr
ME49 <- chr.len$len
names(ME49) <- Chr


## peak expr and peak access of up regulated genes

sc.rna.atac.mu.scale.peak.loc <- readRDS("../Input/toxo_cdc/rds_ME49_59/sc_rna_atac_mu_peak_expr_access_circos_all_genes.rds")

## Circos plot 

tracklist = BioCircosArcTrack("CommGenes", 
                              as.character(sc.rna.atac.mu.scale.peak.loc$V1), 
                              starts = sc.rna.atac.mu.scale.peak.loc$V4, 
                              ends = sc.rna.atac.mu.scale.peak.loc$V5, 
                              colors = "black", 
                              labels = sc.rna.atac.mu.scale.peak.loc$GeneID,
                              maxRadius = 0.98, minRadius = 0.92)


tracklist = tracklist + BioCircosBackgroundTrack("all phase",
                                                 minRadius = 0.3, maxRadius = 0.9,
                                                 borderColors = "black", borderSize = 0.3, 
                                                 fillColors = "white")
tracklist = tracklist + BioCircosBackgroundTrack("C phase",
                                                 minRadius = 0.9, maxRadius = 0.9,
                                                 borderColors = "#b138ee", borderSize = 2.2,
                                                 fillColors = "white")


tracklist = tracklist + BioCircosBackgroundTrack("G phase",
                                                 minRadius = 0.6, maxRadius = 0.6,
                                                 borderColors = "#C2401F", borderSize = 2.2, 
                                                 fillColors = "white")

tracklist = tracklist + BioCircosBackgroundTrack("S phase",
                                                 minRadius = 0.77, maxRadius = 0.77,
                                                 borderColors = "#caae05", borderSize = 2.2, 
                                                 fillColors = "white")

tracklist = tracklist + BioCircosBackgroundTrack("M phase",
                                                 minRadius = 0.8, maxRadius = 0.8,
                                                 borderColors = "#6f883a", borderSize = 2.2, 
                                                 fillColors = "white")

## example plots: AP2s, IMCs, BCs

AP2s <- read.xlsx("../Input/toxo_genomics/genes/AP2s.xlsx") ## AP2s
AP2s <- read.xlsx("../Input/toxo_genomics/genes/Cyclic_AP2s_review_paper.xlsx")
colnames(AP2s) <- gsub("GeneID","TGME49", colnames(AP2s))

IMCs <- read.xlsx("../Input/toxo_genomics/genes/IMCs.xlsx") ## IMCs
## you may want to get the ones that are DEGs then left join with the following table
DEG.sig <- readRDS( '../Input/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')


BCs <- read.xlsx("../Input/toxo_genomics/genes/BCs.xlsx") ## BCs
## you may want to get the ones that are DEGs then left join with the following table
DEG.sig <- readRDS( '../Input/toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')


genes <- AP2s$TGME49
genes <- IMCs$TGME49
genes <- BCs$TGME49

# expression peak 
point.chr <- sc.rna.atac.mu.scale.peak.loc$V1[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
point.loc.start <- sc.rna.atac.mu.scale.peak.loc$V4[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
point.vals.rna <- sc.rna.atac.mu.scale.peak.loc$peak.ord.rna[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
gene.labels <- sc.rna.atac.mu.scale.peak.loc$GeneID[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]

tracklist = tracklist + BioCircosSNPTrack('peak.order.rna',
                                          point.chr, point.loc.start, point.vals.rna,
                                          labels = gene.labels,
                                          colors = c("red"), borderColors = "#AAAAAA",
                                          borderSize = 0.6, size  = 3,
                                          minRadius = 0.3, maxRadius = 0.9)

point.chr <- sc.rna.atac.mu.scale.peak.loc$V1[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
point.loc.start <- sc.rna.atac.mu.scale.peak.loc$V4[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
point.vals.atac <- sc.rna.atac.mu.scale.peak.loc$peak.ord.atac[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]
gene.labels <- sc.rna.atac.mu.scale.peak.loc$GeneID[sc.rna.atac.mu.scale.peak.loc$GeneID %in% genes]

tracklist = tracklist + BioCircosSNPTrack('peak.order.atac',
                                          point.chr, point.loc.start, point.vals.atac,
                                          labels = gene.labels,
                                          colors = c("darkblue"), borderColors = "#AAAAAA",
                                          borderSize = 0.6, size  = 3,
                                          minRadius = 0.3, maxRadius = 0.9)



BioCircos(tracklist,genomeFillColor = "YlOrBr",
          #genomeFillColor = c("gray", "darkgray"),
          #genomeFillColor = "Greys", 
          genome =  as.list(ME49), 
          genomeTicksLen = 3, genomeTicksTextSize = 0, 
          genomeTicksScale = 50000000, 
          genomeLabelTextSize = 18, genomeLabelDy = 0)



## circos plot for all genes
## all genes
# atac peak
point.chr <- sc.rna.atac.mu.scale.peak.loc$V1
point.loc.start <- sc.rna.atac.mu.scale.peak.loc$V4
point.vals.atac <- sc.rna.atac.mu.scale.peak.loc$peak.ord.atac

tracklist = tracklist + BioCircosSNPTrack('peak.order.atac', 
                                          point.chr, point.loc.start, point.vals.atac, 
                                          labels = sc.rna.atac.mu.scale.peak.loc$GeneID,
                                          colors = c("darkblue"),borderColors = "#AAAAAA",
                                          borderSize = 0.6, size = 1.2,
                                          minRadius = 0.3, maxRadius = 0.9)

# expression 
# expression peak 
point.chr <- sc.rna.atac.mu.scale.peak.loc$V1
point.loc.start <- sc.rna.atac.mu.scale.peak.loc$V4
point.vals.rna <- sc.rna.atac.mu.scale.peak.loc$peak.ord.rna
gene.labels <- sc.rna.atac.mu.scale.peak.loc$GeneID

tracklist = tracklist + BioCircosSNPTrack('peak.order.rna',
                                          point.chr, point.loc.start, point.vals.rna,
                                          labels = gene.labels,
                                          colors = c("red"), borderColors = "#AAAAAA",
                                          borderSize = 0.6, size  = 1.2,
                                          minRadius = 0.3, maxRadius = 0.9)


BioCircos(tracklist,genomeFillColor = "YlOrBr",
          #genomeFillColor = c("gray", "darkgray"),
          #genomeFillColor = "Greys", 
          genome =  as.list(ME49), 
          genomeTicksLen = 3, genomeTicksTextSize = 0, 
          genomeTicksScale = 50000000, 
          genomeLabelTextSize = 18, genomeLabelDy = 0)



