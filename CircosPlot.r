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

source('./util_funcs.R')
prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1_ME49.xlsx') 
getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = prob))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}


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



sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

# # ## Filter to include markers only
# sc.rna.spline.fits <- sc.rna.spline.fits %>% dplyr::filter(GeneID %in% DEG.sig$gene)
# sc.atac.spline.fits <- sc.atac.spline.fits %>% dplyr::filter(GeneID %in% c(DEG.sig$gene))
# 

rna_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O_intra_atac_lables_pt.rds')


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


na.ind <- which(apply(sc.rna.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.rna.dtw.wide <- sc.rna.dtw.wide[,-na.ind]
}


sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

na.ind <- which(apply(sc.atac.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.atac.dtw.wide <- sc.atac.dtw.wide[,-na.ind]
}


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


sc.rna.peak.order <- sc.rna.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, sc.rna.peak.order, by = 'GeneID')


sc.atac.peak.order <- sc.atac.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.atac.mu.scale <- left_join(sc.atac.mu.scale, sc.atac.peak.order, by = 'GeneID')


# ## Filter to include markers only
# sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% DEG.sig$gene)
# sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% DEG.sig$gene)

saveRDS(sc.rna.mu.scale, "../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_mu_scale_peak_ord_all_genes.rds")
saveRDS(sc.atac.mu.scale, "../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_mu_scale_peak_ord_all_genes.rds")

#sc.rna.mu.scale <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_mu_scale_phase.rds')
sc.rna.mu.scale.peak.ord <- sc.rna.mu.scale %>% 
  transmute(GeneID = GeneID, peak.ord.rna = peak.ord) %>% distinct()

#sc.atac.mu.scale <- readRDS('../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_mu_scale_phase.rds')
sc.atac.mu.scale.peak.ord <- sc.atac.mu.scale %>% 
  transmute(GeneID = GeneID, peak.ord.atac = peak.ord) %>% distinct()

sc.rna.atac.peak.ord <- left_join(sc.atac.mu.scale.peak.ord,sc.rna.mu.scale.peak.ord, by = "GeneID")
sc.rna.atac.mu.scale.peak.loc <- left_join(sc.rna.atac.peak.ord, ME49.trans, by = "GeneID") 
sc.rna.atac.mu.scale.peak.loc <- sc.rna.atac.mu.scale.peak.loc %>% 
  dplyr::select(GeneID, peak.ord.atac, peak.ord.rna, V1, V4, V5)

saveRDS(sc.rna.atac.mu.scale.peak.loc, "../Input/toxo_cdc/rds_ME49_59/sc_rna_atac_mu_peak_expr_access_circos_all_genes.rds")


#BioCircos()

# tracklist = BioCircosTextTrack('myTextTrack', 'Circos', size = "1em", opacity = 0.5, 
#                                x = -0.67, y = -0.5)

# tracklist = BioCircosArcTrack("peaks", as.character("chrIa"), 
#                                     starts = 800, ends = 100000, labels = 1,
#                                     maxRadius = 0.97, minRadius = 0.03)


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


## gene families

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1_ME49.xlsx') 

## AP2s
AP2s <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx') 
AP2s <- AP2s %>% dplyr::select(GeneID, Ap2Name)
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
AP2s <- left_join(AP2s, TGGT1_ME49, by = c("GeneID" = "TGGT1"))
AP2s$TGME49 <- gsub("_", "-", AP2s$TGME49)
AP2s <- AP2s[grep("AP2", AP2s$Ap2Name),]

write.xlsx(AP2s, "../Input/toxo_genomics/genes/AP2s.xlsx")

## IMCs

## New
IMCs <- read.xlsx("../Input/Toxo_genomics/genes/IMC gene list 090122.xlsx")
IMCs <- IMCs %>% dplyr::select(Gene.ID, new.name)
IMCs <- left_join(IMCs, TGGT1_ME49, by = c("Gene.ID" = "TGGT1"))
IMCs$TGME49 <- gsub('_', '-', IMCs$TGME49)
colnames(IMCs) <- gsub("new.name", "ProductDescription", colnames(IMCs))

write.xlsx(IMCs, "../Input/toxo_genomics/genes/IMCs.xlsx")

## old
# IMCs <- prod.desc[grep('IMC', prod.desc$ProductDescription), ]
# IMCs$ProductDescription <- unlist(lapply(strsplit(IMCs$ProductDescription, split = ' '), function(item){item[length(item)]}))
# IMCs <- left_join(IMCs, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
# IMCs$TGME49 <- gsub('_', '-', IMCs$TGME49)

## BCs


BCs <- read.xlsx("../Input/Toxo_genomics/genes/BC gene list 021723.xlsx")
BCs <- BCs[-1,] %>% dplyr::select(TGGT1, Name )
BCs <- left_join(BCs, TGGT1_ME49, by = "TGGT1")
BCs$TGME49 <- gsub("_", "-", BCs$TGME49)
colnames(BCs) <- gsub("Name", "ProductDescription", colnames(BCs))

write.xlsx(BCs, "../Input/toxo_genomics/genes/BCs.xlsx")

genes <- IMCs$TGME49
genes <- AP2s$TGME49
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








