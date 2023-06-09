library(openxlsx)
library(tidyverse)
library(splines)
library(parallel)
library(ggplot2)
library(tidytext)
library(ggrepel)
library(geomtextpath)
library(ggVennDiagram)




source('./util_funcs.R')

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)



getCurvePeakLoc <- function(t, y, prob = 0.8){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=round(s.0$y, digits = 3), s1=round(s.1$y, digits = 3))
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) >= 1 & any(locs$values== -1)){
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
  if(length(entity.x) == 0){
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}


prod.disc <- read.xlsx('../Input/ID_convert/ProductDescription.xlsx')
GT1.ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.disc <- left_join(prod.disc, GT1.ME49, by = c('GeneID' = 'TGGT1')) %>% na.omit()
prod.disc$TGME49 <- gsub('_', '-', prod.disc$TGME49)

sc.rna.genes.expr.pt <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.genes.expr.pt <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')

## Phase-based DEGs
Intra.markers.sig <- readRDS('../Input_KZ//toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')

## AP2s
AP2s <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx', sheet = 2)
AP2s <- AP2s[grep("AP2", AP2s$Ap2Name),c(1,2)]
AP2s <- left_join(AP2s, GT1.ME49, by = c('GeneID' = 'TGGT1'))
AP2s$GeneID <- gsub('_', '-', AP2s$TGME49)





genes <- unique(sc.rna.genes.expr.pt$GeneID)
fft.stas <- mclapply(1:length(genes), function(i){
  
  rna.sp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = x, y = y) %>% arrange(x)

  atac.sp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = x, y = y) %>% arrange(x)
  
  
  rna.upper.expr <- mean(rna.sp$y[rna.sp$y > quantile(rna.sp$y, p = 0.75)])
  
  rna.transform = fft(rna.sp$y)/length(rna.sp$y)
  rna.magTransform = abs(rna.transform)
  
  rna.amp <- max(rna.magTransform[-1]) * 2
  rna.freq <- which.max(rna.magTransform[-1])
  rna.peak.time <- getCurvePeakLoc(rna.sp$x, rna.sp$y, prob = 0.8)
  
  atac.upper.expr <- mean(atac.sp$y[atac.sp$y > quantile(atac.sp$y, p = 0.75)])
  
  atac.transform = fft(atac.sp$y)/length(atac.sp$y)
  atac.magTransform = abs(atac.transform)
  
  atac.amp <- max(atac.magTransform[-1]) * 2
  atac.freq <- which.max(atac.magTransform[-1])
  atac.peak.time <- getCurvePeakLoc(atac.sp$x, atac.sp$y, prob = 0.8)
  
  L <- list(rna.upper.expr, rna.amp, rna.freq, rna.peak.time,
            atac.upper.expr, atac.amp, atac.freq, atac.peak.time)
}, mc.cores = num.cores)

stats <- data.frame(GeneID = genes, 
                    rna.upper.expr = unlist(lapply(fft.stas, `[[`, 1)),
                    rna.amp = unlist(lapply(fft.stas, `[[`, 2)),
                    rna.freq = unlist(lapply(fft.stas, `[[`, 3)),
                    rna.peak.time = unlist(lapply(fft.stas, `[[`, 4)),
                    atac.upper.expr = unlist(lapply(fft.stas, `[[`, 5)),
                    atac.amp = unlist(lapply(fft.stas, `[[`, 6)),
                    atac.freq = unlist(lapply(fft.stas, `[[`, 7)),
                    atac.peak.time = unlist(lapply(fft.stas, `[[`, 8)))


stats <- left_join(stats, prod.disc, by = c('GeneID' = 'TGME49'))

rna.expr.cutoff <- quantile(stats$rna.upper.expr, prob = 0.01, na.rm = T)
stats$rna.expressed <- ifelse(stats$rna.upper.expr >= rna.expr.cutoff, 1, 0)

atac.expr.cutoff <- quantile(stats$atac.upper.expr, prob = 0.01, na.rm = T)
stats$atac.expressed <- ifelse(stats$atac.upper.expr >= atac.expr.cutoff, 1, 0)

rna.amp.cutoff <- quantile(stats$rna.amp, prob = 0.5)
stats$rna.cyclic<- ifelse(stats$rna.amp >= rna.amp.cutoff , 1, 0)

atac.amp.cutoff <- quantile(stats$atac.amp, prob = 0.5)
stats$atac.cyclic<- ifelse(stats$atac.amp >= atac.amp.cutoff , 1, 0)

stats.expressed <- stats %>% dplyr::filter(rna.expressed == 1)
nrow(stats.expressed) / nrow(stats)
stats.cyclic.rna <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1)
nrow(stats.cyclic.rna) / nrow(stats)
stats.cyclic.atac <- stats %>% dplyr::filter(rna.expressed == 1, atac.cyclic == 1)
stats.cyclic.both <- stats %>% dplyr::filter(rna.expressed == 1, rna.cyclic == 1 & atac.cyclic == 1)

nrow(stats.cyclic.both) / nrow(stats)
#stats$filt <- stats$expressed * stats$cyclic

sum(stats$rna.cyclic) / nrow(stats)




## DEGs
all(Intra.markers.sig$gene %in% stats$GeneID[stats$rna.cyclic == 1])
length(Intra.markers.sig$gene)


stats$rna.constitutive <- ifelse(stats$rna.expressed == 1 & stats$rna.cyclic == 0, 1, 0)
stats$atac.constitutive <- ifelse(stats$atac.expressed == 1 & stats$atac.cyclic == 0, 1, 0)

write.xlsx(stats, '../Output_KZ/tables/all_genes_cyclic_timing.xlsx')

saveRDS(stats, '../Input_KZ/toxo_cdc/rds_ME49_59/all_genes_cyclic_timing.rds')


### Venn Diagram
library(ggVennDiagram)
library(VennDiagram)


my.data <- stats %>% dplyr::filter(rna.expressed == 1) %>% dplyr::select(GeneID, rna.cyclic, atac.cyclic) 

data.list <- list(scRNA = my.data$GeneID[my.data$rna.cyclic == 1],
                  scATAC = my.data$GeneID[my.data$atac.cyclic == 1])

library(VennDiagram) 
p <- venn.diagram(data.list, fill = c("lightblue", "lightgreen"), 
                  alpha = c(0.5, 0.5),  # Numbers
                  category.names = c("scRNA" , "scATAC"),
                  cex = 1,
                  fontface = "bold",
                  fontfamily = "sans",
                  lwd =2, 
                  height = 1400 , 
                  width = 1400 , 
                  lty = 'blank',
                  resolution = 300,
                  "venn_diagram.tiff")


###
##
stats.cyclic <- stats %>% dplyr::filter(atac.cyclic == 1 & rna.cyclic == 1)
stats.cyclic.AP2 <- stats %>% dplyr::filter(atac.cyclic == 1 & rna.cyclic == 1 & grepl("AP2 domain transcription factor", ProductDescription))
plot(stats.cyclic.AP2$rna.peak.time, stats.cyclic.AP2$atac.peak.time)
lm(stats.cyclic.AP2$atac.peak.time~stats.cyclic.AP2$rna.peak.time)
abline(2.8606,  0.1497)

stats.cyclic.AP2$label <- gsub("AP2 domain transcription factor ", "", stats.cyclic.AP2$ProductDescription)
p <- ggplot(data = stats.cyclic.AP2, aes(x = rna.peak.time, y = atac.peak.time)) + 
  geom_point(color = 'black', size = 2)+ 
  geom_smooth(method = "lm", color = 'blue') + 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  ylab('peak atac time') + xlab('peak expr time') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) + 
  geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
                  box.padding = unit(0.6, "lines"),
                  max.overlaps = 300,
                  #segment.angle = 180,
                  nudge_x = 0.25, 
                  nudge_y = 0.25,
                  hjust=0.25,
                  #nudge_x=0.25, 
                  segment.size = 0.1,
                  na.rm = TRUE)+ 
  
  theme(#legend.position = c(0.15, 0.85),
    legend.position = 'none',
    legend.title = element_text(colour="black", size=12, 
                                face="bold"),
    legend.text = element_text(colour="black", size=12, 
                               face="bold"))


p

ggsave(filename='../Output_KZ/figures/AP2_RNA_ATAC.pdf',
       plot=p,
       width = 6, height = 6,
       units = "in")
