library(tidyverse)
library(openxlsx)
library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(tidytext)




## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

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


sc.rna.spline.fits <- readRDS('../Input_YR//toxo_cdc/rds_ME49_59/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input_YR//toxo_cdc/rds_ME49_59/sc_atac_spline_fits_all_genes.rds')


marker.genes <- readRDS('../Input_YR//toxo_cdc/rds_ME49_59/Intra_markers_sig.rds')
rna.trans.marker.genes <- marker.genes

#rna.trans.marker.genes <- readRDS('../Input/toxo_cdc/rds_ME49_59/rna_markers_rna_trns_sig_v2.rds')
rna.trans.marker.genes <- rna.trans.marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()
rna.trans.marker.genes %>% group_by(phase) %>% summarise(n())

## Filter to include transition markers only
sc.rna.spline.fits <- sc.rna.spline.fits %>% dplyr::filter(GeneID %in% rna.trans.marker.genes$GeneID)
sc.atac.spline.fits <- sc.atac.spline.fits %>% dplyr::filter(GeneID %in% rna.trans.marker.genes$GeneID)


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


sc.rna.peak.order <- sc.rna.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, sc.rna.peak.order, by = 'GeneID')
sc.rna.mu.scale <- left_join(sc.rna.mu.scale, rna.trans.marker.genes, by = 'GeneID')

sc.rna.mu.scale$GeneID <- factor(sc.rna.mu.scale$GeneID, 
                                 levels = unique(sc.rna.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))


sc.atac.peak.order <- sc.atac.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(x, expr))
sc.atac.mu.scale <- left_join(sc.atac.mu.scale, sc.atac.peak.order, by = 'GeneID')
sc.atac.mu.scale <- left_join(sc.atac.mu.scale, rna.trans.marker.genes, by = 'GeneID')

## RNA peak order
sc.atac.mu.scale$peak.ord.rna <- sc.rna.mu.scale$peak.ord[match(sc.atac.mu.scale$GeneID, sc.rna.mu.scale$GeneID)]

# sc.atac.mu.scale <- sc.atac.mu.scale  %>% mutate(GeneID = as.factor(GeneID),
#                                                name = reorder_within(GeneID, by = -peak.ord.rna, within = phase))


sc.atac.mu.scale$GeneID <- factor(sc.atac.mu.scale$GeneID,
                                  levels = unique(sc.atac.mu.scale$GeneID[order(-sc.rna.mu.scale$peak.ord)]))


saveRDS(sc.rna.mu.scale, '../Input/toxo_cdc/rds_ME49_59/sc_rna_spline_mu_scale_rna_trans.rds')
saveRDS(sc.atac.mu.scale, '../Input/toxo_cdc/rds_ME49_59/sc_atac_spline_mu_scale_atac_trans.rds')

## Look at Marker genes only (no down genes)
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% rna.trans.marker.genes$GeneID)

p1 <- ggplot(sc.rna.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  # facet_grid(phase~., scales = "free",  space='free',
  #            labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('scRNA') + 
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(size = 20, face = "bold"), 
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=20, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    legend.position = "none") 
#theme(plot.margin=unit(c(2.5,2.5,2.5,2.5),"cm"))

plot(p1)



## Look at Marker genes only (no down genes)
sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% rna.trans.marker.genes$GeneID)

p2 <- ggplot(sc.atac.mu.scale, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  # facet_grid(phase~., scales = "free",  space='free',
  #            labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") + ggtitle('scATAC') +
  #scale_fill_gradientn(colours = hm.palette(10)) +
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.1, "lines")) + 
  theme(
    strip.background = element_rect(fill = "white", color = "black"), 
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=20, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    legend.position = "none") + ylab('')


plot(p2)


pp <- p1 | p2
pp
ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/cyclic_genes_heatmaps_rna_atac_ord_by_peak_expr.png",
       plot=pp,
       width = 6, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


ggsave(filename="../Output_KZ/figures/cyclic_genes_heatmaps_rna_ord_by_peak_expr.png",
       plot=p1,
       width = 3, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

ggsave(filename="../Output_KZ/figures/cyclic_genes_heatmaps_atac_ord_by_peak_expr.png",
       plot=p2,
       width = 3, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## new version
sc.rna.mu.scale$data <- "scRNA"
sc.rna.mu.scale$peak.ord.rna <- sc.rna.mu.scale$peak.ord
sc.atac.mu.scale$data <- "scATAC"
df <- rbind(sc.rna.mu.scale, sc.atac.mu.scale)

df$data <- factor(df$data, levels = c("scRNA", "scATAC"))
p2 <- ggplot(df, aes(x = x, y = GeneID, fill = expr)) + 
  geom_tile() + 
  # facet_grid(phase~., scales = "free",  space='free',
  #            labeller=label_wrap_gen(multi_line = TRUE))+
  ylab("Genes") + xlab("time/cells") +
  #scale_fill_gradientn(colours = hm.palette(10)) +
  facet_grid(. ~ data, scales = "free", space='free', labeller=label_wrap_gen(multi_line = TRUE))+
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(panel.spacing = unit(0.01, "lines")) + 
  theme(
    strip.background = element_rect(fill = "white", color = "white"), 
    panel.spacing = unit(0.01, "lines"),
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y  = element_blank(),
    plot.title = element_text(size=20, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    axis.title.y = element_text(size=20, face="bold"),
    legend.position = "none") + ylab('Genes')
  

plot(p2)

ggsave(filename="../Output/toxo_cdc/ME49_59/figures_paper/cyclic_genes_heatmaps_rna_ord_facet.png",
       plot=p2,
       width = 6, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

