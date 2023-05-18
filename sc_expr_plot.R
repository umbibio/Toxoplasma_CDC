
## plot expression and vln plots 


plot_expr_vln <- function(S.O.wt, S.O.kd, geneID, geneNAME, ident1 = "phase", ident2 = "seurat_clusters", redc ){
  
  geneID <- gsub("_", "-", geneID)
  Idents(S.O.wt) <- ident1
  expr.plt.wt <- FeaturePlot(S.O.wt, features = geneID, reduction = redc) + 
    ggtitle(paste("scRNA", geneID, sep = "-"))
  
  Idents(S.O.kd) <- ident2
  expr.plt.kd <- FeaturePlot(S.O.kd, features = geneID, reduction = redc) + ggtitle("scRNA_KD")
  p1 <- expr.plt.wt | expr.plt.kd
  
  vln.plt.wt <- VlnPlot(S.O.wt, features = geneID) + ggtitle(paste("scRNA", geneID, sep = "-"))
  vln.plt.kd <- VlnPlot(S.O.kd, features = geneID) + ggtitle("scRNA_KD")
  p2 <-   vln.plt.wt | vln.plt.kd
  
  pp <- p1 / p2 
  plot(pp, main = 'title2')
  
}


CutRun.modulated.lfS <- readRDS("../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v4.rds")
S.O.rna.WT <- readRDS("../Input/toxo_cdc/rds_ME49_59/S.O.rna.WT_labels.rds")
S.O.rna.KD <- readRDS('../Input/toxo_cdc/rds_ME49_59/S.O.rna.AP2XII8.KD.new_transferred_lables_bootroyed.rds')

colnames(CutRun.modulated.lfS) <- c("chr", "start_peak", "end_peak")


# KD:S1S2 vs WT:S
gene <- "TGME49_280460" # AP2VIIa-2
gene <- "TGME49_205650" # AP2VIIa-3
gene <- "TGME49_249190" # AP2XII-6
gene <- "TGME49_249190" # AP2VIII-7
gene <- "TGME49_263090" # HDAC3
gene <- "TGME49_305340" # morc

# KD vs WT (phase based)
gene <- "TGME49_250800" # AP2XII-8
gene <- "TGME49_215895" # AP2IX-10

# low level expression - recruiting HDAC3 and morc - switch to merozoizte
gene <- "TGME49_310900" # AP2XI-2
gene <- "TGME49_218960" # AP2XII-1
gene <- "TGME49_315760" #AP2XI-4

# supported both by cut and run and modulated phase based KD vs WT
gene <- "TGME49_280460" # AP2VIIa-2
gene <- "TGME49_269010" # AP2VIII-7
gene <- "TGME49_214840" # AP2X-7

gene <- "TGME49_227600"

gene <- "TGME49_204530"  
plot_expr_vln(S.O.rna.WT, S.O.rna.KD, geneID = gene, redc = "pca")



