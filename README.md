# Toxoplasma_CDC

# Toxoplasma_CDC

prepScOutPut.R - process the output of cellranger

1. readProcScRNA.R - read and process raw rna expression (counts) 

2. readProcScATAC.R - read and process raw atac profile (counts) 

3. peak_gene_assignment_scATAC.R 

4. marker_analysis.R

  perform marker analysis on the following: 
  
  1. rna markers using phase as Idents
  2. rna markers using rna transition points as Idents
  3. rna markers using atac transition points as Idesnts
  
  note: Intersection of 2 & 3 is in this code
  
  
7. fitPseudoTime_scATAC_scRNA_V2.R - fit pseudo-time to scRNA and scATAC    

5. sc_rna_sc_atac_fit_smoothing_splines.R - fitting smoothing splines to scRNA and scATAC 

6. scRNA_scATAC_correlation_analysis.R - correlation analysis between expr and atac curves 


8. map_transition_final.R - Infer transition in peak expression and peak accessibility during  


AP2XII-8 KD 

9. readProcScRNA_AP2XII8_transfer_labels.R - Read and process AP2XII-8 scRNA KD

10. integration_V2.R - integrates WT and KD seurat objects


12. all types of DEGs between KD and WT  - integratedMarkerAnalysis.R
  
  note: Here I used the integrated objects because for finding markers seurat function only 
  accepts one single object. However, when we set the assay to RNA it ignores integration. 

11. Cut&RUN on AP2XII-8 KD - CutRun_MArch_automated_No_frag_filt.R   

  peak_gene assignment on individual narrow peaks 
  peak_gene assignment on union of 4 data 
  moif occurence 
  public chip information
  
12. Cut_run_March_exploratory_final.R - Cut&run and DEGs (KD_vs_WT) exploratory analysis 

  note: all figures and tables for downstream analysis of cut&run and DEGs(KD_vs_WT) has been    done in this script. PowerPoint CutRUN_DEGs_V4_05_05_23 analysis 


13. Dynamic time Warping clustering of gene sets - dtw_clustering.R
  This script takes genes and cluster them and plot the rna and atac profiles according to the   cluster. 
  
14. all_expression_atac_profiles.R 

  If we dont want to cluster the expression or accessibility curves, we can use this script to   plot the expression and accessibility of a gene set. 

15. GO_enrichment_new.R

  takes the output of toxodb enrichment analsys process the tables and plot the top rankes GO    terms. 
  
16. sc_expr_plot.R - plots rna expression of KD and WT + violin plots for each each gene 

17. heatmaps_trans.R - 
  generates heatmap of expression of cell cycle regulated genes ordered by rna peak time

18. heatmaps_trans_order_by_atac.R 
  generates heatmap of expression of cell cycle regulated genes ordered by atac peak time


19. AP2_clustering.R

  Clusters cyclic AP2s (new from the review paper)
  
20. IMCs_clusters.R, BCs_clustersR 
  
  cluster new list of BC and IMC genes 
  

21. circosPlot_app.r
  plots the circos for gene families
  the rds file has been generated in circosPlot.r

