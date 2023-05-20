

prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c("GeneID" = "TGGT1")) %>% na.omit()
MJ_annot <- read.xlsx("../Input/Toxo_genomics/genes/MJ_annotation.xlsx")
MJ_annot <- MJ_annot %>% dplyr::select(!Product.Description)
prod.desc <- left_join(prod.desc, MJ_annot, by= "TGME49" )
prod.desc <- prod.desc %>% 
  mutate(Category = case_when(str_detect(pattern = "ribosomal" ,ProductDescription) ~ "ribosomal",
                              TRUE ~ "others"))
prod.desc <- prod.desc %>% dplyr::select(-c(TGGT1,source_id ))


CutRun.all.info <- readRDS("../Input/toxo_cdc/rds_ME49_59/Union_all_new_peaks_0.05_qval_NO_frag_filt_intersect_info_motif_info_motif_1_motif2.rds" )


CutRun.all.info <- CutRun.all.info %>% 
  dplyr::select(V1.x, start_peak, end_peak, V4, V5, V11,
                gene_name,cutRun.peaks, intersection,has.motif,  motif) %>% distinct()

colnames(CutRun.all.info) <- gsub("gene_name", "TGME49", colnames(CutRun.all.info))

colnames(CutRun.all.info) <- c("chr", "start_peak", "end_peak", "V4", "V5", "V6", "TGME49", 
                               "assigned_to_CutRun_peaks", "intersection_CutRun_dataSets","has.motif", "motif")

modulated <- readRDS("../Input/toxo_cdc/rds_ME49_59/AP2XII8_KD_modulated_genes_all_comparisons.rds")
modulated <- modulated %>% dplyr::select(GeneID, Global_KD_vs_WT, KD_vs_WT_phase_based) %>% distinct()
colnames(modulated)[1] <- "TGME49"


tmp <- full_join(CutRun.all.info, modulated, by ="TGME49")
tmp2 <- left_join(tmp, prod.desc, by ="TGME49")

saveRDS(tmp2, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")
write.xlsx(tmp2, "../OutPut/toxo_cdc/ME49_59/tables/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.xlsx")
write.xlsx(tmp2, "../OutPut/toxo_cdc/ME49_59/tables/Supplements/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.xlsx")
saveRDS(tmp2, "../Input/toxo_cdc/rds_ME49_59/cut_run_union_new_peaks_march_motif_modulated_genes_lfs_v5.rds")

