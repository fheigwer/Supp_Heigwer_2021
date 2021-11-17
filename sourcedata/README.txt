###################################
# Figshare README
###################################

filtered processed single cell RNAseq data are deposit in

  R-readable (requires R/Seurat v>3.6.2 ):
    data_object_filtered_together-PC1.rds
  machine readable csv (filtered genes x filtered cells)
    scRNA_raw_count_matrix.csv.gz (raw count integer)
    scRNA_normalized_count_matrix.csv.gz (normalized and scaled expression double)
  machine readable meta data (cells x attributes)
    meta_data.csv.gz

filtered processed well-average interaction and morphologcal screening data_object_filtered_together
  machine readable csv in gzipped format of all bias correct interactions for each replicate for each gene for each selected feature
    interactions_stat_tested_bias_corrected.csv.gz
  machine readabl csv gzipped normed and corrected feature values
    morphological_feature_dataframe_normed_corrected.csv.gz
