
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)

data_object_filtered@meta.data$type <- factor(data_object_filtered@meta.data$type,levels = c("CTRL","CSN","CDK2","DOUBLE"))
a <- DimPlot(data_object_filtered, reduction = "umap", group.by = "sub_type",split.by = "type",raster=TRUE) + theme_b110()

ggsave("replicate_treatment_umaps.png",a,width = 12,height = 4)

