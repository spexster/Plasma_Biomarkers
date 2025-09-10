
# title: "seurat preProcessing"
# date: "2025-08-19"

# need to apply filtering for low expressed genes and for cells that express low number of genes, too many cells AP: sure too many zeros. e.g. cells <x genes not to keep. select threshold for cells and variants of genes. Then start pre-processing. 
# Adapted from Bioinformatics lectures


library(ggplot2)



set.seed(15)

library(dplyr)
all.genes <- rownames(seu_tcells)

# seu_tcells <- seu_tcells %>%
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
#   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%  # these are the default
#   ScaleData(features = all.genes) %>%  
#   RunPCA()


seu_tcells <- NormalizeData(seu_tcells, normalization.method = "LogNormalize", scale.factor = 10000)

seu_tcells <- FindVariableFeatures(seu_tcells, selection.method = "vst", nfeatures = 2000)

seu_tcells <- ScaleData(seu_tcells, features = all.genes)

seu_tcells <- RunPCA(seu_tcells)   


seu_full <- NormalizeData(seu_full, normalization.method = "LogNormalize", scale.factor = 10000)

seu_full <- FindVariableFeatures(seu_full, selection.method = "vst", nfeatures = 2000)

seu_full <- ScaleData(seu_full, features = all.genes)

seu_full <- RunPCA(seu_full) 



set.seed(15)

seu_tcells <- seu_tcells %>%
  FindNeighbors(dims = 1:20) 

set.seed(15)
seu_tcells <- seu_tcells %>%
  FindClusters(resolution = 1.5) 

set.seed(15)
seu_tcells <- seu_tcells %>%
  RunUMAP(dims = 1:20)





# plots


# dimplot with argument group.by=agegroups

DimPlot(seu_tcells, reduction = "umap")

ggsave("results/DimPlot_ageGroups.pdf")






# save your final seurat object

saveRDS(seu_tcells, file = "../path_to_your_folder/seu_tcells.rds")






# Compare the gene expression for each cluster vs all others, then combines results
# Instead of the cluster labels, in your case you will use the age group levels





all_markers <- FindAllMarkers(seu_tcells,  
                              assay = "originalexp", # put the name of the assay in your seuray object
                              # select the column you want to group by
                              group.by = "AgeGroups",
                              #  min.pct = 0, Include all genes regardless of expression percentage
                              min.pct = 0, 
                              # Statistical test used for finding markers ("wilcox" by default)
                              test.use = "wilcox",           
                              # Include all genes that pass the following log fold change threshold
                              logfc.threshold = 1, 
                              # Keep ONLY the UPREGULATED genes!
                              only.pos = TRUE,      
                              max.cells.per.ident = Inf,
                              random.seed = 15)

dim(all_markers) 

# save this result as a .csv file and then you can right-click on it and open it on an .xslx file to inspect

write.csv(all_markers, "path_to_your_folder/all_markers.csv")







# Dimensions of datasets = number of rows and number of columns (1st is nos rows and 2nd nos is number columns)
dim(seu_tcells@assays[["originalexp"]]@counts)

dim(seu_full@assays[["originalexp"]]@counts)

