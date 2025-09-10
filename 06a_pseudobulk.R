### Pseudobulking tcells
## not used
## pseudobulked with script 26


## pseudobulk by single cell-type
library(Matrix)
cell_type_interest="naive thymus-derived CD4-positive, alpha-beta T cell"
subset.tmp=subset(seu_full_withmeta_latest, subset=cell_type==cell_type_interest)

counts = GetAssayData(seu_full_withmeta_latest, slot="counts")
length(which(counts > 0))
sample_ids = subset.tmp$donor_id

# Create pseudobulk matrix
pseudobulk_list <- split(seq_len(ncol(counts)), sample_ids)

pseudobulk_counts <- sapply(pseudobulk_list, function(cells) {
  rowSums(counts[, cells, drop = FALSE])
})

# Convert to matrix
pseudobulk_matrix <- as.matrix(pseudobulk_counts)

# Final check of the pseudobulk matrix dimensions
print(dim(pseudobulk_matrix))

# Save results
saveRDS(pseudobulk_matrix, file = "results/pseudobulk_naive_thymus_CD4.rds")


##### in separate R session on HPC using powershell
library(edgeR)
pseudobulk_matrix=readRDS(file="pseudobulk_naive_thymus_CD4.rds")



## Automated pseudobulking - all tcells

# Ensure necessary libraries are loaded
library(Seurat)
library(Matrix)

# 1. Start from your full T-cell Seurat object (e.g., 'seu_tcells')
# If you need to subset from the full object first, uncomment the next line:
# seu_tcells <- subset(seu_full, predicted.celltype.l2 %in% c("CD4 Naive", "CD8 Naive", "Treg", ...)) # Add all your T-cell types

# 2. Get all unique T-cell types in your dataset
cell_types <- unique(seu_tcells$cell_type) # Replace 'cell_type' with your actual metadata column name
print("Found these cell types:")
print(cell_types)

# 3. Create a directory to save all the pseudobulk results
output_dir <- "results/pseudobulk_tcells/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 4. Loop through each cell type
for (cell_type_interest in cell_types) {
  
  message("\nProcessing: ", cell_type_interest)
  
  # Subset to the current cell type
  subset.tmp <- subset(seu_tcells, subset = cell_type == cell_type_interest) # Adjust column name if needed
  
  # Check if we have enough cells and donors to proceed
  if (ncol(subset.tmp) < 10) { # Arbitrary minimum cell count threshold
    message("  -> Skipping. Only ", ncol(subset.tmp), " cells found.")
    next
  }
  
  num_donors <- length(unique(subset.tmp$donor_id))
  if (num_donors < 3) { # Arbitrary minimum donor threshold for DE
    message("  -> Skipping. Only ", num_donors, " donors found.")
    next
  }
  
  # Proceed with pseudobulk creation
  counts <- GetAssayData(subset.tmp, slot = "counts")
  sample_ids <- subset.tmp$donor_id
  
  pseudobulk_list <- split(seq_len(ncol(counts)), sample_ids)
  pseudobulk_counts <- sapply(pseudobulk_list, function(cells) {
    rowSums(counts[, cells, drop = FALSE])
  })
  
  pseudobulk_matrix <- as.matrix(pseudobulk_counts)
  
  # 5. Create a clean filename and save the result
  # Remove problematic characters from filename
  clean_name <- gsub("[^[:alnum:]]", "_", cell_type_interest)
  output_file <- file.path(output_dir, paste0("pseudobulk_", clean_name, ".rds"))
  
  saveRDS(pseudobulk_matrix, file = output_file)
  message("  -> Saved pseudobulk matrix for ", cell_type_interest, " to:\n    ", output_file)
}

message("\nPseudobulk analysis complete for all T-cell types!")