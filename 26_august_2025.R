###############################################################################
# Pseudobulk (T-cell type) – working version
# What you change each run: only `cell_type_interest`
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

## ---- Config / paths ---------------------------------------------------------
PROJ   <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
OUTDIR <- "~/scripts"  # where you want the .txt file written
dir.create(path.expand(OUTDIR), showWarnings = FALSE, recursive = TRUE)

## counts matrix (sparse) – try a few sensible locations
count_paths <- c(
  file.path(PROJ, "results/count_matrix_allcells_.rds"),
  "~/results/count_matrix_allcells_.rds",
  "results/count_matrix_allcells_.rds"
)
counts_path <- count_paths[file.exists(path.expand(count_paths))][1]
if (is.na(counts_path)) stop("Could not find counts RDS at any known path.")
counts <- readRDS(path.expand(counts_path))

## Seurat object with metadata
seu_path <- file.path(PROJ, "results/seurat_snapshots/seu_full_withmeta_latest.rds")
if (!file.exists(path.expand(seu_path))) {
  stop("Missing Seurat snapshot: ", seu_path)
}
seu_full_withmeta_latest <- readRDS(path.expand(seu_path))

## ---- 1) CHOOSE THE CELL TYPE OF INTEREST -----------------------------------
# Change this line each time; run the whole file after editing this string:
cell_type_interest <- "double negative thymocyte"
# Examples you can paste over:
# "CD4-positive, alpha-beta T cell"
# "CD8-positive, alpha-beta T cell"
# "gamma-delta T cell"
# "regulatory T cell"
# "naive thymus-derived CD4-positive, alpha-beta T cell"
# "naive thymus-derived CD8-positive, alpha-beta T cell"
# "central memory CD4-positive, alpha-beta T cell"
# "effector memory CD4-positive, alpha-beta T cell"
# "central memory CD8-positive, alpha-beta T cell"
# "effector memory CD8-positive, alpha-beta T cell"
# "CD4-positive, alpha-beta cytotoxic T cell"
# "mucosal invariant T cell"   # (MAIT)
# "double negative thymocyte"

## ---- 2) Robust selection of target cells -----------------------------------
# MAIT and γδT can be labeled differently across columns; handle those flexibly.
is_target <- if (grepl("mucosal.*invariant.*T.*cell|^MAIT$", cell_type_interest, ignore.case = TRUE)) {
  grepl("mucosal.*invariant.*T.*cell", seu_full_withmeta_latest$cell_type, ignore.case = TRUE) |
    grepl("^MAIT$", seu_full_withmeta_latest$predicted.celltype.l2, ignore.case = TRUE)
} else if (grepl("gamma.*delta|γδ|^gdT$", cell_type_interest, ignore.case = TRUE)) {
  grepl("gamma.*delta|γδ", seu_full_withmeta_latest$cell_type, ignore.case = TRUE) |
    grepl("^gdT$", seu_full_withmeta_latest$predicted.celltype.l2, ignore.case = TRUE)
} else {
  seu_full_withmeta_latest$cell_type == cell_type_interest
}

barcodes_target <- colnames(seu_full_withmeta_latest)[is_target]
message("Selected ", length(barcodes_target), " cells for: ", cell_type_interest)
if (length(barcodes_target) == 0L) stop("No cells matched; check the spelling of cell_type_interest.")

subset.tmp <- subset(seu_full_withmeta_latest, cells = barcodes_target)

## ---- 3) Pseudobulk by donor -------------------------------------------------
sample_ids  <- subset.tmp$donor_id
sample_uniq <- unique(as.character(sample_ids))

# Fast sparse-sum per donor
# Fast sparse-sum per donor (drop-in replacement)
tmp <- sapply(sample_uniq, function(s) {
  cs <- names(sample_ids)[sample_ids == s]          # barcodes for donor s
  Matrix::rowSums(counts[, cs, drop = FALSE])       # fast on sparse matrices
})
colnames(tmp) <- sample_uniq


## ---- 4) Save to your requested folder (filename without spaces) -------------
safe_name <- gsub("[^A-Za-z0-9]+", "_", cell_type_interest)
outfile   <- file.path(OUTDIR, paste0(safe_name, "_count.txt"))

write.table(
  tmp,                # Changed from pseudobulk_counts to tmp
  file      = outfile,
  sep       = "\t",
  quote     = FALSE,
  col.names = NA
)

message("Wrote: ", outfile)
message(nrow(tmp), " genes x ", ncol(tmp), " donors")  # Changed from pseudobulk_counts to tmp


###############################################################################
## === Lightweight saves (replace save.image) ===
# targeted saves. Drop at the end of each run (replaces save.image(...)):
###############################################################################

# outdir <- "/rds/general/user/fj524/home/scripts"
# dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# 
# safe <- gsub("[^A-Za-z0-9]+","_", cell_type_interest)
# 
# rdsfile <- file.path(outdir, paste0("pseudobulk_", safe, ".rds"))
# tsvfile <- file.path(outdir, paste0("pseudobulk_", safe, ".tsv"))
# sumfile <- file.path(outdir, paste0("pseudobulk_", safe, "_summary.txt"))
# donfile <- file.path(outdir, paste0("pseudobulk_", safe, "_donors.txt"))

## 1) Save the matrix (compact, reproducible)
# saveRDS(tmp, rdsfile)

## 2) Also save a human-readable TSV
# write.table(tmp, tsvfile, sep = "\t", quote = FALSE, col.names = NA)

## 3) Tiny manifest for provenance
# writeLines(c(
#   paste("cell_type_interest:", cell_type_interest),
#   paste("n_genes:", nrow(tmp)),
#   paste("n_donors:", ncol(tmp)),
#   paste("n_cells_selected:", ncol(subset.tmp)),
#   paste("timestamp:", as.character(Sys.time())),
#   "counts_source: results/count_matrix_allcells_.rds"
# ), sumfile)

## 4) Donor IDs used
# writeLines(colnames(tmp), donfile)
# 
# message("Wrote:\n  ", rdsfile, "\n  ", tsvfile, "\n  ", sumfile, "\n  ", donfile)




###############################################################################
# ------------------------- ORIGINAL CODE (kept for record) --------------------
# Everything below is earlier version, commented out.
###############################################################################
# counts=readRDS(file = "results/count_matrix_allcells_.rds") 
# #Load object: seu_full_withmeta_latest
#
# cell_type_interest="naive thymus-derived CD4-positive, alpha-beta T cell"
#
# # Subset the object
# subset.tmp=subset(seu_full_withmeta_latest, subset=cell_type==cell_type_interest)
#
# sample_ids = subset.tmp$donor_id
#
# sample_uniq=unique(sample_ids)
#
# #tmp=apply(counts[,names(sample_ids[which(sample_ids=="2_2")])],1,function(x){sum(x)})
#
# #save.image(file="26_august__2025.rds")
#
# pseudobulk=function(sample){
#   apply(counts[,names(sample_ids[which(sample_ids==sample)])],1,function(x){sum(x)})
# }
#
# tmp=apply(as.matrix(sample_uniq),1,pseudobulk)
# colnames(tmp)=sample_uniq
# cell_file_name="naive_thymus-derived_CD4-positive_alpha-beta_T_cell_count.txt"
# write.table(file=cell_file_name,tmp,sep="\t")
# #save.image(file="26_august__2025.rds")
###############################################################################



###############################################################################
# ------------------------- ORIGINAL ADJUSTED CODE (kept for record) --------------------
# Everything below is earlier version, commented out, with first set adjustments
###############################################################################



# # counts=readRDS(file = "results/count_matrix_allcells_.rds") 
# #Load object: seu_full_withmeta_latest
# 
# 
# # 1. DEFINE THE CELL TYPE OF INTEREST (ONLY CHANGE THIS LINE)
# # cell_type_interest <- "mucosal invariant T cell"  # <<< REPLACE THIS STRING FOR EACH NEW TYPE
# ### cell_type_interest="naive thymus-derived CD4-positive, alpha-beta T cell"
# 
# 
# 
# # 2. ROBUST CELL SELECTION LOGIC (RUN THIS EVERY TIME)
# is_target <- if (grepl("mucosal.*invariant.*T.*cell|^MAIT$", cell_type_interest, ignore.case = TRUE)) {
#   # Logic for MAIT cells
#   grepl("mucosal.*invariant.*T.*cell", seu_full_withmeta_latest$cell_type, ignore.case = TRUE) |
#     grepl("^MAIT$", seu_full_withmeta_latest$predicted.celltype.l2, ignore.case = TRUE)
# } else if (grepl("gamma.*delta|γδ|^gdT$", cell_type_interest, ignore.case = TRUE)) {
#   # Logic for gamma-delta T cells
#   grepl("gamma.*delta|γδ", seu_full_withmeta_latest$cell_type, ignore.case = TRUE) |
#     grepl("^gdT$", seu_full_withmeta_latest$predicted.celltype.l2, ignore.case = TRUE)
# } else {
#   # Logic for all other T-cell types: exact match
#   seu_full_withmeta_latest$cell_type == cell_type_interest
# }
# 
# 
# # 3. Subset the object
# barcodes_target <- colnames(seu_full_withmeta_latest)[is_target]
# subset.tmp <- subset(seu_full_withmeta_latest, cells = barcodes_target)
# message("Selected ", length(barcodes_target), " cells for: ", cell_type_interest)
# 
# 
# 
# ### subset.tmp=subset(seu_full_withmeta_latest, subset=cell_type==cell_type_interest)
# 
# ### sample_ids = subset.tmp$donor_id
# 
# ### sample_uniq=unique(sample_ids)
# 
# 
# 
# #tmp=apply(counts[,names(sample_ids[which(sample_ids=="2_2")])],1,function(x){sum(x)})
# 
# #save.image(file="26_august__2025.rds")
# 
# pseudobulk=function(sample){
#   apply(counts[,names(sample_ids[which(sample_ids==sample)])],1,function(x){sum(x)})
# }
# 
# tmp=apply(as.matrix(sample_uniq),1,pseudobulk)
# colnames(tmp)=sample_uniq
# 
# #filename
# ### cell_file_name="naive_thymus-derived_CD4-positive_alpha-beta_T_cell_count.txt"
# cell_file_name <- paste0(gsub("[^A-Za-z0-9]+","_", cell_type_interest), "_count.txt")
# 
# write.table(file=cell_file_name,tmp,sep="\t")
# 
# #save.image(file="26_august__2025.rds")
