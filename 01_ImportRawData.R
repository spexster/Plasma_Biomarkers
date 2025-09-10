# ======================================================================
# Project: PlasmaAgeing — Single Dataset Loader & Metadata Extractor
# Platform: Imperial College London HPC (Open OnDemand) — Conda Renv
# Purpose: On-disk H5AD read, robust pathing, metadata harmonisation,
#          optional Seurat QC, session logging.
# ======================================================================

suppressPackageStartupMessages({
  library(reticulate)
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Seurat)
  library(hdf5r)
})

# Try binding reticulate to Conda Python (not required for zellkonverter)
try({ use_condaenv("Renv", conda = "~/miniforge3/bin/conda", required = FALSE) }, silent = TRUE)

# Source neutral helpers if present
if (file.exists("scripts/utils_qc.R")) source("scripts/utils_qc.R")

# ---- Package sanity (no SeuratDisk) ----------------------------------
need <- c("kernlab","Seurat","zellkonverter","hdf5r","SummarizedExperiment")
missing <- need[!vapply(need, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (length(missing)) {
  bioc <- c("zellkonverter","SummarizedExperiment")
  cf  <- setdiff(missing, bioc)
  bc  <- intersect(missing, bioc)
  msg <- paste(
    if (length(cf)) paste0("  conda install -n Renv -c conda-forge ", paste0("r-", tolower(cf), collapse=" ")),
    if (length(bc)) paste0("  conda install -n Renv -c bioconda ", paste0("bioconductor-", tolower(bc), collapse=" ")),
    sep = "\n"
  )
  stop("Missing R packages: ", paste(missing, collapse=", "), "\nInstall via:\n", msg)
}

# ---- Paths (no setwd) ------------------------------------------------
data_dir  <- "data"
h5ad_file <- "scRNA-seq_dataset.h5ad"
primary   <- file.path(data_dir, h5ad_file)
fallbacks <- c(primary, file.path("MScProject_data", h5ad_file), h5ad_file,
               file.path("/rds/general/user/fj524/home","MScProject_data", h5ad_file))
found <- NULL; for (p in unique(fallbacks)) if (file.exists(p)) { found <- normalizePath(p); break }
if (is.null(found)) {
  auto <- tryCatch(list.files(pattern=paste0("^", gsub("([.()+^$|])","\\\\\\1", h5ad_file), "$"),
                              recursive=TRUE, full.names=TRUE), error=function(e) character(0))
  if (length(auto)) found <- normalizePath(auto[1])
}
stopifnot("H5AD file not found" = !is.null(found) && file.exists(found))
h5ad_path <- found
cat("=== PATH OK ===\n",
    "File: ", h5ad_path, "\n",
    sprintf("Size: %.1f MB\n", file.size(h5ad_path)/1024^2),
    "Mod : ", format(file.mtime(h5ad_path), usetz=TRUE), "\n", sep="")
try(zellkonverter::readH5ADSummary(h5ad_path), silent=TRUE)

# ---- Load on-disk (HDF5) ---------------------------------------------
Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
DelayedArray::setAutoBlockSize(1e8)
t0 <- Sys.time()
sce <- readH5AD(h5ad_path, use_hdf5 = TRUE)
cat("Loaded SCE: ", nrow(sce), " genes x ", ncol(sce), " cells in ",
    round(difftime(Sys.time(), t0, units="mins"), 2), " min\n", sep="")

# ---- Metadata standardisation ----------------------------------------
meta <- as.data.frame(colData(sce))
col_map <- list(age="^(age|donor_age)$|^Age$", sex="^(sex|gender)$|^Sex$", ethnicity="ethnic|race", smoking_status="smok|smoking")
standardise_cols <- function(df, map){
  for(nm in names(map)){
    hits <- grep(map[[nm]], names(df), ignore.case=TRUE, value=TRUE)
    if(length(hits)) names(df)[match(hits[1], names(df))] <- nm
  }
  df
}
meta_std <- standardise_cols(meta, col_map)
if (exists("metadata_qc")) meta_std <- metadata_qc(meta_std)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
utils::write.table(head(meta_std, 20), file=file.path("results","metadata_preview.tsv"),
                   sep="\t", quote=FALSE, row.names=TRUE)
det <- intersect(c("age","sex","ethnicity","smoking_status"), names(meta_std))
cat("Metadata detected: ", if(length(det)) paste(det, collapse=", ") else "none", "\n", sep="")

# ---- Optional Seurat QC ----------------------------------------------
DO_SEURAT_QC <- FALSE
if (DO_SEURAT_QC) {
  layer <- if ("counts" %in% assayNames(sce)) "counts" else "X"
  seu <- as.Seurat(sce, counts = layer)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
  if (exists("compute_dimred_qc")) seu <- compute_dimred_qc(seu)
  saveRDS(seu, file.path("results","seurat_object.rds"))
  pdf(file.path("results","umap_by_cluster.pdf")); print(DimPlot(seu, group.by="seurat_clusters")); dev.off()
  pdf(file.path("results","qc_violins.pdf")); print(VlnPlot(seu, c("nFeature_RNA","nCount_RNA","percent.mt"))); dev.off()
  cat("Seurat object & QC plots saved under results/.\n")
}

# ---- Session log ------------------------------------------------------
logfile <- file.path("results","session_info.txt")
writeLines(c(
  "=== PlasmaAgeing run ===",
  paste("Working dir:", getwd()),
  paste("H5AD path:", h5ad_path),
  paste("SCE dims:", paste(dim(sce), collapse=" x ")),
  paste("R version:", R.version.string),
  paste("Conda env:", Sys.getenv("CONDA_PREFIX")),
  "\n--- sessionInfo() ---",
  capture.output(print(sessionInfo()))
), con=logfile)
cat("Session log written to ", logfile, "\n", sep="")

