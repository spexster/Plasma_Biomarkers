# 37b_celltype_specificity_from_csv.R
suppressPackageStartupMessages({ library(readr); library(dplyr); library(pheatmap) })

PROJ <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
EDIR <- file.path(PROJ, "results/edgeR_ISsexadjusted")
TDIR <- file.path(PROJ, "results/tables")
FDIR <- file.path(PROJ, "results/figures")
dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FDIR, recursive = TRUE, showWarnings = FALSE)

cell_types <- c(
  "CD4-positive, alpha-beta T cell",
  "CD8-positive, alpha-beta T cell",
  "gamma-delta T cell",
  "regulatory T cell",
  "naive thymus-derived CD4-positive, alpha-beta T cell",
  "naive thymus-derived CD8-positive, alpha-beta T cell",
  "central memory CD4-positive, alpha-beta T cell",
  "effector memory CD4-positive, alpha-beta T cell",
  "central memory CD8-positive, alpha-beta T cell",
  "effector memory CD8-positive, alpha-beta T cell",
  "CD4-positive, alpha-beta cytotoxic T cell",
  "mucosal invariant T cell",
  "double negative thymocyte"
)

stubify <- function(x) gsub("[^A-Za-z0-9]+","_", x)

sig_lists <- list()
for (ct in cell_types) {
  st <- stubify(ct)
  csv <- file.path(EDIR, paste0("DE_age_effect_", st, "_edgeR.csv"))
  
  # Try alternative file patterns if main file doesn't exist
  if (!file.exists(csv)) {
    alt <- Sys.glob(file.path(EDIR, paste0("*", st, "*_edgeR*.csv")))
    if (length(alt)) {
      csv <- alt[1]
      message("Using alternative file: ", basename(csv))
    } else {
      message("No DE CSV for: ", ct)
      sig_lists[[ct]] <- character(0)
      next
    }
  }
  
  # Read the CSV file - FIXED FOR NO HEADER IN FIRST COLUMN
  tab <- tryCatch({
    # First read without header to check structure
    test_lines <- readLines(csv, n = 3)
    has_header <- length(grep(",", test_lines[1])) > 0 && 
      !grepl("^ENS", test_lines[1])  # If first line doesn't start with ENS, it's probably header
    
    if (has_header) {
      # File has header but first column is empty
      df <- read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE, row.names = NULL)
      # The first column without header will get a default name like "X" or be empty
      # Let's find which column contains ENSEMBL IDs
      ensembl_col <- which(sapply(df, function(x) all(grepl("^ENS", x[1:min(10, nrow(df))]))))[1]
      if (is.na(ensembl_col)) ensembl_col <- 1  # Default to first column
      
      # Find FDR column
      fdr_col <- grep("FDR|fdr|adj.P.Val|padj", names(df), ignore.case = TRUE, value = TRUE)[1]
      
      list(df = df, ensembl_col = ensembl_col, fdr_col = fdr_col)
    } else {
      # No header at all - read and assign column names
      df <- read.csv(csv, header = FALSE, stringsAsFactors = FALSE)
      # Assume first column is ENSEMBL, look for FDR column
      fdr_col_index <- which(sapply(df, function(x) any(grepl("FDR|fdr|adj.P.Val|padj", x[1:min(3, nrow(df))], ignore.case = TRUE))))[1]
      if (is.na(fdr_col_index)) fdr_col_index <- which(sapply(df, function(x) all(x <= 1 & x >= 0, na.rm = TRUE)))[1] # heuristic for p-value column
      
      list(df = df, ensembl_col = 1, fdr_col = fdr_col_index)
    }
  }, error = function(e) {
    message("Read failed: ", csv, " - Error: ", e$message)
    NULL
  })
  
  if (is.null(tab)) {
    sig_lists[[ct]] <- character(0)
    next
  }
  
  df <- tab$df
  ensembl_col <- tab$ensembl_col
  fdr_col <- tab$fdr_col
  
  message("File: ", basename(csv))
  message("  Ensembl column: ", ifelse(is.numeric(ensembl_col), ensembl_col, "unknown"))
  message("  FDR column: ", ifelse(is.numeric(fdr_col), fdr_col, fdr_col))
  
  # Extract Ensembl IDs and FDR values
  ens <- if (is.numeric(ensembl_col)) df[[ensembl_col]] else df[[1]]  # default to first column
  fdr_vals <- if (is.numeric(fdr_col)) {
    as.numeric(df[[fdr_col]])
  } else if (!is.null(fdr_col) && fdr_col %in% names(df)) {
    as.numeric(df[[fdr_col]])
  } else {
    message("  Could not find FDR column, using default column 2")
    as.numeric(df[[2]])  # default to second column
  }
  
  # Get significant genes (FDR < 0.05)
  sig <- ens[!is.na(fdr_vals) & fdr_vals < 0.05 & !is.na(ens) & grepl("^ENS", ens)]
  sig_lists[[ct]] <- unique(na.omit(sig))
  
  message("  Found ", length(sig_lists[[ct]]), " significant genes for ", ct)
}

# union of all genes that were significant somewhere
universe <- sort(unique(unlist(sig_lists)))
mat <- matrix(0L, nrow = length(universe), ncol = length(cell_types),
              dimnames = list(universe, cell_types))

for (j in seq_along(cell_types)) {
  g <- sig_lists[[cell_types[j]]]
  if (length(g)) {
    common_genes <- intersect(g, universe)
    if (length(common_genes) > 0) {
      mat[common_genes, j] <- 1L
    }
  }
}

# write tables + figure
write.csv(cbind(gene = rownames(mat), as.data.frame(mat)),
          file.path(TDIR, "age_sig_specificity_matrix.csv"),
          row.names = FALSE)

present_counts <- data.frame(cell_type = cell_types,
                             n_sig = colSums(mat), row.names = NULL)
write.csv(present_counts, file.path(TDIR, "age_sig_present_counts.csv"), row.names = FALSE)

unique_counts <- data.frame(cell_type = cell_types,
                            n_unique = colSums(mat == 1 & rowSums(mat) == 1), row.names = NULL)
write.csv(unique_counts, file.path(TDIR, "age_sig_unique_counts.csv"), row.names = FALSE)

png(file.path(FDIR, "age_sig_specificity_heatmap_adjSex.png"), width = 1800, height = 1400, res = 200)
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Genes significant for age (FDR<0.05) across T-cell subsets",
         show_rownames = FALSE)
dev.off()

message("Done. Outputs:\n - ", file.path(TDIR, "age_sig_specificity_matrix.csv"),
        "\n - ", file.path(TDIR, "age_sig_present_counts.csv"),
        "\n - ", file.path(TDIR, "age_sig_unique_counts.csv"),
        "\n - ", file.path(FDIR, "age_sig_specificity_heatmap_adjSex.png"))


###########################RUN AFTER ABOVE WHEN NEEDED##########################
# # Unique_genes_proportion_bar.png
# library(ggplot2)
# library(dplyr)
# library(tidyr)
# 
# # Read data
# counts <- read.csv("~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/tables/age_sig_present_counts.csv")
# unique_counts <- read.csv("~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/tables/age_sig_unique_counts.csv")
# 
# # Calculate proportions
# prop_data <- counts %>%
#   left_join(unique_counts, by = "cell_type") %>%
#   rename(total_genes = n_sig, unique_genes = n_unique) %>%
#   mutate(shared_genes = total_genes - unique_genes,
#          prop_unique = unique_genes / total_genes * 100,
#          prop_shared = shared_genes / total_genes * 100) %>%
#   arrange(desc(prop_unique))
# 
# # Create ordered bar plot
# p1 <- ggplot(prop_data, aes(x = reorder(cell_type, prop_unique), y = prop_unique)) +
#   geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Proportion of cell-type specific age-associated genes",
#        x = "Cell Type", 
#        y = "Percentage of unique genes (%)") +
#   theme_minimal() +
#   theme(axis.text.y = element_text(size = 10),
#         plot.title = element_text(hjust = 0.5))
# 
# ggsave("~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/figures/unique_genes_proportion_bar.png", 
#        p1, width = 10, height = 8, dpi = 300)

################################################################################





#old version####################################################################
# # 37b_celltype_specificity_from_csv.R
# suppressPackageStartupMessages({ library(readr); library(dplyr); library(pheatmap) })
# 
# PROJ <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
# # EDIR <- file.path(PROJ, "results/edgeR")
# EDIR <- file.path(PROJ, "results/edgeR_ISsexadjusted")
# TDIR <- file.path(PROJ, "results/tables")
# FDIR <- file.path(PROJ, "results/figures")
# dir.create(TDIR, recursive = TRUE, showWarnings = FALSE)
# dir.create(FDIR, recursive = TRUE, showWarnings = FALSE)
# 
# cell_types <- c(
#   "CD4-positive, alpha-beta T cell",
#   "CD8-positive, alpha-beta T cell",
#   "gamma-delta T cell",
#   "regulatory T cell",
#   "naive thymus-derived CD4-positive, alpha-beta T cell",
#   "naive thymus-derived CD8-positive, alpha-beta T cell",
#   "central memory CD4-positive, alpha-beta T cell",
#   "effector memory CD4-positive, alpha-beta T cell",
#   "central memory CD8-positive, alpha-beta T cell",
#   "effector memory CD8-positive, alpha-beta T cell",
#   "CD4-positive, alpha-beta cytotoxic T cell",
#   "mucosal invariant T cell",
#   "double negative thymocyte"
# )
# 
# stubify <- function(x) gsub("[^A-Za-z0-9]+","_", x)
# 
# sig_lists <- list()
# for (ct in cell_types) {
#   st  <- stubify(ct)
#   # csv <- file.path(EDIR, paste0(st, "_edgeR.csv"))
#   csv <- file.path(EDIR, paste0("DE_age_effect_", st, "_edgeR.csv"))
#   if (!file.exists(csv)) {
#     alt <- Sys.glob(file.path(EDIR, paste0("*", st, "*_edgeR.csv")))
#     if (length(alt)) csv <- alt[1] else { message("No DE CSV for: ", ct); sig_lists[[ct]] <- character(0); next }
#   }
#   # tab <- tryCatch(read.csv(csv, check.names = FALSE), error = function(e) NULL)
#   # tab <- tryCatch(read.csv(csv, header = TRUE, row.names = 1, check.names = FALSE), error = function(e) NULL)
#   tab <- tryCatch(read.csv(csv, header = TRUE, row.names = 1, check.names = FALSE), error = function(e) NULL)
#   if (is.null(tab)) { message("Read failed: ", csv); sig_lists[[ct]] <- character(0); next }
#   
#   ens_col <- grep("^ensembl$", names(tab), ignore.case = TRUE, value = TRUE)
#   if (!length(ens_col)) { message("No Ensembl column in: ", csv); sig_lists[[ct]] <- character(0); next }
#   ens <- tab[[ens_col[1]]]
#   
#   sig <- ens[ tab$FDR < 0.05 ]            # presence regardless of sign
#   sig_lists[[ct]] <- unique(na.omit(sig))
# }
# 
# # union of all genes that were significant somewhere
# universe <- sort(unique(unlist(sig_lists)))
# mat <- matrix(0L, nrow = length(universe), ncol = length(cell_types),
#               dimnames = list(universe, cell_types))
# for (j in seq_along(cell_types)) {
#   g <- sig_lists[[ cell_types[j] ]]
#   if (length(g)) mat[g, j] <- 1L
# }
# 
# # write tables + figure
# write.csv(cbind(gene = rownames(mat), as.data.frame(mat)),
#           file.path(TDIR, "age_sig_specificity_matrix.csv"),
#           row.names = FALSE)
# 
# present_counts <- data.frame(cell_type = cell_types,
#                              n_sig = colSums(mat), row.names = NULL)
# write.csv(present_counts, file.path(TDIR, "age_sig_present_counts.csv"), row.names = FALSE)
# 
# unique_counts <- data.frame(cell_type = cell_types,
#                             n_unique = colSums(mat == 1 & rowSums(mat) == 1), row.names = NULL)
# write.csv(unique_counts, file.path(TDIR, "age_sig_unique_counts.csv"), row.names = FALSE)
# 
# png(file.path(FDIR, "age_sig_specificity_heatmap_adjSex.png"), width = 1800, height = 1400, res = 200)
# pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE,
#          main = "Genes significant for age (FDR<0.05) across T-cell subsets",
#          show_rownames = FALSE)
# dev.off()
# 
# message("Done. Outputs:\n - ", file.path(TDIR,"age_sig_specificity_matrix.csv"),
#         "\n - ", file.path(TDIR,"age_sig_present_counts.csv"),
#         "\n - ", file.path(TDIR,"age_sig_unique_counts.csv"),
#         "\n - ", file.path(FDIR,"age_sig_specificity_heatmap.png"))
