## 33_visuals_qc_and_composition.R
## Quick visuals (UMAP), QC tables, and composition~age analysis for T cells

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
})

## --------- Config ---------
OUTDIR <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results"
FIGDIR <- file.path(OUTDIR, "figures")
TABDIR <- file.path(OUTDIR, "tables")
EDGER_DIR <- file.path(OUTDIR, "edgeR")  # where your DE csvs live (optional volcano)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

## --------- Load Seurat object (if needed) ---------
## FJ: ONLY DO IF NOT IN MEMORY: 'seu_full_withmeta_latest'
## Otherwise, point to your snapshot path:
# seu_full_withmeta_latest <- readRDS(
#   "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/seurat_snapshots/seu_full_withmeta_latest.rds"
# )

stopifnot(exists("seu_full_withmeta_latest"))

## --------- Define T-cell subsets (13) ---------
t_cell_types <- c(
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

## Subset the big object to T cells
seu_t <- subset(seu_full_withmeta_latest, subset = cell_type %in% t_cell_types)
cat("T-cell subset:", ncol(seu_t), "cells across", length(unique(seu_t$donor_id)), "donors\n")

## --- Ensure AgeGroups exists on this object ---
if (!"AgeGroups" %in% colnames(seu_t@meta.data)) {
  if (!"age" %in% colnames(seu_t@meta.data)) {
    stop("No 'age' column in meta.data – cannot bin into AgeGroups.")
  }
  age_num <- suppressWarnings(as.numeric(seu_t$age))
  if (all(is.na(age_num))) stop("'age' is not numeric; please clean/coerce it.")
  
  qs <- quantile(age_num, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, names = FALSE)
  ## guard against identical breakpoints
  qs[1] <- floor(qs[1]); qs[length(qs)] <- ceiling(qs[length(qs)])
  qs <- unique(qs)  # collapse duplicates if any
  if (length(qs) < 2L) stop("Age quantiles collapsed; not enough variation to bin.")
  labs <- paste0("[", floor(qs[-length(qs)]), "-", ceiling(qs[-1]), "]")
  
  seu_t[["AgeGroups"]] <- cut(
    age_num,
    breaks = qs,
    include.lowest = TRUE,
    labels = labs,
    right = TRUE
  )
  rm(age_num, qs, labs)
}

## --------- Create AgeGroups (quartiles, reproducible labels) ---------
# if (is.null(seu_t$AgeGroups)) {
#   qs <- quantile(seu_t$age, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
#   labs <- paste0("[", floor(qs[-length(qs)]), "-", ceiling(qs[-1]), "]")
#   AgeGroups <- cut(seu_t$age, breaks = qs, include.lowest = TRUE, labels = labs)
#   seu_t$AgeGroups <- AgeGroups
# }
cat("AgeGroups levels:", paste(levels(seu_t$AgeGroups), collapse = " | "), "\n")

## --------- Helper: downsample cells for plotting (keeps figures light) ---------
sample_cells <- function(obj, max_n = 200000) {
  cells <- colnames(obj)
  if (length(cells) > max_n) {
    cells <- sample(cells, max_n, replace = FALSE)
    subset(obj, cells = cells)
  } else obj
}

## --------- UMAP plots (use existing reductions if present) ---------
## Prefer an existing embedding to avoid heavy compute
umap_name <- if ("X_umap" %in% names(seu_t@reductions)) "X_umap" else
  if ("umap"   %in% names(seu_t@reductions)) "umap"   else NA

if (is.na(umap_name)) {
  message("No UMAP found. Computing a light UMAP on a 100k downsample for illustration.")
  set.seed(15)
  small <- sample_cells(seu_t, 100000)
  small <- NormalizeData(small) |> FindVariableFeatures(nfeatures = 2000) |>
    ScaleData() |> RunPCA(npca = 30) |> RunUMAP(dims = 1:30)
  umap_obj <- small
  umap_red <- "umap"
} else {
  umap_obj <- sample_cells(seu_t, 200000)
  umap_red <- umap_name
}

## Plot 1: UMAP coloured by T-cell subset
p1 <- DimPlot(umap_obj, reduction = umap_red, group.by = "cell_type", raster = TRUE) +
  ggtitle("T-cell landscape (UMAP) by subset")
ggsave(file.path(FIGDIR, "UMAP_Tcells_by_subset.png"), p1, width = 10, height = 8, dpi = 300)

## Plot 2: UMAP coloured by AgeGroups (qualitative context for DE)
p2 <- DimPlot(umap_obj, reduction = umap_red, group.by = "AgeGroups", raster = TRUE) +
  ggtitle("T-cell landscape (UMAP) by AgeGroups (quartiles)")
ggsave(file.path(FIGDIR, "UMAP_Tcells_by_AgeGroups.png"), p2, width = 10, height = 8, dpi = 300)

## Optional: donor mixing (can be busy with ~1k donors)
# p3 <- DimPlot(umap_obj, reduction = umap_red, group.by = "donor_id", raster = TRUE, label = FALSE) +
#   ggtitle("UMAP coloured by donor (overview)")
# ggsave(file.path(FIGDIR, "UMAP_Tcells_by_donor.png"), p3, width = 10, height = 8, dpi = 300)

## --------- QC tables: counts per cell type & per donor ---------
meta <- seu_t@meta.data[, c("donor_id", "age", "cell_type")]
meta$donor_id <- droplevels(meta$donor_id)

## Table 1: total cells per T-cell subset
tcell_totals <- as.data.frame(meta |> count(cell_type, name = "n_cells")) |>
  arrange(desc(n_cells))
write.csv(tcell_totals, file.path(TABDIR, "tcell_subset_totals.csv"), row.names = FALSE)

## Table 2: counts per donor x T-cell subset (wide)
tab_long <- as.data.frame(meta |> count(donor_id, cell_type, name = "n_cells"))
tab_wide <- tidyr::pivot_wider(tab_long, names_from = cell_type, values_from = n_cells, values_fill = 0)
write.csv(tab_wide, file.path(TABDIR, "celltype_by_donor_counts.csv"), row.names = FALSE)

## Table 3: total T cells per donor
tcell_by_donor <- aggregate(n_cells ~ donor_id, tab_long, sum)
colnames(tcell_by_donor) <- c("donor_id", "tcell_total")
write.csv(tcell_by_donor, file.path(TABDIR, "tcells_by_donor_counts.csv"), row.names = FALSE)

## --------- Composition vs age (fractions per donor, per subset) ---------
## Fractions are within the T-cell compartment per donor
denom <- tcell_by_donor
comp <- merge(tab_long, denom, by = "donor_id")
comp$fraction <- comp$n_cells / comp$tcell_total

## Add donor age
donor_age <- meta |> distinct(donor_id, age)
comp <- merge(comp, donor_age, by = "donor_id")

## Fit fraction ~ age for each subset
summ_list <- lapply(split(comp, comp$cell_type), function(df) {
  if (nrow(df) < 5) return(data.frame(cell_type = unique(df$cell_type),
                                      slope = NA, p = NA, n_donors = length(unique(df$donor_id))))
  fit <- lm(fraction ~ age, data = df)
  co <- summary(fit)$coefficients
  data.frame(cell_type = unique(df$cell_type),
             slope = unname(co["age", "Estimate"]),
             p     = unname(co["age", "Pr(>|t|)"]),
             n_donors = length(unique(df$donor_id)))
})
summ <- do.call(rbind, summ_list)
summ$p_adj <- p.adjust(summ$p, method = "BH")
write.csv(summ, file.path(TABDIR, "composition_vs_age_summary.csv"), row.names = FALSE)

## Faceted scatter with linear fit per subset
gg_comp <- ggplot(comp, aes(x = age, y = fraction)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ cell_type, scales = "free_y") +
  labs(title = "Fraction of T-cell subset per donor vs age",
       x = "Age (years)", y = "Fraction within T cells") +
  theme_bw(base_size = 11)
ggsave(file.path(FIGDIR, "composition_vs_age_faceted.png"), gg_comp, width = 14, height = 10, dpi = 300)

## --------- Optional context: simple volcano for a chosen subset (edgeR output) ---------
## Set one of the 13 names below (or comment this block out)
chosen_subset <- NULL
# chosen_subset <- "mucosal invariant T cell"

if (!is.null(chosen_subset)) {
  clean <- gsub("[^A-Za-z0-9]+", "_", chosen_subset)
  fcsv  <- file.path(EDGER_DIR, paste0("DE_age_effect_", clean, ".csv"))
  if (file.exists(fcsv)) {
    de <- read.csv(fcsv, stringsAsFactors = FALSE)
    req <- c("logFC", "PValue", "FDR")
    if (all(req %in% colnames(de))) {
      de$negLog10P <- -log10(pmax(de$PValue, .Machine$double.xmin))
      pvol <- ggplot(de, aes(x = logFC, y = negLog10P)) +
        geom_point(alpha = 0.5, size = 0.6) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.3) +
        geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3) +
        labs(title = paste0("Volcano: age effect in ", chosen_subset),
             x = "logFC (per year of age)", y = "-log10(P)") +
        theme_bw(base_size = 11)
      ggsave(file.path(FIGDIR, paste0("volcano_", clean, ".png")), pvol, width = 7, height = 6, dpi = 300)
    } else {
      message("EdgeR file found but missing columns: ", fcsv)
    }
  } else {
    message("No edgeR file for subset: ", chosen_subset, " at ", fcsv)
  }
}

cat("\nDone. Files written to:\n - ", FIGDIR, "\n - ", TABDIR, "\n")


#RUN WITH:
# source("~/33_visuals_qc_and_composition.R")


