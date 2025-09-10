# Rscript ~/41e_correlate_adjsex_logFC.R   
# writes sex_adj_logFC_correlations.csv + PNGs
# will allow complete cor_logFC on sex_adj_concordance_summary.csv


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

PROJ   <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
edir   <- file.path(PROJ, "results/edgeR_NOTsexadjusted")
edirS  <- file.path(PROJ, "results/edgeR_ISsexadjusted")
tbldir <- file.path(PROJ, "results/tables")
figdir <- file.path(PROJ, "results/figures")
dir.create(tbldir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# ---- SAME stubify as in 35d ----
stubify <- function(x){
  x |>
    gsub(pattern = "[^A-Za-z0-9]+", replacement = "_") |>
    gsub(pattern = "_+", replacement = "_") |>
    sub(pattern = "^_", replacement = "") |>
    sub(pattern = "_$", replacement = "")
}

# 13 display names (must match your metadata / earlier scripts)
ctypes <- c(
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

# robust CSV reader: normalise 'gene' and 'logFC' columns
read_de <- function(path){
  df <- readr::read_csv(path, show_col_types = FALSE)
  # first column might be gene id
  if (!"gene" %in% names(df)) {
    names(df)[1] <- "gene"
  }
  # standardise case of columns
  nms <- tolower(names(df))
  names(df) <- nms
  # require gene + logfc
  stopifnot("gene" %in% names(df))
  if (!"logfc" %in% names(df)) {
    stop("No logFC column in: ", path)
  }
  df %>% select(gene, logfc)
}

rows <- list()

for (ct in ctypes) {
  stub <- stubify(ct)
  prim <- file.path(edir,  paste0("notsexadjusted_DE_age_effect_", stub, "_edgeR.csv"))
  adj <- file.path(edirS, paste0("DE_age_effect_", stub, "_edgeR.csv"))
  
  if (!file.exists(prim)) {
    message("SKIP (missing unadjusted): ", ct, " -> ", prim)
    next
  }
  if (!file.exists(adj)) {
    message("SKIP (missing adjSex):    ", ct, " -> ", adj)
    next
  }
  
  a <- read_de(prim) %>% rename(logFC_unadj = logfc)
  b <- read_de(adj)  %>% rename(logFC_adj   = logfc)
  
  ab <- inner_join(a, b, by = "gene")
  n_common <- nrow(ab)
  if (n_common < 10) {
    message("⚠ Very few common genes for ", ct, " (", n_common, "). Skipping plot.")
  }
  
  # correlations
  pear <- suppressWarnings(cor(ab$logFC_unadj, ab$logFC_adj, method = "pearson", use = "complete.obs"))
  spear <- suppressWarnings(cor(ab$logFC_unadj, ab$logFC_adj, method = "spearman", use = "complete.obs"))
  
  # scatter
  p <- ggplot(ab, aes(x = logFC_unadj, y = logFC_adj)) +
    geom_point(alpha = 0.3, size = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste0(ct, " — logFC (unadjusted) vs (sex-adjusted)"),
      x = "logFC (unadjusted age model)",
      y = "logFC (sex-adjusted age model)",
      subtitle = paste0("Pearson r = ", round(pear, 3), " | Spearman ρ = ", round(spear, 3),
                        " | N = ", n_common, " genes")
    )
  png(file.path(figdir, paste0(stub, "_sex_adj_logFC_scatter.png")), width = 1400, height = 1000, res = 150)
  print(p)
  dev.off()
  
  rows[[length(rows) + 1]] <- tibble::tibble(
    subset = ct,
    pearson_cor  = pear,
    spearman_cor = spear,
    n_common_genes = n_common,
    unadjusted_csv = prim,
    adjsex_csv     = adj
  )
}

out <- dplyr::bind_rows(rows)
readr::write_csv(out, file.path(tbldir, "sex_adj_logFC_correlations.csv"))
message("Wrote: ", file.path(tbldir, "sex_adj_logFC_correlations.csv"))
message("Figures: *_sex_adj_logFC_scatter.png in ", figdir)
