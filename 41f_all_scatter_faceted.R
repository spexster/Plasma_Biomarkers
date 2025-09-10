# ~/41f_all_scatter_faceted.R
#Rscript ~/41f_all_scatter_faceted.R


# Created to re-reads the unadjusted + sex-adjusted CSVs for all 13 subsets, 
# merge, and draw a single faceted scatter (one panel per subset). 
# Also annotate each facet with Pearson r and N


suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2)})

PROJ   <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
edir   <- file.path(PROJ, "results/edgeR_NOTsexadjusted")
edirS  <- file.path(PROJ, "results/edgeR_ISsexadjusted")
figdir <- file.path(PROJ, "results/figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

stubify <- function(x){
  x |>
    gsub(pattern = "[^A-Za-z0-9]+", replacement = "_") |>
    gsub(pattern = "_+", replacement = "_") |>
    sub(pattern = "^_", replacement = "") |>
    sub(pattern = "_$", replacement = "")
}

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

read_de <- function(path){
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!"gene" %in% names(df)) names(df)[1] <- "gene"
  names(df) <- tolower(names(df))
  stopifnot("gene" %in% names(df))
  if (!"logfc" %in% names(df)) stop("No logFC column in: ", path)
  dplyr::select(df, gene, logfc)
}

get_ab <- function(ct){
  stub <- stubify(ct)
  prim <- file.path(edir,  paste0("notsexadjusted_DE_age_effect_", stub, "_edgeR.csv"))
  adj <- file.path(edirS, paste0("DE_age_effect_", stub, "_edgeR.csv"))
  if (!file.exists(prim) || !file.exists(adj)) return(NULL)
  a <- read_de(prim) %>% dplyr::rename(logFC_unadj = logfc)
  b <- read_de(adj)  %>% dplyr::rename(logFC_adj   = logfc)
  ab <- dplyr::inner_join(a, b, by = "gene")
  ab$subset <- ct
  ab
}

dfs <- lapply(ctypes, get_ab)
dfs <- dfs[!vapply(dfs, is.null, logical(1))]
stopifnot(length(dfs) > 0)

DF <- dplyr::bind_rows(dfs)

# Per-subset correlations (computed on full data)
corr <- DF %>%
  dplyr::group_by(subset) %>%
  dplyr::summarise(
    pearson = suppressWarnings(cor(logFC_unadj, logFC_adj, method = "pearson", use = "complete.obs")),
    n = dplyr::n(),
    .groups = "drop"
  )

# Downsample per subset for plotting (keeps file size reasonable)
set.seed(15)
DF_plot <- DF %>%
  dplyr::group_by(subset) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), 30000))) %>%
  dplyr::ungroup()

lab_map <- setNames(
  paste0(levels(factor(corr$subset, levels = unique(DF$subset))),
         "\n(r=", sprintf("%.2f", corr$pearson), ", N=", corr$n, ")"),
  corr$subset
)

p <- ggplot(DF_plot, aes(logFC_unadj, logFC_adj)) +
  geom_point(alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~subset, ncol = 4, labeller = labeller(subset = lab_map)) +
  labs(
    title = "logFC (unadjusted) vs (sex-adjusted) per T-cell subset",
    subtitle = "Each panel shows Pearson r and gene count (computed on all genes before sampling)",
    x = "logFC (unadjusted age model)",
    y = "logFC (sex-adjusted age model)"
  ) +
  theme_bw()

outfile_png <- file.path(figdir, "sex_adj_logFC_scatter_faceted.png")
outfile_pdf <- file.path(figdir, "sex_adj_logFC_scatter_faceted.pdf")

png(outfile_png, width = 2400, height = 1800, res = 200)
print(p)
dev.off()
ggsave(outfile_pdf, p, width = 12, height = 9)

message("Wrote:\n - ", outfile_png, "\n - ", outfile_pdf)


