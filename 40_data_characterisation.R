# 40_data_characterisation.R
# produce summary tables for data characterisation


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(pheatmap)
})

PROJ <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
FIGDIR <- file.path(PROJ, "results/figures"); dir.create(FIGDIR, showWarnings=FALSE)
TABDIR <- file.path(PROJ, "results/tables")

# Cell counts
sizes <- read.csv(file.path(TABDIR, "tcells_subset_sizes.csv"))
# Plot: total cells per subset
ggplot(sizes, aes(x=reorder(cell_type, -n_cells), y=log10(n_cells))) +
  geom_col(fill="blue") + coord_flip() +
  labs(x="Subset", y="log10(Cell count)") +
  theme_minimal(base_size=12)
ggsave(file.path(FIGDIR, "subset_cell_totals_log10.png"), width=6, height=6)

# UP/DOWN counts from edgeR results
# (assumes *_edgeR.csv files exist in results/edgeR/)
files <- list.files(file.path(PROJ, "results/edgeR"), pattern="_edgeR.csv", full.names=TRUE)
summ <- lapply(files, function(f) {
  tb <- read_csv(f, show_col_types=FALSE)
  sub <- gsub("_edgeR.csv","",basename(f))
  tibble(
    subset=sub,
    n_up=sum(tb$FDR<0.05 & tb$logFC>0, na.rm=TRUE),
    n_down=sum(tb$FDR<0.05 & tb$logFC<0, na.rm=TRUE)
  )
}) |> bind_rows()
write_csv(summ, file.path(TABDIR, "subset_up_down_counts.csv"))

# Top10 subsets by total sig
summ2 <- summ |> mutate(total=n_up+n_down) |> arrange(desc(total)) |> slice_head(n=10)
ggplot(summ2, aes(x=reorder(subset,total), y=-log10(total))) +
  geom_col(fill="red") + coord_flip() +
  labs(x="Subset", y="-log10(# sig genes)") +
  theme_minimal(base_size=12)
ggsave(file.path(FIGDIR, "top10_subsets_siggenes_bar.png"), width=6, height=6)

# Specificity heatmap
spec <- read.csv(file.path(TABDIR, "age_sig_specificity_matrix.csv"), check.names=FALSE)
mat <- as.matrix(spec[,-1]); rownames(mat) <- spec[[1]]
pheatmap(mat, show_rownames=FALSE, fontsize_col=10,
         filename=file.path(FIGDIR,"age_sig_specificity_heatmap.png"))
