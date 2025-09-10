# Rscript ~/41d_concordance_jaccard_bars.R



# Create Jaccard bar chart for all 13 subsets

suppressPackageStartupMessages({library(readr); library(dplyr); library(tidyr); library(ggplot2)})

PROJ   <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
tbldir <- file.path(PROJ, "results/tables")
figdir <- file.path(PROJ, "results/figures")
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

x <- readr::read_csv(file.path(tbldir, "sex_adj_concordance_summary.csv"), show_col_types = FALSE)

plotdat <- x %>%
  select(subset, jaccard_up, jaccard_down) %>%
  pivot_longer(cols = starts_with("jaccard_"),
               names_to = "direction",
               values_to = "jaccard") %>%
  mutate(direction = ifelse(direction == "jaccard_up", "UP", "DOWN"))

p <- ggplot(plotdat, aes(x = reorder(subset, jaccard), y = jaccard, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  coord_flip() +
  labs(title = "Jaccard overlap (sex-adjusted vs primary) by subset",
       x = "T-cell subset", y = "Jaccard overlap") +
  theme_bw() +
  theme(legend.position = "top")

outfile <- file.path(figdir, "sex_adj_concordance_bars.png")
ggsave(outfile, p, width = 9, height = 6, dpi = 300)
message("Wrote: ", outfile)
