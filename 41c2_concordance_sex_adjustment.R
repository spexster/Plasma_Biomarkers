# Rscript ~/41c2_concordance_sex_adjustment.R
# Compare DE - age-only vs sexadj across all 13 tcell subtypes. Includes Jaccard below, commanded out

# 41c2_concordance_sex_adjustment.R
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
})

PROJ    <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025"
UNADJ   <- file.path(PROJ, "results/edgeR_NOTsexadjusted")
ADJ     <- file.path(PROJ, "results/edgeR_ISsexadjusted")
TABDIR  <- file.path(PROJ, "results/tables")
FIGDIR  <- file.path(PROJ, "results/figures")
dir.create(TABDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

# Display names used everywhere else
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

stubify <- function(x){
  x |>
    gsub(pattern = "[^A-Za-z0-9]+", replacement = "_") |>
    gsub(pattern = "_+", replacement = "_") |>
    sub(pattern = "^_", replacement = "") |>
    sub(pattern = "_$", replacement = "")
}

# Robust CSV reader: standardise column names and ensure 'gene' + 'logFC' + 'FDR'
read_de <- function(path){
  df <- readr::read_csv(path, show_col_types = FALSE)
  nms <- tolower(names(df))
  names(df) <- nms
  if (!"gene" %in% names(df)) names(df)[1] <- "gene"
  stopifnot("gene" %in% names(df))
  # logFC
  if (!"logfc" %in% names(df)) stop("No logFC column in: ", path)
  # FDR (try several common names)
  fdr_col <- c("fdr","adj.p.val","padj","qvalue")[c("fdr","adj.p.val","padj","qvalue") %in% names(df)][1]
  if (is.na(fdr_col)) stop("No FDR/adjusted-P column in: ", path)
  df %>% transmute(gene = gene, logFC = as.numeric(logfc), FDR = as.numeric(.data[[fdr_col]]))
}

jaccard <- function(A, B){
  A <- unique(A); B <- unique(B)
  u <- union(A,B); if (length(u) == 0) return(NA_real_)
  length(intersect(A,B)) / length(u)
}

# Corrected 'for' loop for script 41c2

rows <- list()

for (ct in ctypes){
  stub <- stubify(ct)
  
  unadj_csv <- file.path(UNADJ, paste0("notsexadjusted_DE_age_effect_", stub, "_edgeR.csv"))
  adj_csv   <- file.path(ADJ,   paste0("DE_age_effect_", stub, "_edgeR.csv"))
  
  # Check for file existence first and skip if a file is missing
  if (!file.exists(unadj_csv) || !file.exists(adj_csv)) {
    message("SKIP (missing): ", ct)
    next
  }
  
  U <- tryCatch(read_de(unadj_csv), error = function(e) NULL)
  A <- tryCatch(read_de(adj_csv),   error = function(e) NULL)
  
  if (is.null(U) || is.null(A)) {
    message("SKIP (bad columns): ", ct)
    next
  }
  
  # UP/DOWN sets
  up_U  <- U$gene[U$FDR < 0.05 & U$logFC > 0]
  down_U <- U$gene[U$FDR < 0.05 & U$logFC < 0]
  up_A  <- A$gene[A$FDR < 0.05 & A$logFC > 0]
  down_A <- A$gene[A$FDR < 0.05 & A$logFC < 0]
  
  jac_up   <- jaccard(up_U, up_A)
  jac_down <- jaccard(down_U, down_A)
  
  
  # logFC correlation + sign-concordance
  AB <- inner_join(select(U, gene, logFC), select(A, gene, logFC), by = "gene", suffix = c("_unadj","_adj"))
  pear <- suppressWarnings(cor(AB$logFC_unadj, AB$logFC_adj, method = "pearson"))
  spear <- suppressWarnings(cor(AB$logFC_unadj, AB$logFC_adj, method = "spearman"))
  sign_pct <- mean(sign(AB$logFC_unadj) == sign(AB$logFC_adj)) * 100
  
  rows[[length(rows)+1]] <- tibble::tibble(
    subset = ct,
    jaccard_up = jac_up,
    jaccard_down = jac_down,
    cor_logFC = pear,
    spearman_cor = spear,
    sign_concord_pct = sign_pct,
    n_unadj_up = length(up_U),
    n_unadj_down = length(down_U),
    n_adj_up = length(up_A),
    n_adj_down = length(down_A),
    n_common_genes = nrow(AB),
    unadjusted_csv = unadj_csv,
    adjsex_csv = adj_csv
  )
}
  
  # # logFC correlation + sign-concordance - added spearman in chunk above
  # AB <- inner_join(select(U, gene, logFC), select(A, gene, logFC), by = "gene", suffix = c("_unadj","_adj"))
  # pear <- suppressWarnings(cor(AB$logFC_unadj, AB$logFC_adj, method = "pearson"))
  # sign_pct <- mean(sign(AB$logFC_unadj) == sign(AB$logFC_adj)) * 100
  # 
  # rows[[length(rows)+1]] <- tibble::tibble(
  #   subset = ct,
  #   jaccard_up = jac_up,
  #   jaccard_down = jac_down,
  #   cor_logFC = pear,
  #   sign_concord_pct = sign_pct,
  #   n_unadj_up = length(up_U),
  #   n_unadj_down = length(down_U),
  #   n_adj_up = length(up_A),
  #   n_adj_down = length(down_A),
  #   n_common_genes = nrow(AB),
  #   unadjusted_csv = unadj_csv,
  #   adjsex_csv = adj_csv
  # )

  
  # # UP/DOWN sets - DUPLICATE & EXTRA LINE
  # up_U   <- U$gene[U$FDR < 0.05 & U$
  #                    logFC > 0]
  # down_U <- U$gene[U$FDR < 0.05 & U$logFC < 0]
  # up_A   <- A$gene[A$FDR < 0.05 & A$logFC > 0]
  # down_A <- A$gene[A$FDR < 0.05 & A$logFC < 0]
  # 
  # jac_up   <- jaccard(up_U, up_A)
  # jac_down <- jaccard(down_U, down_A)
  # 
  # # logFC correlation + sign-concordance
  # AB <- inner_join(select(U, gene, logFC), select(A, gene, logFC), by = "gene", suffix = c("_unadj","_adj"))
  # pear <- suppressWarnings(cor(AB$logFC_unadj, AB$logFC_adj, method = "pearson"))
  # sign_pct <- mean(sign(AB$logFC_unadj) == sign(AB$logFC_adj)) * 100
  # 
  # rows[[length(rows)+1]] <- tibble::tibble(
  #   subset = ct,
  #   jaccard_up = jac_up,
  #   jaccard_down = jac_down,
  #   cor_logFC = pear,
  #   sign_concord_pct = sign_pct,
  #   n_unadj_up = length(up_U),
  #   n_unadj_down = length(down_U),
  #   n_adj_up = length(up_A),
  #   n_adj_down = length(down_A),
  #   n_common_genes = nrow(AB),
  #   unadjusted_csv = unadj_csv,
  #   adjsex_csv = adj_csv
  # )


res <- bind_rows(rows)
out_csv <- file.path(TABDIR, "sex_adj_concordance_from_csv_summary.csv")
write_csv(res, out_csv)
message("Wrote: ", out_csv)

# --- Make the UP/DOWN Jaccard bar chart
if (nrow(res)){
  long <- res %>%
    select(subset, jaccard_up, jaccard_down) %>%
    pivot_longer(-subset, names_to = "set", values_to = "jaccard") %>%
    mutate(set = recode(set, jaccard_up = "UP", jaccard_down = "DOWN")) %>%
    filter(!is.na(jaccard))
  
  if (nrow(long)){
    p <- ggplot(long, aes(x = reorder(subset, jaccard, na.rm = TRUE), y = jaccard, fill = set)) +
      geom_col(position = "dodge") +
      geom_text(aes(label = sprintf("%.2f", jaccard)),
                position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
      coord_flip() +
      labs(title = "Sex-adjusted vs Primary:\nJaccard Overlap of Significant Sets (FDR<0.05)",
           x = "T-cell subset", y = "Jaccard overlap", fill = "") +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
      theme_minimal(base_size = 12)
    
    ggsave(file.path(FIGDIR, "sex_adj_concordance_bar.png"), p, width = 8, height = 4.5, dpi = 300, bg = "white")
    message("Wrote: ", file.path(FIGDIR, "sex_adj_concordance_bar.png"))
  } else {
    message("No non-NA Jaccard values to plot.")
  }
} else {
  message("No rows in summary; check file paths and naming.")
}




################################################################################
# Create UP/DOWN Jaccard bar chart - run it after above code - latest v.
################################################################################

# tbldir <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/tables"
# figdir <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/figures"
# dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
# 
# suppressPackageStartupMessages({
#   library(readr); library(dplyr); library(tidyr); library(ggplot2)
# })
# 
# res <- readr::read_csv(file.path(tbldir, "sex_adj_concordance_from_csv_summary.csv"), show_col_types = FALSE)
# 
# # Long format for UP/DOWN bars
# long <- res %>%
#   select(subset, jaccard_up, jaccard_down) %>%
#   pivot_longer(-subset, names_to = "set", values_to = "jaccard") %>%
#   mutate(set = recode(set, jaccard_up = "UP", jaccard_down = "DOWN"))
# 
# p <- ggplot(long, aes(x = reorder(subset, jaccard, na.rm = TRUE), y = jaccard, fill = set)) +
#   geom_col(position = "dodge") +
#   geom_text(aes(label = sprintf("%.2f", jaccard)),
#             position = position_dodge(width = 0.9), vjust = -0.3, size = 3) +
#   coord_flip() +
#   labs(title = "Sex-adjusted vs Primary:\nJaccard Overlap of Significant Sets",
#        x = "T-cell subset", y = "Jaccard overlap", fill = "") +
#   scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
#   theme_minimal(base_size = 12)
# 
# ggsave(file.path(figdir, "sex_adj_concordance_bar.png"), p, width = 8, height = 4.5, dpi = 300, bg="white")
# 
# message("Wrote: ", file.path(figdir, "sex_adj_concordance_bar.png"))



