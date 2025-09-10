#31a_advanced_tcell_ageing_analysis_NOTadjSex.R
# WebGestalt - NOTsexadjusted
# GSEA


## === Configure once ===
PROJ <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/"
indir  <- file.path(PROJ, "results/edgeR_NOTsexadjusted")
outdir <- file.path(PROJ, "results/webgestalt_NOTsexadjusted")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## Choose cell type tag (matches your DE file suffix)
cell <- "double_negative_thymocyte_edgeR"   # <-- change for each T-cell type

# File Names - substitute into cell code above
# CD4_positive_alpha_beta_T_cell_edgeR
# CD8_positive_alpha_beta_T_cell_edgeR
# gamma_delta_T_cell_edgeR
# regulatory_T_cell_edgeR
# naive_thymus_derived_CD4_positive_alpha_beta_T_cell_edgeR
# naive_thymus_derived_CD8_positive_alpha_beta_T_cell_edgeR
# central_memory_CD4_positive_alpha_beta_T_cell_edgeR
# effector_memory_CD4_positive_alpha_beta_T_cell_edgeR
# central_memory_CD8_positive_alpha_beta_T_cell_edgeR
# effector_memory_CD8_positive_alpha_beta_T_cell_edgeR
# CD4_positive_alpha_beta_cytotoxic_T_cell_edgeR
# mucosal_invariant_T_cell_edgeR
# double_negative_thymocyte_edgeR


# Tcell Names - reference only
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


de_csv <- file.path(indir, paste0("notsexadjusted_DE_age_effect_", cell, ".csv"))
stopifnot(file.exists(de_csv))

tab <- read.csv(de_csv, check.names = FALSE)

## --- Identify the gene ID column & clean Ensembl IDs ---
cand_gene_cols <- c("gene","gene_id","GeneID","ensembl","Ensembl","EnsemblID","Ensembl_ID")
idcol <- cand_gene_cols[cand_gene_cols %in% names(tab)][1]
if (is.na(idcol)) {
  if (!is.null(rownames(tab))) {
    tab$gene <- rownames(tab); idcol <- "gene"
  } else stop("No gene/Ensembl column or rownames found in DE table.")
}
tab$Ensembl <- sub("\\..*$", "", as.character(tab[[idcol]]))  # strip version suffixes

## --- Locate p-value & FDR columns  ---
p_candidates  <- c("PValue","P.Value","pvalue","p.value","PValue.edgeR","P.Value.edgeR")
q_candidates  <- c("FDR","adj.P.Val","adj.PVal","padj","qvalue","FDR.edgeR")
pname <- p_candidates[p_candidates %in% names(tab)][1]
qname <- q_candidates[q_candidates %in% names(tab)][1]
if (is.na(pname)) stop("No P-value column found.")
if (is.na(qname)) stop("No FDR/adjusted P column found.")
stopifnot("logFC" %in% names(tab))

## --- UNIVERSE: all tested genes (Ensembl, no version) ---
#universe <- unique(tab$Ensembl[!is.na(tab$Ensembl) & nzchar(tab$Ensembl)])

universe <- unique(tab[,1])
writeLines(universe, file.path(outdir, paste0(cell, "_UNIVERSE_ensembl.txt")))

## --- UP & DOWN: FDR < 0.05 split by logFC sign ---
#up   <- unique(tab$Ensembl[tab[[qname]] < 0.05 & tab$logFC > 0])
#down <- unique(tab$Ensembl[tab[[qname]] < 0.05 & tab$logFC < 0])


up=tab[which(tab$FDR<0.05 & tab$logFC > 0),1]
down=tab[which(tab$FDR<0.05 & tab$logFC < 0),1]

writeLines(up,   file.path(outdir, paste0(cell, "_UP_ensembl.txt")))
writeLines(down, file.path(outdir, paste0(cell, "_DOWN_ensembl.txt")))

cat(sprintf("%s: %d UP, %d DOWN, %d in UNIVERSE\n",
            cell, length(up), length(down), length(universe)))

## --- RANKED list for pre-ranked GSEA (signed score) ---
score <- sign(tab$logFC) * -log10(pmax(tab[[pname]], .Machine$double.xmin))
rk <- data.frame(Ensembl = tab$Ensembl, Score = score, stringsAsFactors = FALSE)
rk <- rk[!is.na(rk$Ensembl) & nzchar(rk$Ensembl), ]
rk <- rk[!duplicated(rk$Ensembl), ]
rk <- rk[order(rk$Score, decreasing = TRUE), ]

write.table(rk,
            file = file.path(outdir, paste0(cell, "_RANKED_ensembl.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################################################
# HOW TO RUN
# source("~/31a_advanced_tcell_ageing_analysis_NOTadjSex.R")
################################################################################