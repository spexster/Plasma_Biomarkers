# 30_edgeR_age_sex.R
# tcell regression analysis
# sessionInfo() - gives detail packages and products running

suppressPackageStartupMessages({library(edgeR); library(limma)})

# --- inputs that can be change per run ---
cell_type <- "double negative thymocyte"  # <- set to the T-cell type for this run
pb_dir    <- "~/scripts"                 # where *_count.txt files were written


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

# --- paths / Seurat metadata source ---
seu_path <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/seurat_snapshots/seu_full_withmeta_latest.rds"
stopifnot(file.exists(seu_path))
seu_full_withmeta_latest <- readRDS(seu_path)

# pseudobulk counts for this cell type
cell_file <- file.path(pb_dir, paste0(gsub("[^A-Za-z0-9]+","_", cell_type), "_count.txt"))
stopifnot(file.exists(cell_file))
cts <- read.table(cell_file, sep="\t", header=TRUE, row.names=1, check.names=FALSE)  # genes x donors

# align donor-level metadata
donors <- colnames(cts)
meta_all <- as.data.frame(seu_full_withmeta_latest@meta.data[, c("donor_id","age","sex")])
meta <- unique(meta_all)                    # one row per donor
rownames(meta) <- as.character(meta$donor_id)
meta <- meta[donors, , drop=FALSE]
stopifnot(all(donors == rownames(meta)))

# edgeR pipeline: filter, TMM, design, QL-F test for age effect (adjusting sex)
y <- DGEList(counts = cts)
keep <- rowSums(cpm(y) > 1) >= ceiling(0.2 * ncol(y))  # CPM>1 in ≥20% of donors
y <- y[keep,, keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")

meta$age_sc <- as.numeric(scale(meta$age))
meta$sex    <- relevel(factor(meta$sex), ref="female")
design <- model.matrix(~ age_sc + sex, data=meta)

y  <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef="age_sc")

tab <- topTags(res, n=Inf)$table

# write results
outdir <- "~/projects/Jafari_PlasmaAgeing_Srivastava-Papadaki_2025/results/edgeR"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, paste0("_DE_age_effect_", gsub("[^A-Za-z0-9]+","_", cell_type), "_edgeR.csv"))
write.csv(tab, outfile, row.names=TRUE)

message("Done: ", outfile, "  (", nrow(tab), " genes)")


## HOW TO RUN SCRIPT 
# source("~/30_edgeR_age_sex.R")

