suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

counts_file <- getArg("--counts")
samples_file <- getArg("--samples")
outdir <- getArg("--outdir", "results/06_limma")

if (is.null(counts_file) || is.null(samples_file)) {
  stop("Usage: Rscript R/02_limma_voom_DE_heatmap.R --counts featureCounts.txt --samples samples.tsv --outdir outdir")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 1) Read featureCounts output
# -------------------------
fc <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

# featureCounts columns: Geneid Chr Start End Strand Length <sample1.bam> <sample2.bam> ...
counts <- fc[, 7:ncol(fc)]
rownames(counts) <- fc$Geneid

# Clean column names to sample IDs (basename without .bam)
colnames(counts) <- gsub("\\.bam$", "", basename(colnames(counts)))

# -------------------------
# 2) Read sample sheet (samples.tsv)
# -------------------------
meta <- read.delim(samples_file, stringsAsFactors = FALSE)
stopifnot(all(c("sample","group","replicate") %in% colnames(meta)))

missing <- setdiff(meta$sample, colnames(counts))
if (length(missing) > 0) {
  stop("These samples are in samples.tsv but missing in counts: ", paste(missing, collapse=", "))
}

# Order columns by sample sheet
counts <- counts[, meta$sample, drop = FALSE]

meta$group <- factor(meta$group)
meta$replicate <- factor(meta$replicate)

write.csv(meta, file.path(outdir, "sample_metadata.csv"), row.names = FALSE)

# -------------------------
# 3) edgeR filtering + normalization
# -------------------------
dge <- DGEList(counts = counts, samples = meta)
keep <- filterByExpr(dge, group = meta$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# -------------------------
# 4) limma-voom
# -------------------------
design <- model.matrix(~0 + group, data = meta)
colnames(design) <- levels(meta$group)

v <- voom(dge, design)

pdf(file.path(outdir, "voom_meanvariance.pdf"))
plot(v)
dev.off()

fit <- lmFit(v, design)

# -------------------------
# 5) Contrasts (default: each group vs NT)
#    Edit this section if your paper uses different comparisons.
# -------------------------
if (!("NT" %in% colnames(design))) {
  stop("Group 'NT' not found. Available groups: ", paste(colnames(design), collapse=", "))
}

contr_list <- lapply(setdiff(colnames(design), "NT"), function(g) paste0(g, " - NT"))
names(contr_list) <- paste0(setdiff(colnames(design), "NT"), "_vs_NT")

contrasts_mat <- makeContrasts(contrasts = unlist(contr_list), levels = design)

fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2)

# -------------------------
# 6) Save DE results per contrast
# -------------------------
all_results <- list()

for (coef_name in colnames(contrasts_mat)) {
  res <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "P")
  res$gene_id <- rownames(res)
  res <- res[, c("gene_id", setdiff(colnames(res), "gene_id"))]

  write.csv(res, file.path(outdir, paste0("DE_", coef_name, ".csv")), row.names = FALSE)
  all_results[[coef_name]] <- res
}

saveRDS(all_results, file.path(outdir, "DE_all_contrasts.rds"))

# -------------------------
# 7) QC plots
# -------------------------
pdf(file.path(outdir, "MDS_plot.pdf"))
plotMDS(v, labels = meta$group, col = as.numeric(meta$group))
dev.off()

logCPM <- cpm(dge, log = TRUE, prior.count = 1)

cors <- cor(logCPM, method = "pearson")
ann <- data.frame(group = meta$group)
rownames(ann) <- meta$sample

pdf(file.path(outdir, "SampleCorrelation_heatmap.pdf"), width = 7, height = 6)
pheatmap(cors, annotation_col = ann, annotation_row = ann)
dev.off()

# -------------------------
# 8) Heatmap of top DE genes per contrast
# -------------------------
topN <- 50

for (coef_name in names(all_results)) {
  res <- all_results[[coef_name]]
  res_sig <- res[!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05, , drop = FALSE]
  if (nrow(res_sig) < 10) res_sig <- res

  top_genes <- head(res_sig$gene_id[order(res_sig$adj.P.Val, res_sig$P.Value)], topN)
  mat <- logCPM[top_genes, , drop = FALSE]

  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0

  pdf(file.path(outdir, paste0("Heatmap_top", topN, "_", coef_name, ".pdf")), width = 8, height = 8)
  pheatmap(
    mat_z,
    annotation_col = ann,
    show_colnames = TRUE,
    show_rownames = TRUE,
    main = paste0("Top ", topN, " genes: ", coef_name)
  )
  dev.off()
}

# -------------------------
# 9) Volcano plots
# -------------------------
for (coef_name in names(all_results)) {
  res <- all_results[[coef_name]]
  res$negLog10FDR <- -log10(pmax(res$adj.P.Val, 1e-300))
  res$Sig <- (res$adj.P.Val < 0.05) & (abs(res$logFC) >= 1)

  p <- ggplot(res, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(alpha = Sig), size = 1.2) +
    labs(title = paste0("Volcano: ", coef_name),
         x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_bw()

  ggsave(file.path(outdir, paste0("Volcano_", coef_name, ".png")),
         plot = p, width = 6, height = 5, dpi = 300)
}

message("DONE. Results in: ", outdir)
