#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}
library(BiocManager)

packages <- c("DESeq2", "AnnotationDbi", "GenomicFeatures", "Rsamtools",
              "GenomicAlignments", "vsn", "ggplot2", 
              "pheatmap", "argparse", "ggrepel")

parser <- ArgumentParser(description='Análisis de DESeq2 para una comparación específica.')
parser$add_argument('--organism_db', type="character", required=TRUE)
parser$add_argument('--counting_method', type="character", required=TRUE)
parser$add_argument('--counts_file', type="character")
parser$add_argument('--bam_dir', type="character")
parser$add_argument('--gtf_file', type="character")
parser$add_argument('--fasta_file', type="character")
parser$add_argument('--metadata_file', type="character", required=TRUE)
parser$add_argument('--design_formula', type="character", required=TRUE)
parser$add_argument('--control_group', type="character", required=TRUE)
parser$add_argument('--key_type', type="character", required=TRUE)
parser$add_argument('--output_dir', type="character", required=TRUE)
parser$add_argument('--padj_threshold', type="double", default=0.05)
parser$add_argument('--log2fc_threshold', type="double", default=1.0)
parser$add_argument('--treatment_group', type="character", required=TRUE)
parser$add_argument('--control_group_contrast', type="character", required=TRUE)
args <- parser$parse_args()

packages <- c(packages, args$organism_db)

cat("📦 Verificando e instalando paquetes necesarios...\n")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("   -> Instalando", pkg, "...\n"))
    BiocManager::install(pkg, ask = FALSE)
  }
}

suppressPackageStartupMessages({
  invisible(lapply(packages, library, character.only = TRUE))
})

contrast_name <- paste0(args$treatment_group, "_vs_", args$control_group_contrast)
cat(paste("✅ 1. Preparando análisis para el contraste:", contrast_name, "\n"))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

if (tolower(args$counting_method) == "summarizeoverlaps") {
  cat("✅ 2. Generando matriz de conteo con summarizeOverlaps...\n")
  txdb <- makeTxDbFromGFF(args$gtf_file, format="gtf")
  exons_by_gene <- exonsBy(txdb, by="gene")
  metadata <- read.csv(args$metadata_file, row.names = 1)
  bam_files <- file.path(args$bam_dir, paste0(rownames(metadata), "_Aligned.sortedByCoord.out.bam"))
  bam_list <- BamFileList(bam_files, yieldSize=2000000)
  se <- summarizeOverlaps(features=exons_by_gene, reads=bam_list,
                          mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
  matriz_conteos <- assay(se)
} else {
  cat("✅ 2. Cargando matriz de conteo desde featureCounts...\n")
  matriz_completa <- read.table(args$counts_file, header = TRUE, sep = "\t", comment.char = "#")
  matriz_conteos <- as.matrix(matriz_completa[, 7:ncol(matriz_completa)])
  rownames(matriz_conteos) <- matriz_completa$Geneid
  
  # Limpieza de nombres de columna más robusta
  col_names_raw <- colnames(matriz_conteos)
  col_names_clean <- gsub("_Aligned.sortedByCoord.out.bam$", "", basename(col_names_raw))
  colnames(matriz_conteos) <- col_names_clean
  
  metadata <- read.csv(args$metadata_file, row.names = 1)
  matriz_conteos <- matriz_conteos[, rownames(metadata)]
}

cat("✅ 3. Ejecutando DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = matriz_conteos,
                              colData = metadata,
                              design = as.formula(args$design_formula))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = args$control_group)
dds <- DESeq(dds)

cat(paste("📊 Generando resultados para:", contrast_name, "\n"))
res <- results(dds, contrast=c("condition", args$treatment_group, args$control_group_contrast), alpha = args$padj_threshold)
res_ordered <- res[order(res$padj),]

cat("✅ 4. Anotando genes y guardando resultados...\n")
res_ordered_df <- as.data.frame(res_ordered)
res_ordered_df$gene_id <- rownames(res_ordered_df)

if (args$key_type == "ENSEMBL") {
  ids_mapeo <- sub("\\..*$", "", res_ordered_df$gene_id)
} else {
  ids_mapeo <- res_ordered_df$gene_id
}

res_ordered_df$symbol <- mapIds(get(args$organism_db), keys=ids_mapeo, column="SYMBOL", keytype=args$key_type, multiVals="first")
res_ordered_df <- res_ordered_df[, c("gene_id", "symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(res_ordered_df,
            file = file.path(args$output_dir, paste0("resultados_completos_", contrast_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

res_sig_annotated <- res_ordered_df[!is.na(res_ordered_df$padj) & res_ordered_df$padj < args$padj_threshold, ]
res_sig_annotated <- res_sig_annotated[order(res_sig_annotated$padj), ]

write.table(res_sig_annotated,
            file = file.path(args$output_dir, paste0("genes_significativos_", contrast_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat("✅ 5. Generando gráficos...\n")
vsd <- vst(dds, blind=FALSE)

if (!file.exists(file.path(args$output_dir, "grafico_PCA.pdf"))) {
    pca_plot <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("PCA Plot") + theme_bw()
    ggsave(file.path(args$output_dir, "grafico_PCA.pdf"), plot = pca_plot)
}

pdf(file.path(args$output_dir, paste0("grafico_MA_", contrast_name, ".pdf")))
plotMA(res, main=paste("MA Plot -", contrast_name))
dev.off()

res_df_volcano <- res_ordered_df
res_df_volcano$diffexpressed <- "NO"
res_df_volcano$diffexpressed[res_df_volcano$log2FoldChange > args$log2fc_threshold & !is.na(res_df_volcano$padj) & res_df_volcano$padj < args$padj_threshold] <- "UP"
res_df_volcano$diffexpressed[res_df_volcano$log2FoldChange < -args$log2fc_threshold & !is.na(res_df_volcano$padj) & res_df_volcano$padj < args$padj_threshold] <- "DOWN"
res_df_labeled <- head(res_df_volcano[order(res_df_volcano$padj), ], 10)

volcano_plot <- ggplot(data=res_df_volcano, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
  geom_point() + theme_minimal() +
  geom_text_repel(data=res_df_labeled, aes(label=symbol), max.overlaps = Inf) +
  scale_color_manual(values=c(DOWN="blue", NO="black", UP="red")) +
  geom_vline(xintercept=c(-args$log2fc_threshold, args$log2fc_threshold), col="red", linetype="dashed") +
  geom_hline(yintercept=-log10(args$padj_threshold), col="red", linetype="dashed") +
  ggtitle(paste("Gráfico de Volcán -", contrast_name))
ggsave(file.path(args$output_dir, paste0("grafico_volcan_", contrast_name, ".pdf")), plot = volcano_plot)

if (nrow(res_sig_annotated) > 1) {
  top_genes <- head(res_sig_annotated$gene_id, 30)
  mat <- assay(vsd)[top_genes,]
  mat <- mat - rowMeans(mat)
  annotation_col = data.frame(condition = colData(dds)$condition)
  rownames(annotation_col) = colnames(mat)
  pheatmap_file <- file.path(args$output_dir, paste0("mapa_calor_top30_", contrast_name, ".pdf"))
  pheatmap(mat, annotation_col = annotation_col, filename = pheatmap_file, main = paste("Top 30 DE Genes -", contrast_name))
} else {
  cat(paste("ℹ️ No se generó mapa de calor para", contrast_name, " (menos de 2 genes significativos).\n"))
}

cat(paste("🎉 Análisis para", contrast_name, "finalizado.\n"))