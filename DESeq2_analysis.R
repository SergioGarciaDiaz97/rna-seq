#!/usr/bin/env Rscript

# ==============================
# Script DESeq2 Completo con QC
# ==============================

# --- 1. INSTALACIÓN Y CARGA DE ARGPARSE ---
cat("📦 Verificando 'argparse'...\n")
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse", repos = "http://cran.us.r-project.org")
}
library(argparse)

# --- 2. DEFINICIÓN Y LECTURA DE ARGUMENTOS ---
parser <- ArgumentParser(description='Análisis de DESeq2 con QC y enriquecimiento opcional.')
parser$add_argument('--counting_method', type="character", required=TRUE, help="Método de conteo: featureCounts o summarizeOverlaps")
parser$add_argument('--counts_file', type="character", help="Ruta al archivo de conteos (para featureCounts)")
parser$add_argument('--bam_dir', type="character", help="Directorio con archivos BAM (para summarizeOverlaps)")
parser$add_argument('--gtf_file', type="character", help="Ruta al archivo GTF (para summarizeOverlaps)")
parser$add_argument('--metadata_file', type="character", required=TRUE, help="Archivo de metadatos (CSV)")
parser$add_argument('--design_formula', type="character", required=TRUE, help="Fórmula de diseño DESeq2 (ej. '~ condition')")
parser$add_argument('--control_group', type="character", required=TRUE, help="Grupo control de referencia")
parser$add_argument('--treatment_group', type="character", required=TRUE, help="Grupo de tratamiento para el contraste")
parser$add_argument('--control_group_contrast', type="character", required=TRUE, help="Grupo control para el contraste")
parser$add_argument('--organism_db', type="character", required=TRUE, help="Paquete de anotación (ej. 'org.Hs.eg.db')")
parser$add_argument('--key_type', type="character", required=TRUE, help="Tipo de ID de gen (ej. 'ENSEMBL')")
parser$add_argument('--strip_gene_version', type="logical", default=FALSE, help="TRUE para quitar la versión de los IDs de gen (ej. '.1')")
parser$add_argument('--output_dir', type="character", required=TRUE, help="Directorio de salida")
parser$add_argument('--padj_threshold', type="double", default=0.05, help="Umbral de p-valor ajustado")
parser$add_argument('--log2fc_threshold', type="double", default=1.0, help="Umbral de log2 Fold Change")
parser$add_argument('--run_kegg', type="logical", default=FALSE, help="Activar análisis de enriquecimiento")
parser$add_argument('--kegg_organism_code', type="character", help="Código de organismo para gprofiler (ej. 'athaliana')")
parser$add_argument('--kegg_padj_threshold', type="double", default=0.05, help="Umbral de p-valor ajustado en enriquecimiento")
args <- parser$parse_args()

# --- 3. INSTALACIÓN Y CARGA DE PAQUETES ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}
packages_to_install <- c("DESeq2", "AnnotationDbi", "GenomicFeatures", "Rsamtools", "GenomicAlignments",
                         "vsn", "ggplot2", "pheatmap", "ggrepel", "dplyr", "tibble",
                         args$organism_db, "gprofiler2")
cat("📦 Verificando dependencias...\n")
for (pkg in packages_to_install) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste("   -> Instalando", pkg, "...\n"))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({ invisible(lapply(packages_to_install, library, character.only = TRUE)) })

# --- 4. INICIO DEL ANÁLISIS ---
contrast_name <- paste0(args$treatment_group, "_vs_", args$control_group_contrast)
cat(paste("\n✅ Análisis iniciado para el contraste:", contrast_name, "\n"))
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 5. MATRIZ DE CONTEOS ---
if (tolower(args$counting_method) == "summarizeoverlaps") {
  cat("✅ Usando summarizeOverlaps...\n")
  txdb <- makeTxDbFromGFF(args$gtf_file, format="gtf")
  exons_by_gene <- exonsBy(txdb, by="gene")
  metadata_temp <- read.csv(args$metadata_file, row.names = 1, header=TRUE)
  bam_files <- file.path(args$bam_dir, paste0(rownames(metadata_temp), "_Aligned.sortedByCoord.out.bam"))
  bam_list <- BamFileList(bam_files, yieldSize=2000000)
  se <- summarizeOverlaps(features=exons_by_gene, reads=bam_list, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)
  matriz_conteos <- assay(se)
} else {
  cat("✅ Usando featureCounts...\n")
  matriz_completa <- read.table(args$counts_file, header = TRUE, sep = "\t", comment.char = "#", check.names = FALSE)
  matriz_conteos <- as.matrix(matriz_completa[, 7:ncol(matriz_completa)])
  rownames(matriz_conteos) <- matriz_completa$Geneid
  colnames(matriz_conteos) <- gsub("_Aligned.sortedByCoord.out.bam$", "", basename(colnames(matriz_conteos)))
}

# --- 6. CONTROL DE CALIDAD ---
cat("🔍 Realizando QC inicial...\n")
metadata <- read.csv(args$metadata_file, header = TRUE)
library_sizes <- colSums(matriz_conteos)

grouping_variable <- tail(all.vars(as.formula(args$design_formula)), n = 1)
stats_df <- data.frame(Muestra = names(library_sizes)) %>% 
  left_join(metadata, by = setNames(colnames(metadata)[1], "Muestra")) %>%
  mutate(Conteos_Totales = library_sizes[Muestra])
full_stats_table <- stats_df %>%
  group_by(.data[[grouping_variable]]) %>%
  mutate(Media_del_Grupo = round(mean(Conteos_Totales),0),
         Mediana_del_Grupo = round(median(Conteos_Totales),0),
         Desv_Estandar_Grupo = round(sd(Conteos_Totales),0)) %>%
  ungroup()

write.table(full_stats_table,
            file = file.path(args$output_dir, "estadisticas_conteos_crudos.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pdf(file.path(args$output_dir, "histograma_conteos_crudos_log10.pdf"))
hist(log10(as.vector(matriz_conteos) + 1), breaks=50,
     main="Histograma de Conteos Crudos (log10)", xlab="log10(Conteos + 1)", col="lightblue")
dev.off()

pdf(file.path(args$output_dir, "boxplot_conteos_crudos_log10.pdf"), width = 12, height = 8)
par(mar = c(10,4,4,2))
boxplot(log10(matriz_conteos + 1), main="Boxplot de Conteos Crudos por Muestra", las=2, col="lightblue")
par(mar = c(5,4,4,2))
dev.off()
cat("✅ QC inicial completado.\n")

# --- 7. DESeq2 ---
rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1, drop=FALSE]
matriz_conteos <- matriz_conteos[, rownames(metadata)]

dds <- DESeqDataSetFromMatrix(countData = matriz_conteos, colData = metadata, design = as.formula(args$design_formula))
dds <- dds[rowSums(counts(dds)) >= 10,]
dds[[grouping_variable]] <- relevel(dds[[grouping_variable]], ref = args$control_group)
dds <- DESeq(dds)

res <- results(dds, contrast=c(grouping_variable, args$treatment_group, args$control_group_contrast), alpha = args$padj_threshold)
res_ordered <- res[order(res$padj),]
res_ordered_df <- as.data.frame(res_ordered)
res_ordered_df$gene_id <- rownames(res_ordered_df)

if (args$strip_gene_version) {
  ids_mapeo <- sub("\\..*$", "", res_ordered_df$gene_id)
} else {
  ids_mapeo <- res_ordered_df$gene_id
}
res_ordered_df$symbol <- mapIds(get(args$organism_db),
                                keys=ids_mapeo,
                                column="SYMBOL",
                                keytype=args$key_type,
                                multiVals="first")

write.table(res_ordered_df,
            file = file.path(args$output_dir, paste0("resultados_completos_", contrast_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# --- 8. GENES SIGNIFICATIVOS ---
res_sig <- subset(res_ordered_df, padj < args$padj_threshold & abs(log2FoldChange) >= args$log2fc_threshold)
write.table(res_sig,
            file = file.path(args$output_dir, paste0("genes_significativos_", contrast_name, ".txt")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# --- 9. PCA ---
rld <- vst(dds, blind=FALSE)
pdf(file.path(args$output_dir, paste0("PCA_", contrast_name, ".pdf")))
plotPCA(rld, intgroup=grouping_variable)
dev.off()

# --- 10. MA PLOT ---
pdf(file.path(args$output_dir, paste0("MAplot_", contrast_name, ".pdf")))
plotMA(res, ylim=c(-5,5))
dev.off()

# --- 11. VOLCANO PLOT ---
pdf(file.path(args$output_dir, paste0("Volcano_", contrast_name, ".pdf")))
ggplot(res_ordered_df, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=(padj < args$padj_threshold & abs(log2FoldChange) >= args$log2fc_threshold)), alpha=0.5) +
  scale_color_manual(values=c("grey", "red")) +
  theme_minimal() +
  labs(title=paste("Volcano Plot -", contrast_name),
       x="Log2 Fold Change", y="-log10(padj)")
dev.off()

# --- 12. HEATMAP ---
if (nrow(res_sig) >= 2) {
  topgenes <- head(order(res$padj), 30)
  pdf(file.path(args$output_dir, paste0("Heatmap_", contrast_name, ".pdf")))
  pheatmap(assay(rld)[topgenes,], cluster_rows=TRUE, show_rownames=TRUE,
           cluster_cols=TRUE, annotation_col=as.data.frame(colData(rld)[, grouping_variable, drop=FALSE]))
  dev.off()
}

# --- 13. ENRIQUECIMIENTO FUNCIONAL ---
if (args$run_kegg) {
  cat("🔬 Ejecutando análisis de enriquecimiento con gprofiler2...\n")
  
  significant_genes <- res_sig$gene_id
  
  if (length(significant_genes) > 5) {
    if (args$strip_gene_version) {
        significant_genes <- sub("\\..*$", "", significant_genes)
    }
    
    gost_result <- gost(query = significant_genes,
                        organism = args$kegg_organism_code,
                        ordered_query = FALSE,
                        multi_query = FALSE,
                        significant = TRUE,
                        exclude_iea = FALSE,
                        measure_underrepresentation = FALSE,
                        evcodes = TRUE,
                        user_threshold = args$kegg_padj_threshold,
                        correction_method = "g_SCS",
                        sources = c("GO:BP", "GO:MF", "KEGG", "REAC"))

    if (!is.null(gost_result) && nrow(gost_result$result) > 0) {
      
      results_df <- gost_result$result
      is_list_col <- sapply(results_df, is.list)
      
      for (col in names(results_df)[is_list_col]) {
        results_df[[col]] <- sapply(results_df[[col]], function(x) paste(x, collapse = ","))
      }
      
      output_file_path <- file.path(args$output_dir, paste0("analisis_funcional_", contrast_name, ".txt"))
      write.table(results_df, output_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
      
      cat(paste("✅ Análisis de enriquecimiento guardado en:", args$output_dir, "\n"))

      # Copiar el fichero de resultados a la carpeta existente ../../R_CODES
      r_codes_dir <- file.path(args$output_dir, "../../R_CODES")
      file.copy(from = output_file_path, to = r_codes_dir, overwrite = TRUE)
      
      cat(paste("✅ Copia del resultado de enriquecimiento guardada en:", r_codes_dir, "\n"))

    } else {
      cat("ℹ️ No se encontraron términos significativos en GO, KEGG o Reactome.\n")
    }
  } else {
    cat("ℹ️ No hay suficientes genes significativos para realizar el análisis de enriquecimiento.\n")
  }
}

cat(paste("🎉 Análisis completo para", contrast_name, "finalizado.\n"))
