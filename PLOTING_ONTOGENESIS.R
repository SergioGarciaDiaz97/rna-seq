# ===================================================================
# --- R Script for Batch Enrichment Analysis Visualization ---
# This script finds ALL enrichment results files in its directory
# and generates a corresponding summary plot for each one.
# ===================================================================

# --- 1. Package Installation ---
required_packages <- c("ggplot2", "dplyr", "stringr")
packages_to_install <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(packages_to_install)) {
  cat("Installing required packages:", paste(packages_to_install, collapse=", "), "\n")
  install.packages(packages_to_install)
}

# Load the necessary libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
})


# --- 2. Automatic File Finder ---
cat("INFO: Searching for functional analysis results files...\n")
results_filenames <- list.files(pattern = "analisis_funcional_.*\\.txt$")

# Check if at least one file was found
if (length(results_filenames) == 0) {
  stop("ERROR: No 'analisis_funcional...' text files found in this directory.")
}

cat(paste0("INFO: Found ", length(results_filenames), " result file(s) to process.\n\n"))


# --- 3. Loop Through Each File and Create a Plot ---
for (current_file in results_filenames) {
  
  cat(paste0("--- Processing file: ", current_file, " ---\n"))
  
  # --- a. Read and Process Data ---
  results_df <- read.table(current_file, header = TRUE, sep = "\t", quote = "")
  
  plot_data <- results_df %>%
    group_by(source) %>%
    slice_min(order_by = p_value, n = 10, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      term_name_short = str_trunc(term_name, 50),
      term_name_ordered = reorder(paste0(term_name_short, " (", term_id, ")"), p_value)
    )
  
  # --- b. Create the Plot ---
  enrichment_plot <- ggplot(plot_data, aes(x = intersection_size, y = term_name_ordered)) +
    geom_point(aes(size = intersection_size, color = p_value)) +
    facet_wrap(~ source, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "red", high = "blue", name = "P-value") +
    scale_size_continuous(name = "Gene Count") +
    theme_bw(base_size = 12) +
    labs(
      x = "Number of Genes in Term (Intersection Size)",
      y = "Term / Biological Pathway",
      title = "Functional Enrichment Analysis",
      subtitle = paste("Showing top 10 terms for:", gsub("analisis_funcional_|.txt", "", current_file))
    ) +
    theme(
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 10)
    )
  
  # --- c. Display and Save the Plot ---
  print(enrichment_plot)
  
  # Create a unique name for the output PDF based on the input filename
  output_filename <- gsub("\\.txt$", "_plot.pdf", current_file)
  ggsave(output_filename, plot = enrichment_plot, width = 12, height = 10)
  
  cat(paste0("✅ Plot saved as '", output_filename, "'\n\n"))
}

cat("--- All files processed. ---\n")
