#FOR GENE ONTOLOGY, we filtered only immune pathways with the following keywords:
  immune_keywords <- c("leukocyte", "immune", "antigen", "presentation",
                       "T cell", "B cell", "macrophage", "neutrophil",
                       "lymphocyte", "cytokine", "chemokine", "migration", "adhesion")
  
  
  
library(ggplot2)
library(ggrepel)
library(dplyr)
  
plot_volcano_control_labels_tumor_vs_blood <- function(df, max_labels = 15) {
    
  # Create group column if it doesn't exist
    if (!"group" %in% colnames(df)) {
      df <- df %>%
        mutate(group = case_when(
          avg_log2FC > log2(3) & p_val_adj < 0.05 ~ "Up in Tumor",
          avg_log2FC < -log2(3) & p_val_adj < 0.05 ~ "Up in Blood",
          TRUE ~ "Not significant"
        ))
    }
    
  # Filter genes for labeling
    sig_genes <- df %>%
      filter(abs(avg_log2FC) > log2(3), p_val_adj < 0.05) %>%
      arrange(p_val_adj)
    
    # Limit to max_labels
    sig_genes <- head(sig_genes, max_labels)
    
    # Define colors
    colors <- c("Up in Tumor" = "blue", "Up in Blood" = "red", "Not significant" = "gray80")
    
    # Plot
    p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
      geom_point(alpha = 0.8, size = 1.5) +
      scale_color_manual(values = colors) +
      theme_minimal(base_size = 12) +
      theme(panel.grid = element_blank(), legend.position = "top") +
      labs(title = "Volcano plot: Tumor vs Blood (B cells)",
           x = "Log2 Fold Change",
           y = "-Log10 Adjusted P-value",
           color = "Expression")
    
    # Add labels
    p + geom_text_repel(
      data = sig_genes,
      aes(label = Gene),
      size = 3.5,
      max.overlaps = max_labels,
      box.padding = 0.5,
      point.padding = 0.4,
      segment.color = 'gray50'
    )
  }
  
# Usage example:
# plot_volcano_control_labels_tumor_vs_blood(DEGS_p005_tumor_vs_blood, max_labels = 15)

## for the binder heatmap with tumor and dln-like genes
# Gradient white -> firebrick
  heatmap_colors <- colorRampPalette(c("#FEF3E2", "#F3C623", "#FA812F","firebrick", "#4E1F00"))(100)
  
  # Colors for categories
  ann_colors <- list(Category = c("tumor-like" = "#2973B2", "dLNs-like" = "#99BC85"))
  
  pheatmap(expr_merged,
           annotation_row = annotation_genes,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 9,
           fontsize_col = 9,
           main = "Heatmap of tumor-like and dLNs-like DEGs",
           color = heatmap_colors,
           border_color = "black")
  