# First, it's important that the format of the matrix is correct. Here, we use columns as the Abs, and rows as the targets of each ELISA. example:
#             Ab1     Ab2     Ab3
# PD-L1    3          1          1
# Spike     1          1          1 
# CD73     1          1          1
# Isotype IgA1   IgG2     IgM

AUC_matrix <- as.data.frame(AUC_matrix)

# Set first column as rownames
rownames(AUC_matrix) <- AUC_matrix[,1]
# Remove the first column
AUC_matrix <- AUC_matrix[,-1]

library(pheatmap)

# Extract the last row (Isotype)
isotype_row <- as.character(AUC_matrix[nrow(AUC_matrix), ]) 
AUC_matrix_num <- AUC_matrix[-nrow(AUC_matrix), ]  # Remove last row to retain only numeric values

# Convert main matrix to numeric
AUC_matrix_num <- apply(AUC_matrix_num, 2, as.numeric)  
rownames(AUC_matrix_num) <- rownames(AUC_matrix)[-nrow(AUC_matrix)]

# Create annotation dataframe for Isotype
annotation_df <- data.frame(Isotype = isotype_row)
rownames(annotation_df) <- colnames(AUC_matrix_num)

# Define custom color palette for isotypes
isotype_colors <- c("IgA1" = "darkgreen", "IgA2" = "lightgreen", 
                    "IgG1" = "lightcoral", "IgG2" = "red", 
                    "IgG3" = "darkred", "IgM" = "orange", "IgD" = "orange")

# Convert Isotype column to factor with correct levels
annotation_df$Isotype <- factor(annotation_df$Isotype, levels = names(isotype_colors))

# Create the list of colors for annotations
annotation_colors <- list(Isotype = isotype_colors)

# Reorder matrix columns based on Isotypes
col_order <- order(annotation_df$Isotype)  # Order by Isotype factor levels
AUC_matrix_num <- AUC_matrix_num[, col_order] # Reorder the matrix
annotation_df <- annotation_df[col_order, , drop = FALSE]  # Reorder the annotation

# Define color palette for numeric values
breaks <- seq(0, 4, length.out = 100)  
color_palette <- colorRampPalette(c("beige", "lightgrey", "grey", "lightblue", "blue", "darkblue"))(length(breaks) - 1)

# Create heatmap with Isotype annotation and grouped columns
pheatmap(AUC_matrix_num, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = color_palette,  
         breaks = breaks,  
         annotation_col = annotation_df,  
         annotation_colors = annotation_colors,  
         main = "Monoclonals specificity (AUC)",
         angle_col = 0,  
         fontsize_col = 14, 
         fontsize_row = 16, border_color = "black")

# Create a new color sequence: blue -> white -> red
color_palette <- colorRampPalette(c("blue", "white","pink", "red"))(length(breaks) - 1)


###updated heatmap with all binders using a df called heatmap_binders
library(dplyr)
library(pheatmap)
library(RColorBrewer)

# Select the numeric variable
heatmap_data <- heatmap_binders %>%
  select(cell_id, PD_L1_AUC) %>%
  column_to_rownames("cell_id")

# Prepare categorical annotations
annotation_data <- heatmap_binders %>%
  select(cell_id, Sample_Name, Isotype) %>%
  column_to_rownames("cell_id")

# Generate color palettes for annotations
# Sample_Name: defined previously
sample_colors <- c(
  "Pul117" = "#66C2A5",
  "Pul119" = "#FC8D62",
  "Pul122" = "#8DA0CB",
  "Pul123" = "#E78AC3"
)

# Full annotations
annotation_colors <- list(
  Sample_Name = sample_colors,
  Isotype = isotype_colors
)

# Plot the heatmap

pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_data,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("white", "firebrick3"))(50),
  show_rownames = FALSE,
  fontsize = 12,
  border_color = "black",  
  linewidths = 0.2         
)

#or if we want to plot it in horizontal
pheatmap(
  t(heatmap_data),  
  cluster_rows = FALSE,  
  cluster_cols = TRUE,   
  annotation_col = annotation_data,  
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("white", "firebrick3"))(50),
  show_colnames = FALSE,  
  fontsize = 12,
  border_color = "black",
  linewidths = 0.2
)

#FOR THE POLYCLONAL HEATMAP
library(pheatmap)
library(dplyr)
library(tibble)

# Filter patients with at least one positive AUC
heatmap_filtered <- heatmap %>%
  filter(`aPD-1` > 0 | `aPD-L1` > 0 | `aCTLA-4` > 0)

# Prepare numeric values for the heatmap
heatmap_data <- heatmap_filtered %>%
  select(Patient_ID = 1, `aPD-1`, `aPD-L1`, `aCTLA-4`) %>%
  column_to_rownames("Patient_ID")

# Prepare categorical annotation (Stage, Autoimmunity, Death)
annotation_data <- heatmap_filtered %>%
  select(Patient_ID = 1, Stage, Autoimmunity, Death) %>%
  mutate(Stage = as.character(Stage),
         Autoimmunity = as.character(Autoimmunity),
         Death = as.character(Death)) %>%
  mutate(Stage = ifelse(is.na(Stage), "NA", Stage)) %>%
  mutate(Stage = ifelse(Stage == "IA", "IA1", Stage)) %>%
  column_to_rownames("Patient_ID")

# Define color palettes for categories
stage_colors <- c(
  "IA1" = "#7CFC00",
  "IA2" = "#32CD32",
  "IA3" = "#228B22",
  "IB"  = "#006400",
  "IIA" = "#FFFACD",
  "IIB" = "#FFD700",
  "IIIA" = "#FFA500",
  "IIIB" = "#FF8C00",
  "IIIC" = "#FF7F50",
  "IVA" = "firebrick",
  "IVB" = "#8B0000",
  "NA"  = "gray"
)

autoimmunity_colors <- c(
  "Yes" = "purple",
  "No" = "gray"
)

death_colors <- c(
  "Not Deceased" = "lightyellow",
  "Deceased" = "black",
  "Deceased (not from cancer)" = "gray"
)

annotation_colors <- list(
  Stage = stage_colors,
  Autoimmunity = autoimmunity_colors,
  Death = death_colors
)

# Verify that all values are mapped
check_missing_colors <- function(annotation_data, annotation_colors) {
  for(col in colnames(annotation_data)) {
    valores <- unique(annotation_data[[col]])
    cores <- names(annotation_colors[[col]])
    faltando <- setdiff(valores, cores)
    if(length(faltando) > 0) {
      cat("❗ Coluna:", col, "- valores sem cor definida:\n")
      print(faltando)
      stop("Adicione esses valores ao vetor de cores antes de prosseguir.")
    }
  }
}
check_missing_colors(annotation_data, annotation_colors)

# Generate the heatmap
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_data,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#91C8E4", "white", "pink" , "red", "firebrick", "darkred"))(100),
  breaks = seq(min(heatmap_data, na.rm = TRUE), max(heatmap_data, na.rm = TRUE), length.out = 101),
  show_rownames = TRUE,
  fontsize = 9,
  border_color = "black",
  linewidths = 0.2,
  main = "Heatmap de AUCs por paciente com anotações categóricas"
)

# HEATMAP VERSION 2

# Define colors for categories
stage_colors <- c(
  "I" = "#228B22",
  "II" = "#FFD700",
  "III" = "#FF8C00",
  "IV" = "firebrick",
  "NA"  = "gray"
)

autoimmunity_colors <- c(
  "Yes" = "purple",
  "No" = "gray"
)

death_colors <- c(
  "Not Deceased" = "lightyellow",
  "Deceased" = "black",
  "Deceased (not from cancer)" = "gray"
)

annotation_colors <- list(
  Stage = stage_colors,
  Autoimmunity = autoimmunity_colors,
  Death = death_colors
)

# Verify that all values are mapped
check_missing_colors <- function(annotation_data, annotation_colors) {
  for(col in colnames(annotation_data)) {
    valores <- unique(annotation_data[[col]])
    cores <- names(annotation_colors[[col]])
    faltando <- setdiff(valores, cores)
    if(length(faltando) > 0) {
      cat("❗ Coluna:", col, "- valores sem cor definida:\n")
      print(faltando)
      stop("Adicione esses valores ao vetor de cores antes de prosseguir.")
    }
  }
}
check_missing_colors(annotation_data, annotation_colors)

# Generate the heatmap
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_data,
  annotation_colors = annotation_colors,
  color = colorRampPalette(c("#91C8E4", "white", "pink" , "red", "firebrick", "darkred"))(100),
  breaks = seq(min(heatmap_data, na.rm = TRUE), max(heatmap_data, na.rm = TRUE), length.out = 101),
  show_rownames = TRUE,
  fontsize = 9,
  border_color = "black",
  linewidths = 0.2,
  main = "Heatmap de AUCs por paciente com anotações categóricas"
)

## SET ZERO VALUES TO GRAY


# Define limits and number of colors
max_val <- max(heatmap_data_positive, na.rm = TRUE)

# Number of colors for the continuous part
n_breaks <- 100

# Define breaks: start with 0, then a sequence above 0
breaks <- c(
  0,                                           
  seq(1, max_val, length.out = n_breaks)     
)

# Define colors: one for zero, the rest for the gradient
colors <- c(
  "darkgrey",  
  colorRampPalette(c("white", "pink", "red", "firebrick", "darkred"))(n_breaks)
)

# Create the heatmap
pheatmap(
  heatmap_data_positive,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_positive,
  annotation_colors = annotation_colors,
  color = colors,
  breaks = breaks,
  show_rownames = TRUE,
  fontsize = 9,
  border_color = "black",
  linewidths = 0.2,
  main = "Heatmap com cor especial para valores zero"
)


#THE ONLY METHOD THAT WORKED: USING COMPLEXHEATMAP

library(ComplexHeatmap)

# Create annotation with two categorical variables
ha <- HeatmapAnnotation(
  Stage = annotation_positive$Stage,
  Autoimmunity = annotation_positive$Autoimmunity,
  col = list(
    Stage = stage_colors,
    Autoimmunity = autoimmunity_colors
  )
)

# Plot heatmap with multiple annotations
Heatmap(
  t(heatmap_data_positive),
  col = colors,
  row_names_side = "left",
  top_annotation = ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = "black", lwd = 0.5)
)
