#All single cell files are currently being downloaded from seven bridges after alignment of the fastqs with the reference genome and sample multiplexing. 
#Seurat files are .rds and this script works for R.

# Load packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(readxl)
library(openxlsx)
library(scales)
library(gridExtra)
library(presto)
library(sctransform)
library(glmGamPoi)
library(tidyr)

# Set updated file path
file_path <- "C:/caminho/do/arquivo.rds"

# Load Seurat object
seurat <- readRDS(file_path)

#Before preprocessing, let's take a look in the overall quality of this seurat object**

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
summary(seurat$percent.mt)

hist(seurat$percent.mt, xlab = "percent.mt", main = "Histogram of MT genes")

# Plot distributions of nFeature_RNA, nCount_RNA, and percent.mt
plot.cellQC = VlnPlot(seurat,
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                      pt.size = 0,
                      ncol = 3)

# Plot scatter plots of different features
QCdotplot1 = FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt") + 
  theme(legend.position = "none")
QCdotplot2 = FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")
QCdotplot3 = FeatureScatter(seurat, feature1 = "percent.mt", feature2 = "nFeature_RNA") + 
  theme(legend.position = "none")

# Adjust violin plot parameters
nFeature.per.library = VlnPlot(seurat, features = "nFeature_RNA", pt.size = 0) + 
  theme(legend.position = "none")
nCount.per.library = VlnPlot(seurat, features = "nCount_RNA", pt.size = 0) + 
  theme(legend.position = "none")

##To save in PDF all the plots above**

# Set PDF output path
pdf("C:/Users/Usuario/Documents/Rh5_singlecell_data/QC_plots.pdf", width = 12, height = 8)

# Plot violin and scatter plots
print(plot.cellQC)
print(QCdotplot1)
print(QCdotplot2)
print(QCdotplot3)
print(nFeature.per.library)
print(nCount.per.library)

# Close PDF file
dev.off()


# Now, let's clean this data! We are using a function that excludes cells that do 
# not fit in these parameters: number of genes and mRNA transcripts 5x the SD of the median, 
# and cells with a percentage of mitochondrial genes 3x SD above the median

# Function to identify outliers based on a specific metric
is_outlier <- function(SeuratObject, metric, nmads) {
  
  eval(parse(text = paste0("M <- SeuratObject$", metric)))
  
  # Identify cells that meet outlier criteria
  outlier <- (M < median(M) - nmads * mad(M)) | (M > median(M) + nmads * mad(M))
  
  return(outlier)
}

# Identify statistical outliers across specific quality control metrics
check_outliers_nFeature <- is_outlier(seurat, 'nFeature_RNA', 5)
check_outliers_nCount <- is_outlier(seurat, 'nCount_RNA', 5)
check_outliers_mito <- is_outlier(seurat, 'percent.mt', 3)

# Extract identifiers for cells meeting the inclusion criteria across all metrics
non_outliers_nFeature <- names(check_outliers_nFeature[!check_outliers_nFeature])
non_outliers_nCount <- names(check_outliers_nCount[!check_outliers_nCount])
non_outliers_mito <- names(check_outliers_mito[!check_outliers_mito])

# Identify cells meeting consensus inclusion criteria across all QC parameters
non_outliers <- intersect(non_outliers_nFeature, non_outliers_nCount) %>%
  intersect(non_outliers_mito)

# Refine the dataset by preserving only the validated, non-outlier cell population
seurat <- subset(seurat, cells = non_outliers)

#Now, let's normalize this data to a log of 10,000**

seurat_filtrado <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#We have to exclude the 'Multiplet' and 'Undetermined' tags from Sample_Name now

# Create a vector of cells to be retained
cels_para_manter <- rownames(seurat_filtrado@meta.data[!seurat_filtrado@meta.data$Sample_Name %in% c("Multiplet", "Undetermined"), ])

# Update the Seurat object to retain only the selected cells
seurat_filtrado <- subset(seurat_filtrado, cells = cels_para_manter)

#Since 'Multiplet' and 'Undetermined' categories now have zero counts, unused factor levels can be dropped to refine the metadata.
table(seurat_filtrado@meta.data$Sample_Name)

seurat_filtrado@meta.data$Sample_Name <- droplevels(seurat_filtrado@meta.data$Sample_Name)

table(seurat_filtrado@meta.data$Sample_Name)

##The data is normalized and clean, we can now look at the markers**

# Specify main marker genes for T cells, B cells, macrophages, myeloid cells, and dendritic cells
genes_t_cells <- c("CD3D", "CD3E", "CD4", "CD8A")
genes_b_cells <- c("CD19", "MS4A1", "CD79A", "CD79B")
genes_macrophages <- c("CD14", "CSF1R", "MARCO", "MRC1")
genes_myeloid_cells <- c("LYZ", "ITGAM", "CD11B", "CD86")
genes_dendritic_cells <- c("CD11C", "ITGAX", "CLEC9A", "BATF3")

# Visualization of T cell marker expression (t-SNE)
FeaturePlot(seurat_filtrado, features = genes_t_cells, reduction = "tsne")

# Visualization of B cell marker expression (t-SNE)
FeaturePlot(seurat_filtrado, features = genes_b_cells, reduction = "tsne")

DimPlot(seurat_filtrado, reduction = "tsne", group.by = "Cell_Type_Experimental") +
  ggtitle("")

# Visualization of macrophage marker expression (t-SNE)
FeaturePlot(seurat_filtrado, features = genes_macrophages, reduction = "tsne")

# Visualization of myeloid cell marker expression (t-SNE)
FeaturePlot(seurat_filtrado, features = genes_myeloid_cells, reduction = "tsne")

# Visualization of dendritic cell marker expression (t-SNE)
FeaturePlot(seurat_filtrado, features = genes_dendritic_cells, reduction = "tsne")

# For UMAP, you can use the same code by replacing reduction = "tsne" with reduction = "umap"

##To plot a bar graph and see the differences of cell types across sample_names, use this**

# Identify the desired cell types for analysis
cell_types <- c("B", "Dendritic", "Monocyte_classical", "Monocyte_nonclassical", 
                "Natural_killer", "T_CD4_memory", "T_CD4_naive", 
                "T_CD8_memory", "T_CD8_naive", "T_gamma_delta")

# Subset the metadata to retain only the target cell populations
filtered_data <- seurat_filtrado@meta.data[seurat_filtrado$Cell_Type_Experimental %in% cell_types, ]

# Quantify cell abundance across cell types and samples
cell_type_counts <- table(filtered_data$Sample_Name, filtered_data$Cell_Type_Experimental)

# Convert counts to a data frame format
cell_type_counts_df <- as.data.frame(cell_type_counts)

# Rename columns for improved interpretability
colnames(cell_type_counts_df) <- c("Sample_Name", "Cell_Type_Experimental", "Count")

# Calculate total cell counts per experimental Sample_Name
total_cells_per_sample <- table(filtered_data$Sample_Name)

# Calculate the percentage of each cell type per Sample_Name
cell_type_percentages <- prop.table(table(filtered_data$Sample_Name, filtered_data$Cell_Type_Experimental), 1) * 100

# Convert percentages to a data frame format
cell_type_percentages_df <- as.data.frame(cell_type_percentages)

# Rename columns for the percentage data frame
colnames(cell_type_percentages_df) <- c("Sample_Name", "Cell_Type_Experimental", "Percentage")

library(ggplot2)

# Create a bar plot for percentages faceted by Cell_Type_Experimental
ggplot(cell_type_percentages_df, aes(x = Sample_Name, y = Percentage, fill = Sample_Name)) +
  geom_bar(stat = "identity", position = "dodge") +  
  theme_minimal() +
  labs(title = "Porcentagem de Populações de Células por 'Sample_Name'",
       x = "Sample Name",
       y = "Porcentagem (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  
  scale_fill_brewer(palette = "Set3") +  
  facet_wrap(~ Cell_Type_Experimental, scales = "free_y") 


#We have to integrate our datasets now. first, we create a list of the objects that will be integrated
lista_objetos <- list(rh3_rh8, rh5, rh9)

#and then we select 3000 features to group these datasets
features_integracao <- SelectIntegrationFeatures(object.list = lista_objetos, nfeatures = 3000)
anchors <- FindIntegrationAnchors(object.list = lista_objetos, anchor.features = features_integracao, dims = 1:30)
seurat_integrado <- IntegrateData(anchorset = anchors, dims = 1:30)

DefaultAssay(seurat_integrado) <- "integrated"

seurat_integrado <- ScaleData(seurat_integrado, verbose = FALSE)
seurat_integrado <- RunPCA(seurat_integrado, npcs = 30, verbose = FALSE)
seurat_integrado <- RunUMAP(seurat_integrado, reduction = "pca", dims = 1:30)
seurat_integrado <- FindNeighbors(seurat_integrado, dims = 1:30)
seurat_integrado <- FindClusters(seurat_integrado, resolution = 0.5)
seurat_integrado <- RunTSNE(seurat_integrado, dims = 1:30)

DimPlot(seurat_integrado, reduction = "umap", group.by = "Cell_Type_Experimental")


##now I will filter only B cells and other subtypes that have a VH molecule count to a new seurat so I can analyze B cell subtypes.

b_or_bcr_cells <- subset(
  seurat_integrado,
  subset = Cell_Type_Experimental == "B" | BCR_Heavy_Molecule_Count > 0
)
table(b_or_bcr_cells$Cell_Type_Experimental)
summary(b_or_bcr_cells$BCR_Heavy_Molecule_Count)

saveRDS(b_or_bcr_cells, file = "C:/Users/Usuario/Documents/Rh5_singlecell_data/B_or_BCR_cells.rds")

# 1. Set assay to RNA (if it is still set to 'integrated')
DefaultAssay(b_cells_filtered) <- "RNA"

# 2. Re-normalize if needed (optional)
# b_cells_filtered <- NormalizeData(b_cells_filtered)

# 3. Identify variable genes
b_cells_filtered <- FindVariableFeatures(b_cells_filtered)

# 4. Scale the data
b_cells_filtered <- ScaleData(b_cells_filtered)

# 5. PCA
b_cells_filtered <- RunPCA(b_cells_filtered)

# 6. Determine number of PCs with ElbowPlot (optional)
ElbowPlot(b_cells_filtered)

# 7. Clustering
b_cells_filtered <- FindNeighbors(b_cells_filtered, dims = 1:20)
b_cells_filtered <- FindClusters(b_cells_filtered, resolution = 0.6)

# 8. t-SNE
b_cells_filtered <- RunTSNE(b_cells_filtered, dims = 1:20)

# 9. UMAP
b_cells_filtered <- RunUMAP(b_cells_filtered, dims = 1:20)

DimPlot(b_cells_filtered, reduction = "tsne", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("B cells - Clusters via t-SNE")

saveRDS(b_cells_filtered, "C:/Users/Usuario/Documents/Rh5_singlecell_data/B_cells_reclustered.rds")

genes_b <- c("CXCR5", "BCL6", "MS4A1", "JCHAIN", "GZMB", "CD274",
             "IGHD", "FCRL1", "IGHG1", "IGHA1", "IGHA2", "IGHM", "IL2RG")

DotPlot(b_cells_filtered, features = genes_b) +
  RotatedAxis() +
  labs(title = "Expressão de genes em subtipos de células B")

FeaturePlot(b_cells_filtered, features = c("JCHAIN", "GZMB", "CD274"), reduction = "tsne")
FeaturePlot(b_cells_filtered, features = c("IGHD", "IGHM", "FCRL1"), reduction = "tsne")
FeaturePlot(b_cells_filtered, features = c("IGHG1", "IGHA1", "IGHA2"), reduction = "tsne")

# Define gene lists organized by signature
gene_lists <- list(
  Plasmablast = c("IGHG1","IGHG2","IGHG3", "IGHA1","IGHA2", "PRDM1", "CD38", "FNDC3B", "XBP1", "MZB1"),
  Plasma_cell_reg = c("JCHAIN", "GZMB", "XBP1", "CXCR3"),
  Naive = c("IGHD", "IGHM", "FCRL1"),
  Germinal_center = c("IGHM", "CXCR5", "CD40", "HLA-DRB1", "CD80", "CD86", "FAS"),
  Immature = c("CD2", "IL7R", "FOXP1"),
  Memory = c("IGHG1","IGHG2","IGHG3", "IGHA1","IGHA2", "CD27", "CD24"),
  Marginal_Zone = c("IGHM", "CD24", "CXCR5", "CD40"),
  Atypical_memory = c("ZEB2", "TOX", "CD22", "BCL2", "CXCR5", "CCR7")
)

# Flatten the list and create a vector of unique genes
all_genes <- unique(unlist(gene_lists))

# Create a DotPlot with all genes
DotPlot(b_cells_filtered, features = all_genes, group.by = "Cell_Type",scale.max = 75) +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "")

# Filter only "Blood" cells
blood_cells <- subset(b_cells_filtered, subset = Sample == "Blood")


# Group by Specificity and Cell_Type and calculate proportions
prop_blood <- blood_cells@meta.data %>%
  group_by(Specificity, Cell_Type) %>%
  summarise(Freq = n(), .groups = "drop") %>%
  group_by(Specificity) %>%
  mutate(Proportion = Freq / sum(Freq))

ggplot(prop_blood, aes(x = Specificity, y = Proportion, fill = Cell_Type)) +
  geom_col(position = "stack", color = "black", size = 0.2, alpha = 0.6) +
  labs(
    title = "Proporção de Cell_Type nas células Blood",
    x = "Specificity", y = "Proporção"
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 16),
    plot.title = element_text(color = "black", face = "bold", size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12)
  )

#PACKAGES I USED TO INVESTIGATE DEGs BETWEEN GROUPS:
#findmarkers from Seurat
#I put the groups as Identities to do this analysis



## Plot composition differences between groups (integrated Seurat)
library(ggplot2)

# Filter only Tumor Pre Neo and Tumor Pos Neo
df_plot <- seurat_integrado@meta.data[
  seurat_integrado$Treatment %in% c("Tumor Pre Neo", "Tumor Pos Neo"), 
]

# Stacked histogram in proportions (%)
ggplot(df_plot, aes(x = Treatment, fill = Cell_Type_Experimental)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.3) +  
  scale_y_continuous(labels = scales::percent) +                  
  scale_fill_brewer(palette = "Pastel1") +                        
  theme_minimal(base_family = "Arial") +                          
  theme(
    panel.grid = element_blank(),       
    axis.title = element_text(size = 16, face = "bold", color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "Treatment",
    y = "Proporção de células",
    fill = "Cell Type (Experimental)"
  )


#### HIGHLIGHT A SPECIFIC CELL GROUP

# Define cells of interest
cells_blood_pos_neo <- WhichCells(seurat_obj, expression = Treatment == "Blood Pos Neo")

# Generate FeaturePlot highlighting only selected cells
FeaturePlot(
  seurat_obj,
  features = "PDCD1",
  cells.highlight = cells_blood_pos_neo,
  cols = c("lightgrey", "red"), # grey for background, red for the group
  pt.size = 1
)