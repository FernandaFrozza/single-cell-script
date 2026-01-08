#description of colors and code to plot BCR isotypes according to each patient

#consider a dataframe already preprocessed only with the heavy chain of your BCRs. Sample_Name is the column with each patient name, and Isotype is the column with your Igs.

# Define custom color palette for isotypes
isotype_colors <- c("IgA1" = "darkgreen", "IgA2" = "lightgreen", 
                    "IgG1" = "lightcoral", "IgG2" = "red", 
                    "IgG3" = "darkred", "IgG4" = "#67000D", "IgM" = "orange")

# Count cells per Sample_Name and Isotype
isotype_counts <- igh %>%
  count(Sample_Name, Isotype)

# Create the plot with the same visual settings
ggplot(isotype_counts, aes(x = Sample_Name, y = n, fill = Isotype)) +
  geom_bar(stat = "identity", alpha = 0.4, color = "black", size = 0.1) + 
  scale_fill_manual(values = isotype_colors) +
  labs(x = "", y = "Number of cells", fill = "Isotype") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black")
  )