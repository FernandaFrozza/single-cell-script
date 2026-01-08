
# Assuming the dataframe is named 'df' and contains the following columns: 
# log2FoldChange, padj, and group (where group is based on log2FoldChange: 
# if >0 = Positive, if <0 = Negative)
library(deseq2)

df$group <- case_when(
  df$padj < 0.05 & df$log2FoldChange > 1.25  ~ "Positive",
  df$padj < 0.05 & df$log2FoldChange < -1.25 ~ "Negative",
  TRUE ~ "NS"
)

# Simplified volcano plot WITHOUT labels
ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = group)) +
  geom_point(size = 2.5, alpha = 1, shape = 16) +
  scale_color_manual(values = c("Positive" = "#DC143C", 
                                "Negative" = "#DDDAD0", 
                                "NS" = "darkgrey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1.25, 1.25), linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    x = "log2 Fold Change (Positive vs Negative)",
    y = "-log10(p-adjusted)",
    color = ""
  )