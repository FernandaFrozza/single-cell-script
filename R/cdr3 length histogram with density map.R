ggplot() +
  # Bars
  geom_col(
    data = cdr3_distribution,
    aes(x = cdr3_aa_length, y = percentage, fill = Group),
    position = "identity",
    alpha = 0.6,
    color = "black"
  ) +
  # Density curves
  geom_density(
    data = vdj_cdr3_plot,
    aes(x = cdr3_aa_length, y = ..scaled.. * 40, color = Group),
    size = 1
  ) +
  scale_fill_manual(values = c("PD-L1 Binders" = "#578FCA", "Total Peripheral Memory" = "firebrick")) +
  scale_color_manual(values = c("PD-L1 Binders" = "#578FCA", "Total Peripheral Memory" = "firebrick")) +
  labs(
    x = "CDR3 length (aa)",
    y = "Percentage of BCRs (%)",
    fill = "",
    color = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.x = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 15),
    legend.text = element_text(color = "black", size = 13),
    legend.title = element_blank()
  )