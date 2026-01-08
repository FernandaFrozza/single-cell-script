#manual plotting of the region of mutations across different groups. 
#we basically compute a value of 1 nucleotide specific mutation when comparing the sequence_alignment column 
#with the germline_alignment column and calculate the mean frequency of mutations in each nucleotide for each group.


# First, according to Immcantation, the position of FR and CDR regions is around this (https://shazam.readthedocs.io/en/stable/topics/setRegionBoundaries/):
# obs: we are considering a mean junction length of 40, but you can change it according to your data:
#Região | Início | Fim
#-- | -- | --
#FR1 | 1 | 78
#CDR1 | 79 | 114
#FR2 | 115 | 165
#CDR2 | 166 | 195
#FR3 | 196 | 312
#CDR3 | 313 | 313 + 40 - 6 = 347
#FR4 | 348 | Fim da sequência

#we are now going to create a df with this info so we can color the graph later
region_rects <- data.frame(
  xmin = c(1, 79, 115, 166, 196, 313, 348),
  xmax = c(78, 114, 165, 195, 312, 347, 400),
  RegionType = c("FR1", "CDR1", "FR2", "CDR2", "FR3", "CDR3", "FR4")
)

region_colors <- c(
  "FR1" = "#D3D3D3",
  "FR2" = "#D3D3D3",
  "FR3" = "#D3D3D3",
  "FR4" = "#D3D3D3",
  "CDR1" = "#FFD700",
  "CDR2" = "#FFA500",
  "CDR3" = "#FF6347"  # tomato red to highlight CDR3
)


#to add the junction_length column:
vdj <- vdj %>%
  mutate(junction_length = nchar(junction))

#you can also do it for the aminoacid alignment in case you want to look only at replacement mutations
#in the example below, I'm filtering only non-PD-L1 Memory B cells, but you can filter any group you prefer

vdj_memory_ns <- vdj %>%
  filter(`Cell Type` == "Peripheral Memory")

# Identify mutations ignoring "N" in the germline
vdj_memory_ns <- vdj_memory_ns %>%
  mutate(mut_flags = map2(sequence_alignment, germline_alignment, ~ {
    seq_vec <- strsplit(.x, "")[[1]]
    germ_vec <- strsplit(.y, "")[[1]]
    if (length(seq_vec) != length(germ_vec)) return(rep(NA, length(seq_vec)))
    sapply(seq_along(seq_vec), function(i) {
      if (germ_vec[i] == "N") NA else as.integer(seq_vec[i] != germ_vec[i])
    })
  }))

# Expand to one row per position
vdj_long_ns <- vdj_memory_ns %>%
  mutate(seq_length = map_int(mut_flags, length)) %>%
  unnest_longer(mut_flags, indices_to = "Position") %>%
  filter(!is.na(mut_flags)) %>%
  mutate(Position = as.integer(Position))

# Calculate mean mutations per position by Isotype
plot_df_ns <- vdj_long_ns %>%
  group_by(Position, Isotype) %>%
  summarise(mean_mut = mean(mut_flags), .groups = "drop")


# Create the plot
ggplot() +
  # Background bands by region
  geom_rect(data = region_rects, 
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = RegionType),
            alpha = 0.3) +
  scale_fill_manual(values = region_colors, guide = "none") +
  
  # Add mutation bars
  geom_col(data = plot_df_ns, aes(x = Position, y = mean_mut), fill = "brown") +
  
  # Facet by Isotype with fixed y-axis
  facet_wrap(~Isotype, scales = "fixed", ncol = 1) +
  ylim(0, 1) +
  
  # Themes
  theme_minimal(base_size = 14) +
  labs(
    x = "Nucleotide Position",
    y = "Mean Mutation Frequency",
    title = "Somatic Mutation Frequency (Total Memory)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
