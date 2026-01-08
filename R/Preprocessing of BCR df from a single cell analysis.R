#in case your df also contains TCR sequences, we are first going to remove them and keep only BCRs:
library(dplyr)
#our raw dataframe is called 'vdj'

bcr <- vdj %>% filter(locus %in% c("IGH", "IGK", "IGL"))
library(dplyr)

#Now, we are going to remove the double light chains from each cell_id:

# Separate light and heavy chains
bcr_light <- bcr %>% filter(locus %in% c("IGK", "IGL"))
bcr_heavy <- bcr %>% filter(locus == "IGH")

# Identify cell_id with more than one light chain
bcr_light_filtered <- bcr_light %>%
  group_by(cell_id) %>%
  filter(n() == 2) %>%  # only cell_ids with 2 light chains
  slice_max(order_by = umi_count, n = 1, with_ties = FALSE) %>%  # Keep the one with the most UMIs
  ungroup()

# Now also get cell_ids with a single light chain (no conflict)
bcr_light_single <- bcr_light %>%
  group_by(cell_id) %>%
  filter(n() == 1) %>%
  ungroup()

# Merge everything: heavy chain + single light chain per cell
bcr_final <- bind_rows(bcr_heavy, bcr_light_filtered, bcr_light_single)

#now, let's remove those lines where the quality of the sequence is low:

bcr_final <- bcr_final %>% 
  filter(high_quality_cell_tcr_bcr == TRUE)