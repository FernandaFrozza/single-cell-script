#Right after we process our single cell or bulk data, we will have a tsv file with all the information regarding 
#our BCRs/TCRs. Before using Immcantation, we have to clean this table.

# First thing, we have to remove double light chains from the same cell_id. it's pretty common that some cells 
#will have counts/umis to more than one light chain, probably because of gene editing during maturation. 
#We consider the real light chain that one thas has more counts and more UMIs. It is usually clear the difference; 
#one light chain will have a lot more UMIs and counts than the other. we can use this code to filter these chains for us.

#let's call our df "bcr"

library(dplyr)

# Identify cell_ids with more than one light chain
cell_ids_com_duplicadas <- bcr %>%
  filter(locus %in% c("IGK", "IGL")) %>%
  count(cell_id) %>%
  filter(n > 1) %>%
  pull(cell_id)

# Conflicting light chains -- keep only the one with the highest umi_count
cadeias_leves_filtradas <- bcr %>%
  filter(cell_id %in% cell_ids_com_duplicadas, locus %in% c("IGK", "IGL")) %>%
  group_by(cell_id) %>%
  slice_max(order_by = umi_count, n = 1, with_ties = FALSE) %>%
  ungroup()

# Light chains without conflict -- keep all
cadeias_leves_unicas <- bcr %>%
  filter(!(cell_id %in% cell_ids_com_duplicadas), locus %in% c("IGK", "IGL"))

# Other chains (IGH or other non-light chains)
outras_cadeias <- bcr %>%
  filter(!(locus %in% c("IGK", "IGL")))

# 4. Merge all parts
bcr_filtrado <- bind_rows(cadeias_leves_filtradas, cadeias_leves_unicas, outras_cadeias)