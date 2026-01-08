#A pallete of colors for each VDJ gene.

#here, the df I'm using is called primeira_leva_monoclonais
#what I'm doing is to create new Dfs filtering them according to the graph I'm plotting. #for example, for the d_call donut, I created a df called d_call_high only with lines that contain any value in d_call:

library(dplyr)
library(ggplot2)
library(patchwork)

# Filter IGL/IGK and remove NAs in j_call
j_call_igl <- vdj_igl %>%
  filter(!is.na(j_call))

# Adapted function for j_call donut chart using j_call_light_colors palette
donut_plot_j_light <- function(data, title) {
  data <- data %>%
    count(j_call) %>%
    mutate(percent = n / sum(n) * 100,
           label = paste0(j_call, " (", round(percent, 1), "%)"),
           ymax = cumsum(percent),
           ymin = c(0, head(ymax, n = -1)))
  
  ggplot(data, aes(ymax = ymax, ymin = ymin, xmax = 3, xmin = 2, fill = j_call)) +
    geom_rect(color = "black", size = 0.3, alpha = 0.85) +
    coord_polar(theta = "y") +
    xlim(c(1, 3)) +
    theme_void() +
    scale_fill_manual(values = j_call_light_colors) +
    labs(title = title) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(color = "black")
    )
}

# Create the three plots
donut_j_igl_binder <- donut_plot_j_light(j_call_igl %>% filter(bind_assay == "Binder"), "Binder")
donut_j_igl_non_binder <- donut_plot_j_light(j_call_igl %>% filter(bind_assay == "Non Binder"), "Non Binder")
donut_j_igl_non_specific <- donut_plot_j_light(j_call_igl %>% filter(bind_assay == "non specific"), "Non specific")

# Arrange plots side-by-side
donut_j_igl_binder + donut_j_igl_non_binder + donut_j_igl_non_specific


v_call_colors <- c(
  # IGHV1 - green
  "IGHV1-2*02"   = "#ADDD8E",
  "IGHV1-18*01"  = "#238443",
  "IGHV1-18*04"  = "#41AB5D",
  "IGHV1-24*01"  = "#78C679",
  "IGHV1-3*01"   = "darkgreen",
  "IGHV1-3*02"   = "#D9F0A3",
  "IGHV1-46*01"  = "#A6DBA0",
  "IGHV1-46*03"  = "#66C2A4",
  "IGHV1-58*02"  = "#2CA25F",
  "IGHV1-69*01"  = "#006D2C",
  "IGHV1-69*02"  = "#00441B",
  "IGHV1-69*06"  = "#005A32",
  "IGHV1-69*09"  = "#006837",
  "IGHV1-69*10"  = "#007F3B",
  "IGHV1-8*01"   = "#00994C",
  
  # IGHV2 - orange
  "IGHV2-26*01"  = "#FDAE6B",
  "IGHV2-5*01"   = "#F16913",
  "IGHV2-5*02"   = "#D94801",
  "IGHV2-70*04" = "#FE7743",
  
  # IGHV3 - blue
  "IGHV3-11*01"  = "#08306B",
  "IGHV3-11*04"  = "lightblue",
  "IGHV3-11*05"  = "#2171B5",
  "IGHV3-13*01"  = "#4292C6",
  "IGHV3-15*01"  = "#1C6BA0",
  "IGHV3-15*07"  = "#9ECAE1",
  "IGHV3-20*04"  = "#C6DBEF",
  "IGHV3-21*01"  = "#DEEBF7",
  "IGHV3-21*02"  = "#BDD7E7",
  "IGHV3-21*06"  = "#9ECAE1",
  "IGHV3-23*01"  = "#0B3C5D",
  "IGHV3-23*04"  = "#4292C6",
  "IGHV3-30*02"  = "#2171B5",
  "IGHV3-30*03"  = "#BDC7C6",
  "IGHV3-30*07"  = "#08306B",
  "IGHV3-30*09"  = "#0C3B66",
  "IGHV3-30*14"  = "#115E99",
  "IGHV3-30*18"  = "#1D7DB2",
  "IGHV3-30-3*01"= "#F0E5B2",
  "IGHV3-30-5*03"= "#41B6C4",
  "IGHV3-33*01"  = "#4B87BF",
  "IGHV3-33*06"  = "#99D8C9",
  "IGHV3-43*02"  = "lightblue",
  "IGHV3-48*02"  = "blue",
  "IGHV3-48*03"  = "darkblue",
  "IGHV3-49*02"  = "cyan",
  "IGHV3-49*03"  = "#D1D2B3",
  "IGHV3-49*04"  = "#3479B0",
  "IGHV3-53*01"  = "#104E8B",
  "IGHV3-53*02"  = "#08306B",
  "IGHV3-53*04"  = "#2171B5",
  "IGHV3-64*07"  = "#6BAED6",
  "IGHV3-64D*06" = "#9ECAE1",
  "IGHV3-64D*08" = "#C6DBEF",
  "IGHV3-66*02"  = "#DEEBF7",
  "IGHV3-7*01"   = "#08306B",
  "IGHV3-7*03"   = "#2171B5",
  "IGHV3-7*05"   = "#6BAED6",
  "IGHV3-73*01"  = "#9ECAE1",
  "IGHV3-73*02"  = "#C6DBEF",
  "IGHV3-74*01"  = "#DEEBF7",
  "IGHV3-74*03"  = "#F7FBFF",
  "IGHV3-9*01"   = "#7BA3DD",
  "IGHV3-30*20" = "#333446",
  "IGHV3-66*04" = "#7F8CAA",
  "IGHV3-21*03" = "#B8CFCE",
  "IGHV3-30*19" = "#EAEFEF",
  
  # IGHV4 - red
  "IGHV4-28*01"  = "#99000D",
  "IGHV4-30-2*01"= "#CB181D",
  "IGHV4-30-4*01"= "#EF3B2C",
  "IGHV4-30-4*09"= "#FB6A4A",
  "IGHV4-34*01"  = "#FC9272",
  "IGHV4-34*02"  = "#FCBBA1",
  "IGHV4-38-2*01"= "#FEE0D2",
  "IGHV4-39*01"  = "#67000D",
  "IGHV4-39*02"  = "#A50F15",
  "IGHV4-39*07"  = "#D7301F",
  "IGHV4-4*02"   = "#EF3B2C",
  "IGHV4-4*07"   = "#FB6A4A",
  "IGHV4-4*10"   = "#FC9272",
  "IGHV4-59*01"  = "#FCBBA1",
  "IGHV4-59*12"  = "#FEE0D2",
  "IGHV4-61*01"  = "#67000D",
  "IGHV4-61*02"  = "#A50F15",
  "IGHV4-61*12"  = "#D7301F",
  "IGHV4-59*02"  = "#E55050",
  "IGHV4-31*02"  = "#732255",
  "IGHV4-61*08"   = "#E7F2E4",
  "IGHV4-38-2*02" = "#E69DB8",
  "IGHV4-30-2*04" = "#FFD0C7",
  "IGHV4-30-4*08" = "#7D0A0A",
  
  
  # IGHV5 - golden/yellow
  "IGHV5-10-1*03"= "#FFB000",
  "IGHV5-51*01"  = "#FFC300",
  "IGHV5-51*03"  = "#FFD700",
  
  # IGHV6 - lilac
  # (no IGHV6 listed, but prepared in case you want to add it)
  
  # IGHV7 - fuchsia
  "IGHV7-4-1*02" = "#E377C2"
)

# Color palette with distinct shades for each family
d_call_colors <- c(
  # IGHD1 (red shades)
  "IGHD1-1*01"       = "#E48080",  # light red
  "IGHD1-14*01"      = "#E58D8D",  # softer shade
  "IGHD1-20*01"      = "#D56B6B",  # deeper shade
  "IGHD1-26*01"      = "#F5C8C8",  # pinkish beige
  "IGHD1/OR15-1a*01" = "#E3B3B3",  # soft pink
  "IGHD1-7*01" = "darkred",
  
  # IGHD2 (blue shades)
  "IGHD2-15*01"      = "#A6C8E3",  # Soft blue
  "IGHD2-2*01"       = "#7CA9D5",  # Medium blue
  "IGHD2-2*02"       = "#6097C0",  # Dark blue
  "IGHD2-2*03"       = "#4D83A6",  # Navy blue
  "IGHD2-21*01"      = "#2C6F8A",  # Intense blue
  "IGHD2-21*02"      = "#396B8A",  # Soft blue
  "IGHD2-8*01"       = "#50728F",  # Greenish blue
  "IGHD2-8*02"       = "#3E6C85",  # Deep blue
  "IGHD2/OR15-2a*01" = "#517A9E",  # Light blue
  
  # IGHD3 (gray shades)
  "IGHD3-10*01"      = "#B0B0B0",  # Light gray
  "IGHD3-10*02"      = "#A5B68D",  # Medium gray
  "IGHD3-10*03"      = "#A08963",
  "IGHD3-16*01"      = "#B17F59",  # Gray
  "IGHD3-16*02"      = "#6E6E6E",  # Dark gray
  "IGHD3-16*03"      = "#EDE8DC",  # Soft gray
  "IGHD3-22*01"      = "#D1D1D1",  # Light gray
  "IGHD3-3*01"       = "#8A8A8A",  # Dark gray
  "IGHD3-3*02"       = "#A1A1A1",  # Medium gray
  "IGHD3-9*01"       = "#F2E2B1",  # Soft gray
  "IGHD3/OR15-3a*01" = "#E3D2C3",  # Very light gray
  
  # IGHD4 (green shades)
  "IGHD4-11*01"      = "#A8D8A8",  # Light green
  "IGHD4-17*01"      = "#88C688",  # Medium green
  "IGHD4-23*01"      = "#6DAA6D",  # Dark green
  "IGHD4/OR15-4a*01" = "#76A076",  # Moss green
  
  # IGHD5 (yellow shades)
  "IGHD5-12*01"      = "#E1D8A6",  # Light yellow
  "IGHD5-18*01"      = "#D0C687",  # Golden yellow
  "IGHD5-18*02"      = "#B8B35A",  # Dark yellow
  "IGHD5-24*01"      = "#C8C481",  # Soft yellow
  "IGHD5/OR15-5a*01" = "#C8B85F",  # Mustard yellow
  
  # IGHD6 (purple shades)
  "IGHD6-13*01"      = "#A84B94",  # Dark purple
  "IGHD6-19*01"      = "#D366B2",  # Medium purple
  "IGHD6-25*01"      = "#C186C2",  # Light purple
  "IGHD6-6*01"       = "#BC58A9",  # Soft purple
  
  # IGHD7 (orange shades)
  "IGHD7-27*01"      = "#F4A300"   # Orange
)


j_call_colors <- c(
  # IGHJ4 (Beige and greenish tones)
  "IGHJ4*02" = "#D9DFC6",  # Pale greenish beige
  "IGHJ4*01" = "beige",    # Beige
  
  # IGHJ5 (Gray tones)
  "IGHJ5*02" = "#A8A8A8",  # Medium gray
  "IGHJ5*01" = "lightgrey", # Light gray
  
  # IGHJ6 (Brown and bronze tones)
  "IGHJ6*02" = "brown",    # Standard brown
  "IGHJ6*03" = "#6E4B3A",  # Dark earthy brown
  "IGHJ6*04" = "#9C6F30",  # Bronze / Strong gold
  
  # IGHJ1, IGHJ2 and IGHJ3 (Golden tones)
  "IGHJ1*01" = "#D9A83A",  # Golden yellow
  "IGHJ3*02" = "#E6C08F",  # Pale gold / Sand
  "IGHJ2*01" = "#D0A354"   # Golden bronze
)

v_call_light_colors <- c(
  v_call_light_colors <- c(
    # IGKV1 - green shades
    "IGKV1-16*01" = "#a6dba0",
    "IGKV1-12*01" = "#b8e186",
    "IGKV1-5*03" = "#7fbf7b",
    "IGKV1-5*01" = "#a1d99b",
    "IGKV1-39*01" = "#74c476",
    "IGKV1D-8*01" = "#66c2a4",
    "IGKV1-8*01" = "#4daf4a",
    "IGKV1-33*01" = "#238b45",
    "IGKV1-9*01" = "#006d2c",
    "IGKV1-8*03" = "#41ab5d",
    "IGKV1-6*01" = "#78c679",
    "IGKV1-13*02" = "#31a354",
    "IGKV1-17*01" = "#2ca25f",
    "IGKV1-16*02" = "#006d2c",
    "IGKV1D-13*01" = "#4d9c3b",
    "IGKV1D-12*01" = "#62b246",
    "IGKV1D-8*02" = "#7ec850",
    
    # IGKV2 - yellow/golden tones
    "IGKV2-24*01" = "#ffe082",
    "IGKV2-30*01" = "#ffd54f",
    "IGKV2-28*01" = "#ffca28",
    "IGKV2D-29*02" = "#ffc107",
    "IGKV2-30*02" = "#ffb300",
    
    # IGKV3 - light blue shades
    "IGKV3-15*01" = "#a6cee3",
    "IGKV3-11*01" = "#1f78b4",
    "IGKV3-20*01" = "#6baed6",
    "IGKV3D-15*01" = "#4292c6",
    
    # IGKV4 - dark blue/gray-blue shades
    "IGKV4-1*01" = "#9ecae1",
    "IGKV4-1*03" = "#6baed6",
    "IGKV4-1*02" = "#4292c6",
    
    # IGLV1 - orange shades
    "IGLV1-51*01" = "#fdae6b",
    "IGLV1-40*01" = "#fd8d3c",
    "IGLV1-44*01" = "#f16913",
    "IGLV1-47*01" = "#d94801",
    "IGLV1-51*02" = "#fdd0a2",
    
    # IGLV2 - pinkish/salmon shades
    "IGLV2-14*01" = "#fcbba1",
    "IGLV2-8*01" = "#fc9272",
    "IGLV2-23*04" = "#fb6a4a",
    "IGLV2-23*03" = "#ef3b2c",
    "IGLV2-14*03" = "#cb181d",
    "IGLV2-11*01" = "#fb9a99",
    "IGLV2-23*01" = "#e41a1c",
    "IGLV2-14*02" = "#de2d26",
    "IGLV2-18*02" = "#f768a1",
    "IGLV2-18*01" = "#dd3497",
    "IGLV2-23*02" = "#ae017e",
    "IGLV2-5*01" = "#7a0177",
    
    # IGLV3 - purple/purplish-blue shades
    "IGLV3-9*01" = "#dadaeb",
    "IGLV3-21*02" = "#bcbddc",
    "IGLV3-1*01" = "#9e9ac8",
    "IGLV3-21*04" = "#807dba",
    "IGLV3-19*01" = "#6a51a3",
    "IGLV3-21*03" = "#54278f",
    "IGLV3-10*01" = "#756bb1",
    "IGLV3-21*01" = "#9f9fd3",
    "IGLV3-25*03" = "#b39ddb",
    
    # IGLV4 - soft neutral shade
    "IGLV4-69*01" = "#bdbdbd",
    
    # IGLV5 - soft violet shade
    "IGLV5-45*02" = "#cab2d6",
    
    # IGLV6 - olive green shades
    "IGLV6-57*02" = "#c7e9c0",
    "IGLV6-57*01" = "#a1d99b",
    
    # IGLV7 - bronze/light brown shades
    "IGLV7-46*01" = "#d8b365",
    "IGLV7-43*01" = "#b69e55",
    "IGLV7-46*04" = "#997d3b",
    
    # IGLV8 - light blue-green shade
    "IGLV8-61*01" = "#ccebc5",
    
    # IGLV9 - neutral pastel shade
    "IGLV9-49*01" = "#f0f0f0",
    
    # IGLV10 - sand shades
    "IGLV10-54*04" = "#d9d0c9",
    "IGLV10-54*01" = "#cbc4bc"
  )
) 
  
  j_call_light_colors <- c(
    "IGLJ1*01" = "#A89D91",  # Yellowish-gray for Lambda
    "IGLJ3*02" = "#B8AA9E",
    "IGLJ2*01" = "#C8B8AB",
    
    "IGKJ2*01" = "#D4A017",  # Strong gold for Kappa
    "IGKJ4*01" = "#E3B23C",
    "IGKJ1*01" = "#F0C75E",
    "IGKJ5*01" = "#F7D891",
    "IGKJ3*01" = "#FAE5B6",
    "IGKJ4*02" = "beige",
    "IGKJ2*03" = "#ffca28",
    "IGKJ2*02" = "#ffb300",
    "IGLJ6*01" = "#ffe082"
  )
  
  
  ##This code can be used to plot column graphs comparing the gene usage of BCRs across different groups. 
  
  library(ggplot2)
  library(dplyr)
  library(ggalluvial)
  library(forcats)
  
  # Prepare data
  v_usage <- merged_df %>%
    filter(!is.na(v_call)) %>%
    group_by(bind_assay, v_call) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(bind_assay) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    ungroup()
  
  # Order bind_assay so Binder appears in the middle
  v_usage$bind_assay <- factor(v_usage$bind_assay, levels = c("Non Binder", "Binder", "Non specific"))
  
  # Order v_call for prettier plotting
  v_usage$v_call <- fct_reorder(v_usage$v_call, v_usage$percent, .fun = sum, .desc = TRUE)
  
  # Prepare data for alluvial
  alluvial_data <- v_usage %>%
    mutate(axis = bind_assay)
  
  # Alluvial plot with Binder in the middle
  ggplot(alluvial_data,
         aes(x = axis, stratum = v_call, alluvium = v_call,
             y = percent, fill = v_call, label = v_call)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", color = "darkgray", alpha = 0.7) +
    geom_stratum(width = 0.3, color = "black") +
    scale_fill_manual(values = v_call_colors) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = "V gene usage comparison: Non Binder → Binder → Non specific")
  
  ## OR ALTERNATIVELY WITH NON SPECIFIC in the middle instead of Binders
  # Filter 'aPD-L1' group and NAs in v_call
  vdj_igl_filtrado_v <- vdj_igl %>%
    filter(bind_assay != "aPD-L1", !is.na(v_call))
  
  # Calculate V gene usage
  v_usage <- vdj_igl_filtrado_v %>%
    group_by(bind_assay, v_call) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(bind_assay) %>%
    mutate(percent = 100 * count / sum(count)) %>%
    ungroup()
  
  # Define bind_assay order: Binder -> Non specific -> Non Binder
  v_usage$bind_assay <- factor(v_usage$bind_assay,
                               levels = c("Binder", "non specific", "Non Binder"))
  
  # Order v_call by usage
  v_usage$v_call <- fct_reorder(v_usage$v_call, v_usage$percent, .fun = sum, .desc = TRUE)
  
  # Prepare for plotting
  alluvial_data_v <- v_usage %>%
    mutate(axis = bind_assay)
  
  # Generate plot
  ggplot(alluvial_data_v,
         aes(x = axis, stratum = v_call, alluvium = v_call,
             y = percent, fill = v_call, label = v_call)) +
    geom_flow(stat = "alluvium", lode.guidance = "forward", color = "darkgray", alpha = 0.7) +
    geom_stratum(width = 0.3, color = "black") +
    scale_fill_manual(values = v_call_light_colors) +
    theme_minimal(base_size = 14) +
    theme(
      axis.title = element_blank(),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    ) +
    labs(title = "VL gene usage comparison: Binder → Non specific → Non Binder")
  
  
  
  ### Palette for families
  IGHV_colors <- c(
    "IGHV1" = "#E3B23C",
    "IGHV2" = "#ffb300",
    "IGHV3" = "#997d3b",
    "IGHV4" = "brown",
    "IGHV5" = "#6E4B3A",
    "IGHV6" = "#E6C08F",
    "IGHV7" = "#ffe082"
  )
  
  IGHD_colors <- c(
    "IGHD1" = "#a6cee3",
    "IGHD2" = "#4292c6",
    "IGHD3" = "#E1D8A6",
    "IGHD4" = "#E48080",
    "IGHD5" = "#E3B3B3",
    "IGHD6" = "#08306B",
    "IGHD7" = "#DEEBF7"
  )
  
  IGHJ_colors <- c(
    "IGHJ1" = "#fdd0a2",
    "IGHJ2" = "#fd8d3c",
    "IGHJ3" = "#d94801",
    "IGHJ4" = "#a6dba0",
    "IGHJ5" = "beige",
    "IGHJ6" = "#76A076"
  )
  