#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Geographic distribution of the datasets to be analysed ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source config file ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
phylo_div <- 
  readr::read_rds("Inputs/Data/phylogenetic_diversity_estimated_191223.rds")

#--------------------------------------------------------#
# 3. Merge the similar climate zones ----
#--------------------------------------------------------#
dat_full1 <- 
  dat_full %>% 
  dplyr::mutate(
    Climate_zone = dplyr::case_when(
      ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
      ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
      ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
      ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
      ecozone_koppen_15 == "Polar_Frost" ~ "Polar frost",
      ecozone_koppen_15 == "Cold_Without_dry_season" ~ "Cold without dry season",
      ecozone_koppen_15 == "Tropical_Rainforest" ~ "Tropical rainforest",
      ecozone_koppen_15 == "Tropical_Monsoon" ~ "Tropical monsoon",
      ecozone_koppen_15 == "Tropical_Savannah" ~ "Tropical savannah",
      TRUE ~ ecozone_koppen_15
    )
  )  %>% 
  dplyr::filter(!long < 0) 
  

#--------------------------------------------------------#
# 3. Base map with 5 modified Köppen-Geiger climate zones (in Beck et al. 2018) ----
#--------------------------------------------------------#
# Read the raster points from the geo-tiff file published in Beck et al. 2018
raster_file <- 
  raster::raster("Inputs/Data/Biomes/Beck_tif/Beck_KG_V1_present_0p0083.tif")
  
# Read the raster value-climatic zone tranlation table
koppen_tranlation_table <-
  readr::read_csv("Inputs/Data/Biomes/Beck_tif/koppen_link.csv")
    
# Extract the required raster points
# Reduce the raster dimension (by factor of 10), otherwise would be too big file
raster_file1 <- raster::aggregate(raster_file, fact = 10) 

raster_df <- 
  # Convert raster points into a dataframe
  as.data.frame(raster_file1, xy = TRUE) %>% 
  
  # Extract the rater points of only required area
  dplyr::filter(x > 20 & y > 0) %>%  # for required lat, long 
  dplyr::rename(raster_values = Beck_KG_V1_present_0p0083) %>% 
  dplyr::mutate(raster_values = round(raster_values, digits = 0)) %>% 
  
  # Assign the names of climate zone to the raster values
  dplyr::left_join(koppen_tranlation_table, by = c("raster_values")) %>% 
  dplyr::filter(!raster_values == 0) %>% 
  dplyr::rename(
    ecozone_koppen_30 = genzone,
    ecozone_koppen_15 = genzone_cluster,
    ecozone_koppen_5 = broadbiome
    ) %>% 
  
  # Modify the climate zones as suggested by John
  dplyr::mutate(
    Climate_zone = dplyr::case_when(
      ecozone_koppen_15 == "Arid_Desert" ~ "Arid",
      ecozone_koppen_15 == "Arid_Steppe" ~ "Arid",
      ecozone_koppen_15 == "Cold_Dry_Summer" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Cold_Dry_Winter" ~ "Cold-seasonally dry",
      ecozone_koppen_15 == "Temperate_Dry_Summer" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Dry_Winter" ~ "Temperate",
      ecozone_koppen_15 == "Temperate_Without_dry_season" ~ "Temperate",
      ecozone_koppen_15 == "Polar_Tundra" ~ "Polar",
      ecozone_koppen_15 == "Polar_Frost" ~ "Polar frost",
      ecozone_koppen_15 == "Cold_Without_dry_season" ~ "Cold without dry season",
      ecozone_koppen_15 == "Tropical_Rainforest" ~ "Tropical rainforest",
      ecozone_koppen_15 == "Tropical_Monsoon" ~ "Tropical monsoon",
      ecozone_koppen_15 == "Tropical_Savannah" ~ "Tropical savannah",
      TRUE ~ ecozone_koppen_15
    )
  )       

# Assign unique colour to each climate zone
cbPalette1 <-
  c("#FFCC99",
    "#993300",
    "#FF6600",
    "#3399FF",
    "#999999",
    "#00CCCC",
    "#99CC00",
    "#006600",
    "#996600"
  ) %>% 
  rlang::set_names(
    nm = c(
      "Arid",
      "Cold-seasonally dry",
      "Cold without dry season",
      "Polar",
      "Polar frost",
      "Temperate",
      "Tropical monsoon",
      "Tropical rainforest",
      "Tropical savannah"
    )
  )

# Base map (modified climate zones)
base_map <-
  dat_full1 %>%
  ggplot2::ggplot(
    aes(
      x = long, 
      y = lat
      )
    ) +
  ggplot2::coord_fixed(
    ylim = c(20.00, 80.00),
    xlim = c(65.00, 135.00)
    ) +
  ggplot2::geom_tile(
    data = raster_df,
    aes(
      x = x, 
      y = y, 
      fill = Climate_zone
      ),
    inherit.aes = FALSE,
    alpha = 0.75
    ) +
  ggplot2::scale_fill_manual(values = cbPalette1) +
  labs(
    x = expression(
      paste('Longitude ', (degree ~ E))
      ),
    y = expression(
      paste('Latitude ', (degree ~ N))
      ),
    fill = "Climate zones"
    ) +
  ggplot2::theme_classic() +
  borders(
    colour = "black",
    linewidth = 0.2
    ) +
  guides(
    fill = guide_legend(
      nrow = 3,
      byrow = TRUE,
      title.position = "top"
    ),
    size = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  ) +
  
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.key.size = unit(0.60, "cm"),
    legend.title = element_text(
      size = 14
      ),
    legend.text = element_text(
      size = 10
      ),
    axis.title = element_text(
      color = "black", 
      size = 16
      ),
    axis.text = element_text(
      colour = "black", 
      size = 12
      ),
    plot.margin = margin(0, 0, 0.1, 0, "cm")
    ) # t, r, b, l

# Distribution of sequences across ecozones 
main_fig <- 
  base_map +
  ggplot2::geom_point(
    color = "#000099", 
    size = 3
    ) 

# Sequence density 
data_density <- 
  dat_full1 %>%
  dplyr::filter(!lat < 20 & !lat > 80) %>% 
  dplyr::filter(!long < 65 & !long > 135)

xbp <-
  ggpubr::gghistogram(
    data_density$long,
    fill = "#2CA388",
    color = "#2CA388",
    binwidth = 2.5,
    size = 0.1,
    alpha = 0.5) +
  ggpubr::theme_transparent() 

ybp <-
  ggpubr::gghistogram(
    data_density$lat,
    fill = "#2CA388",
    color = "#2CA388",
    binwidth = 2.5,
    size = 0.1,
    alpha = 0.5) +
  ggpubr::rotate() +
  ggpubr::theme_transparent() 

xbp_grob <-  ggplot2::ggplotGrob(xbp)
ybp_grob <-  ggplot2::ggplotGrob(ybp)

# combine figures
combined_fig <- 
  main_fig +
  ggplot2::annotation_custom(
    grob = xbp_grob,
    xmin = 55.5,
    xmax = 142,
    ymin = 14.5,
    ymax = 24.5) +
  ggplot2::annotation_custom(
    grob = ybp_grob,
    xmin = 59.1,
    xmax = 69,
    ymin = 15,
    ymax = 80)

# Add bounding box for filtering the records
filtered_data <-
  data_density %>%
  dplyr::filter(!phylogenetic_diversity == "NA") %>%
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) # 99 records

insetrect <- 
  data.frame(
    xmin = min(filtered_data$long) - 0.1, 
    xmax = max(filtered_data$long) + 0.1,
    ymin = min(filtered_data$lat) - 0.1,  
    ymax = max(filtered_data$lat) + 0.1
    ) 

final_fig <-
  combined_fig +
  ggplot2::geom_rect(
    data = insetrect,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax),
    alpha = 0,
    colour = "red",
    linewidth = 1.5,
    inherit.aes = FALSE
    )

ggplot2::ggsave(final_fig,
       file = "Outputs/Figure/Fileterd_datasets_191223.tiff", 
       height = 15, 
       width = 15, 
       units = "cm", 
       dpi = 400,
       compression = "lzw")

