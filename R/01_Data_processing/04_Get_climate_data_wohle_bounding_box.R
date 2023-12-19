#----------------------------------------------------------#
#
#         Latitudinal analysis of phylogenetic dispersion
#
#           Extract Chelsa climate for the study area ----
#                          
#----------------------------------------------------------#
# Source config file ----
source("R/00_Config_file.R")

# Get raw Chelsa data for entire study area at 1 degree resolution
min_age <- 0
max_age <- 12e3
timestep <- 100

data_meta <-
 base::expand.grid(
    lat = seq(25, 66, by = 1),
    long = seq(75, 125,by = 1)
    )

time_reference_table <- 
  readr::read_rds(
    "Inputs/Data/Chelsa_climate/time_reference_table.rds"
    ) %>% 
  dplyr::filter(endyear > -12001)


# 'tasmin', and the 'bio_var_selected' are additional here, will reduce time to 
#  download and extract a lot of data not needed.

raw_data <- 
  get_climate_data(
    variables_selected = c("bio", "tasmin"),
    bio_var_selected = c(1, 6, 12, 15, 18, 19),
    time_var_selected = c(20:-200),
    month_var_selected = c(1:12),
    xy = data_meta
    )

data_climate <-
  get_climate_indices(
    data_source = raw_data,
    time_ref = time_reference_table) %>%
  tidyr::unnest(climate_data)

readr::write_rds(
  data_climate,
  file = "Inputs/Data/Chelsa_climate/Chelsa_climate_whole_area.rds"
  )
