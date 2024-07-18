#---------------------------------------------------------------------------------#
#
# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia during 
#                                the Holocene
#
#                             Bhatta et al. 2024
#
#---------------------------------------------------------------------------------#

#----------------------------------------------------------#
#
# Make a matrix of data with 1 degree latitude and 200-yr temporal resolution ----
#                          
#----------------------------------------------------------#

# Source config file ----
source("R/00_Config_file.R")

# Climate ----
data_climate <-
  readr::read_rds("Inputs/Data/Chelsa_climate/Chelsa_climate_whole_area.rds") %>%
  dplyr::group_by(data_id) %>%
  tidyr::nest() %>% 
  dplyr::ungroup()

data_meta <-
  base::expand.grid(
    lat = seq(25, 66, by = 1),
    long = seq(75, 125, by = 1)
    ) %>%
  rownames_to_column("data_id") %>%
  tibble::as_tibble()

climate_matrix <- 
  dplyr::inner_join(
    data_meta,
    data_climate, 
    by = "data_id"
    ) %>% 
  dplyr::select(
    -c(data_id, long)
    ) %>% 
  dplyr::mutate(
    data = purrr::map(
      data, 
      ~ .x %>% 
        dplyr::select(-time_id) %>%
        dplyr::arrange(age) %>%
        dplyr::filter(!age > 12000) %>% 
        dplyr::mutate(
          age = ceiling(
            age / 200)*200
          ) %>%
        dplyr::group_by(age) %>%
        dplyr::summarise_all(., mean) %>% 
        dplyr::ungroup()
      )
    ) %>% 
  tidyr::unnest(data) %>% 
  dplyr::group_by(lat, age) %>% 
  dplyr::summarise_all(., mean) %>% 
  dplyr::ungroup()

  
# Number of samples ----
source_data <- 
  read_rds("Inputs/Data/source_data_191223.rds")

blank_dat <- 
  base::expand.grid(
    lat = seq(25, 66, by = 1),
    age = seq(200, 12000, by = 200)
    )

n_samples <-
  source_data %>%
  dplyr::select(lat, age) %>%
  dplyr::arrange(lat) %>%
  dplyr::mutate(
    lat = ceiling(lat - 1 / 1) * 1
    ) %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(
    age = ceiling(
      age / 200) * 200
    ) %>%
  dplyr::group_by(lat, age) %>%
  dplyr::summarise(n_samples = n()) %>% 
  dplyr::ungroup() 

n_samples_matrix <-
  dplyr::left_join(
    blank_dat,
    n_samples,
    by = c("lat", "age")
    ) 

# Age uncertainty ----
age_uncertainty <- 
  source_data %>% 
  dplyr::select(
    lat, 
    age,
    upper,
    lower,
    age_uncertainty_index
    ) %>% 
   dplyr::mutate(
    age_error = abs(
      lower - upper
      )
    ) %>% 
  dplyr::arrange(lat) %>%
  dplyr::mutate(
    lat = ceiling(lat - 1 / 1) * 1
    ) %>%
  dplyr::arrange(age) %>%
  dplyr::mutate(
    age = ceiling(age / 200) * 200
    ) %>% 
  dplyr::group_by(lat, age) %>% 
  dplyr::summarise_all(., mean) %>% 
  dplyr::ungroup()

age_uncertainty_matrix <- 
  dplyr::left_join(
    blank_dat, 
    age_uncertainty,
    by = c("lat", "age")
    ) %>%
  dplyr::select(
    lat,
    age,
    age_error, 
    age_uncertainty_index
    ) 
  
#----------------------------------------------------------#
# Import Holocene-wide phylogenetic metrics ----
#----------------------------------------------------------#
output_gam_pd <- 
  readr::read_rds("Outputs/Data/Overall_gam_lat_PD_201223.rds") 

mpd_matrix <- 
  output_gam_pd[1,] %>% 
  dplyr::select(predicted_gam) %>% 
  tidyr::unnest(predicted_gam) %>% 
  dplyr::select(-dataset_id) %>% 
  dplyr::rename(ses_MPD = var) %>% 
  dplyr::group_by(lat, age) %>% 
  dplyr::summarise_all(., mean) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    lat, 
    age, 
    ses_MPD
    )

mntd_matrix <- 
  output_gam_pd[2,] %>% 
  dplyr::select(predicted_gam) %>% 
  tidyr::unnest(predicted_gam) %>% 
  dplyr::select(-dataset_id) %>% 
  dplyr::rename(ses_MNTD = var) %>% 
  dplyr::group_by(lat, age) %>% 
  dplyr::summarise_all(., mean) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    lat, 
    age, 
    ses_MNTD
    )

phylo_matrix <- 
  dplyr::inner_join(
    mpd_matrix, 
    mntd_matrix,
    by = c("lat", "age")
    )

full_matrix <-
  dplyr::inner_join(
    climate_matrix, 
    n_samples_matrix,
    by = c("lat", "age")
    ) %>%
  dplyr::inner_join(
    age_uncertainty_matrix,
    by = c("lat", "age")
    ) %>%
  dplyr::inner_join(
    phylo_matrix,
    by = c("lat", "age")
    )

#----------------------------------------------------------#
# Export full matrix ----
#----------------------------------------------------------#
readr::write_rds(
  full_matrix,
  file = "Outputs/Data/Space_time_matrix_201223.rds",
  compress = "gz"
  )
