#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Temporal variation in latitudinal pattern of phylogenetic dispersion ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Load configuration file ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
data_full <- 
  readr::read_rds(
    paste(
      "Inputs/Data/",
      "data_for_main_analysis_121023.rds",
      sep = ""
    )
  ) 

#--------------------------------------------------------#
# 3. Filter the datasets ----
#--------------------------------------------------------#
source_data_filtered <-
  data_full %>%
  dplyr::select(
    dataset_id,
    lat,
    long,
    phylodiversity_age_combined
    ) %>%
   tidyr::unnest(phylodiversity_age_combined) %>%
  dplyr::filter(age > 0) %>%
  dplyr::mutate(
    age_uncertainty_index = 
      mean(
        abs(lower - upper)
        ) / abs(lower - upper)
    ) 

#--------------------------------------------------------#
# 4. Load the climate data predicted for each sample of the datasets ----
#--------------------------------------------------------#
climate_data <- 
  readr::read_rds("Inputs/Data/Chelsa_climate/data_climate_pred-2021-12-16.rds") %>% 
  dplyr::filter(dataset_id %in% source_data_filtered$dataset_id) %>% 
  dplyr::select(
    dataset_id, 
    clim_data_pred
    ) %>% 
  tidyr::unnest(clim_data_pred)

#--------------------------------------------------------#
# 5. Combine the climate data and filtered phylodiversity data ----
#--------------------------------------------------------#
combined_data <-
  source_data_filtered %>%
  dplyr::select(
    dataset_id,
    lat,
    sample_id,
    age_uncertainty_index,
    ses_mpd = mpd_phylogeny_pool_abundance_wt,
    ses_mntd = mntd_phylogeny_pool_abundance_wt
    ) %>%
  dplyr::inner_join(
    climate_data,
    by = c(
      "dataset_id", 
      "sample_id"
      )
    ) %>%
  dplyr::filter(age > 0) %>%
  arrange(age) %>%
  dplyr::mutate(period = ceiling(age / 1000)) %>%
  dplyr::mutate(period = period * 1000) %>%
  dplyr::mutate_at("dataset_id", as_factor) %>%
  dplyr::mutate_at("period", as_factor)

#--------------------------------------------------------#
# 6. Fit and plot the GAM models ----
#--------------------------------------------------------#
# A. MPD/MNTD ----
data_gam_phylodiv <-
  combined_data %>%
  dplyr::select(dataset_id,
                sample_id,
                lat,
                age,
                age_uncertainty_index,
                period,
                ses_mpd,
                ses_mntd
                ) %>% 
  tidyr::gather(
    c(
      ses_mpd,
      ses_mntd),
    key = "vars",
    value = "estimate"
    ) %>%
  dplyr::group_by(vars) %>%
  tidyr::nest() %>%
  dplyr::ungroup() 

set.seed(2330)

gam_mod_phylodiv <-
  data_gam_phylodiv %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
        data <- .x
        mod <-
          mgcv::gam(
            estimate ~
              lat + 
              s(lat, 
                by = period, 
                bs = 'tp',
                m = 1) +
              s(age, 
                k = 10, 
                bs = "tp") +
              s(age,            
               by = period, 
               bs = 'tp',      
               m = 1) +
              s(dataset_id, 
                k = 99, 
                bs = 're') +
              s(period, 
                k = 12, 
                bs = 'fs') +
            ti(lat, age, 
               by = period,
               bs = c("tp", "tp")
               ),
            data = data,
            method = "REML",
            family = "gaussian",
            weights = age_uncertainty_index,
            control = gam.control(trace = TRUE, maxit = 200)
          )
        }
      )
    )

# B. Climate ----
data_gam_climate <-
  combined_data %>%
  dplyr::select(
    dataset_id,
    lat,
    age,
    age_uncertainty_index,
    period,
    temp_cold,
    prec_summer,
    prec_winter
    ) %>% 
  tidyr::gather(
    c(
      temp_cold,
      prec_summer,
      prec_winter
      ),
    key = "vars",
    value = "estimate"
    ) %>%
  dplyr::group_by(vars) %>%
  tidyr::nest() %>%
  dplyr::ungroup() 

set.seed(2330)

gam_mod_climate <-
  data_gam_climate %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
        data <- .x
        mod <-
          mgcv::gam(
            estimate ~
              lat + 
              s(lat, 
                by = period, 
                bs = 'tp',
                m = 1 ) +
              s(age, 
                k = 10, 
                bs = "tp") +
              s(age,            
                by = period, 
                bs = 'tp',      
                m = 1) +
              s(dataset_id, 
                k = 99, 
                bs = 're') +
              s(period, 
                k = 12, 
                bs = 'fs') +
              ti(lat, age, 
                 by = period,
                 bs = c('tp', 'tp')
              ),
            data = data,
            method = "REML",
            family = "gaussian",
            weights = age_uncertainty_index,
            control = gam.control(trace = TRUE, maxit = 200)
          )
        }
        )
    )

#--------------------------------------------------------#
# 7. Test of differences in slopes of models of different periods ----
#--------------------------------------------------------#
# Temporal variation in MPD and MNTD ----
set.seed(2330)

diff_pattern_phylodiv <-
  gam_mod_phylodiv %>%
  dplyr::mutate(
    new_data =
      purrr::map(
        data,
        ~ with(
          .x,
          base::expand.grid(
            lat = seq(
              min(lat),
              max(lat),
              by = 0.25
              ),
            period = seq(
              1000,
              12000,
              by = 1000
              )
            )
          )
        ),
    temporal_variation_phylodiv =
      purrr::map2(
        .x = gam_model,
        .y = new_data,
        .f = ~ {
          diff_pattern <-
            gratia::difference_smooths(
              .x,
              smooth = "s(lat)",
              newdata = .y,
              ci_level = 0.95,
              n = 300,
              method = "REML"
            )
          return(diff_pattern)
        }
      )
    )

# Temporal variation in climate ----
set.seed(2330)

diff_pattern_climate <-
  gam_mod_climate %>%
  dplyr::mutate(
    new_data =
      purrr::map(
        data, 
        ~ with(
          .x,
          base::expand.grid(
            lat = seq(
              min(lat),
              max(lat),
              by = 0.25
              ),
            period = seq(
              1000,
              12000,
              by = 1000
              )
            )
          )
        ),
    temporal_variation_climate =
      purrr::map2(
        .x = gam_model,
        .y = new_data,
        .f = ~ {
          diff_pattern <-
            gratia::difference_smooths(
              .x,
              smooth = "s(lat)",
              newdata = .y,
              ci_level = 0.95,
              n = 300,
              method = "REML"
            )
          return(diff_pattern)
        }
      )
    )

mpd_diff <- 
  diff_pattern_phylodiv[1,]$temporal_variation_phylodiv[[1]] %>% 
  dplyr::select(
    lat, 
    level_1, 
    level_2, 
    mpd_diff = diff
    )

mntd_diff <- 
  diff_pattern_phylodiv[2,]$temporal_variation_phylodiv[[1]] %>% 
  dplyr::select(
    lat, 
    level_1, 
    level_2, 
    mntd_diff = diff
    )

temp_cold <-
  diff_pattern_climate[1,]$temporal_variation_climate[[1]] %>% 
  dplyr::select(
    lat, 
    level_1, 
    level_2, 
    temp_cold_diff = diff
    )

prec_summer <-
  diff_pattern_climate[2,]$temporal_variation_climate[[1]] %>% 
  dplyr::select(
    lat, 
    level_1, 
    level_2, 
    prec_summer_diff = diff
    )

prec_winter <-
  diff_pattern_climate[3,]$temporal_variation_climate[[1]] %>% 
  dplyr::select(
    lat, 
    level_1, 
    level_2, 
    prec_winter_diff = diff
    )

diff_matrix <- 
  dplyr::inner_join(
    mpd_diff, 
    mntd_diff, 
    by = c(
      "lat", 
      "level_1", 
      "level_2"
      )
    ) %>% 
  dplyr::inner_join(
    temp_cold, 
    by = c(
      "lat", 
      "level_1", 
      "level_2"
      )
    ) %>% 
  dplyr::inner_join(
    prec_summer, 
    by = c(
      "lat", 
      "level_1", 
      "level_2"
      )
    ) %>%
  dplyr::inner_join(
    prec_winter, 
    by = c(
      "lat", 
      "level_1", 
      "level_2"
      )
    ) %>% 
  dplyr::mutate_at(
    c(
      "level_1",
      "level_2"), 
    as.numeric
    ) %>% 
  dplyr::group_by(level_1, level_2) %>% 
  dplyr::summarise_all(., mean) %>% 
  dplyr::ungroup()

#--------------------------------------------------------#
# 8. Procrustes test between the temporal variation in latitudinal pattern of ----
#  climate and that in latitudinal pattern of NRI/NTI ---- 
#--------------------------------------------------------#
# Make species by samples dataframe of NRI/NTI and each of climate variable
mpd_dif_df <-
  diff_matrix %>%
  dplyr::select(
    lat,
    level_1,
    level_2, 
    mpd_diff
    ) %>%
  dplyr::select(-lat) %>% 
  dplyr::mutate(
    level_2 = paste(
      "s_",
      row.names(.), 
      sep = ""
      )
    ) %>% 
  tidyr::spread(
    key = level_1,
    value = mpd_diff,
    fill = 0
    ) %>%
  column_to_rownames("level_2")

mntd_dif_df <-
  diff_matrix %>%
  dplyr::select(
    lat,
    level_1, 
    level_2,
    mntd_diff
    ) %>%
  dplyr::select(-lat) %>% 
  dplyr::mutate(
    level_2 = paste(
      "s_",
      row.names(.),
      sep = ""
      )
    ) %>% 
  tidyr::spread(
    key = level_1,
    value = mntd_diff,
    fill = 0) %>%
  column_to_rownames("level_2")
  

temp_cold_diff_df <-
  diff_matrix %>%
  dplyr::select(
    level_1,
    level_2, 
    temp_cold_diff
    ) %>%
  dplyr::mutate(
    level_2 = paste(
      "s_",
      row.names(.), 
      sep = ""
      )
    ) %>% 
  tidyr::spread(
    key = level_1,
    value = temp_cold_diff,
    fill = 0
    ) %>%
  column_to_rownames("level_2")
  

prec_summer_diff_df <-
  diff_matrix %>%
  dplyr::select(
    level_1,
    level_2,
    prec_summer_diff
    ) %>%
  dplyr::mutate(
    level_2 = paste(
      "s_", 
      row.names(.), 
      sep = ""
      )
    ) %>% 
  tidyr::spread(
    key = level_1,
    value = prec_summer_diff,
    fill = 0
    ) %>%
  column_to_rownames("level_2")

  
prec_winter_diff_df <-
  diff_matrix %>%
  dplyr::select(
    level_1, 
    level_2, 
    prec_winter_diff
    ) %>%
  ungroup() %>%
  dplyr::mutate(
    level_2 = paste(
      "sample_", 
      row.names(.), 
      sep = ""
      )
    ) %>% 
  tidyr::spread(
    key = level_1, 
    value = prec_winter_diff,
    fill = 0
    ) %>%
  column_to_rownames("level_2")
  
# Standardise the data
mpd_dif_std <- 
  decostand(
    mpd_dif_df, 
    method = "standardize"
    )
mntd_dif_std <- 
  decostand(
    mntd_dif_df, 
    method = "standardize"
    )
temp_cold_std <- 
  decostand(
    temp_cold_diff_df, 
    method = "standardize"
    )
prec_summer_std <- 
  decostand(
    prec_summer_diff_df, 
    method = "standardize"
    )
prec_winter_std <- 
  decostand(
    prec_winter_diff_df, 
    method = "standardize"
    )

# Perform PCA of NRI/NTI and each of climate variable
pca_mpd <- rda(mpd_dif_std)
pca_mntd <- rda(mntd_dif_std)
pca_temp_cold <- rda(temp_cold_std)
pca_prec_summer <- rda(prec_summer_std)
pca_prec_winter <- rda(prec_winter_std)

# Test of significance
set.seed(2330)
protest(
  X = pca_mpd,
  Y = pca_temp_cold,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  ) 
# Procrustes Sum of Squares (m12 squared) 0.8508, Correlation in a symmetric Procrustes rotation: 0.3863, Significance: 0.001, Number of permutations: 999

set.seed(2330)
protest(
  X = pca_mpd,
  Y = pca_prec_summer,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  )
# Procrustes Sum of Squares (m12 squared):0.6723, Correlation in a symmetric Procrustes rotation: 0.5725, Significance: 0.001, Number of permutations: 999

set.seed(2330)
protest(
  X = pca_mpd,
  Y = pca_prec_winter,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  )
# Procrustes Sum of Squares (m12 squared): 0.9922, Correlation in a symmetric Procrustes rotation: 0.08808, Significance: 0.68, Number of permutations: 999


set.seed(2330)
protest(
  X = pca_mntd,
  Y = pca_temp_cold,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  ) 
# Procrustes Sum of Squares (m12 squared): 0.7148, Correlation in a symmetric Procrustes rotation: 0.5340, Significance: 0.001, Number of permutations: 999

set.seed(2330)
protest(
  X = pca_mntd,
  Y = pca_prec_summer,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  )
#Procrustes Sum of Squares (m12 squared): 0.6121, Correlation in a symmetric Procrustes rotation: 0.6229, Significance: 0.001, Number of permutations: 999

set.seed(2330)
protest(
  X = pca_mntd,
  Y = pca_prec_winter,
  scores = "sites",
  permutations = 999,
  choises = c(1, 2)
  )
# Procrustes Sum of Squares  0.9989, Correlation in a symmetric Procrustes rotation: 0.0330, Significance: 0.982, Number of permutations: 999

procrustes_summary <- 
  tibble("Matrix_1 (PCA)" = c(
    rep("ses_MPD", 3),
    rep("ses_MNTD", 3)
    ), 
    "Matrix_2 (PCA)" = c(
      "temp_cold",
      "prec_summer",
      "prec_winter",
      "temp_cold",
      "prec_summer",
      "prec_winter"
      ),
    "Procrustes Sum of Squares (m12 squared)" = c(
      0.8508,
      0.6723,
      0.9922,
      0.7148,
      0.6121,
      0.9989
      ),
    "Correlation in a symmetric Procrustes rotation" = c(
      0.3863,
      0.5725,
      0.0880,
      0.5340,
      0.6229,
      0.0330
      ),
    "Significance" = c(
      0.001,
      0.001,
      0.68,
      0.001,
      0.001,
      0.982
      ),
    "Number of permutations" = c(
      rep(999, 6)
      )
    )

write_csv(
  procrustes_summary,
  file = "Outputs/Table/v2_121023/Procrustes_test_metrics_climate_271023.csv"
  )
