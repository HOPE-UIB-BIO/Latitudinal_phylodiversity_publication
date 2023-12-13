#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Important data properties ----

#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
data_filtered <- 
  read_rds("Inputs/Data/data_for_main_analysis_121023.rds") 

#--------------------------------------------------------#
# 3. Filter the data to avoid marginal datasets ----
#--------------------------------------------------------#
data_filtered_1 <-
  data_filtered %>%
  dplyr::select(
    dataset_id,
    long,
    lat,
    phylodiversity_age_combined
    ) %>%
  tidyr::unnest(phylodiversity_age_combined ) %>%
  dplyr::filter(age > 0) %>% 
  dplyr::select(
    dataset_id,
    lat,
    long,
    age,
    mpd = mpd_phylogeny_pool_abundance_wt,
    mntd = mntd_phylogeny_pool_abundance_wt
    )
  

#--------------------------------------------------------#
# 4. Extract data summary ----
#--------------------------------------------------------#
summary_1 <- 
  data_filtered_1 %>% 
  dplyr::summarise(
    lat_min = min(lat),
    lat_max = max(lat),
    age_min = min(age),
    age_max = max(age),
    mpd_min = min(mpd),
    mpd_max = max(mpd),
    mpd_mean = mean(mpd),
    mntd_min = min(mntd),
    mntd_max = max(mntd),
    mntd_mean = mean(mntd)
    ) 

summary_2 <- 
  data_filtered %>% 
  dplyr::select(
    raw_counts,
    harmonised_fam_angiosperms
    ) %>% 
  dplyr::mutate(
    raw_samples = purrr::map_dbl(
      raw_counts,
      ~ nrow(.x)
      ),
    filtered_samples = purrr::map_dbl(
      harmonised_fam_angiosperms,
      ~ nrow(.x)
    ),
    raw_taxa = purrr::map_dbl(
      raw_counts,
      ~ .x %>% 
        dplyr::select(-sample_id) %>% 
        ncol(.)
      ),
    harmonised_taxa = purrr::map_dbl(
      harmonised_fam_angiosperms ,
      ~ .x %>% 
        dplyr::select(-sample_id) %>% 
        ncol(.)
      ),
    
    raw_taxa_per_sample = purrr::map(
      raw_counts,
      ~ .x %>% 
        dplyr::mutate(
          dplyr::across(
            where(is.numeric),
            ~ tidyr::replace_na(., 0)
          )
        ) %>% 
        dplyr::group_by(sample_id) %>% 
        tidyr::nest(data = -group_cols()) %>% 
        ungroup() %>% 
        dplyr::mutate(ntaxa_raw_per_sample = 
                        purrr::map_dbl(
                          data,
                          ~ .x %>% 
                            dplyr::select_if(colSums(.) != 0) %>% 
                            ncol(.)
                          ) 
                      ) %>%
        dplyr::select(ntaxa_raw_per_sample)
      ),
    harmonised_taxa_per_sample = purrr::map(
      harmonised_fam_angiosperms ,
      ~ .x %>% 
        dplyr::mutate(dplyr::across(
          where(is.numeric),
          ~ tidyr::replace_na(., 0)
          )
          ) %>%
        dplyr::group_by(sample_id) %>%
        tidyr::nest(data = -group_cols()) %>%
        ungroup() %>%
        dplyr::mutate(ntaxa_harmonised_per_sample = 
                        purrr::map_dbl(
                          data,
                          ~ .x %>%
                            dplyr::select_if(
                              colSums(.) != 0
                              ) %>%
                            ncol(.)
                          ) 
                      ) %>% 
        dplyr::select(ntaxa_harmonised_per_sample)
      )
    )



min_raw_samples_per_record <- min(summary_2$raw_samples) #10
max_raw_samples_per_record <- max(summary_2$raw_samples) #623
mean_raw_samples_per_record  <- mean(summary_2$raw_samples) #71

min_filtered_samples_per_record  <- min(summary_2$filtered_samples) #9
max_filtered_samples_per_record  <- max(summary_2$filtered_samples) #611
mean_filtered_samples_per_record  <- mean(summary_2$filtered_samples) #67


  
min_raw_taxa_per_record <- min(summary_2$raw_taxa) #7
max_raw_taxa_per_record <- max(summary_2$raw_taxa) #296
mean_raw_taxa_per_record  <- mean(summary_2$raw_taxa) #53

min_harmonised_taxa_per_record  <- min(summary_2$harmonised_taxa) #4
max_harmonised_taxa_per_record  <- max(summary_2$harmonised_taxa) #73
mean_harmonised_taxa_per_record  <- mean(summary_2$harmonised_taxa) #20

per_sample_taxa_raw <- 
  summary_2 %>% 
  dplyr::select(raw_taxa_per_sample) %>% 
  tidyr::unnest(cols = c(raw_taxa_per_sample))
  
per_sample_taxa_harmonised <- 
  summary_2 %>% 
  dplyr::select(harmonised_taxa_per_sample) %>% 
  tidyr::unnest(cols = c(harmonised_taxa_per_sample))

min_raw_taxa_per_sample <- min(per_sample_taxa_raw$ntaxa_raw_per_sample) #1
max_raw_taxa_per_sample <- max(per_sample_taxa_raw$ntaxa_raw_per_sample) #128
mean_raw_taxa_per_sample <- mean(per_sample_taxa_raw$ntaxa_raw_per_sample) #19


min_harmonised_taxa_per_sample <- 
  min(per_sample_taxa_harmonised$ntaxa_harmonised_per_sample) #2
max_harmonised_taxa_per_sample <- 
  max(per_sample_taxa_harmonised$ntaxa_harmonised_per_sample) #41
mean_harmonised_taxa_per_sample <- 
  mean(per_sample_taxa_harmonised$ntaxa_harmonised_per_sample)#11