#---------------------------------------------------------------------------------#
#
# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia during 
#                                the Holocene
#
#                             Bhatta et al. 2024
#
#---------------------------------------------------------------------------------#

#----------------------------------------------------------#

# Prepare data for estimation of MPD and MNTD ----
#                          
#----------------------------------------------------------#

#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load data the filtered data for Asia ----
#-----------------------------------------------#
# Filtered data for Asia [207 datasets], which was used in Asian synthesis
asia_dat <-
  read_rds("Inputs/Data/Data_processed_2022-08-30.rds") %>%
  dplyr::rename(
    harmonised_counts = pollen_counts
    ) %>% 
  dplyr::filter(
    harmonisation_region == "Asia_Main" |
    harmonisation_region == "Asia_Siberia"
    ) # 193 datasets

#-----------------------------------------------#
# Harmonise taxa at family level ----
#-----------------------------------------------#
dat_harmonised <-
  asia_dat %>%
  dplyr::mutate(
    harmonised_fam_angiosperms =
                  purrr::pmap(
                    .l = list(
                      dataset_id, # ..1
                      harmonisation_region, # ..2
                      raw_counts # ..3
                      ),
                    .f = ~ {
                      harm_table <-
                        readr::read_csv(
                          paste(
                            "Inputs/Tables/Phylo_div_harmonisation/", 
                            ..2,
                            "_280923",
                            ".csv", 
                            sep = ""
                            )
                          )
                      message(
                        msg = paste0(
                          ..1,
                          ": ",
                          ..2,
                          ", ",
                          paste(
                            "Inputs/Tables/Phylo_div_harmonisation/", 
                            ..2,
                            "_280923",
                            ".csv",
                            sep = ""
                            ),
                          "\n"
                          )
                        )
                      # transform all taxa to numeric
                      col_names <- names(..3)
                      taxa_names <- col_names[!col_names %in% "sample_id"]
                      
                      suppressWarnings(dat <-
                                         ..3 %>%
                                         dplyr::mutate_at(
                                           taxa_names,
                                           as.numeric
                                           )
                                       )
                      res <- 
                        dat %>%
                        as.data.frame() %>%
                        dplyr::mutate(
                          dplyr::across(
                            where(is.numeric),
                            ~ tidyr::replace_na(., 0)
                            )
                          ) %>%
                        tidyr::gather(
                          key = taxon_name, 
                          value = counts,
                          -sample_id
                          ) %>%
                        dplyr::left_join(
                          harm_table, 
                          by = "taxon_name"
                          ) %>%
                        dplyr::filter(!level_3_phylo_div == "delete") %>%
                        dplyr::filter(!gymnosperm == "TRUE") %>% 
                        dplyr::group_by(
                          sample_id,
                          level_3_phylo_div
                          ) %>%
                        dplyr::summarise(
                          .groups = "keep",
                          counts = sum(counts)
                          ) %>%
                        dplyr::rename(harmonised_counts = `level_3_phylo_div`) %>%
                        tidyr::spread(
                          harmonised_counts, 
                          counts
                          ) %>%
                        dplyr::ungroup() %>%
                        tibble::column_to_rownames("sample_id") %>%
                        dplyr::select_if(colSums(.) != 0) %>%
                        tibble::rownames_to_column("sample_id") %>%
                        tibble::as_tibble()
                      return(res)
                      }
                    )
    )

#-----------------------------------------------#
# Prepare the data for phylogenetic diversity analysis ----
#-----------------------------------------------#
# Keep only the samples in the 'harmonised_fam_angiosperms' that are present in
# 'harmonised_counts' because harmonised counts are filtered based on HOPE 
# criteria

dat_harmonised_filtered <- 
  dat_harmonised %>% 
  dplyr::mutate(
    harmonised_fam_angiosperms = 
      purrr::map2(
        .x = harmonised_counts,
        .y = harmonised_fam_angiosperms,
        .f = ~ .x %>% 
        dplyr::select(sample_id) %>% 
          inner_join(
            .y,
            by = "sample_id"
            )
        ),
    n_sample = purrr::map_dbl(
      harmonised_fam_angiosperms,
      ~ nrow(.x)
      ),
    levels_filtered =
      purrr::map(
        levels, 
        ~ .x %>%
          dplyr::filter(age < 12000)
        ),
    # Filter out the levels (rows) with < 3 taxa (otherwise such levels would 
    #  produce NA/NaN during phylogenetic diversity analyses).
    harmonised_fam_angiosperms = 
      purrr::map(
        harmonised_fam_angiosperms,
        ~ .x[rowSums(.x > 0) > 2,] %>%
    # Remove columns (taxa) with zero colsums      
          tibble::column_to_rownames("sample_id") %>%
          dplyr::select_if(colSums(.) != 0) %>%
          tibble::rownames_to_column("sample_id") %>%
          tibble::as_tibble()
        ),
    ) %>% 
  dplyr::filter(n_sample >= min_n_levels) %>%  #[config_criteria]
  dplyr::select(-n_sample)

suppressWarnings(
  dat_harmonised_filtered %>%
  dplyr::mutate(test = purrr::map_lgl(
    harmonised_fam_angiosperms,
    ~ ifelse(
      is.null(
        any(
          colSums(.x %>%
                    column_to_rownames("sample_id")
                  )
              )
         ), TRUE, FALSE
       )
     )
   ) 
  ) %>%
  dplyr::select(test) %>%
  dplyr::filter(test == "TRUE") # 0, OK!

#-----------------------------------------------#
# Convert harmonised families of counts data to percentages ----
#-----------------------------------------------#
dat_harmonised_filtered_1 <- 
  dat_harmonised_filtered %>% 
  dplyr::mutate(
    harmonised_fam_angiosperms_percentages = 
      # Convert counts to percentages (no need to convert pre-existing percentages)    
      ifelse(orig_pollen_percentage == TRUE, 
             purrr::map(
               .x = harmonised_fam_angiosperms, 
               ~ .x %>% 
                 column_to_rownames("sample_id") %>% 
                 round(., digits = 3) %>% 
                 rownames_to_column("sample_id")
               ),
             
             purrr::map(
               .x = harmonised_fam_angiosperms, 
               .f = ~ {
                 counts <- 
                   .x %>% 
                   column_to_rownames("sample_id")
                 percentages <-
                 ((counts/rowSums(counts)) * 100) %>%
                 round(., digits = 3) %>% 
                   rownames_to_column("sample_id")
                   
                 return(percentages)
                 }
               )
             )
    ) %>% 
  dplyr::select(
    dataset_id,
    sitename,
    long,
    lat,
    depositionalenvironment,
    harmonisation_region,
    ecozone_koppen_5,
    ecozone_koppen_15,
    levels,
    levels_filtered,
    harmonised_fam_angiosperms,
    harmonised_fam_angiosperms_percentages,
    orig_pollen_percentage,
    source_of_data
  )

write_rds(
  dat_harmonised_filtered_1, 
  file = paste(
    "Inputs/Data/",
    "data_processed_for_",
    "phylodiversity_estimation_191223.rds",
    sep = ""
    ),
  compress = "gz"
  )
