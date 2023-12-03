#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Estimation of MPD and MNTD ----
                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Load configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load data ----
#-----------------------------------------------#
full_dat <- 
  readr::read_rds(
    paste(
      "Inputs/Data/",
      "data_processed_for_phylodiversity_estimation_101023.rds",
      sep = ""
      )
  )

full_dat %>% 
  dplyr::filter(
    orig_pollen_percentage == "FALSE"
    ) # 48 datasets; 145 datasets are counts
# Of the 99 final datasets, 46 datasets are percentages
# Therefore, all counts are converted to percentages for consistency

#-----------------------------------------------#
# Load full backbone tree ----
#-----------------------------------------------#
# NOTE: run this only for the first time!
# After we create a regional tree, no need to re-run this code

mytree <-
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Ramirez_Barahona_etal_2020_raxml_processed.tre",
      sep = ""
      ) 
  )

# To prune the backbone tree to keep only keep the taxa from Asian data, we need 
# the list of taxa present in the data.
req_taxa <-
  full_dat %>%
  dplyr::select(harmonised_fam_angiosperms) %>%
  dplyr::mutate(taxa = purrr::map(
    harmonised_fam_angiosperms,
    ~ colnames(.x) %>%
      enframe(name = NULL,
              value = "taxa")
    )
  ) %>%
  dplyr::select(taxa) %>%
  unnest(col = c(taxa)) %>%
  distinct() %>% 
  arrange(taxa) # 149 taxa in the datasets

req_taxa %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list <- 
 mytree$tip.label[!mytree$tip.label %in% req_taxa$taxa] # 295

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree <- ape::drop.tip(mytree, drop_list) # 148 families     

#-----------------------------------------------#
# Save the regional backbone tree ----
#-----------------------------------------------#
#ape::write.tree(
#  final_tree, 
#  paste(
#    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
#    "pruned_Ramirez_Barahona_etal_2020_raxml_processed.tre",
#    sep = ""
#    )
#  )

#-----------------------------------------------#
# Load the regional tree for further analyses ----
# No need to re-run the code above everytime!
#-----------------------------------------------#
cl <- parallel::detectCores()
cluster <- multidplyr::new_cluster(cl - 1)
multidplyr::cluster_library(cluster, "tidyverse")
multidplyr::cluster_library(cluster, "ape")
multidplyr::cluster_library(cluster, "picante")

# copy custom function to each cluster
multidplyr::cluster_copy(cluster, "get_phylogenetic_diversity") 

set.seed(1234)

start_time <- Sys.time()  

phylo_div <-
  full_dat %>%
  multidplyr::partition(cluster) %>% 
  dplyr::mutate(phylogenetic_diversity =
                  purrr::map2(
                    .x = dataset_id,
                    .y = harmonised_fam_angiosperms_percentages,
                    .f = get_phylogenetic_diversity
                    )
                ) %>% 
  collect() 

finished <- Sys.time()
total_time <- start_time - finished
total_time # -1.111461 days


#-----------------------------------------------#
# Test if the estimates are produced correctly for all the levels ----
#-----------------------------------------------#
"mpd_taxa_labels_abundance_wt" 
"mpd_taxa_labels_no_wt"
"mpd_phylogeny_pool_abundance_wt"
"mpd_phylogeny_pool_no_wt"  
"mntd_taxa_labels_abundance_wt"
"mntd_taxa_labels_no_wt"
"mntd_phylogeny_pool_abundance_wt"
"mntd_phylogeny_pool_no_wt"

test_1 <-
  phylo_div %>%
  dplyr::mutate(
    test = purrr::map_chr(
      phylogenetic_diversity,
      ~ ifelse(
        any(
          is.na(
            .x$mpd_taxa_labels_abundance_wt
            )
          ),
        "FAIL", 
        "PASS"
        )
      )
    ) %>%
  dplyr::filter(test == "FAIL") # 0

test_2 <-
  phylo_div %>%
  dplyr::mutate(
    test = purrr::map_chr(
      phylogenetic_diversity,
      ~ ifelse(
        any(
          is.na(
            .x$mpd_taxa_labels_no_wt
          )
        ),
        "FAIL", 
        "PASS"
      )
    )
  ) %>%
  dplyr::filter(test == "FAIL") # 0

test_3 <-
  phylo_div %>%
  dplyr::mutate(
    test = purrr::map_chr(
      phylogenetic_diversity,
      ~ ifelse(
        any(
          is.na(
            .x$mntd_taxa_labels_abundance_wt
          )
        ),
        "FAIL", 
        "PASS"
      )
    )
  ) %>%
  dplyr::filter(test == "FAIL") # 0

test_4 <-
  phylo_div %>%
  dplyr::mutate(
    test = purrr::map_chr(
      phylogenetic_diversity,
      ~ ifelse(
        any(
          is.na(
            .x$mntd_taxa_labels_no_wt
          )
        ),
        "FAIL", 
        "PASS"
      )
    )
  ) %>%
  dplyr::filter(test == "FAIL") # 0

#-------------------------------------------------#
# Load full data to join additional information ----
#-------------------------------------------------#
asia_dat <- 
  readr::read_rds(
    "Inputs/Data/Data_processed_2022-08-30.rds"
  ) %>% 
  dplyr::select(
    dataset_id,
    sitename,
    depositionalenvironment,
    ecozone_koppen_5,
    ecozone_koppen_15,
    source_of_data,
    data_publicity
  )

full_data <- 
  dplyr::inner_join(
    phylo_div,
    asia_dat,
    by = "dataset_id"
  ) %>% 
  dplyr::select(
    dataset_id,
    sitename,
    depositionalenvironment,
    lat,
    long,
    harmonisation_region,
    ecozone_koppen_5,
    ecozone_koppen_15,
    raw_counts,
    everything() 
  )

#-------------------------------------------------#
# Save the data of all sites ----
#-------------------------------------------------#
readr::write_rds(
  full_data, 
  file = "Inputs/Data/phylogenetic_diversity_estimated_121023.rds",
  compress = "gz"
  )

#-------------------------------------------------#
# Prepare the data ready for further analyses ----
#-------------------------------------------------#
filtered_dat <- 
  full_data %>% 
  dplyr::filter(!phylogenetic_diversity == "NA") %>%
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>%  # 99 records
  dplyr::mutate(
    phylodiversity_age_combined = purrr::map2(
      .x = phylogenetic_diversity,
      .y = levels_filtered,
      .f = ~ {
        dat <- .x
        ages <- 
          .y %>% 
          dplyr::select(
            sample_id,
            age,
            upper,
            lower
            )
        pd <- 
          dat$mpd_taxa_labels_abundance_wt[[1]] %>% 
          dplyr::select(
            sample_id,
            mpd_taxa_labels_abundance_wt = mpd
            ) %>% 
          dplyr::inner_join(
            dat$mpd_taxa_labels_no_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mpd_taxa_labels_no_wt = mpd
              ),
            by = "sample_id"
            ) %>%
          dplyr::inner_join(
            dat$mpd_phylogeny_pool_abundance_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mpd_phylogeny_pool_abundance_wt = mpd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            dat$mpd_phylogeny_pool_no_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mpd_phylogeny_pool_no_wt = mpd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            dat$mntd_taxa_labels_abundance_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mntd_taxa_labels_abundance_wt = mntd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            dat$mntd_taxa_labels_no_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mntd_taxa_labels_no_wt = mntd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            dat$mntd_phylogeny_pool_abundance_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mntd_phylogeny_pool_abundance_wt = mntd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            dat$mntd_phylogeny_pool_no_wt[[1]] %>% 
              dplyr::select(
                sample_id,
                mntd_phylogeny_pool_no_wt = mntd
              ),
            by = "sample_id"
            ) %>% 
          dplyr::inner_join(
            ages,
            by = "sample_id"
            ) %>% 
          dplyr::select(
            sample_id,
            age,
            upper,
            lower,
            everything()
          )
        }
      )
    )

readr::write_rds(
  filtered_dat,
  file = "Inputs/Data/data_for_main_analysis_121023.rds",
  compress = "gz"
  )
