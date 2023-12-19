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
dat_harmonised_filtered_1 <- 
  readr::read_rds(
    paste(
      "Inputs/Data/",
      "data_processed_for_",
      "phylodiversity_estimation_191223.rds",
      sep = ""
    )
  )

dat_harmonised_filtered_1 %>% 
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

pruned_tree_2 <-
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
  dat_harmonised_filtered_1 %>%
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
  dplyr::filter(!taxa %in% pruned_tree_2$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list <- 
  pruned_tree_2$tip.label[!pruned_tree_2$tip.label %in% req_taxa$taxa] # 295

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree <- ape::drop.tip(pruned_tree_2, drop_list) # 148 families     

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
# No need to re-run the code above every time!
#-----------------------------------------------#
set.seed(1234)
phylo_div <-
  dat_harmonised_filtered_1 %>%
  dplyr::mutate(
    phylogenetic_diversity_mpd =
      purrr::map2(
        .x = harmonised_fam_angiosperms_percentages,
        .y = dataset_id,
        .f = ~ get_phylogenetic_diversity(
          .x,
          .y,
          type = "mpd",
          null.model = "phylogeny.pool", 
          abundance.weighted = TRUE,
          runs = 9999
          )
        )
    ) %>%
  dplyr::mutate(
    phylogenetic_diversity_mntd =
      purrr::map2(
        .x = harmonised_fam_angiosperms_percentages,
        .y = dataset_id,
        .f = ~ get_phylogenetic_diversity(
          .x,
          .y,
          type = "mntd",
          null.model = "phylogeny.pool", 
          abundance.weighted = TRUE,
          runs = 9999
          )
        )
    ) 

#-------------------------------------------------#
# Save the data of all sites ----
#-------------------------------------------------#
readr::write_rds(
  phylo_div, 
  file = "Inputs/Data/phylogenetic_diversity_estimated_191223.rds",
  compress = "gz"
  )

#-------------------------------------------------#
# Prepare the data for further analyses ----
#-------------------------------------------------#
data_filtered_phylodiversity <- 
  phylo_div %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>%  # 99 records
  dplyr::mutate(
    phylodiversity_combined = purrr::pmap(
      list(
        levels_filtered,
        phylogenetic_diversity_mpd,
        phylogenetic_diversity_mntd
        ),
      .f = ~ {
        ..1 %>%
          dplyr::select(
            sample_id, 
            age, 
            upper, 
            lower
            ) %>%
          dplyr::left_join(..2 %>% 
                             dplyr::select(
                               sample_id, 
                               mpd.obs.z
                               ),
                           by = "sample_id") %>%
          dplyr::left_join(..3 %>% 
                             dplyr::select(
                               sample_id, 
                               mntd.obs.z
                               ),
                           by = "sample_id") %>% 
          dplyr::rename(
            mpd = mpd.obs.z,
            mntd = mntd.obs.z
            ) %>%
          return()
      }
    )
  ) 
        
#-------------------------------------------------#
# Save the filtered data  ----
#-------------------------------------------------#
readr::write_rds(
  data_filtered_phylodiversity,
  file = "Inputs/Data/data_for_main_analysis_191223.rds",
  compress = "gz"
  )
