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
# cl <- parallel::detectCores()
# cluster <- multidplyr::new_cluster(cl - 1)
# multidplyr::cluster_library(cluster, "tidyverse")
# multidplyr::cluster_library(cluster, "ape")
# multidplyr::cluster_library(cluster, "picante")
# 
# # copy custom function to each cluster
# multidplyr::cluster_copy(cluster, "get_phylogenetic_diversity") 
# 
# set.seed(1234)
# 
# start_time <- Sys.time()  

phylo_div <-
  full_dat %>%
 # multidplyr::partition(cluster) %>% 
  dplyr::mutate(phylogenetic_diversity_mpd =
                  purrr::map2(
                    .x = harmonised_fam_angiosperms_percentages,
                    .y = dataset_id,
                    .f = ~ get_phylogenetic_diversity(
                      .x, 
                      .y,
                      type = "mpd",
                      null.model = "phylogeny.pool", 
                      abundance.weighted = TRUE,
                      runs = 99)
                    )
  ) %>%
  dplyr::mutate(phylogenetic_diversity_mntd =
                  purrr::map2(
                    .x = harmonised_fam_angiosperms_percentages,
                    .y = dataset_id,
                    .f = ~ get_phylogenetic_diversity(
                      .x, 
                      .y,
                      type = "mntd",
                      null.model = "phylogeny.pool", 
                      abundance.weighted = TRUE,
                      runs = 99)
                    )
  ) 
                
#%>% 
#  collect() 

#finished <- Sys.time()
#total_time <- start_time - finished
#total_time # -1.111461 days



#-------------------------------------------------#
# Save the data of all sites ----
#-------------------------------------------------#
readr::write_rds(
  phylo_div, 
  file = "Inputs/Data/phylogenetic_diversity_estimated_121023.rds",
  compress = "gz"
  )

#-------------------------------------------------#
# Prepare the data ready for further analyses ----
#-------------------------------------------------#
get_combined_data <- function(.x, 
                              .y, 
                              select_var = NULL
                              ){
  dat <- .x %>%
    dplyr::select(sample_id, age,upper, lower) %>%
    dplyr::left_join(
      .y %>% dplyr::select(all_of(select_var)) ,
      by = "sample_id"
    ) 
  
  return(dat)
}

data_filtered <- 
  phylo_div %>% 
  dplyr::filter(!long < 75 & !long > 125) %>%
  dplyr::filter(!lat < 25 & !lat > 66) %>%  # 99 records
  dplyr::mutate(
    phylodiversity_combined = purrr::pmap(
      list(levels_filtered,
           phylogenetic_diversity_mpd,
           phylogenetic_diversity_mntd),
      .f = ~ {
        ..1 %>%
          dplyr::select(sample_id, age,upper, lower) %>%
          dplyr::left_join(..2 %>% 
                             dplyr::select(sample_id, mpd.obs.z),
                           by = "sample_id") %>%
          dplyr::left_join(..3 %>% 
                             dplyr::select(sample_id, mntd.obs.z),
                           by = "sample_id") %>% 
          rename(mpd = mpd.obs.z,
                 mntd = mntd.obs.z) %>%
          return()
      })
    ) 
        
       

readr::write_rds(
  data_filtered,
  file = "Inputs/Data/data_for_main_analysis_121023.rds",
  compress = "gz"
  )
