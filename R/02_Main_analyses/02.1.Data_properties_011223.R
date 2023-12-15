#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Important data properties ----

#----------------------------------------------------------#


### Question: Is this script needed? I see the data_filtered is used in 04. with the same code and some additional step, thinking this script 2 is redundant.



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
data_filtered_long <-
  data_filtered %>%
  dplyr::select(
    dataset_id,
    long,
    lat,
    phylodiversity_combined
    ) %>%
  tidyr::unnest(phylodiversity_combined) %>%
  dplyr::filter(age > 0) # why filter and not calculate to zero?



# Is the following data summary relevant anymore ?
#--------------------------------------------------------#
# 4. Extract data summary ----
#--------------------------------------------------------#


summary_1 <- 
  data_filtered_long %>% 
  dplyr::summarise(
    dplyr::across(where(is.numeric),
                  list(min = ~ min(.x, na.rm = TRUE),
                       max = ~ max(.x, na.rm = TRUE),
                       mean = ~ mean(.x, na.rm = TRUE)
                       )))
    

