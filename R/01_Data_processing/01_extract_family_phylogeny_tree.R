#------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Process family level phylogeny of angiosperms from "Ramírez-Barahona et al. 
# 2020. The delayed and geographically heterogeneous diversification of 
# flowering plant families. Nature Ecology and Evolution. 4: 1232-1238"

#------------------------------------------#
library(tidyverse)
library(ape)

# We used maximum likelihood angiosperm phylogeny (NEWICK format) estimated with 
# RAxML: Data1_RAxMLTree.tre by Ramírez-Barahona et al. 2020 ----

tree_2 <-
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Data1_RAxMLTree.tre",
      sep = ""
    )
  ) # 1216 tip labels, each label corresponds to authors' samples

tips <- tree_2$tip.label 

families_raw <- 
  tips %>% 
  strsplit(., "_") %>% 
  purrr::map(
    .,
    ~pluck(.x[1])
  ) %>% 
  unlist()

pruned_tree_2 <- 
  ape::drop.tip(
    tree_2, 
    tip = tips[
      which(
        duplicated(
          families_raw)
      )
    ]
  )

pruned_tree_2_tip_label <- 
  pruned_tree_2$tip.label %>% 
  strsplit(., "_") %>% 
  purrr::map(
    .,
    ~pluck(.x[1])
  ) %>% 
  unlist()

pruned_tree_2$tip.label <- pruned_tree_2_tip_label
ape::is.rooted(pruned_tree_2)       # TRUE
ape::is.ultrametric(pruned_tree_2)  # FALSE

#write.tree(
#  pruned_tree_2,
#  file = paste(
#    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
#    "Ramirez_Barahona_etal_2020_raxml_processed.tre",
#    sep = ""
#  )
#)
