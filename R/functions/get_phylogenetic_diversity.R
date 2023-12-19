#-----------------------------------------------#
# Function for estimating phylogenetic diversity ----
#-----------------------------------------------#

get_phylogenetic_diversity <- 
  function(counts,
           dataset_id,
           type = "mpd",
           null.model = "phylogeny.pool", 
           abundance.weighted = TRUE,
           runs = 99,
           ...) {
    message(
    msg = dataset_id
    )
  
  dat <- 
    counts %>%
    as.data.frame() %>% 
    column_to_rownames("sample_id")
  
  final_tree <- 
    ape::read.tree(
      paste(
        "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
        "pruned_Ramirez_Barahona_etal_2020_raxml_processed.tre",
        sep = ""
      )
    )
  
  # Make a vector of families from the tree that are missing from the counts data
  drop_list <- 
    final_tree$tip.label[!final_tree$tip.label %in% colnames(dat)]
  
  # Remove all the families that do not occur in our sample using drop.tip() 
  #  function in the 'ape' package.
  pruned_tree <- ape::drop.tip(final_tree, drop_list) 
  
  # Arrange taxa in the same order as the pruned tree
  data_ordered <- dat[,c(pruned_tree$tip.label)]
  
  # Calculate MPD and MNTD
  # Create cophenetic distance from the pruned tree.
  phy_dist <- cophenetic(pruned_tree) 
  
  # MPD using ses.mpd() function in the 'picante' package
  
  if (type == "mpd") {
    mpd_phylogeny <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = null.model, 
        abundance.weighted = abundance.weighted,
        runs
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() 
    
    return(mpd_phylogeny) 
    
  } else if (type == "mntd") {  
    
    mntd_phylogeny <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = null.model, 
        abundance.weighted = abundance.weighted,
        runs = runs
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() 
    
    return(mntd_phylogeny) 
  } else {
    stop("type must be either 'mpd' or 'mntd'")
  }
}
