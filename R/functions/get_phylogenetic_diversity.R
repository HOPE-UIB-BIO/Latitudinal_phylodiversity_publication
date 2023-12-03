#-----------------------------------------------#
# Function for estimating phylogenetic diversity ----
#-----------------------------------------------#

get_phylogenetic_diversity <- 
  
  function(dataset_id, counts) {
    
    message(
      msg = dataset_id
    )
    dat <- 
      counts %>%
      as.data.frame() %>% 
      column_to_rownames("sample_id")
    
    tree <- 
      ape::read.tree(
        paste(
          "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
          "pruned_Ramirez_Barahona_etal_2020_raxml_processed.tre",
          sep = ""
          )
        )
    
    # Make a vector of families from the tree that are missing from the counts data
    drop_list <- 
      tree$tip.label[!tree$tip.label %in% colnames(dat)]
    
    # Remove all the families that do not occur in our sample using drop.tip() 
    #  function in the 'ape' package.
    pruned_tree <- ape::drop.tip(tree, drop_list) 
    
    # Arrange taxa in the same order as the pruned tree
    data_ordered <- dat[,c(pruned_tree$tip.label)]
    
    # Calculate MPD and MNTD
    # Create cophenetic distance from the pruned tree.
    phy_dist <- cophenetic(pruned_tree) 
    
    # MPD using ses.mpd() function in the 'picante' package
    mpd_taxa_labels_abundance_wt <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = "taxa.labels", 
        abundance.weighted = TRUE,
        runs = 9999
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mpd = mpd.obs.z) 
    
    mpd_taxa_labels_no_wt <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = "taxa.labels",
        abundance.weighted = FALSE,
        runs = 9999
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mpd = mpd.obs.z) 
    
    mpd_phylogeny_pool_abundance_wt <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = "phylogeny.pool", 
        abundance.weighted = TRUE,
        runs = 9999
        ) %>% 
    # phylogeny.pool:	Randomizes community data matrix by drawing species from 
    # pool of phylogeny with the species occurring in at least one community 
    # (sample pool) with equal probability
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mpd = mpd.obs.z) 
    
    mpd_phylogeny_pool_no_wt <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = "phylogeny.pool", 
        abundance.weighted = FALSE,
        runs = 9999
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mpd = mpd.obs.z) 
    
    # The output contains 8 columns. The sixth column 'mpd.obs.z' is the 
    #  Standardized effect size of MPD vs. null communities (equivalent to -NRI).
    
    # "abundance.weighted": Should mean nearest taxon distances for each species 
    #  be weighted by species abundance? (default = FALSE).
    
    
    # MNTD/NTI using ses.mntd() function in the 'picante' package
    mntd_taxa_labels_abundance_wt <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = "taxa.labels", 
        abundance.weighted = TRUE,
        runs = 9999
        ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mntd = mntd.obs.z) 
    # The output contains 8 columns. The sixth column 'mntd.obs.z' is the 
    #  Standardized effect size of MNTD vs. null communities (equivalent to -NTI).
    
    mntd_taxa_labels_no_wt <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = "taxa.labels", 
        abundance.weighted = FALSE,
        runs = 9999
      ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mntd = mntd.obs.z) 
    
    mntd_phylogeny_pool_abundance_wt <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = "phylogeny.pool", 
        abundance.weighted = TRUE,
        runs = 9999
      ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mntd = mntd.obs.z) 
    
    mntd_phylogeny_pool_no_wt <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = "phylogeny.pool", 
        abundance.weighted = FALSE,
        runs = 9999
      ) %>% 
      rownames_to_column("sample_id") %>% 
      as_tibble() %>% 
      dplyr::rename(mntd = mntd.obs.z) 
    
    result <- 
      tibble(
        mpd_taxa_labels_abundance_wt = list(mpd_taxa_labels_abundance_wt),
        mpd_taxa_labels_no_wt = list(mpd_taxa_labels_no_wt),
        mpd_phylogeny_pool_abundance_wt = list(mpd_phylogeny_pool_abundance_wt),
        mpd_phylogeny_pool_no_wt = list(mpd_phylogeny_pool_no_wt),
        mntd_taxa_labels_abundance_wt = list(mntd_taxa_labels_abundance_wt),
        mntd_taxa_labels_no_wt = list(mntd_taxa_labels_no_wt),
        mntd_phylogeny_pool_abundance_wt = list(mntd_phylogeny_pool_abundance_wt),
        mntd_phylogeny_pool_no_wt = list(mntd_phylogeny_pool_no_wt)
      )
     
    return(result) 
  }
