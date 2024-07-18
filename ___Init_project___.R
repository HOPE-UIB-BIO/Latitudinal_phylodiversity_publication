#---------------------------------------------------------------------------------#
#
# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia during 
#                                the Holocene
#
#                             Bhatta et al. 2024
#
#---------------------------------------------------------------------------------#

#----------------------------------------------------------#
#
#                     Project setup
#                 
#----------------------------------------------------------#

# Script to prepare all necessary components of environment to run the Project.
#   Needs to be run only once

#----------------------------------------------------------#
# Step 1: Install 'renv' package -----
#----------------------------------------------------------#

utils::install.packages("renv")

#----------------------------------------------------------#
# Step 2: Deactivate 'renv' package -----
#----------------------------------------------------------#

# deactivate to make sure that packages are updated on the machine
renv::deactivate()

#----------------------------------------------------------#
# Step 3: Create a list of packages
#----------------------------------------------------------#

package_list <- 
  c(
    "ape",
    "assertthat",
    "car",
    "corrplot", 
    "cowplot",
    "devtools",
    "effects", 
    "flextable",
    "ggpubr",
    "ggmap",
    "ggspatial",
    "ggtext",
    "glmm",
    "glmmTMB",
    "gratia",
    "grid",
    "gridExtra",
    "here",   
    "leaflet",
    "lemon",
    "lme4", 
    "maps",
    "mgcv", 
    "multidplyr",
    "MuMIn", 
    "parallel",
    "performance",
    "picante",
    "raster",
    "RColorBrewer",
    "remotes",
    "renv", 
    "rpart",
    "roxygen2", 
    "sf",
    "sjPlot", 
    "sjmisc",
    "sjstats", 
    "targets",
    "tidymv", 
    "tidyverse",  
    "usethis",
    "vegan",
    "viridis",
    "xaringan" 
  )


#----------------------------------------------------------#
# Step 4: Install packages to the machine
#----------------------------------------------------------#

sapply(package_list, utils::install.packages, character.only = TRUE)
remotes::install_github("HOPE-UIB-BIO/R-Fossilpol-package")
remotes::install_github("HOPE-UIB-BIO/R-Ecopol-package")

#----------------------------------------------------------#
# Step 5: Activate 'renv' project
#----------------------------------------------------------#

renv::activate()

#----------------------------------------------------------#
# Step 6: Install packages to the project
#----------------------------------------------------------#
package_list <- 
  c(
    "ape",
    "assertthat",
    "car",
    "corrplot", 
    "cowplot",
    "devtools",
    "effects", 
    "flextable",
    "ggpubr",
    "ggmap",
    "ggspatial",
    "ggtext",
    "glmm",
    "glmmTMB",
    "gratia",
    "grid",
    "gridExtra",
    "here",   
    "leaflet",
    "lemon",
    "lme4", 
    "maps",
    "mgcv", 
    "multidplyr",
    "MuMIn", 
    "parallel",
    "performance",
    "picante",
    "raster",
    "RColorBrewer",
    "remotes",
    "renv", 
    "rpart",
    "roxygen2", 
    "sf",
    "sjPlot", 
    "sjmisc",
    "sjstats", 
    "targets",
    "tidymv", 
    "tidyverse",  
    "usethis",
    "vegan",
    "viridis",
    "xaringan" 
  )


sapply(package_list, utils::install.packages, character.only = TRUE)
remotes::install_github("HOPE-UIB-BIO/R-Fossilpol-package")
remotes::install_github("HOPE-UIB-BIO/R-Ecopol-package")

#----------------------------------------------------------#
# Step 7: Synchronize package versions with the project 
#----------------------------------------------------------#
renv::snapshot(lockfile = "renv/library_list.lock")

library(here)
renv::restore(lockfile = here::here( "renv/library_list.lock"))
renv::snapshot()
#----------------------------------------------------------#
# Step 8: Run the project 
#----------------------------------------------------------#
