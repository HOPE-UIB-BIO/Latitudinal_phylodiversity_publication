#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Holocene-wide latitudinal pattern ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
source_data <- 
  read_rds("Inputs/Data/source_data_191223.rds") 

output_gam_pd <- 
  read_rds("Outputs/Data/Overall_gam_lat_PD_201223.rds") 

#--------------------------------------------------------#
# 3. Plot GAM models ----
#--------------------------------------------------------#
# 3.1 Plot the model spatially ----
plot_gam <-
  purrr::map2(
    .x = output_gam_pd[1,]$predicted_gam,
    .y = output_gam_pd[1,]$vars,
    .f = ~ {
      dat <-
        .x %>%
        dplyr::select(-dataset_id)
      
      dat1 <-
        .x %>%
        dplyr::select(-dataset_id) %>%
        dplyr::group_by(lat) %>%
        dplyr::summarise_all(mean)
      
      ggplot2::ggplot(
        aes(
          x = lat,
          y = var,
          group = age,
          colour = age
          ),
        data = dat
        ) +
        
        ggplot2::theme_classic() +
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
            ),
          linewidth = 1.5,
          colour = color_pd_curve,
          data = dat1
          ) +
        ggplot2::labs(
          x = expression(
            paste(
            'Latitude ', (degree ~ N)
            )
            ),
          y = "Phylogenetic dispersion"
          ) +
        
        ggplot2::theme(
          axis.title = element_text(
            size = 19,
            color = color_common
            ),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = margin(
            0.1, 0.5, 1, 1, 
            "cm"
          )  #t, r, b, l
          
          ) 
      }
    )

gam_curve_mpd <- 
  plot_gam[[1]] + 
  ggplot2::labs(
    y = "Phylogenetic dispersion"
    )


# 3.2 Plot the model temporally ----
plot_gam_temporal <-
  purrr::map2(
    .x = output_gam_pd[1,]$predicted_gam,
    .y = output_gam_pd[1,]$vars,
    
    .f = ~ {
      dat <-
        .x %>%
        dplyr::select(-dataset_id)
      
      dat1 <-
        .x %>%
        dplyr::select(-dataset_id) %>%
        dplyr::group_by(age) %>%
        dplyr::summarise_all(., mean) 
        
      ggplot2::ggplot(
        aes(
          x = age,
          y = var,
          group = lat,
          colour = lat
          ),
        data = dat
        ) +
        ggplot2::geom_line(
          aes(
            x = age,
            y = var
            ),
          linewidth = 1.5,
          colour = color_pd_curve,
          data = dat1
          ) +
       
        ggplot2::theme_classic() +
        
        ggplot2::labs(
          x = "Time",
          y = "Phylogenetic dispersion"
          ) +
        
        ggplot2::theme(
          axis.title = element_text(
            size = 19,
            color = color_common
            ),
          axis.text = element_blank(),
          legend.position = "none",
          plot.margin = margin(
            0.1, 0.5, 1, 1, 
            "cm"
          )  #t, r, b, l
          ) 
      }
    )

gam_curve_mpd_temporal <- 
  plot_gam_temporal[[1]] 
gam_curve_mpd_temporal_rev_x <- 
  gam_curve_mpd_temporal + 
  ggplot2::scale_x_reverse(
    limits = c(12000, 0),
    breaks = seq(12000, 0, by = -2000)
  ) 

# Merge figures into one
final_fig_mpd_rev_x <-
  ggpubr::ggarrange(
    gam_curve_mpd,
    gam_curve_mpd_temporal_rev_x,
    ncol = 2,
    nrow = 1,
    labels = c("(A)", "(B)"),
    font.label = list(size = 20),
    hjust = -2.5,
    vjust = 13
    ) 

ggplot2::ggsave(
  final_fig_mpd_rev_x,
  filename = paste(
    "Outputs/Figure/",
    "Diagrammatic_spatiotemporal_PD_170424.tiff",
    sep = ""
    ),
  height = 9,
  width = 25,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )

