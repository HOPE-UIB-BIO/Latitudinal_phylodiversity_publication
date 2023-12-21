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
# 2. Load the source data ----
#--------------------------------------------------------#
source_data <- 
  read_rds("Inputs/Data/source_data_191223.rds") 

#--------------------------------------------------------#
# 3. Fit the GAM models ----
#--------------------------------------------------------#
data_gam <-
  source_data %>%
  tidyr::pivot_longer(
    mpd:mntd, 
    names_to = "vars",
    values_to = "estimate"
    ) %>%
  dplyr::select(
    dataset_id,
    lat, 
    age, 
    age_uncertainty_index,
    vars, 
    estimate
    ) %>%
  tidyr::nest(data = -vars)  

output_gam_pd <-
  data_gam %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
        data <- .x
        set.seed(2330)
        mod <-
          mgcv::gam(
            estimate ~
              s(lat, 
                k = 10, 
                bs = 'tp'
                ) +
              s(age, 
                k = 10, 
                bs = 'tp'
                ) +
              s(age, 
               by = dataset_id, 
               bs = 'tp',      
               m = 1
               ) +
              s(dataset_id, 
                k = 99, 
                bs = 're'
                ) +
            ti(lat, age, 
               bs = c('tp', 'tp')
               ),
            weights = age_uncertainty_index,
            data = data,
            method = "REML",
            family = "gaussian",
            control = gam.control(
              trace = TRUE, 
              maxit = 200
              )
          )
      }
    ),
    
    predicted_gam =
      purrr::map2(
        .x = data, 
        .y = gam_model, 
        .f = ~ {
        data <- .x
        new_data_gam <-
          with(
            data,
            base::expand.grid(
              lat = seq(25, 66, by = 1),
              dataset_id = dataset_id[1],
              age = seq(200, 12000, by = 200)
              )
            )
        not_inlude <-
          gratia::smooths(.y) %>%
          str_subset(., 
                     c("dataset_id")
                     )
        
        crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
        
        set.seed(2330)
        predicted_mod <-
          new_data_gam %>%
          dplyr::bind_cols(
            predict(
              .y,
              newdata = new_data_gam,
              type = "response",
              se.fit = TRUE,
              exclude = not_inlude
              )
            ) %>%
          
          dplyr::mutate(
            var = fit,
            lwr = fit - (crit * se.fit),
            upr = fit + (crit * se.fit)
            ) %>%
          dplyr::select(
            !dplyr::any_of(
              c("fit", "se.fit")
              )
            )
        
        return(predicted_mod)
        }
      )
    )

#--------------------------------------------------------#
# Save GAM models ----
#--------------------------------------------------------#
write_rds(output_gam_pd,
          file = "Outputs/Data/Overall_gam_lat_PD_201223.rds",
          compress = "gz")


#--------------------------------------------------------#
# Extract and save summary of GAM models  ----
#--------------------------------------------------------#

output_gam_pd <- 
  read_rds("Outputs/Data/Overall_gam_lat_PD_201223.rds") 

gam_summary_mpd <-
  flextable::as_flextable(output_gam_pd[1, ]$gam_model[[1]])

save_as_docx(
  gam_summary_mpd,
  path = "Outputs/Table/overall_gam_mpd_201223.docx"
  )

gam_summary_mntd <-
  flextable::as_flextable(output_gam_pd[2, ]$gam_model[[1]])

save_as_docx(
  gam_summary_mntd,
  path = "Outputs/Table/overall_gam_mntd_201223.docx"
  )

#--------------------------------------------------------#
# 5. Plot GAM models ----
#--------------------------------------------------------#
# 5.1 Plot the model spatially ----
plot_gam <-
  purrr::map2(
    .x = output_gam_pd$predicted_gam,
    .y = output_gam_pd$vars,
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
        ggplot2::geom_line(
          linewidth = 0.3,
          alpha = 1
          ) +
        ggplot2::scale_color_gradient(
          high = color_high_age,
          low = color_low_age
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
          y = paste(.y)
          ) +
        
        ggplot2::theme(
          axis.title = element_text(
            size = 15,
            color = color_common
            ),
          axis.text = element_text(
            size = 13,
            color = color_common
            ),
          legend.position = c(0.38, 0.3),
          legend.background = element_rect(
            fill = "transparent"
            ),
          legend.title = element_text(
            size = 12,
            color = color_common
            ),
          legend.text = element_text(
            size = 11,
            color = color_common
            ),
          legend.spacing.y = unit(0.35, 'cm'),
          plot.margin = margin(
            1, 0, 0, 0.25, 
            "cm"
            )  #t, r, b, l
          ) +
        ggplot2::guides(
          colour = guide_colorbar(
          barheight = 3,
          barwidth = 0.75
          )
        )
      }
    )

gam_curve_mpd <- 
  plot_gam[[1]] + 
  ggplot2::labs(
    y = "ses_MPD",
    colour = "Time (cal yr BP)"
    )

gam_curve_mntd <- 
  plot_gam[[2]] + 
  ggplot2::labs(
    y = "ses_MNTD",
    colour = "Time (cal yr BP)"
    )


# 5.2 Plot the model temporally ----
plot_gam_temporal <-
  purrr::map2(
    .x = output_gam_pd$predicted_gam,
    .y = output_gam_pd$vars,
    
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
        ggplot2::geom_line(
          linewidth = 0.3,
          alpha = 1
          ) +
        ggplot2::scale_color_gradient(
          high = color_high_lat,
          low = color_low_lat
          ) +
        ggplot2::theme_classic() +
        
        ggplot2::labs(
          x = "Time (cal yr BP)",
          y = paste(.y)
          ) +
        ggplot2::scale_x_continuous(
          limits = c(0, 12000),
          breaks = seq(0, 12000, 2000)
          ) +
        ggplot2::theme(
          axis.title = element_text(
            size = 15,
            color = color_common
            ),
          axis.text = element_text(
            size = 13,
            color = color_common
            ),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
            ),
          legend.position = c(0.8, 0.8),
          legend.background = element_rect(
            fill = "transparent"
            ),
          legend.title = element_text(
            size = 12,
            color = color_common
            ),
          legend.text = element_text(
            size = 11,
            color = color_common
            ),
          plot.margin = margin(
            1, 0, 0, 0.25, 
            "cm"
            )  #t, r, b, l
          ) +
        ggplot2::guides(
          colour = guide_colorbar(
            barheight = 5,
            barwidth = 0.75
            )
        )
      }
    )

gam_curve_mpd_temporal <- 
  plot_gam_temporal[[1]] + 
  ggplot2::labs(
    y = "ses_MPD",
    colour = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      )
    ) 
gam_curve_mpd_temporal_rev_x <- 
  gam_curve_mpd_temporal + 
  ggplot2::scale_x_reverse(
    limits = c(12000, 0),
    breaks = seq(12000, 0, by = -2000)
  ) 

gam_curve_mntd_temporal <- 
  plot_gam_temporal[[2]] + 
  ggplot2::labs(
    y = "ses_MNTD",
    colour = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      )
    ) 

gam_curve_mntd_temporal_rev_x <- 
  gam_curve_mntd_temporal + 
  ggplot2::scale_x_reverse(
    limits = c(12000, 0),
    breaks = seq(12000, 0, by = -2000)
  ) 

# Merge all the figures into one
final_fig_mpd <-
  ggpubr::ggarrange(
    gam_curve_mpd,
    gam_curve_mpd_temporal,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    hjust = -0.4,
    vjust = 2.5
    )
final_fig_mpd_rev_x <-
  ggpubr::ggarrange(
    gam_curve_mpd,
    gam_curve_mpd_temporal_rev_x,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    hjust = -0.4,
    vjust = 2.5
    )


final_fig_mntd <-
 ggpubr::ggarrange(
    gam_curve_mntd,
    gam_curve_mntd_temporal,
    ncol = 1,
    nrow = 2,
    labels = c("(c)", "(d)"),
    hjust = -0.4,
    vjust = 2.5
    )

final_fig_mntd_rev_x <-
  ggpubr::ggarrange(
    gam_curve_mntd,
    gam_curve_mntd_temporal_rev_x,
    ncol = 1,
    nrow = 2,
    labels = c("(c)", "(d)"),
    hjust = -0.4,
    vjust = 2.5
    )

final_composite <- 
  ggpubr::ggarrange(
    final_fig_mpd,
    final_fig_mntd,
    ncol = 2,
    nrow = 1
    )

final_composite_rev_x <- 
  ggpubr::ggarrange(
    final_fig_mpd_rev_x,
    final_fig_mntd_rev_x,
    ncol = 2,
    nrow = 1
    )

ggplot2::ggsave(
  final_composite,
  filename = paste(
    "Outputs/Figure/",
    "Holocene_wide_spatiotemporal_pattern_201223.tiff",
    sep = ""
    ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )

ggplot2::ggsave(
  final_composite_rev_x,
  filename = paste(
    "Outputs/Figure/",
    "Holocene_wide_spatiotemporal_pattern_rev_x_201223.tiff",
    sep = ""
    ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )

#---------------------------------------------------#
# 5.3 Plot figure with actual values (points in the background) ----
#---------------------------------------------------#
# Curve
plot_gam <-
  purrr::pmap(
    .l = list(
      output_gam_pd$data, # ..1
      output_gam_pd$predicted_gam, # ..2
      output_gam_pd$vars # ..3
      ),
    .f = ~ {
      
      dat_1 <- 
        ..2 %>%
        dplyr::select(-dataset_id) %>%
        dplyr::group_by(lat) %>%
        dplyr::summarise_all(mean)
      
      ggplot2::ggplot(
        aes(
          x = lat,
          y = var
          ),
        data = ..2
        ) +
        ggplot2::geom_point(
          aes(
            x = lat, 
            y = estimate,
            colour = age
            ),
          size = 0.75,
          alpha = 0.15,
          data = ..1
          ) + 
        ggplot2::scale_color_gradient(
          high = color_high_age,
          low = color_low_age
          ) +
        ggplot2::theme(
          legend.position = "none"
          ) +
        ggplot2::geom_line(
          aes(
            group = age, 
            colour = age
            ),
          linewidth = 0.15,
          alpha = 1,
          data = ..2
          ) +
        ggplot2::geom_line(
          aes(
            x = lat,
            y = var
          ),
          linewidth = 1,
          colour = color_pd_curve,
          data = dat_1
          ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = expression(
            paste(
            'Latitude ', (degree ~ N)
            )
          ),
          y = paste(..3)
          ) +
        ggplot2::theme(
          axis.title = element_text(
            size = 15,
            color = color_common
            ),
          axis.text = element_text(
            size = 13,
            color = color_common
            ),
          legend.position = c(0.8, 0.9),
          legend.background = element_rect(
            fill = "transparent"
            ),
          legend.title = element_text(
            size = 11,
            color = color_common
            ),
          legend.text = element_text(
            size = 10,
            color = color_common
            ),
          plot.margin = margin(
            1, 0, 0, 0.25, 
            "cm"
            )  #t, r, b, l
          ) +
        ggplot2::guides(
          colour = guide_colorbar(
            barheight = 3,
            barwidth = 0.75,
            reverse = TRUE,
            label.vjust = -0.25
            )
        )
      }
    )
    
gam_curve_mpd <- 
  plot_gam[[1]] + 
  labs(
    y = "ses_MPD",
    colour = "Time (cal yr BP)"
    )


gam_curve_mntd <- 
  plot_gam[[2]] + 
  labs(
    y = "ses_MNTD",
    colour = "Time (cal yr BP)"
    )

# 5.4 Plot the model temporally ----
plot_gam_temporal <-
  purrr::pmap(
    .l = list(
      output_gam_pd$data, # ..1
      output_gam_pd$predicted_gam, # ..2
      output_gam_pd$vars # ..3
      ),
    
    .f = ~ {
      
      dat_1 <-
        ..2 %>%
        dplyr::select(-dataset_id) %>%
        dplyr::group_by(age) %>%
        dplyr::summarise_all(., mean)
      
      ggplot2::ggplot(
        aes(
          x = age,
          y = var
          ),
        data = ..2
        ) +
        ggplot2::geom_point(
          aes(
            x = age, 
            y = estimate,
            colour = lat
            ),
          size = 0.5,
          alpha = 0.2,
          data = ..1
          ) + 
        ggplot2::scale_color_gradient(
          high = color_high_lat,
          low = color_low_lat
          ) +
        ggplot2::theme(
          legend.position = "none"
          ) +
        ggplot2::geom_line(
          aes(
            group = lat,
            colour = lat
            ),
          linewidth = 0.4,
          alpha = 1,
          data = ..2
          ) +
        ggplot2::geom_line(
          aes(
              x = age,
              y = var,
              colour = lat
              ),
          linewidth = 1,
          colour = color_pd_curve,
          data = dat_1
          ) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          x = "Time (cal yr BP)",
          y = paste(..3)
          ) +
        ggplot2::scale_x_continuous(
          limits = c(0, 12000),
          breaks = seq(0, 12000, 2000)
          ) +
        ggplot2::theme(
          axis.title = element_text(
            size = 15,
            color = color_common
            ),
          axis.text = element_text(
            size = 13,
            color = color_common
            ),
          axis.text.x = element_text(
            angle = 45,
            hjust = 1
            ),
          legend.position = c(0.80, 0.85),
          legend.background = element_rect(
            fill = "transparent"
            ),
          legend.title = element_text(
            size = 11,
            color = color_common
            ),
          legend.text = element_text(
            size = 10,
            color = color_common
            ),
          plot.margin = margin(
            1, 0, 0, 0.25, 
            "cm"
            )  #t, r, b, l
          ) +
        ggplot2::guides(
          colour = guide_colorbar(
            barheight = 3,
            barwidth = 0.75,
            reverse = TRUE
          )
        )
    }
  )

gam_curve_mpd_temporal <- 
  plot_gam_temporal[[1]] + 
  ggplot2::labs(
    y = "ses_MPD",
    colour = expression(
      paste(
        'Latitude ', (degree ~ N)
      )
    )
  ) 
gam_curve_mpd_temporal_rev_x <- 
  gam_curve_mpd_temporal + 
  ggplot2::scale_x_reverse(
    limits = c(12000, 0),
    breaks = seq(12000, 0, by = -2000)
  ) 


gam_curve_mntd_temporal <- 
  plot_gam_temporal[[2]] + 
  ggplot2::labs(
    y = "ses_MNTD",
    colour = expression(
      paste(
        'Latitude ', (degree ~ N)
      )
    )
  ) 
gam_curve_mntd_temporal_rev_x <- 
  gam_curve_mntd_temporal + 
  ggplot2::scale_x_reverse(
    limits = c(12000, 0),
    breaks = seq(12000, 0, by = -2000)
  ) 


# Merge all the figures into one
final_fig_mpd <-
  ggpubr::ggarrange(
    gam_curve_mpd,
    gam_curve_mpd_temporal,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    hjust = -0.4,
    vjust = 2.5
  )

final_fig_mpd_rev_x <-
  ggpubr::ggarrange(
    gam_curve_mpd,
    gam_curve_mpd_temporal_rev_x,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    hjust = -0.4,
    vjust = 2.5
  )

final_fig_mntd <-
  ggpubr::ggarrange(
    gam_curve_mntd,
    gam_curve_mntd_temporal,
    ncol = 1,
    nrow = 2,
    labels = c("(c)", "(d)"),
    hjust = -0.4,
    vjust = 2.5
  )
final_fig_mntd_rev_x <-
  ggpubr::ggarrange(
    gam_curve_mntd,
    gam_curve_mntd_temporal_rev_x,
    ncol = 1,
    nrow = 2,
    labels = c("(c)", "(d)"),
    hjust = -0.4,
    vjust = 2.5
  )

final_composite <- 
  ggpubr::ggarrange(
    final_fig_mpd,
    final_fig_mntd,
    ncol = 2,
    nrow = 1
  )

final_composite_rev_x <- 
  ggpubr::ggarrange(
    final_fig_mpd_rev_x,
    final_fig_mntd_rev_x,
    ncol = 2,
    nrow = 1
  )

ggplot2::ggsave(
  final_composite,
  filename = paste(
    "Outputs/Figure/",
    "Holocene_wide_spatiotemporal_pattern_201223_suppl.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
)

ggplot2::ggsave(
  final_composite_rev_x,
  filename = paste(
    "Outputs/Figure/",
    "Holocene_wide_spatiotemporal_pattern_201223_suppl_rev_x.tiff",
    sep = ""
  ),
  height = 20,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
)
