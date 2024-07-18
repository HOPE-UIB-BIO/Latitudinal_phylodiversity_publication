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
#  Plot the matrix of variables ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Source configuration ----
#--------------------------------------------------------#

source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
full_matrix <- 
  readr::read_rds("Outputs/Data/Space_time_matrix_201223.rds")

#--------------------------------------------------------#
# 3. Plot line and heat plots ----
#--------------------------------------------------------#

# 3.1 Number of samples ----
n_sample_curve <-
  full_matrix %>%
  ggplot2::ggplot() +
  ggplot2::geom_smooth(
    aes(
      x = lat,
      y = log(n_samples)
      ),
    method = 'gam',
    colour = color_common,
    linewidth = 1.5,
    se = FALSE
    ) +
  ggplot2::geom_line(
    aes(
    x = lat,
    y = log(n_samples),
    group = age,
    colour = age),
  linewidth = 0.3,
  alpha = 0.3
  ) +
  ggplot2::scale_color_gradient(
    high = color_high_age,
    low = color_low_age
    ) +
  ggplot2::theme_classic() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = expression(
      paste('Latitude ', (degree ~ N)
            )
      ),
    y = "Number of samples",
    colour = "Time (cal yr BP)"
    ) +
  ggplot2::scale_y_continuous(
    breaks = c(1, 2, 3, 4),
    labels = round(
      c(exp(1), exp(2), exp(3), exp(4)), 0
      )
    ) +
  ggplot2::theme(
    axis.text = element_text(
      color = color_common,
      size = 12
      ),
    axis.title = element_text(
      color = color_common,
      size = 14
      ),
    legend.spacing.y = unit(0.35, 'cm'),
    legend.position = "bottom"
    )

plot_n_samples <-  
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = lat,
      fill = n_samples
      )
    ) +
  ggplot2::geom_tile() +
  ggpubr::theme_classic2() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0,0),# to coincide the heat plot and ggplot background
    breaks = seq(25, 65, by = 5)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),# to coincide the heat plot and ggplot background
    breaks = seq(0, 12000, by = 1000)
    ) +
  ggplot2::labs(
    x = "Time (cal yr BP)",
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "n_samples"
    ) +
  ggplot2::theme(
    axis.text.x = 
      element_text(
        color = color_common,
        size = 12,
        angle = 45,
        hjust = 1
        ),
    axis.text.y = 
      element_text(
        color = color_common,
        size = 12
        ),
    axis.title = 
      element_text(
        color = color_common,
        size = 14),
    legend.position = "bottom"
    )

final_n_samples <- 
  ggpubr::ggarrange(
    n_sample_curve,
    plot_n_samples,
    ncol = 2,
    nrow = 1,
    labels = c("(A)", "(B)")
    )


# Save final plot of samples ----
ggplot2::ggsave(
  final_n_samples,
  filename = "Outputs/Figure/N_samples_201223.tiff",
  width = 20,
  height = 10, 
  units = "cm",
  dpi = 400,
  compress = "lzw"
  )


# 3.2 Age uncertainty ----
age_error_curve <- 
  full_matrix %>%
  ggplot2::ggplot() +
  ggplot2::geom_smooth(
    aes(
      x = lat, 
      y = log(age_error)
      ),
    method = 'gam',
    colour = color_common,
    linewidth = 1.5,
    se = FALSE
    ) +
  ggplot2::geom_line(
    aes(
      x = lat,
      y = log(age_error),
      group = age,
      colour = age
      ),
    linewidth = 0.3,
    alpha = 0.3
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_age,
    low = color_low_age
    ) +
  ggplot2::theme_classic() +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = expression(
    paste(
      'Latitude ', (degree ~ N)
      )
    ),
    y = "Age uncertainty (years)",
    colour = "Age_uncertainty"
    ) + 
  ggplot2::scale_y_continuous(
    breaks = c(5, 6, 7, 8, 9), 
    labels = round(
      c(exp(5), exp(6), exp(7), exp(8), exp(9)), 0)
    ) + 
  ggplot2::theme(
    axis.text = element_text(
      color = color_common,
      size = 12
      ),
    axis.title = element_text(
      color = color_common,
      size = 14
      ),
    legend.position = "bottom"
    )

plot_age_error <-  
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = lat,
      fill = age_error
      )
    ) +
  ggplot2::geom_tile() +
  ggpubr::theme_classic2() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0,0),# to coincide the heat plot and ggplot background
    breaks = seq(25, 65, by = 5)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),# to coincide the heat plot and ggplot background
    breaks = seq(0, 12000, by = 1000)
    ) +
  ggplot2::labs(
    x = "Time (yr BP)",
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "Age_uncertainty"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 12,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 12
      ),
    axis.title = element_text(
      color = color_common,
      size = 14
      ),
    legend.position = "bottom"
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(barwidth = 7)
    )

final_age_error <- 
  ggpubr::ggarrange(
    age_error_curve,
    plot_age_error,
    ncol = 2,
    nrow = 1,
    labels = c("(A)", "(B)")
    )

# Save final plot of age uncertainty ----
ggplot2::ggsave(
  final_age_error,
  filename = "Outputs/Figure/Age_uncertainty_201223.tiff",
  width = 20,
  height = 10, 
  units = "cm",
  dpi = 400,
  compress = "lzw"
  )

# 3.3 Climate variables ----
dat1 <- 
  full_matrix %>%
  dplyr::group_by(lat) %>%
  dplyr::summarise_all(mean)

dat2 <- 
  full_matrix %>%
  dplyr::group_by(age) %>%
  dplyr::summarise_all(mean)

# Annual mean temperature ----
annual_temp_curve <-
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = lat,
      y = temp_annual,
      group = age,
      colour = age
      )
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
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
      y = temp_annual
      ),
    linewidth = 1,
    colour = color_common,
    data = dat1
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'ann_temp  ', (degree ~ C)
        )
      ),
    colour = "Time \n(cal yr BP)"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )

annual_temp_temporal <-
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = temp_annual,
      group = lat,
      colour = lat
      )
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_lat,
    low = color_low_lat
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) +
  ggplot2::geom_line(
    aes(x = age,
        y = temp_annual
        ),
    linewidth = 1,
    colour = color_common,
    data = dat2
    ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = element_blank(),
    y = element_blank(),
    colour = expression(
      paste(
        'Lat ', (degree ~ N)
        )
      )
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    breaks = seq(0, 12000, 2000)
    ) +
  ggplot2::theme(
    axis.title = element_text(
      size = 16,
      color = color_common
      ),
    axis.text = element_text(
      size = 14,
      color = color_common
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    legend.position = "right",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )   

annual_temp_heatplot <-
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = age,
      y = lat, 
      fill = temp_annual
      )
    ) +
  ggplot2::geom_tile() +
  ggpubr::theme_classic2() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0,0),
    breaks = seq(25, 65, by = 10)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),
    breaks = seq(0, 12000, by = 2000)
    ) + 
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "ann_temp"
    ) +
  ggplot2::theme(
    axis.text.x = 
      element_text(
        color = color_common,
        size = 14,
        angle = 45,
        hjust = 1
        ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(
      reverse = TRUE
      )
    )

final_annual_temp <- 
  ggpubr::ggarrange(
    annual_temp_curve,
    annual_temp_temporal,
    annual_temp_heatplot,
    ncol = 3,
    nrow = 1
    ) +
  ggplot2::theme(
    plot.margin = margin(1, 0, 0, 0.1, "cm")
    ) 


# Minimum temperature of the coldest month ----
temp_cold_curve <- 
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = lat,
      y = temp_cold,
      group = age,
      colour = age)
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) + 
  ggplot2::scale_color_gradient(
    high = color_high_age,
    low = color_low_age
    ) +
  ggplot2::theme_classic() +
  ggplot2::geom_line(
    aes(x = lat,
        y = temp_cold
        ),
    linewidth = 1,
    colour = color_common,
    data = dat1
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'temp_cold ', (degree ~ C)
        )
      ),
    colour = "Time \n(cal yr BP)"
    ) + 
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )

temp_cold_heatplot <-
  full_matrix %>% 
  ggplot2::ggplot(
    aes(x = age, 
        y = lat, 
        fill = temp_cold)
    ) +
  ggplot2::geom_tile() +
  ggpubr::theme_classic2() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66),
    expand = c(0,0),
    breaks = seq(25, 65, by = 10)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),
    breaks = seq(0, 12000, by = 2000)
    ) + 
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "temp_cold"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(
      reverse = TRUE
    )
  )

temp_cold_temporal <-
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = temp_cold,
      group = lat,
      colour = lat
      )
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_lat,
    low = color_low_lat
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) +
  ggplot2::geom_line(
    aes(
      x = age,
      y = temp_cold
      ),
    linewidth = 1,
    colour = color_common,
    data = dat2
    ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = element_blank(),
    y = element_blank(),
    colour = expression(
      paste(
        'Lat ', (degree ~ N)
        )
      )
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    breaks = seq(0, 12000, 2000)
    ) +
  ggplot2::theme(
    axis.title = element_text(
      size = 16,
      color = color_common
      ),
    axis.text = element_text(
      size = 14,
      color = color_common
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    legend.position = "right",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )   

final_temp_cold <- 
  ggpubr::ggarrange(
    temp_cold_curve,
    temp_cold_temporal,
    temp_cold_heatplot,
    ncol = 3,
    nrow = 1
    ) +
  ggplot2::theme(
    plot.margin = margin(1, 0, 0, 0.1, "cm")
    ) 


# Annual precipitation ----
annual_precip_curve <- 
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = lat,
      y = prec_annual,
      group = age,
      colour = age
      )
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
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
      y = prec_annual
      ),
    linewidth = 1,
    colour = color_common,
    data = dat1
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = bquote(
      'ann_prec (kg '*m^-2~ year^-1*')'
      ),
    colour = "Time \n(cal yr BP)"
    ) + 
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )

annual_precip_heatplot <-
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = age, 
      y = lat, 
      fill = prec_annual
      )
    ) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0, 0),
    breaks = seq(25, 65, by = 10)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),
    breaks = seq(0, 12000, by = 2000)
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
            )
      ),
    fill = "ann_prec"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(
      reverse = TRUE
    )
  )

annual_precip_temporal <-
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = prec_annual,
      group = lat,
      colour = lat
      )
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_lat,
    low = color_low_lat
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) +
  ggplot2::geom_line(
    aes(
      x = age,
      y = prec_annual
      ),
    linewidth = 1,
    colour = color_common,
    data = dat2
    ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = element_blank(),
    y = element_blank(),
    colour = expression(
      paste(
        'Lat ', (degree ~ N)
        )
      )
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    breaks = seq(0, 12000, 2000)
    ) +
  ggplot2::theme(
    axis.title = element_text(
      size = 16,
      color = color_common
      ),
    axis.text = element_text(
      size = 14,
      color = color_common
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    legend.position = "right",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )   

final_ann_prec <- 
  ggpubr::ggarrange(
    annual_precip_curve,
    annual_precip_temporal,
    annual_precip_heatplot,
    ncol = 3,
    nrow = 1
    ) +
  ggplot2::theme(
    plot.margin = margin(1, 0, 0, 0.1, "cm")
    ) 


# Summer precipitation ----
summer_precip_curve <- 
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = lat,
      y = prec_summer,
      group = age,
      colour = age
      )
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
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
      y = prec_summer
      ),
    linewidth = 1,
    colour = color_common,
    data = dat1
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = bquote(
      'sum_prec (kg '*m^-2~ quarter^-1*')'
      ),
    colour = "Time \n(cal yr BP)"
    ) + 
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )

summer_precip_heatplot <-
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = age,
      y = lat, 
      fill = prec_summer
      )
    ) + 
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0,0),
    breaks = seq(25, 65, by = 10)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),
    breaks = seq(0, 12000, by = 2000)
    ) +
  ggplot2::labs(
    x = element_blank(),
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "sum_prec"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(
      reverse = TRUE
    )
  )

summer_precip_temporal <- 
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = prec_summer,
      group = lat,
      colour = lat
      )
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_lat,
    low = color_low_lat
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) +
  ggplot2::geom_line(
    aes(
      x = age,
      y = prec_summer
      ),
    linewidth = 1,
    colour = color_common,
    data = dat2
    ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = element_blank(),
    y = element_blank(),
    colour = expression(
      paste(
        'Lat ', (degree ~ N)
        )
      )
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    breaks = seq(0, 12000, 2000)
    ) +
  ggplot2::theme(
    axis.title = element_text(
      size = 16,
      color = color_common
      ),
    axis.text = element_text(
      size = 14,
      color = color_common
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    legend.position = "right",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )   

final_prec_summer <- 
  ggpubr::ggarrange(
    summer_precip_curve,
    summer_precip_temporal,
    summer_precip_heatplot,
    ncol = 3,
    nrow = 1
    ) +
  ggplot2::theme(
    plot.margin = margin(1, 0, 0, 0.1, "cm")
    ) 

# Winter precipitation ----
winter_precip_curve <- 
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = lat,
      y = prec_win,
      group = age,
      colour = age
      )
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
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
      y = prec_win
      ),
    linewidth = 1,
    colour = color_common,
    data = dat1
    ) +
  ggplot2::labs(
    x = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    y = bquote(
         'wint_prec (kg '*m^-2~ quarter^-1*')'
         ),
    colour = "Time \n(cal yr BP)"
    ) + 
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )

winter_precip_heatplot <-
  full_matrix %>% 
  ggplot2::ggplot(
    aes(
      x = age, 
      y = lat, 
      fill = prec_win
      )
    ) + 
  ggplot2::geom_tile() +
  ggplot2::scale_fill_viridis_c(
    direction = -1,
    option = "magma"
    ) + 
  ggplot2::scale_y_continuous(
    limits = c(25, 66), 
    expand = c(0,0),
    breaks = seq(25, 65, by = 10)
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    expand = c(0,0),
    breaks = seq(0, 12000, by = 2000)
    ) +
  ggplot2::labs(
    x = "Time (cal yr BP)",
    y = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    fill = "wint_prec"
    ) +
  ggplot2::theme(
    axis.text.x = element_text(
      color = color_common,
      size = 14,
      angle = 45,
      hjust = 1
      ),
    axis.text.y = element_text(
      color = color_common,
      size = 14
      ),
    axis.title = element_text(
      color = color_common,
      size = 16
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    ) + 
  ggplot2::guides(
    fill = guide_colorbar(
      reverse = TRUE
    )
  )

prec_winter_temporal <- 
  full_matrix %>%
  ggplot2::ggplot(
    aes(
      x = age,
      y = prec_win,
      group = lat,
      colour = lat
      )
    ) +
  ggplot2::scale_color_gradient(
    high = color_high_lat,
    low = color_low_lat
    ) +
  ggplot2::geom_line(
    linewidth = 0.5,
    alpha = 1
    ) +
  ggplot2::geom_line(
    aes(
      x = age,
      y = prec_win
      ),
    linewidth = 1,
    colour = color_common,
    data = dat2
    ) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Time (cal yr BP)",
    y = element_blank(),
    colour = expression(
      paste(
        'Lat ', (degree ~ N)
        )
      )
    ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 12000),
    breaks = seq(0, 12000, 2000)
    ) +
  ggplot2::theme(
    axis.title = element_text(
    size = 16,
    color = color_common
    ),
    axis.text = element_text(
      size = 14,
      color = color_common
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    legend.position = "right",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.title = element_text(
      size = 14,
      color = color_common
      ),
    legend.text = element_text(
      size = 12,
      color = color_common
      )
    )
    
final_prec_winter <- 
  ggpubr::ggarrange(
    winter_precip_curve,
    prec_winter_temporal,
    winter_precip_heatplot,
    ncol = 3,
    nrow = 1
    ) +
  ggplot2::theme(
    plot.margin = margin(1, 0, 0, 0.1, "cm")
    )  #t, r, b, l

final_figure <- 
  ggpubr::ggarrange(
    final_annual_temp,
    final_temp_cold,
    final_ann_prec,
    final_prec_summer, 
    final_prec_winter, 
    nrow = 5,
    labels = c("(A)", "(B)", "(C)", "(D)", "(E)"),
    hjust = -2,
    vjust = 1)

#-----------------------------------------------#
# Save final climate trends ---
#-----------------------------------------------#
ggplot2::ggsave(
  final_figure,
  filename = "Outputs/Figure/Climate_trends_201223.tiff",
  width = 30,
  height = 35,
  units = "cm",
  dpi = 400,
  compress = "lzw")



