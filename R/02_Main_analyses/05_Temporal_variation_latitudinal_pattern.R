#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Temporal variation in latitudinal pattern of phylogenetic dispersion ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Load configuration file ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load the data ----
#--------------------------------------------------------#
source_data <- 
  read_rds("Inputs/Data/source_data_191223.rds")

#--------------------------------------------------------#
# 3. Fit and plot the GAM models ----
#--------------------------------------------------------#
# 1000-year time bin ----
data_filtered <- 
  source_data %>%
  dplyr::arrange(age) %>% 
  dplyr::mutate(
    period = 
      ceiling(age / 1000)
    ) %>%
  dplyr::mutate(period = period * 1000) %>%
  dplyr::mutate_at("period", as_factor) 

data_gam_period <-
  data_filtered %>%
  tidyr::gather(
    c(
      mpd,
      mntd
      ),
    key = "vars",
    value = "estimate"
    ) %>%
  dplyr::group_by(vars) %>%
  tidyr::nest() %>%
  dplyr::ungroup() 

gam_mod_temporal_1k <-
  data_gam_period %>%
  dplyr::mutate(
    gam_model =
      purrr::map(
        .x = data, 
        .f = ~ {
          data <- .x
          
          set.seed(2468)
          
          mod <-
            mgcv::gam(
              estimate ~
                lat + 
                s(lat, 
                  by = period, 
                  bs = 'tp',
                  m = 1
                  ) +
                s(age, 
                  k = 10, 
                  bs = 'tp'
                  ) +
                s(age,  
                  by = period, 
                  bs = 'tp', 
                  m = 1
                  ) +
                s(dataset_id,
                  k = 99, 
                  bs = 're'
                  ) +
                s(period, 
                  k = 12, 
                  bs = 'fs'
                  ) +
                ti(lat, age, 
                   by = period,
                   bs = c('tp', 'tp')
                   ),
              data = data,
              method = "REML",
              family = "gaussian",
              control = gam.control(
                trace = TRUE, 
                maxit = 200
                ),
              weights = age_uncertainty_index
              
            )
          }
        )
    )


#--------------------------------------------------------#
# Save model ----
#--------------------------------------------------------#
readr::write_rds(
  gam_mod_temporal_1k,
  file = "Outputs/Data/Model_1k_period_201223.rds",
  compress = "gz"
  )

#--------------------------------------------------------#
# 4. Extract summary of GAM models ----
#--------------------------------------------------------#
temporal_gam_summary_mpd <-
  flextable::as_flextable(
    gam_mod_temporal_1k[1,]$gam_model[[1]]
    )

save_as_docx(
  temporal_gam_summary_mpd,
  path = "Outputs/Table/MPD_1k_201223.docx"
  )

temporal_gam_summary_mntd <-
  flextable::as_flextable(
    gam_mod_temporal_1k[2, ]$gam_model[[1]]
    )

save_as_docx(
  temporal_gam_summary_mntd,
  path = "Outputs/Table/MNTD_1k_201223.docx"
  )

#--------------------------------------------------------#
# 5. Plot the models ----
#--------------------------------------------------------#
data_plot_1k <- 
  gam_mod_temporal_1k %>% 
  dplyr::mutate(
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
                lat = seq(
                  min(lat), 
                  max(lat), 
                  by = 0.25
                ),
                dataset_id = dataset_id[1],
                age = mean(age),
                period = seq(1000, 12000, by = 1000)
              )
            )
          not_inlude <-
            gratia::smooths(.y) %>%
            str_subset(., "dataset_id|age")
          crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
          
          set.seed(2468)
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
  ) %>% 
  dplyr::select(
    vars, 
    predicted_gam
  ) %>% 
  tidyr::unnest(predicted_gam) 


data_plot_1k$vars <- 
  factor(data_plot_1k$vars, 
         levels = c(
           "mpd",
           "mntd")
         )
indices <- 
  c(
    `mpd` = "ses_MPD",
    `mntd` = "ses_MNTD"
    )

plot_1k <- 
  data_plot_1k %>% 
  ggplot(
    aes(
      x = lat,
      y = var,
      group = period,
      colour = period
      )
    ) +
  ggplot2::theme_classic() +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_ribbon(
    aes(
      ymin = lwr,
      ymax = upr,
      fill = period
      ),
    alpha = 1/8, 
    colour = NA
    ) + 
  ggplot2::scale_color_gradient(
    high = color_high_age, 
    low = color_low_age,
    trans = "reverse"
    ) + 
  ggplot2::scale_fill_gradient(
    high = color_high_age,
    low = color_low_age
    ) +
  ggplot2::facet_wrap(
    ~ vars, 
    scales = "free_y",
    labeller = as_labeller(indices)
    ) +
  ggplot2::labs(
    fill = "Period (cal yr BP)",
    x = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    y = "Estimate"
    ) +
  ggplot2::theme(
    axis.title = element_text(
      size = 14,
      color = color_common
      ),
    axis.text = element_text(
      size = 12, 
      color = color_common
      ),
    strip.text.x = element_text(
      size = 14
      ),
    legend.position = "bottom",
    legend.background = element_rect(
      fill = "transparent"
      ),
    legend.spacing.x = unit(0.75, "cm"),
    legend.title = element_text(
      size = 12
    )
  ) + 
  ggplot2::guides(
    colour = "none",
    fill = guide_colorbar(
      barheight = 1, 
      barwidth = 7.5
    )
  )

# Save plot ----
ggsave(
  plot_1k,
  filename = "Outputs/Figure/GAM_1k_bin_201223.tiff",
  dpi = 400,
  width = 20,
  height = 10,
  units = "cm",
  compress = "lzw"
)

#--------------------------------------------------------#
# 6. Test of differences in slopes of models of different periods ----
#--------------------------------------------------------#
# Based on "https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/"

# Make a new dataset to be predicted
gam_mod_temporal_1k <- 
  readr::read_rds(
    "Outputs/Data/Model_1k_period_201223.rds"
    )

new_data <-
  with(gam_mod_temporal_1k[1, ]$data[[1]],
       base::expand.grid(
         lat = seq(
           min(lat),
           max(lat), 
           by = 0.25
           ),
         dataset_id = dataset_id[1],
         age = mean(age),
         period = seq(1000, 12000, by = 1000)
         )
       )

# MPD
mod_mpd <- 
  gam_mod_temporal_1k[1,]$gam_model[[1]]


# Estimate the difference in the smooth term (slope) of the model between the groups
# Where the confidence interval excludes zero, we might infer significant differences 
# between a pairs of estimated smooths.
diff_pattern_mpd <-
  gratia::difference_smooths(
    mod_mpd,
    smooth = "s(lat)",
    newdata = new_data,
    ci_level = 0.95,
    n = 500,
    method = "REML"
    ) %>%
  dplyr::mutate_at("level_1", as.factor) %>%
  dplyr::mutate_at("level_2", as.factor) 


# Set desired order of periods
diff_pattern_mpd$level_1 <- 
  factor(diff_pattern_mpd$level_1, 
         levels = c(
           "1000",
           "2000",
           "3000",
           "4000",
           "5000",
           "6000",
           "7000",
           "8000",
           "9000",
           "10000",
           "11000"
           )
         )
diff_pattern_mpd$level_2 <- 
  factor(diff_pattern_mpd$level_2, 
         levels = c(
           "2000",
           "3000",
           "4000",
           "5000",
           "6000",
           "7000",
           "8000",
           "9000",
           "10000",
           "11000",
           "12000"
           )
         )

# Data for shaded bars at the points of significant change in smooths
plus_mpd <- 
  diff_pattern_mpd %>% 
  dplyr::filter(lower > 0.05 & upper > 0.05) %>% 
  dplyr::group_by(level_1, level_2) %>% 
  dplyr::summarise(
    min_lat = min(lat),
    max_lat = max(lat)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    diff_lat = max_lat - min_lat
  ) %>% 
  dplyr::filter(diff_lat > 2)

minus_mpd <- 
  diff_pattern_mpd %>% 
  dplyr::filter(lower < -0.05 & upper < -0.05) %>% 
  group_by(level_1, level_2) %>% 
  summarise(
    min_lat = min(lat),
    max_lat = max(lat)
  ) %>% 
  ungroup() %>% 
  dplyr::mutate(
    diff_lat = max_lat - min_lat
  ) %>% 
  dplyr::filter(diff_lat > 2)

mpd_period_1000 <-
  ggplot2::ggplot(
    diff_pattern_mpd,
    aes(
      x = lat,
      y = diff
      )
    ) +
  ggplot2::geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
      ),
    alpha = 0.15
    ) +
  ggplot2::geom_line(
    color = "#0066CC",
    linewidth = 1
  ) +
  lemon::facet_rep_grid(
    level_2 ~ level_1,
    scales = "free_y"
    ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.75
    ) +
  geom_rect(data = plus_mpd,
            aes(
              xmin = min_lat,
              xmax = max_lat,
              ymin = -Inf,
              ymax = Inf,
              linewidth = 0.1
            ),
            fill = "#E69F00",
            alpha = 0.4,
            inherit.aes = FALSE
            ) +
  geom_rect(data = minus_mpd,
            aes(
              xmin = min_lat,
              xmax = max_lat,
              ymin = -Inf,
              ymax = Inf,
              linewidth = 0.1
              ),
            fill = "#E69F00",
            alpha = 0.4,
            inherit.aes = FALSE
            ) +
  ggplot2::labs(
    x = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    y = 'Difference in latitudinal trends of ses_MPD'
    ) + 
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    strip.text.x = element_text(
      size = 23
      ),
    strip.text.y = element_text(
      size = 23
      ),
    axis.text = element_text(
      color = color_common, 
      size = 20
      ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
      ),
    axis.title = element_text(
      color = color_common,
      size = 30
      ),
    panel.border = element_blank(), 
    axis.line = element_line(),
    panel.spacing = unit(-2, "lines")
    )

ggplot2::ggsave(
  mpd_period_1000,
  filename = "Outputs/Figure/MPD_slope_diff_1000_201223.tiff",
  dpi = 400,
  width = 40,
  height = 35,
  units = "cm",
  compress = "lzw"
  )


# MNTD
mod_mntd <- 
  gam_mod_temporal_1k[2,]$gam_model[[1]]

# Estimate the difference in the smooth term (slope) of the model between the groups
diff_pattern_mntd <-
  gratia::difference_smooths(
    mod_mntd,
    smooth = "s(lat)",
    newdata = new_data,
    ci_level = 0.95,
    n = 500,
    method = "REML"
    ) 

# Set desired order of periods
diff_pattern_mntd$level_1 <- 
  factor(diff_pattern_mntd$level_1,
         levels = c(
           "1000",
           "2000",
           "3000",
           "4000",
           "5000",
           "6000",
           "7000",
           "8000",
           "9000",
           "10000",
           "11000"
           )
         )
diff_pattern_mntd$level_2 <- 
  factor(diff_pattern_mntd$level_2, 
         levels = c(
           "2000",
           "3000",
           "4000",
           "5000",
           "6000",
           "7000",
           "8000",
           "9000",
           "10000",
           "11000",
           "12000"
           )
         )

# Data for shaded bars at the points of significant change in smooths
plus <- 
  diff_pattern_mntd %>% 
  dplyr::filter(lower > 0.05 & upper > 0.05) %>% 
  dplyr::group_by(level_1, level_2) %>% 
  dplyr::summarise(
    min_lat = min(lat),
    max_lat = max(lat)
    ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    diff_lat = max_lat - min_lat
    ) %>% 
  dplyr::filter(diff_lat > 2)

minus <- 
  diff_pattern_mntd %>% 
  dplyr::filter(lower < -0.05 & upper < -0.05) %>% 
  group_by(level_1, level_2) %>% 
  summarise(
    min_lat = min(lat),
    max_lat = max(lat)
    ) %>% 
  ungroup() %>% 
  dplyr::mutate(
    diff_lat = max_lat - min_lat
  ) %>% 
  dplyr::filter(diff_lat > 2)

mntd_period_1000 <-
  ggplot2::ggplot(
    diff_pattern_mntd,
    aes(x = lat,
        y = diff)
    ) +
  ggplot2::geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
      ),
    alpha = 0.15
    ) +
  ggplot2::geom_line(color = "#0066CC",
                     linewidth = 1
                     ) +
  lemon::facet_rep_grid(
    level_2 ~ level_1,
    scales = "free_y"
    ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "red",
    linewidth = 0.75
    ) +
    geom_rect(data = plus,
              aes(
                xmin = min_lat,
                xmax = max_lat,
                ymin = -Inf,
                ymax = Inf,
                linewidth = 0.1
              ),
              fill = "#E69F00",
              alpha = 0.4,
              inherit.aes = FALSE
              ) +
  geom_rect(data = minus,
            aes(
              xmin = min_lat,
              xmax = max_lat,
              ymin = -Inf,
              ymax = Inf,
              linewidth = 0.1
            ),
            fill = "#E69F00",
            alpha = 0.4,
            inherit.aes = FALSE
            ) +
  ggplot2::labs(
    x = expression(
      paste(
        'Latitude ', (degree ~ N)
        )
      ),
    y = 'Difference in latitudinal trends of ses_MNTD'
    ) + 
  ggplot2::theme_classic() +
  ggplot2::theme(
    legend.position = "none",
    strip.text.x = element_text(
      size = 23
      ),
    strip.text.y = element_text(
      size = 23
      ),
    axis.text = element_text(
      color = color_common,
      size = 20
      ),
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1
      ),
    axis.title = element_text(
      color = color_common, 
      size = 30
      ),
    panel.border = element_blank(), 
    axis.line = element_line(),
    panel.spacing = unit(-2, "lines")
    )

#--------------------------------------------------------#
# Save plot ----
#--------------------------------------------------------#
ggplot2::ggsave(
  mntd_period_1000,
  filename = "Outputs/Figure//MNTD_slope_diff_1000_201223.tiff",
  dpi = 400,
  width = 40,
  height = 35,
  units = "cm",
  compress = "lzw"
  )


