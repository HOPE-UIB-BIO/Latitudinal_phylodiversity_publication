#----------------------------------------------------------#

# Latitudinal gradients in the phylogenetic assembly of angiosperms in Asia 
# during the Holocene

# Relationship between phylogenetic dispersion (PD) and climatic variables ----
#                          
#----------------------------------------------------------#

#--------------------------------------------------------#
# 1. Load configuration file ----
#--------------------------------------------------------#
source("R/00_Config_file.R")

#--------------------------------------------------------#
# 2. Load Holocene-wide latitudinal model of PD ----
#--------------------------------------------------------#
output_gam_pd <- 
  read_rds("Outputs/Data/Overall_gam_lat_PD_201223.rds")

#--------------------------------------------------------#
# 3. Load the source data ----
#--------------------------------------------------------#
source_data <- 
  read_rds("Inputs/Data/source_data_191223.rds")

#--------------------------------------------------------#
# 4. Load the climate data predicted for each sample of the datasets ----
# Climatic data was predicted for the global fossil pollen data of the 
# project, and data used here is just sourced from that.
#--------------------------------------------------------#
climate_data <- 
  readr::read_rds(
    paste(
      "Inputs/Data/Chelsa_climate/",
      "data_climate_pred-2021-12-16.rds",
      sep = "")
    ) %>% 
  dplyr::filter(dataset_id %in% unique(source_data$dataset_id)) %>% 
  dplyr::select(
    dataset_id, 
    clim_data_pred
    ) %>% 
  tidyr::unnest(clim_data_pred) %>% 
  dplyr::select(-depth)

#--------------------------------------------------------#
# 5. Combine the climate data and filtered phylodiversity data ----
#--------------------------------------------------------#
combined_data <-
  source_data %>%
  dplyr::inner_join(
    climate_data,
    by = c(
      "dataset_id", 
      "sample_id",
      "age"
      )
    ) 
#--------------------------------------------------------#
# 6. Fit the latitudinal GAM models of climate ----
#--------------------------------------------------------#
data_gam <-
  combined_data %>%
  dplyr::select(
    dataset_id,
    sample_id,
    lat,
    age,
    age_uncertainty_index,
    temp_cold,
    prec_summer,
    prec_winter
    ) %>% 
  dplyr::mutate_at("dataset_id", as.factor) %>% 
  tidyr::gather(
    c(
      temp_cold,
      prec_summer,
      prec_winter
      ),
    key = "vars",
    value = "climate") %>%
  dplyr::group_by(vars) %>%
  tidyr::nest() %>%
  dplyr::ungroup() 

gam_mod_climate <-
  data_gam %>%
  dplyr::mutate(
    gam_mod =
      purrr::map(
        .x = data,
        .f = ~ {
          data <- .x
          set.seed(246)
          mod_mpd_climate <-
            mgcv::gam(
              climate ~
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
              data = data,
              method = "REML",
              family = "gaussian",
              weights = age_uncertainty_index,
              control = gam.control(
                trace = TRUE, 
                maxit = 200
                )
            )
          return(mod_mpd_climate)
        }
      ),
    predicted_mod = purrr::map2(
      .x = data,
      .y = gam_mod,
      .f = ~ {
        data <- .x
        new_data <-
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
          str_subset(.,c("dataset_id"))
        crit <- qnorm((1 - 0.89) / 2, lower.tail = FALSE)
        
        set.seed(246)
        predicted_mod <-
          new_data %>%
          dplyr::bind_cols(
            predict(
              .y,
              newdata = new_data,
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
              c("fit",
                "se.fit")
            )
          )
        return(predicted_mod)
      }
    )
  )

#--------------------------------------------------------#
# 7. Save the model ----
#--------------------------------------------------------#
readr::write_rds(
  gam_mod_climate,
  file = "Outputs/Data/GAM_latitude_climate_201223.rds"
  )

#--------------------------------------------------------#
# 8. Plot the models (predicted values of PD vs predicted values of climate) ----
#--------------------------------------------------------#
temp_cold_predicted <- 
  gam_mod_climate[1,]$predicted_mod[[1]] %>% 
  dplyr::select(
    lat, 
    age, 
    temp_cold = var
  )

prec_summer_predicted <- 
  gam_mod_climate[2,]$predicted_mod[[1]] %>% 
  dplyr::select(
    lat, 
    age, 
    prec_summer = var
  )
prec_winter_predicted <- 
  gam_mod_climate[3,]$predicted_mod[[1]] %>% 
  dplyr::select(
    lat, 
    age, 
    prec_winter = var
  )        

mpd_predicted <-    
  output_gam_pd[1,]$predicted_gam[[1]] %>% 
  dplyr::select(
    lat, 
    age, 
    ses_MPD = var
  )  

mntd_predicted <-    
  output_gam_pd[2,]$predicted_gam[[1]] %>% 
  dplyr::select(
    lat, 
    age, 
    ses_MNTD = var
  )  

predicted_climate_pd <- 
  dplyr::inner_join(
    mpd_predicted,
    mntd_predicted,
    by = c("lat", "age")
  ) %>% 
  dplyr::inner_join(
    temp_cold_predicted,
    by = c("lat", "age")
  ) %>% 
  dplyr::inner_join(
    prec_summer_predicted,
    by = c("lat", "age")
  ) %>% 
  dplyr::inner_join(
    prec_winter_predicted,
    by = c("lat", "age")
  ) %>% 
  tidyr::gather(
    c(
      temp_cold,
      prec_summer,
      prec_winter
    ),
    key = "climate",
    value = "estimate"
  ) %>% 
  dplyr::mutate_at("climate", as.factor) 

predicted_climate_pd$climate <- 
  factor(
    predicted_climate_pd$climate,
    levels = (
      c(
        "temp_cold",
        "prec_summer",
        "prec_winter"
      )
    )
  )


# plots_mpd
dat_label <- 
  data.frame(
    label = c("(A)", "(B)", "(C)"),
    climate   = c("temp_cold", "prec_summer", "prec_winter")
  ) %>% 
  dplyr::mutate_at("label", as.factor)
dat_label$climate <- 
  factor(
    dat_label$climate,
    levels = c("temp_cold", "prec_summer", "prec_winter")
  )

plots_mpd <-   
  predicted_climate_pd %>%
  dplyr::group_by(climate) %>%
  ggplot(aes(
    x = estimate,
    y = ses_MPD,
    group = as_factor(climate)
  )
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    colour = color_common,
    linewidth = 1
  ) +
  ggplot2::facet_wrap( 
    ~ climate,
    scales = "free",
    ncol = 1
  ) +
  ggplot2::geom_text(
    aes(
      x = -Inf, 
      y = Inf, 
      hjust = -0.1,
      vjust = 1,
      label = label
    ),
    data = dat_label,
    size = 5
  ) + 
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Climate"
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
      size = 14,
      color = color_common
    )
  ) 


#plots_mntd
dat_label <- 
  data.frame(
    label = c("(D)", "(E)", "(F)"),
    climate   = c("temp_cold", "prec_summer", "prec_winter")
  ) %>% 
  dplyr::mutate_at("label", as.factor)
dat_label$climate <- 
  factor(
    dat_label$climate,
    levels = c("temp_cold", "prec_summer", "prec_winter")
  )

plots_mntd <-   
  predicted_climate_pd %>%
  dplyr::group_by(climate) %>%
  ggplot(aes(
    x = estimate,
    y = ses_MNTD,
    group = as_factor(climate)
  )
  ) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    colour = color_common,
    linewidth = 1
  ) +
  ggplot2::facet_wrap( 
    ~ climate,
    scales = "free",
    ncol = 1
  ) +
  ggplot2::geom_text(
    aes(
      x = -Inf, 
      y = Inf, 
      hjust = -0.1,
      vjust = +1,
      label = label
    ),
    data = dat_label,
    size = 5
  ) + 
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Climate"
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
      size = 14,
      color = color_common
    )
  ) 


final_plot <- 
  ggpubr::ggarrange(
    plots_mpd, 
    plots_mntd, 
    ncol = 2
  )
ggplot2::ggsave(
  final_plot,filename = "Outputs/Figure/MPD_MNTD_vs_climate_211223.tiff",
  height = 15,
  width = 20,
  units = "cm",
  dpi = 400,
  compression = "lzw"
  )



