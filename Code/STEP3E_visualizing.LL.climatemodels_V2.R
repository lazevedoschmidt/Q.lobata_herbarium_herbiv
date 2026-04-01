#Purpose: Function to plot forest plots for all LEAF LEVEL models
#date started: 08/22/2024
#Last updated: 03/20/2025
  #Code was rerun following a cleaning of the original data
  #models with high Rhat (>1.00) values were removed from the dataset. 
    #this was overlooked last time but is now fixed. 
#Previous updated: 02/26/2025
  #Code was updated to include all model updates as well as create "simplified" figures
  #for the manuscript. These figures only plot the variables that do not cross the 0 line
  #with 80% CI
#Author: LAS
#R version: 4.3.3.

library(tidybayes)
library(ggplot2)
library(dplyr)
library(purrr)
library(stringr)
library(ggeffects)
library(cowplot)
library(stringr)
library(tidyverse)
library(wesanderson)

#NOTES: 
  #Update function to work with .csv files instead of models
  #import via path because of all the files. There should be some code for doing this in the TW function code?? 
  #This will also let me to clear the global environment that is crazy high right now

# Functions for Forest Plots
# Updated for combined year_doy models

# Full Forest Plot Function ----
plot_forest_models <- function(csv_paths, model_names = NULL, show_all = TRUE) {
  # Read all CSV files
  models_data <- map(csv_paths, read_csv)
  
  # If model names aren't provided, generate default names
  if (is.null(model_names)) {
    model_names <- paste("Model", seq_along(models_data))
  }
  
  # Function to process a single model's data
  process_model <- function(model_data, model_name) {
    # Get parameter names
    params <- colnames(model_data)
    params <- params[str_detect(params, "^b_") & !str_detect(params, "^b_Intercept")]
    
    summary_stats <- model_data %>%
      select(all_of(params)) %>%
      pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
      group_by(parameter) %>%
      summarize(
        mean = mean(value),
        lower_95 = quantile(value, 0.025), # 95% CI lower bound
        upper_95 = quantile(value, 0.975), # 95% CI upper bound
        lower_80 = quantile(value, 0.1),   # 80% CI lower bound
        upper_80 = quantile(value, 0.9),   # 80% CI upper bound
        positive = sum(value > 0) / length(value),
        negative = sum(value < 0) / length(value)
      ) %>%
      mutate(
        parameter = gsub("b_", "", parameter),
        model = model_name,
        # Add a flag to identify estimates that don't cross zero (80% CI)
        significant_80 = (lower_80 > 0) | (upper_80 < 0)
      )
    
    return(summary_stats)
  }
  
  # Process all models
  all_summaries <- map2_df(models_data, model_names, process_model)
  
  # Create ordered factor for parameters to group by color
  # Get all unique parameters
  all_params <- unique(all_summaries$parameter)
  
  # Group parameters by color categories
  time_params <- c("scaledoy.clean", "scaleyear.x", "scaleyear.x.scaledoy.clean")
  temp_params <- c("scaleTD", "scaleMWMT", "scaleEXT", "scaleEMT", "scaleeFFP", 
                   "scaleMAR", "scaleDD5", "scaleDD1040", "scaleDD_18", "scaleDD_0", 
                   "scaleMAT", "scalebFFP", "scaleDD18","scaleNNFD","scaleFFP")
  moisture_params <- c("scaleMAP", "scaleMSP", "scaleCMI", "scaleRH", "scaleCMD", "scaleEref","scalePAS")
  stress_params <- c("scaleAHM", "scaleSHM")
  
  # Create ordered list: Time variables first, then group by color
  param_order <- c(
    time_params[time_params %in% all_params],
    temp_params[temp_params %in% all_params], 
    moisture_params[moisture_params %in% all_params],
    stress_params[stress_params %in% all_params]
  )
  
  # Add any remaining parameters that weren't categorized
  remaining_params <- setdiff(all_params, param_order)
  param_order <- c(param_order, remaining_params)
  
  # Convert to ordered factor (reverse order so top items appear at top of plot)
  all_summaries$parameter <- factor(all_summaries$parameter, 
                                    levels = rev(param_order), 
                                    ordered = TRUE)
  
  # Color definitions for parameters (climate variables)
  parameter_color_def <- c(
    "scaleyear.x.scaledoy.clean" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scaledoy.clean" = "#C1AE7C", 
    "scaleTD" = "#E16036", 
    #"scaleMWMT" = "#E16036", 
    "scaleEXT" = "#E16036",
    #"scaleEMT" = "#E16036", 
    #"scaleeFFP" = "#E16036", 
    "scaleAHM" = "#6D4F92", 
    "scaleMAP" = "#5EB1BF", 
    "scaleMAR" = "#E16036",
    "scaleMSP" = "#5EB1BF", 
    "scaleSHM" = "#6D4F92", 
    #"scaleCMI" = "#5EB1BF", 
    "scaleRH" = "#5EB1BF", 
    #"scaleCMD" = "#5EB1BF",
    #"scaleDD5" = "#E16036", 
    #"scaleDD1040" = "#E16036", 
    #"scaleDD_18" = "#E16036",
    #"scaleDD_0" = "#E16036",
    #"scaleMAT" = "#E16036", 
    #"scalebFFP" = "#E16036",
    #"scaleEref" = "#5EB1BF", 
    #"scaleDD18" = "#E16036",
    #"scaleNFFD" = "#E16036",
    #"scaleFFP" = "#E16036",
    "scalePAS" = "#5EB1BF"
  )
  
  # Filter data based on show_all parameter
  plot_data <- if (show_all) {
    all_summaries
  } else {
    all_summaries %>% filter(significant_80)
  }
  
  # Create line type mapping for models
  # model_linetype_def <- setNames(
  #   rep(c("solid", "dashed", "dotted", "dotdash", "longdash"),
  #       length.out = length(unique(plot_data$model))),
  #   sort(unique(plot_data$model))
  # )
  model_linetype_def <- c(
    "Model 1" = "solid",
    "Model 2" = "dashed", 
    "Model 3" = "dotted",
    "Model 4" = "dotdash",
    "Model 5" = "longdash"
  )
  
  # Create the forest plot
  forest_plot <- plot_data %>%
    ggplot(aes(y = parameter, x = mean, color = parameter, linetype = model)) +
    # 95% CI with thin line
    geom_linerange(aes(xmin = lower_95, xmax = upper_95),
                   position = position_dodge(width = 0.6),
                   linewidth = 0.5) +
    # 80% CI with thick line
    geom_linerange(aes(xmin = lower_80, xmax = upper_80),
                   position = position_dodge(width = 0.6),
                   linewidth = 2.5) +
    # Point estimate
    geom_point(position = position_dodge(width = 0.6), size = 2) +
    # Zero reference line
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    theme_minimal() +
    labs(x = "Estimate", y = "Parameter") +
    # Colors by parameter, line types by model
    scale_color_manual(values = parameter_color_def, name = "Parameter") +
    scale_linetype_manual(values = model_linetype_def, name = "Model") +
    # Improved layout with legend at the right
    theme(
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  return(forest_plot)
}

# parameter_linetype_def <- c(
#   "scaleAHM" = "twodash", "scaleRH" = "solid", "scaleMAR" = "longdash", "scaleEXT" = "dotted", "scaleTD" = "dotdash",
#   "scaleSHM" = "dashed", "scaleMSP" = "dotdash", "scalePAS" = "dashed", "scaleMAP" = "dotted"
# )

# # Convenience function for showing only significant estimates (80% CI) ----
# plot_forest_significant <- function(csv_paths, model_names = NULL) {
#   return(plot_forest_models(csv_paths, model_names, show_all = FALSE))
# }
# 
# # Alternative function with separate plots for comparison ----
# plot_forest_comparison <- function(csv_paths, model_names = NULL) {
#   # Create both plots
#   full_plot <- plot_forest_models(csv_paths, model_names, show_all = TRUE)
#   sig_plot <- plot_forest_models(csv_paths, model_names, show_all = FALSE)
#   
#   # Return a list of both plots
#   return(list(
#     full = full_plot + ggtitle("All Parameter Estimates"),
#     significant = sig_plot + ggtitle("Significant Parameter Estimates (80% CI)")
#   ))
# }
######################CLIMATE MODELS##########################################
# Total plot ----

# Old leaves - frequency
otf_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/",
                       pattern = "^leaftotal.freq",
                       full.names = TRUE)
otf_list <- otf_list[!(grepl(".summary\\.csv$", otf_list) | grepl("_R2\\.csv$", otf_list))]
otf_names <- paste("Model", 1:length(otf_list))

# Young leaves - frequency
ytf_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/",
                       pattern = "^leaftotal.freq",
                       full.names = TRUE)
ytf_list <- ytf_list[!(grepl(".summary\\.csv$", ytf_list) | grepl("_R2\\.csv$", ytf_list))]
ytf_names <- paste("Model", 1:length(ytf_list))

# Old leaves - richness
otd_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/",
                       pattern = "^leaftotal.rich",
                       full.names = TRUE)
otd_list <- otd_list[!(grepl(".summary\\.csv$", otd_list) | grepl("_R2\\.csv$", otd_list))]
otd_names <- paste("Model", 1:length(otd_list))

# Young leaves - richness
ytd_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/",
                       pattern = "^leaftotal.rich",
                       full.names = TRUE)
ytd_list <- ytd_list[!(grepl(".summary\\.csv$", ytd_list) | grepl("_R2\\.csv$", ytd_list))]
ytd_names <- paste("Model", 1:length(ytd_list))

# Create plots using the updated function (show_all = TRUE for full plots)
otdplot <- plot_forest_models(otd_list, otd_names, show_all = TRUE)
otdplot
otdplot <- otdplot +
  labs(title = "Richness of \nAll Damage (late-season)") +
  scale_y_discrete(labels = c("AHM", "RH", "MAR", "EXT", "TD", "Year", "DOY")) +
  theme(legend.position = "none")
otdplot

otfplot <- plot_forest_models(otf_list, otf_names, show_all = TRUE)
otfplot
otfplot <- otfplot +
  labs(title = "Frequency of \nAll Damage (late-season)") +
  scale_y_discrete(labels = c("SHM", "AHM", "MSP", "MAR", "EXT", "Year", "DOY")) +
  theme(legend.position = "none")
otfplot

# Young leaves
ytdplot <- plot_forest_models(ytd_list, ytd_names, show_all = TRUE)
# Extract legend before removing it
modellegend <- get_legend(ytdplot)
ytdplot

ytdplot <- ytdplot +
  labs(title = "Richness of All \nDamage (early-season)") +
  scale_y_discrete(labels = c("AHM", "PAS", "MSP", "MAR", "EXT", "Year", "DOY")) +
  theme(legend.position = "none")
ytdplot

ytfplot <- plot_forest_models(ytf_list, ytf_names, show_all = TRUE)
ytfplot
ytfplot <- ytfplot +
  labs(title = "Frequency of All \nDamage (early-season)") +
  scale_y_discrete(labels = c("AHM", "PAS", "MSP", "MAR", "EXT", "Year", "DOY")) +
  theme(legend.position = "none")

# Combine all plots
totalplot <- plot_grid(ytfplot, otfplot,
                       ytdplot, otdplot,
                       nrow = 2,
                       ncol = 2)

totalplot

# 80% Total plots (significant estimates only) ----
# Create plots using show_all = FALSE for significant estimates only
otdplot80 <- plot_forest_models(otd_list, otd_names, show_all = FALSE)
otdplot80
otdplot80 <- otdplot80 +
  labs(title = "Richness of \nAll Damage (late-season)") +
  scale_y_discrete(labels = c("TD", "DOY")) +
  theme(legend.position = "none")
otdplot80

otfplot80 <- plot_forest_models(otf_list, otf_names, show_all = FALSE)
otfplot80
# otfplot80 <- otfplot80 +
#   labs(title = "Frequency of \nTotal Damage (late-season)") +
#   scale_y_discrete(labels = c("TD")) +
#   theme(legend.position = "none")

# Young leaves
ytdplot80 <- plot_forest_models(ytd_list, ytd_names, show_all = FALSE)
ytdplot80
ytdplot80 <- ytdplot80 +
  labs(title = "Richness of All \nDamage (early-season)") +
  scale_y_discrete(labels = c("AHM","MAR", "EXT", "Year")) +
  theme(legend.position = "none")
ytdplot80

ytfplot80 <- plot_forest_models(ytf_list, ytf_names, show_all = FALSE)
ytfplot80
ytfplot80 <- ytfplot80 +
  labs(title = "Frequency of All \nDamage (early-season)") +
  scale_y_discrete(labels = c("AHM","MAR", "Year")) +
  theme(legend.position = "none")
ytfplot80

# Combine significant plots
totalplot80 <- plot_grid(ytfplot80, NULL,
                         ytdplot80, otdplot80,
                         nrow = 2,
                         ncol = 2)

totalplot80

#Chewing plot----
# Updated code for new folder structure (year_doy)

# Old leaves - percent area chew
opac_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
opac_list <- opac_list[!(grepl(".summary\\.csv$", opac_list) | grepl("_R2\\.csv$", opac_list))]
opac_names <- paste("Model", 1:length(opac_list))

# Young leaves - percent area chew
ypac_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
ypac_list <- ypac_list[!(grepl(".summary\\.csv$", ypac_list) | grepl("_R2\\.csv$", ypac_list))]
ypac_names <- paste("Model", 1:length(ypac_list))

# Old leaves - chewing herbs
ochew_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ochew_list <- ochew_list[!(grepl(".summary\\.csv$", ochew_list) | grepl("_R2\\.csv$", ochew_list))]
ochew_names <- paste("Model", 1:length(ochew_list))

# Young leaves - chewing herbs
ychew_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ychew_list <- ychew_list[!(grepl(".summary\\.csv$", ychew_list) | grepl("_R2\\.csv$", ychew_list))]
ychew_names <- paste("Model", 1:length(ychew_list))

#pac plots
opacplot <- plot_forest_models(opac_list, opac_names, show_all = TRUE)
opacplot
opacplot <- opacplot +
  labs(title = "Percent Leaf \nArea Chewed (late-season)") +
  scale_y_discrete(labels=c("SHM", "RH", "MSP", "MAR", "EXT","Year","DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
opacplot

ypacplot <- plot_forest_models(ypac_list, ypac_names, show_all = TRUE)
ypacplot
ypacplot <- ypacplot +
  labs(title = "Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("SHM", "AHM", "PAS", "MAR", "EXT","Year","DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot

#chewing
ochewplot <- plot_forest_models(ochew_list, ochew_names, show_all = TRUE)
ochewplot
ochewplot <- ochewplot +
  labs(title = "Frequency of Chewing \nDamage (late-season)") +
  scale_y_discrete(labels=c("SHM", "RH", "MSP", "MAR", "EXT","Year","DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot

ychewplot <- plot_forest_models(ychew_list, ychew_names, show_all = TRUE)
ychewplot
ychewplot <- ychewplot +
  labs(title = "Frequency of Chewing Damage \n(early-season)") +
  scale_y_discrete(labels=c("AHM", "PAS", "MSP", "MAR","EXT", "Year","DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ychewplot

#chew plot
chewplot <- plot_grid(ypacplot,opacplot,
                     ychewplot, ochewplot,
                      nrow = 2,
                      ncol = 2)
chewplot

#80 Chew plots----
opacplot80 <- plot_forest_models(opac_list, opac_names, show_all = FALSE)
opacplot80
# opacplot80 <- opacplot80 +
#   labs(title = "Percent Leaf Area \nChewed (late-season)")+
#   scale_y_discrete(labels=c("MAP","TD"))+
#   theme(legend.position = "none")
# opacplot80

ypacplot80 <- plot_forest_models(ypac_list, ypac_names, show_all = FALSE)
ypacplot80
ypacplot80 <- ypacplot80 +
  labs(title = "Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("AHM"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none")
ypacplot80

#chewing
ochewplot80 <- plot_forest_models(ochew_list, ochew_names, show_all = FALSE)
ochewplot80
ochewplot80 <- ochewplot80 +
  labs(title = "Frequency of \nChewing Damage (late-season)")+
  scale_y_discrete(labels=c("RH"))+
  theme(legend.position = "none")
ochewplot80

ychewplot80 <- plot_forest_models(ychew_list, ychew_names, show_all = FALSE)
ychewplot80
ychewplot80 <- ychewplot80 +
  labs(title = "Frequency of \nChewing Damage (early-season)")+
  scale_y_discrete(labels=c("AHM", "MAR","DOY"))+
  theme(legend.position = "none")
ychewplot80

#chew plot
chewplot80 <- plot_grid(ypacplot80,NULL,
                      ychewplot80, ochewplot80,
                      nrow = 2,
                      ncol = 2)
chewplot80

#Generalized plot----
# ychew_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
#                          pattern = "^chewvalue", 
#                          full.names = TRUE)
# ychew_list <- ychew_list[!(grepl(".summary\\.csv$", ychew_list) | grepl("_R2\\.csv$", ychew_list))]
# ychew_names <- paste("Model", 1:length(ychew_list))

ogenf_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                         pattern = "^leafgen.freq", 
                         full.names = TRUE)
ogenf_list <- ogenf_list[!(grepl(".summary\\.csv$", ogenf_list) | grepl("_R2\\.csv$", ogenf_list))]
ogenf_names <- paste("Model", 1:length(ogenf_list))

ygenf_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                         pattern = "^leafgen.freq", 
                         full.names = TRUE)
ygenf_list <- ygenf_list[!(grepl(".summary\\.csv$", ygenf_list) | grepl("_R2\\.csv$", ygenf_list))]
ygenf_names <- paste("Model", 1:length(ygenf_list))

ogend_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                         pattern = "^leafgen.rich", 
                         full.names = TRUE)
ogend_list <- ogend_list[!(grepl(".summary\\.csv$", ogend_list) | grepl("_R2\\.csv$", ogend_list))]
ogend_names <- paste("Model", 1:length(ogend_list))

ygend_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                         pattern = "^leafgen.rich", 
                         full.names = TRUE)
ygend_list <- ygend_list[!(grepl(".summary\\.csv$", ygend_list) | grepl("_R2\\.csv$", ygend_list))]
ygend_names <- paste("Model", 1:length(ygend_list))

ogendplot <- plot_forest_models(ogend_list, ogend_names, show_all = TRUE)
ogendplot
ogendplot <- ogendplot +
  labs(title = "Richness of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "RH", "MSP", "MAR", "Year", "DOY"))+
  theme(legend.position = "none")
ogendplot

ogenfplot <- plot_forest_models(ogenf_list, ogenf_names, show_all = TRUE)
ogenfplot
ogenfplot <- ogenfplot +
  labs(title = "Frequency of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("SHM", "MSP", "MAP", "MAR", "EXT", "Year", "DOY"))+
  theme(legend.position = "none")
ogenfplot

#young leaves
ygendplot <- plot_forest_models(ygend_list, ygend_names, show_all = TRUE)
ygendplot
ygendplot <- ygendplot +
  labs(title = "Richness of Generalized \nDamage (early-season)")+
  scale_y_discrete(labels =c("SHM", "AHM", "MSP", "MAR", "EXT", "Year", "DOY"))+
  theme(legend.position = "none")
ygendplot

ygenfplot <- plot_forest_models(ygenf_list, ygenf_names, show_all = TRUE)
ygenfplot
ygenfplot <- ygenfplot +
  labs(title = "Frequency of Generalized \nDamage (early-season)")+
  scale_y_discrete(labels= c("SHM", "AHM", "MSP", "MAR", "EXT", "Year", "DOY"))+
  theme(legend.position = "none")
ygenfplot

genplot <- plot_grid(ygenfplot,ogenfplot, 
                     ygendplot,ogendplot,
                     nrow = 2,
                     ncol = 2)
genplot

#80 Generalized plot----
ogendplot80 <- plot_forest_models(ogend_list, ogend_names, show_all = FALSE)
ogendplot80
ogendplot80 <- ogendplot80 +
  labs(title = "Richness of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("Year", "DOY"))+
  theme(legend.position = "none")
ogendplot80

ogenfplot80 <- plot_forest_models(ogenf_list, ogenf_names, show_all = FALSE)
ogenfplot80
ogenfplot80 <- ogenfplot80 +
  labs(title = "Frequency of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels=c("MAP"))+
  theme(legend.position = "none")
ogenfplot80

#young leaves
ygendplot80 <- plot_forest_models(ygend_list, ygend_names, show_all = FALSE)
ygendplot80
ygendplot80 <- ygendplot80+
  labs(title = "Richness of \nGeneralized Damage (early-season)")+
  scale_y_discrete(labels= c("AHM"))+
  theme(legend.position = "none")
ygendplot80

ygenfplot80 <- plot_forest_models(ygenf_list, ygenf_names, show_all = FALSE)
ygenfplot80
#all cross 0

genplot80 <- plot_grid(NULL,ogenfplot80, 
                     ygendplot80,ogendplot80, 
                     nrow = 2,
                     ncol = 2)
genplot80


#Specialized plot----
osf_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE)
osf_list <- osf_list[!(grepl(".summary\\.csv$", osf_list) | grepl("_R2\\.csv$", osf_list))]
osf_names <- paste("Model", 1:length(osf_list))


ysf_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE)
ysf_list <- ysf_list[!(grepl(".summary\\.csv$", ysf_list) | grepl("_R2\\.csv$", ysf_list))]
ysf_names <- paste("Model", 1:length(ysf_list))


osd_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE)

osd_list <- osd_list[!(grepl(".summary\\.csv$", osd_list) | grepl("_R2\\.csv$", osd_list))]
osd_names <- paste("Model", 1:length(osd_list))

ysd_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE)
ysd_list <- ysd_list[!(grepl(".summary\\.csv$", ysd_list) | grepl("_R2\\.csv$", ysd_list))]
ysd_names <- paste("Model", 1:length(ysd_list))

osdplot <- plot_forest_models(osd_list, osd_names, show_all = TRUE)
osdplot
osdplot <- osdplot +
  labs(title = "Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "RH", "MAR", "EXT", "Year", "DOY"))+
  theme(legend.position = "none")
osdplot

osfplot <- plot_forest_models(osf_list, osf_names, show_all = TRUE)
osfplot
osfplot <- osfplot +
  labs(title = "Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("SHM","AHM", "MSP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
osfplot

#young leaves
ysdplot <- plot_forest_models(ysd_list, ysd_names, show_all = TRUE)
ysdplot
ysdplot <- ysdplot +
  labs(title = "Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("AHM", "PAS", "MSP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ysdplot

ysfplot <- plot_forest_models(ysf_list, ysf_names, show_all = TRUE)
ysfplot
ysfplot <- ysfplot +
  labs(title = "Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("AHM", "PAS", "MSP",  "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ysfplot

specplot <- plot_grid(ysfplot,osfplot, 
                      ysdplot,osdplot, 
                      nrow = 2,
                      ncol = 2)
specplot

#80 Specialized plots----
osdplot80 <- plot_forest_models(osd_list, osd_names, show_all = FALSE)
osdplot80
osdplot80 <- osdplot80 +
  labs(title = "Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("AHM", "RH"))+
  theme(legend.position = "none")
osdplot80

osfplot80 <- plot_forest_models(osf_list, osf_names, show_all=FALSE)
osfplot80
# osfplot80 <- osfplot80 +
#   labs(title = "Frequency of \nSpecialized Damage (late-season)")+
#   scale_y_discrete(labels = c("CMD", "TD"))+
#   theme(legend.position = "none")
# osfplot80

#young leaves
ysdplot80 <- plot_forest_models(ysd_list, ysd_names, show_all = FALSE)
ysdplot80
ysdplot80 <- ysdplot80 +
  labs(title = "Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("AHM","MAR", "EXT", "Year"))+
  theme(legend.position = "none")
ysdplot80

ysfplot80 <- plot_forest_models(ysf_list, ysf_names, show_all = FALSE)
ysfplot80
ysfplot80 <- ysfplot80 +
  labs(title = "Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("AHM","MAR", "Year"))+
  theme(legend.position = "none")
ysfplot80

specplot80 <- plot_grid(ysfplot80,NULL, 
                      ysdplot80,osdplot80, 
                      nrow = 2,
                      ncol = 2)
specplot80

#Mine plot----
opam_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
opam_list <- opam_list[!(grepl(".summary\\.csv$", opam_list) | grepl("_R2\\.csv$", opam_list))]
opam_names <- paste("Model", 1:length(opam_list))

ypam_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
ypam_list <- ypam_list[!(grepl(".summary\\.csv$", ypam_list) | grepl("_R2\\.csv$", ypam_list))]
ypam_names <- paste("Model", 1:length(ypam_list))

omf_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE)
omf_list <- omf_list[!(grepl(".summary\\.csv$", omf_list) | grepl("_R2\\.csv$", omf_list))]
omf_names <- paste("Model", 1:length(omf_list))

ymf_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE)
ymf_list <- ymf_list[!(grepl(".summary\\.csv$", ymf_list) | grepl("_R2\\.csv$", ymf_list))]
ymf_names <- paste("Model", 1:length(ymf_list))

omd_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE)
omd_list <- omd_list[!(grepl(".summary\\.csv$", omd_list) | grepl("_R2\\.csv$", omd_list))]
omd_names <- paste("Model", 1:length(omd_list))

ymd_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE)
ymd_list <- ymd_list[!(grepl(".summary\\.csv$", ymd_list) | grepl("_R2\\.csv$", ymd_list))]
ymd_names <- paste("Model", 1:length(ymd_list))


opamplot <- plot_forest_models(opam_list, opam_names, show_all = TRUE)
opamplot
opamplot <- opamplot +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("AHM", "RH", "MAP", "MAR", "EXT","Year", "DOY"))+
  theme(legend.position = "none")
opamplot

ypamplot <- plot_forest_models(ypam_list, ypam_names, show_all = TRUE)
ypamplot
ypamplot <- ypamplot +
  labs(title = "Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("SHM", "AHM", "MSP", "MAR", "EXT", "Year","DOY"))+ 
  theme(legend.position = "none")
ypamplot

omdplot <- plot_forest_models(omd_list, omd_names, show_all = TRUE)
omdplot
omdplot <- omdplot +
  labs(title = "Richness of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "MAP", "MAR", "EXT","Year","DOY"))+
  theme(legend.position = "none")
omdplot

omfplot <- plot_forest_models(omf_list, omf_names, show_all = TRUE)
omfplot
omfplot <- omfplot +
  labs(title = "Frequency of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "MAP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
omfplot

#young leaves
ymdplot <- plot_forest_models(ymd_list, ymd_names, show_all = TRUE)
ymdplot
ymdplot <- ymdplot +
  labs(title = "Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "MSP", "MAR", "EXT","Year","DOY"))+
  theme(legend.position = "none")
ymdplot

ymfplot <- plot_forest_models(ymf_list, ymf_names, show_all = TRUE)
ymfplot
ymfplot <- ymfplot +
  labs(title = "Frequency of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("SHM", "AHM","MSP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ymfplot

mineplot <- plot_grid(ypamplot,opamplot, 
                      ymfplot,omfplot, 
                      ymdplot,omdplot, 
                      nrow = 3,
                      ncol = 2)
mineplot

#80 Mine plot----
opamplot80 <- plot_forest_models(opam_list, opam_names, show_all = FALSE)
opamplot80
opamplot80 <- opamplot80 +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("Year"))+
  theme(legend.position = "none")
opamplot80

ypamplot80 <- plot_forest_models(ypam_list, ypam_names, show_all = FALSE)
ypamplot80
# all cross zero

omdplot80 <- plot_forest_models(omd_list, omd_names, show_all = FALSE)
omdplot80
# omdplot80 <- omdplot80 +
#   labs(title = "Richness of \nMine Damage (late-season)")+
#   scale_y_discrete(labels = c("eFFP"))+
#   theme(legend.position = "none")
# omdplot80

omfplot80 <- plot_forest_models(omf_list, omf_names, show_all = FALSE)
omfplot80
# omfplot80 <- omfplot80 +
#   labs(title = "Frequency of \nMine Damage (late-season)")+
#   scale_y_discrete(labels = c("FFP"))+
#   theme(legend.position = "none")
# omfplot80

#young leaves
ymdplot80 <- plot_forest_models(ymd_list, ymd_names, show_all = FALSE)
ymdplot80
ymdplot80 <- ymdplot80 +
  labs(title = "Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("MAR"))+
  theme(legend.position = "none")
ymdplot80

ymfplot80 <- plot_forest_models(ymf_list, ymf_names, show_all = FALSE)
ymfplot80
#all cross zero

mineplot80 <- plot_grid(ymdplot80,opamplot80,
                      nrow = 1,
                      ncol = 2)
mineplot80


# Gall plot----
opag_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
opag_list <- opag_list[!(grepl(".summary\\.csv$", opag_list) | grepl("_R2\\.csv$", opag_list))]
opag_names <- paste("Model", 1:length(opag_list))

ypag_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
ypag_list <- ypag_list[!(grepl(".summary\\.csv$", ypag_list) | grepl("_R2\\.csv$", ypag_list))]
ypag_names <- paste("Model", 1:length(ypag_list))

ogf_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE)
ogf_list <- ogf_list[!(grepl(".summary\\.csv$", ogf_list) | grepl("_R2\\.csv$", ogf_list))]
ogf_names <- paste("Model", 1:length(ogf_list))

ygf_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE)
ygf_list <- ygf_list[!(grepl(".summary\\.csv$", ygf_list) | grepl("_R2\\.csv$", ygf_list))]
ygf_names <- paste("Model", 1:length(ygf_list))

ogd_list <- list.files("Script_output/Brms/climate_models/l.lvl/old/year_doy/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE)
ogd_list <- ogd_list[!(grepl(".summary\\.csv$", ogd_list) | grepl("_R2\\.csv$", ogd_list))]
ogd_names <- paste("Model", 1:length(ogd_list))

ygd_list <- list.files("Script_output/Brms/climate_models/l.lvl/young/year_doy/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE)
ygd_list <- ygd_list[!(grepl(".summary\\.csv$", ygd_list) | grepl("_R2\\.csv$", ygd_list))]
ygd_names <- paste("Model", 1:length(ygd_list))

#plot--
opagplot <- plot_forest_models(opag_list, opag_names,show_all = TRUE)
opagplot
opagplot <- opagplot +
  labs(title = "Percent Area \nDamaged by Galls (late-season)")+
  scale_y_discrete(labels = c("SHM", "MSP", "MAP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
opagplot

#young
ypagplot <- plot_forest_models(ypag_list, ypag_names, show_all = TRUE)
ypagplot
ypagplot <- ypagplot +
  labs(title = "Percent Area Damaged by \nGalls (early-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "RH", "MAR","EXT", "Year","DOY"))+ 
  theme(legend.position = "none")
ypagplot

ogdplot <- plot_forest_models(ogd_list, ogd_names, show_all = TRUE)
ogdplot
ogdplot <- ogdplot +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c("SHM", "MSP", "MAP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ogdplot

ogfplot <- plot_forest_models(ogf_list, ogf_names, show_all = TRUE)
ogfplot
ogfplot <- ogfplot +
  labs(title = "Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("SHM", "MSP", "MAP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ogfplot

#young
ygdplot <- plot_forest_models(ygd_list, ygd_names, show_all = TRUE)
ygdplot
ygdplot <- ygdplot +
  labs(title = "Richness of \nGalls (early-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "MSP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ygdplot

ygfplot <- plot_forest_models(ygf_list, ygf_names,show_all = TRUE)
ygfplot
ygfplot <- ygfplot +
  labs(title = "Frequency \nof Galls (early-season)")+
  scale_y_discrete(labels = c("SHM", "AHM", "MSP", "MAR", "EXT", "Year","DOY"))+
  theme(legend.position = "none")
ygfplot


gallplot <- plot_grid(ypagplot,opagplot, 
                      ygfplot,ogfplot,
                      ygdplot,ogdplot,
                      nrow = 3,
                      ncol = 2)
gallplot

#80 Gall plot----
opagplot80 <- plot_forest_models(opag_list, opag_names,show_all = FALSE)
opagplot80
#all overlap 0

#young
ypagplot80 <- plot_forest_models(ypag_list, ypag_names,show_all = FALSE)
ypagplot80
# all overlap 0

ogdplot80 <- plot_forest_models(ogd_list, ogd_names, show_all = FALSE)
ogdplot80
ogdplot80 <- ogdplot80 +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c( "SHM", "MSP","MAP", "DOY"))+
  theme(legend.position = "none")
ogdplot80

ogfplot80 <- plot_forest_models(ogf_list, ogf_names,show_all = FALSE)
ogfplot80
ogfplot80 <- ogfplot80 +
  labs(title = "Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("SHM","MSP"))+
  theme(legend.position = "none")
ogfplot80

#young
ygdplot80 <- plot_forest_models(ygd_list, ygd_names,show_all = FALSE)
ygdplot80
#all cross zero

ygfplot80 <- plot_forest_models(ygf_list, ygf_names,show_all = FALSE)
ygfplot80
#all overlap 0


gallplot80 <- plot_grid(ogfplot80,ogdplot80,
                      nrow = 1,
                      ncol = 2)
gallplot80


#saving manuscript (80% CI) plots----
ggsave("figures/l.lvl/totalbrms80_plot.pdf", totalplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/chewbrms80_plot.pdf", chewplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/specbrms80_plot.pdf", specplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/minebrms80_plot.pdf", mineplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/gallbrms80_plot.pdf", gallplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/genbrms80_plot.pdf", genplot80, height = 15, width = 13, units = "in")

#saving SI plots----
ggsave("figures/l.lvl/SI/totalbrms_plot.pdf", totalplot, height = 15, width = 13, units = "in") 
ggsave("figures/l.lvl/SI/chewbrms_plot.pdf", chewplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/specbrms_plot.pdf", specplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/minebrms_plot.pdf", mineplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/gallbrms_plot.pdf", gallplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/genbrms_plot.pdf", genplot, height = 15, width = 13, units = "in")



#Previous versions of function----
# plot_forest_models <- function(csv_paths, model_names = NULL) {
#   # Read all CSV files
#   models_data <- map(csv_paths, read_csv)
#   
#   # If model names aren't provided, generate default names
#   if (is.null(model_names)) {
#     model_names <- paste("Model", seq_along(models_data))
#   }
#   
#   # Function to process a single model's data
#   process_model <- function(model_data, model_name) {
#     # Get parameter names
#     params <- colnames(model_data)
#     params <- params[str_detect(params, "^b_") & !str_detect(params, "^b_Intercept")]
#     
#     summary_stats <- model_data %>%
#       select(all_of(params)) %>%
#       pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
#       group_by(parameter) %>%
#       summarize(
#         mean = mean(value),
#         lower_95 = quantile(value, 0.025), #95% CI lower bound
#         upper_95 = quantile(value, 0.975), #95% CI upper bound
#         lower_80 = quantile(value, 0.1),  # 80% CI lower bound
#         upper_80 = quantile(value, 0.9),  # 80% CI upper bound
#         positive = sum(value > 0) / length(value),
#         negative = sum(value < 0) / length(value)
#       ) %>%
#       mutate(
#         parameter = gsub("b_", "", parameter),
#         model = model_name
#       )
#     
#     return(summary_stats)
#   }
#   
#   # Process all models
#   all_summaries <- map2_df(models_data, model_names, process_model)
#   
#   # Determine which models should have solid vs dashed lines based on path or pattern
#   # For the example provided, we can extract this from the model names or paths
#   
#   # Determine model type (year vs doy) from csv_paths
#   model_types <- rep(NA, length(csv_paths))
#   
#   for (i in seq_along(csv_paths)) {
#     if (grepl("/year/", csv_paths[i])) {
#       model_types[i] <- "year"
#     } else if (grepl("/doy/", csv_paths[i])) {
#       model_types[i] <- "doy"
#     } else {
#       # If we can't determine from path, use model number
#       model_num <- as.numeric(gsub("Model ", "", model_names[i]))
#       if (!is.na(model_num) && model_num <= 5) {
#         model_types[i] <- "year"
#       } else {
#         model_types[i] <- "doy"
#       }
#     }
#   }
#   
#   # Create a data frame with model names and their types
#   model_info <- data.frame(
#     model = model_names,
#     type = model_types,
#     label = ifelse(model_types == "year", "+ Year", "+ DOY")
#   )
#   
#   # Define model line types based on type
#   model_linetype_def <- setNames(
#     ifelse(model_types == "year", "solid", "dashed"),
#     model_names
#   )
#   
#   # Create model labels
#   model_labels <- setNames(
#     model_info$label,
#     model_info$model
#   )
#   
#   # Color definitions for parameters (climate variables)
#   parameter_color_def <- c(
#     "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
#     "scaleyear.x" = "#C1AE7C", 
#     "scalestartDayOfYear" = "#C1AE7C", 
#     "scaleTD" = "#E16036", 
#     "scaleMWMT" = "#E16036", 
#     "scaleEXT" = "#E16036",
#     "scaleEMT" = "#E16036", 
#     "scaleeFFP" = "#E16036", 
#     "scaleAHM" = "#6D4F92", 
#     "scaleMAP" = "#5EB1BF", 
#     "scaleMAR" = "#E16036",
#     "scaleMSP" = "#5EB1BF", 
#     "scaleSHM" = "#5EB1BF", 
#     "scaleCMI" = "#5EB1BF", 
#     "scaleRH" = "#5EB1BF", 
#     "scaleCMD" = "#5EB1BF",
#     "scaleDD5" = "#E16036", 
#     "scaleDD1040" = "#E16036", 
#     "scaleDD_18" = "#E16036", 
#     "scaleMAT" = "#E16036", 
#     "scalebFFP" = "#E16036",
#     "scaleEref" = "#5EB1BF", 
#     "scaleDD18" = "#E16036"
#   )
#   
#   # Create the forest plot
#   forest_plot <- all_summaries %>%
#     ggplot(aes(y = parameter, x = mean, color = parameter, linetype = model)) +
#     # 95% CI with thin line
#     geom_linerange(aes(xmin = lower_95, xmax = upper_95),
#                    position = position_dodge(width = 0.5),
#                    linewidth = 0.5) +
#     # 80% CI with thick line
#     geom_linerange(aes(xmin = lower_80, xmax = upper_80),
#                    position = position_dodge(width = 0.5),
#                    linewidth = 2.5) +
#     # Point estimate
#     geom_point(position = position_dodge(width = 0.5), size = 2) +
#     # Zero reference line
#     geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
#     theme_minimal() +
#     labs(x = "Estimate", y = "Parameter") +
#     # Colors by parameter, line types by model
#     scale_color_manual(values = parameter_color_def, name = "Parameter") +
#     scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
#     # Improved layout with legend at the right
#     theme(
#       legend.position = "right",
#       legend.box = "vertical",
#       panel.grid.minor = element_blank()
#     )
#   
#   return(forest_plot)
# }


#previous function. I have trust issues :) 
# plot_forest_models <- function(csv_paths, model_names = NULL) {
#   
#   # Read all CSV files
#   models_data <- map(csv_paths, read_csv)
#   
#   # If model names aren't provided, generate default names
#   if (is.null(model_names)) {
#     model_names <- paste("Model", seq_along(models_data))
#   }
#   
#   # Function to process a single model's data
#   process_model <- function(model_data, model_name) {
#     # Get parameter names
#     params <- colnames(model_data)
#     params <- params[str_detect(params, "^b_") & !str_detect(params, "^b_Intercept")]
#     
#     summary_stats <- model_data %>%
#       select(all_of(params)) %>%
#       pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
#       group_by(parameter) %>%
#       summarize(
#         mean = mean(value),
#         lower_95 = quantile(value, 0.025), #95% CI lower bound
#         upper_95 = quantile(value, 0.975), #95% CI upper bound
#         lower_80 = quantile(value, 0.1),  # 80% CI lower bound
#         upper_80 = quantile(value, 0.9),  # 80% CI upper bound
#         positive = sum(value > 0) / length(value),
#         negative = sum(value < 0) / length(value)
#       ) %>%
#       mutate(
#         parameter = gsub("b_", "", parameter),
#         model = model_name
#       )
#     
#     return(summary_stats)
#   }
#   
#   # Process all models
#   all_summaries <- map2_df(models_data, model_names, process_model)
#   
#   #color definitions: 
#   color_def <- c("scaleyear.x.scalestartDayOfYear" = "#C1AE7C","scaleyear.x" = "#C1AE7C", "scalestartDayOfYear" = "#C1AE7C", "scaleTD" = "#E16036", "scaleMWMT" = "#E16036", "scaleEXT" = "#E16036",
#                  "scaleEMT" = "#E16036", "scaleeFFP" = "#E16036", "scaleAHM" = "#6D4F92", "scaleMAP" = "#5EB1BF", "scaleMAR" = "#E16036", 
#                  "scaleMSP" = "#5EB1BF", "scaleSHM" = "#5EB1BF", "scaleCMI" = "#5EB1BF", "scaleRH" = "#5EB1BF", "scaleCMD" = "#5EB1BF", 
#                  "scaleDD5" = "#E16036", "scaleDD1040" = "#E16036", "scaleDD_18" = "#E16036", "scaleMAT" = "#E16036", "scalebFFP" = "#E16036", 
#                  "scaleEref" = "#5EB1BF", "scaleDD18" = "#E16036")
#   
#   # Create the forest plot
#   forest_plot <- all_summaries %>%
#     ggplot(aes(y = parameter, x = mean, color = parameter, linetype=model)) +
#     # 95% CI with thin line
#     geom_linerange(aes(xmin = lower_95, xmax = upper_95), 
#                    position = position_dodge(width = 0.5),
#                    linewidth = 0.5) +
#     # 80% CI with thick line
#     geom_linerange(aes(xmin = lower_80, xmax = upper_80), 
#                    position = position_dodge(width = 0.5),
#                    linewidth = 2.5) +
#     # Point estimate
#     geom_point(position = position_dodge(width = 0.5), size = 2) +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#   # Add positive values at upper bound
#   geom_text(aes(x = upper_80,
#                 label = sprintf("%.1f%%", positive*100)),  # Only positive values
#             position = position_dodge(width = 0.5),
#             hjust = -0.9,
#             size = 3) +
#   # Add negative values at lower bound
#   geom_text(aes(x = lower_80,
#                 label = sprintf("%.1f%%", negative*100)),  # Only negative values
#             position = position_dodge(width = 0.5),
#             hjust = 2.5,
#             size = 3) +
#     theme_minimal() +
#     labs(x = "Estimate", y = "Parameter") +
#     scale_color_manual(values = color_def)+
#     scale_linetype_manual(values = c("Model 1" = "solid", "Model 2" = "solid", "Model 3" = "solid", "Model 4" = "solid", "Model 5" = "solid",
#                                      "Model 6" = "solid", "Model 7" = "solid", "Model 8" = "solid", "Model 9" = "solid", "Model 10" = "solid"))
#     #scale_color_manual(values = wes_palette("Zissou1Continuous", n = length(unique(all_summaries$parameter)))); #Change this so that the colors are coded to temp vs. precip??
#     # theme_minimal() +
#     # labs(x = "Estimate", y = "Parameter") 
#     #coord_flip()
#     #theme(legend.position = "bottom")
#   
#   return(forest_plot)
#   #return(all_summaries)
# }

# plot80_forest_models <- function(csv_paths, model_names = NULL) {
#   
#   # Read all CSV files
#   models_data <- map(csv_paths, read_csv)
#   
#   # If model names aren't provided, generate default names
#   if (is.null(model_names)) {
#     model_names <- paste("Model", seq_along(models_data))
#   }
#   
#   # Function to process a single model's data
#   process_model <- function(model_data, model_name) {
#     # Get parameter names
#     params <- colnames(model_data)
#     params <- params[str_detect(params, "^b_") & !str_detect(params, "^b_Intercept")]
#     
#     summary_stats <- model_data %>%
#       select(all_of(params)) %>%
#       pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
#       group_by(parameter) %>%
#       summarize(
#         mean = mean(value),
#         lower_95 = quantile(value, 0.025), #95% CI lower bound
#         upper_95 = quantile(value, 0.975), #95% CI upper bound
#         lower_80 = quantile(value, 0.1),  # 80% CI lower bound
#         upper_80 = quantile(value, 0.9),  # 80% CI upper bound
#         positive = sum(value > 0) / length(value),
#         negative = sum(value < 0) / length(value)
#       ) %>%
#       mutate(
#         parameter = gsub("b_", "", parameter),
#         model = model_name,
#         # Add a flag to identify estimates that don't cross zero (80% CI)
#         significant_80 = (lower_80 > 0) | (upper_80 < 0)
#       )
#     
#     return(summary_stats)
#   }
#   
#   # Process all models
#   all_summaries <- map2_df(models_data, model_names, process_model)
#   
#   #color definitions: 
#   color_def <- c("scaleyear.x.scalestartDayOfYear" = "#C1AE7C","scaleyear.x" = "#C1AE7C", "scalestartDayOfYear" = "#C1AE7C", "scaleTD" = "#E16036", "scaleMWMT" = "#E16036", "scaleEXT" = "#E16036",
#                  "scaleEMT" = "#E16036", "scaleeFFP" = "#E16036", "scaleAHM" = "#6D4F92", "scaleMAP" = "#5EB1BF", "scaleMAR" = "#E16036", 
#                  "scaleMSP" = "#5EB1BF", "scaleSHM" = "#5EB1BF", "scaleCMI" = "#5EB1BF", "scaleRH" = "#5EB1BF", "scaleCMD" = "#5EB1BF", 
#                  "scaleDD5" = "#E16036", "scaleDD1040" = "#E16036", "scaleDD_18" = "#E16036", "scaleMAT" = "#E16036", "scalebFFP" = "#E16036", 
#                  "scaleEref" = "#5EB1BF", "scaleDD18" = "#E16036")
#   
#   # Filter data to only keep estimates that don't cross 0 (80% CI)
#   significant_summaries <- all_summaries %>%
#     filter(significant_80)
#   
#   # Create the forest plot with only significant estimates
#   forest_plot <- significant_summaries %>%
#     ggplot(aes(y = parameter, x = mean, color = parameter, linetype=model)) +
#     # 95% CI with thin line
#     geom_linerange(aes(xmin = lower_95, xmax = upper_95), 
#                    position = position_dodge(width = 0.5),
#                    linewidth = 0.5) +
#     # 80% CI with thick line
#     geom_linerange(aes(xmin = lower_80, xmax = upper_80), 
#                    position = position_dodge(width = 0.5),
#                    linewidth = 2.5) +
#     # Point estimate
#     geom_point(position = position_dodge(width = 0.5), size = 2) +
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     # Add positive values at upper bound
#     geom_text(aes(x = upper_80,
#                   label = sprintf("%.1f%%", positive*100)),  # Only positive values
#               position = position_dodge(width = 0.5),
#               hjust = -0.9,
#               size = 3) +
#     # Add negative values at lower bound
#     geom_text(aes(x = lower_80,
#                   label = sprintf("%.1f%%", negative*100)),  # Only negative values
#               position = position_dodge(width = 0.5),
#               hjust = 2.5,
#               size = 3) +
#     theme_minimal() +
#     labs(x = "Estimate", y = "Parameter") +
#     scale_color_manual(values = color_def) +
#     scale_linetype_manual(values = c("Model 1" = "solid", "Model 2" = "solid", "Model 3" = "solid", "Model 4" = "solid", "Model 5" = "solid",
#                                      "Model 6" = "solid", "Model 7" = "solid", "Model 8" = "solid", "Model 9" = "solid", "Model 10" = "solid"))
#   
#   return(forest_plot)
# }



