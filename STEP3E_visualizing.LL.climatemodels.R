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

#Functions
#Forest plot function----
plot_forest_models <- function(csv_paths, model_names = NULL) {
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
        lower_95 = quantile(value, 0.025), #95% CI lower bound
        upper_95 = quantile(value, 0.975), #95% CI upper bound
        lower_80 = quantile(value, 0.1),  # 80% CI lower bound
        upper_80 = quantile(value, 0.9),  # 80% CI upper bound
        positive = sum(value > 0) / length(value),
        negative = sum(value < 0) / length(value)
      ) %>%
      mutate(
        parameter = gsub("b_", "", parameter),
        model = model_name
      )
    
    return(summary_stats)
  }
  
  # Process all models
  all_summaries <- map2_df(models_data, model_names, process_model)
  
  # Determine model type (year vs doy) from csv_paths
  model_types <- rep(NA, length(csv_paths))
  
  for (i in seq_along(csv_paths)) {
    if (grepl("_year_model\\.csv$", csv_paths[i]) || grepl("/year/", csv_paths[i])) {
      model_types[i] <- "year"
    } else if (grepl("_doy_model\\.csv$", csv_paths[i]) || grepl("/doy/", csv_paths[i])) {
      model_types[i] <- "doy"
    } else {
      # If we can't determine from path, use model number
      model_num <- as.numeric(gsub("Model ", "", model_names[i]))
      if (!is.na(model_num) && model_num <= 5) {
        model_types[i] <- "year"
      } else {
        model_types[i] <- "doy"
      }
    }
  }
  
  # Create a data frame with model names and their types
  model_info <- data.frame(
    model = model_names,
    type = model_types,
    label = ifelse(model_types == "year", "+ Year", "+ DOY")
  )
  
  # Define model line types based on type
  model_linetype_def <- setNames(
    ifelse(model_types == "year", "solid", "dashed"),
    model_names
  )
  
  # Create model labels
  model_labels <- setNames(
    model_info$label,
    model_info$model
  )
  
  # Color definitions for parameters (climate variables)
  parameter_color_def <- c(
    "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scalestartDayOfYear" = "#C1AE7C", 
    "scaleTD" = "#E16036", 
    "scaleMWMT" = "#E16036", 
    "scaleEXT" = "#E16036",
    "scaleEMT" = "#E16036", 
    "scaleeFFP" = "#E16036", 
    "scaleAHM" = "#6D4F92", 
    "scaleMAP" = "#5EB1BF", 
    "scaleMAR" = "#E16036",
    "scaleMSP" = "#5EB1BF", 
    "scaleSHM" = "#5EB1BF", 
    "scaleCMI" = "#5EB1BF", 
    "scaleRH" = "#5EB1BF", 
    "scaleCMD" = "#5EB1BF",
    "scaleDD5" = "#E16036", 
    "scaleDD1040" = "#E16036", 
    "scaleDD_18" = "#E16036", 
    "scaleMAT" = "#E16036", 
    "scalebFFP" = "#E16036",
    "scaleEref" = "#5EB1BF", 
    "scaleDD18" = "#E16036"
  )
  
  # Create the forest plot
  forest_plot <- all_summaries %>%
    ggplot(aes(y = parameter, x = mean, color = parameter, linetype = model)) +
    # 95% CI with thin line
    geom_linerange(aes(xmin = lower_95, xmax = upper_95),
                   position = position_dodge(width = 0.5),
                   linewidth = 0.5) +
    # 80% CI with thick line
    geom_linerange(aes(xmin = lower_80, xmax = upper_80),
                   position = position_dodge(width = 0.5),
                   linewidth = 2.5) +
    # Point estimate
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    # Zero reference line
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
    theme_minimal() +
    labs(x = "Estimate", y = "Parameter") +
    # Colors by parameter, line types by model
    scale_color_manual(values = parameter_color_def, name = "Parameter") +
    scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
    # Improved layout with legend at the right
    theme(
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  return(forest_plot)
}

#80% Forest plot function----
plot80_forest_models <- function(csv_paths, model_names = NULL) {
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
        lower_95 = quantile(value, 0.025), #95% CI lower bound
        upper_95 = quantile(value, 0.975), #95% CI upper bound
        lower_80 = quantile(value, 0.1),  # 80% CI lower bound
        upper_80 = quantile(value, 0.9),  # 80% CI upper bound
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
  
  # Determine model types (year vs doy) from file paths or names
  model_types <- rep(NA, length(csv_paths))
  
  for (i in seq_along(csv_paths)) {
    if (grepl("_year_model\\.csv$", csv_paths[i])) {
      model_types[i] <- "year"
    } else if (grepl("_doy_model\\.csv$", csv_paths[i])) {
      model_types[i] <- "doy"
    } else if (grepl("/year/", csv_paths[i])) {
      model_types[i] <- "year"
    } else if (grepl("/doy/", csv_paths[i])) {
      model_types[i] <- "doy"
    } else {
      # If we can't determine from path pattern, use model number as fallback
      model_num <- as.numeric(gsub("Model ", "", model_names[i]))
      if (!is.na(model_num) && model_num <= 5) {
        model_types[i] <- "year"
      } else {
        model_types[i] <- "doy"
      }
    }
  }
  
  # Create a data frame with model names and their types
  model_info <- data.frame(
    model = model_names,
    type = model_types,
    label = ifelse(model_types == "year", "+ Year", "+ DOY")
  )
  
  # Define model line types based on type - IMPORTANT: Keep this and remove the hard-coded version
  model_linetype_def <- setNames(
    ifelse(model_types == "year", "solid", "dashed"),
    model_names
  )
  
  # Create model labels
  model_labels <- setNames(
    model_info$label,
    model_info$model
  )
  
  # Color definitions for parameters (climate variables)
  parameter_color_def <- c(
    "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scalestartDayOfYear" = "#C1AE7C", 
    "scaleTD" = "#E16036", 
    "scaleMWMT" = "#E16036", 
    "scaleEXT" = "#E16036",
    "scaleEMT" = "#E16036", 
    "scaleeFFP" = "#E16036", 
    "scaleAHM" = "#6D4F92", 
    "scaleMAP" = "#5EB1BF", 
    "scaleMAR" = "#E16036", 
    "scaleMSP" = "#5EB1BF", 
    "scaleSHM" = "#5EB1BF", 
    "scaleCMI" = "#5EB1BF", 
    "scaleRH" = "#5EB1BF", 
    "scaleCMD" = "#5EB1BF", 
    "scaleDD5" = "#E16036", 
    "scaleDD1040" = "#E16036", 
    "scaleDD_18" = "#E16036", 
    "scaleMAT" = "#E16036", 
    "scalebFFP" = "#E16036", 
    "scaleEref" = "#5EB1BF", 
    "scaleDD18" = "#E16036"
  )
  
  # REMOVE the hard-coded model_linetype_def that was here before
  
  # Filter data to only keep estimates that don't cross 0 (80% CI)
  significant_summaries <- all_summaries %>%
    filter(significant_80)
  
  # Create the forest plot with only significant estimates
  forest_plot <- significant_summaries %>%
    ggplot(aes(y = parameter, x = mean, color = parameter, linetype = model)) +
    # 95% CI with thin line
    geom_linerange(aes(xmin = lower_95, xmax = upper_95),
                   position = position_dodge(width = 0.5),
                   linewidth = 0.5) +
    # 80% CI with thick line
    geom_linerange(aes(xmin = lower_80, xmax = upper_80),
                   position = position_dodge(width = 0.5),
                   linewidth = 2.5) +
    # Point estimate
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    labs(x = "Estimate", y = "Parameter") +
    # Colors by parameter, line types by model
    scale_color_manual(values = parameter_color_def, name = "Parameter") +
    scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
    # Improved layout with legend at the right
    theme(
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  return(forest_plot)
}

######################CLIMATE MODELS##########################################
#Total plot----
otf_list <- c(
  list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                        pattern = "^leaftotal.freq", 
                        full.names = TRUE),
  list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
           pattern = "^leaftotal.freq", 
           full.names = TRUE)
  )
otf_list <- otf_list[!(grepl(".summary\\.csv$", otf_list) | grepl("_R2\\.csv$", otf_list))]
otf_names <- paste("Model", 1:10)
otf_names <- factor(otf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ytf_list <- c(
  list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
             pattern = "^leaftotal.freq", 
             full.names = TRUE),
  list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
             pattern = "^leaftotal.freq", 
             full.names = TRUE)
)
ytf_list <- ytf_list[!(grepl(".summary\\.csv$", ytf_list) | grepl("_R2\\.csv$", ytf_list))]
ytf_names <- paste("Model", 1:10)
ytf_names <- factor(ytf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

otd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leaftotal.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leaftotal.rich", 
                         full.names = TRUE)
)
otd_list <- otd_list[!(grepl(".summary\\.csv$", otd_list) | grepl("_R2\\.csv$", otd_list))]
otd_names <- paste("Model", 1:10)
otd_names <- factor(otd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ytd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leaftotal.rich", 
                       full.names = TRUE), 
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leaftotal.rich", 
                         full.names = TRUE)
)
ytd_list <- ytd_list[!(grepl(".summary\\.csv$", ytd_list) | grepl("_R2\\.csv$", ytd_list))]
ytd_names <- paste("Model", 1:10)
ytd_names <- factor(ytd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

otdplot <- plot_forest_models(otd_list, otd_names)
otdplot
otdplot <- otdplot +
  labs(title = "Richness of \nTotal Damage (late-season)")+
  scale_y_discrete(labels = c("eFFP", "EMT", "EXT", "MWMT", "DOY","TD", "Year"))+
  theme(legend.position = "none")
otdplot

otfplot <- plot_forest_models(otf_list, otf_names)
otfplot
otfplot <- otfplot +
  labs(title = "Frequency of \nTotal Damage (late-season)")+
  scale_y_discrete (labels = c("eFFP", "EMT", "EXT", "MWMT", "DOY","TD", "Year"))+
  theme(legend.position = "none")
otfplot

#young leaves
ytdplot <- plot_forest_models(ytd_list, ytd_names)
ytdplot
#pulling legend
modellegend <- get_legend(ytdplot)
ytdplot <- ytdplot +
  labs(title = "Richness of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ytdplot

ytfplot <- plot_forest_models(ytf_list, ytf_names)
ytfplot
ytfplot <- ytfplot +
  labs(title = "Frequency of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ytfplot

totalplot <- plot_grid(ytfplot,otfplot,
                       ytdplot,otdplot,
                       nrow = 2,
                       ncol = 2)
totalplot
totalplot <- plot_grid(totalplot,
                       modellegend,
                       nrow = 1,
                       ncol = 2, 
                       rel_widths = c(1, .10))
totalplot

#80 Total plots----
otdplot80 <- plot80_forest_models(otd_list, otd_names)
otdplot80
otdplot80 <- otdplot80 +
  labs(title = "Richness of \nTotal Damage (late-season)")+
  scale_y_discrete(labels = c("eFFP", "EMT", "MWMT", "DOY","TD", "Year"))+
  theme(legend.position = "none")
otdplot80

otfplot80 <- plot80_forest_models(otf_list, otf_names)
otfplot80
otfplot80 <- otfplot80 +
  labs(title = "Frequency of \nTotal Damage (late-season)")+
  scale_y_discrete (labels = c("eFFP", "EMT", "TD", "Year"))+
  theme(legend.position = "none")
otfplot80

#young leaves
ytdplot80 <- plot80_forest_models(ytd_list, ytd_names)
ytdplot80
ytdplot80 <- ytdplot80 +
  labs(title = "Richness of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "Year"))+
  theme(legend.position = "none")
ytdplot80

ytfplot80 <- plot80_forest_models(ytf_list, ytf_names)
ytfplot80
ytfplot80 <- ytfplot80 +
  labs(title = "Frequency of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("AHM"))+
  theme(legend.position = "none")
ytfplot80

totalplot80 <- plot_grid(ytfplot80,otfplot80,
                       ytdplot80,otdplot80,
                       nrow = 2,
                       ncol = 2)
totalplot80
totalplot80 <- plot_grid(totalplot80,
                       modellegend,
                       nrow = 1,
                       ncol = 2, 
                       rel_widths = c(1, .10))
totalplot80


#Chewing plot----
# pac_list <- list.files("Script_output/Brms/", 
#                        pattern = "^leafperc_area_chew", 
#                        full.names = TRUE)
# pac_list <- pac_list[!grepl("summary\\.csv$", pac_list)]
# pac_names <- paste("Model", 1:5)

opac_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                          pattern = "^leafperc_area_chew", 
                          full.names = TRUE)
)
opac_list <- opac_list[!(grepl(".summary\\.csv$", opac_list) | grepl("_R2\\.csv$", opac_list))]
opac_names <- paste("Model", 1:10)
opac_names <- factor(opac_names, 
                     levels = paste("Model", 1:10),
                     ordered = TRUE)

ypac_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                          pattern = "^leafperc_area_chew", 
                          full.names = TRUE)
               )
ypac_list <- ypac_list[!(grepl(".summary\\.csv$", ypac_list) | grepl("_R2\\.csv$", ypac_list))]
ypac_names <- paste("Model", 1:10)
ypac_names <- factor(ypac_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

#chewing herbs "
ochew_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE),
                list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                           pattern = "^chewvalue", 
                           full.names = TRUE)
)
ochew_list <- ochew_list[!(grepl(".summary\\.csv$", ochew_list) | grepl("_R2\\.csv$", ochew_list))]
ochew_names <- paste("Model", 1:10)
ochew_names <- factor(ochew_names, 
                      levels = paste("Model", 1:10),
                      ordered = TRUE)

ychew_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE),
                list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                           pattern = "^chewvalue", 
                           full.names = TRUE)
)
ychew_list <- ychew_list[!(grepl(".summary\\.csv$", ychew_list) | grepl("_R2\\.csv$", ychew_list))]
ychew_names <- paste("Model", 1:10)
ychew_names <- factor(ychew_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

#pac plots
# pacplot <- plot_forest_models(pac_list, pac_names)
# pacplot
# pacplot <- pacplot +
#   labs(title = "Percent Leaf Area Chewed") +
#   #scale_y_discrete(labels=c("CMI", "TD", "Leaf Age", "MAP", "MSP","Year"))+ #we are replacing the y-axis (cord flip in function)
#   theme(legend.position = "none") 
# pacplot

opacplot <- plot_forest_models(opac_list, opac_names)
opacplot
opacplot <- opacplot +
  labs(title = "Percent Leaf \nArea Chewed (late-season)") +
  scale_y_discrete(labels=c("AHM", "EXT", "MAP", "MSP", "SHM","DOY","Year"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
opacplot

ypacplot <- plot_forest_models(ypac_list, ypac_names)
ypacplot
ypacplot <- ypacplot +
  labs(title = "Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("AHM", "CMI", "MAP", "MAR","MSP", "DOY","Year"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot

#chewing
ochewplot <- plot_forest_models(ochew_list, ochew_names)
ochewplot
ochewplot <- ochewplot +
  labs(title = "Chewing \nDamage (proportion; late-season)") +
  scale_y_discrete(labels=c("EXT", "MAP", "MWMT", "RH", "DOY","TD","Year"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot

ychewplot <- plot_forest_models(ychew_list, ychew_names)
ychewplot
ychewplot <- ychewplot +
  labs(title = "Chewing Damage \n(proportion; early-season)") +
  scale_y_discrete(labels=c("AHM", "MAP", "MAR", "MSP","SHM", "DOY","Year"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ychewplot

#chew plot
chewplot <- plot_grid(ypacplot,opacplot,
                     ychewplot, ochewplot,
                      nrow = 2,
                      ncol = 2)
chewplot
chewplot <- plot_grid(chewplot,
                      modellegend,
                      nrow = 1,
                      ncol = 2, 
                      rel_widths = c(1, .10))
chewplot

#80 Chew plots----
opacplot80 <- plot80_forest_models(opac_list, opac_names)
opacplot80
#All variables cross 0

ypacplot80 <- plot80_forest_models(ypac_list, ypac_names)
ypacplot80
ypacplot80 <- ypacplot80 +
  labs(title = "Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("AHM", "MAP", "MSP"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot80

#chewing
ochewplot80 <- plot80_forest_models(ochew_list, ochew_names)
ochewplot80
#all cross 0

ychewplot80 <- plot80_forest_models(ychew_list, ychew_names)
ychewplot80
#all cross 0

#Specialized plot----
osf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafspec.freq", 
                         full.names = TRUE)
)
osf_list <- osf_list[!(grepl(".summary\\.csv$", osf_list) | grepl("_R2\\.csv$", osf_list))]
osf_names <- paste("Model", 1:10)
osf_names <- factor(osf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ysf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafspec.freq", 
                         full.names = TRUE)
)
ysf_list <- ysf_list[!(grepl(".summary\\.csv$", ysf_list) | grepl("_R2\\.csv$", ysf_list))]
ysf_names <- paste("Model", 1:10)
ysf_names <- factor(ysf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

osd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafspec.rich", 
                         full.names = TRUE)
)

osd_list <- osd_list[!(grepl(".summary\\.csv$", osd_list) | grepl("_R2\\.csv$", osd_list))]
osd_names <- paste("Model", 1:10)
osd_names <- factor(osd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ysd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafspec.rich", 
                         full.names = TRUE)
)
ysd_list <- ysd_list[!(grepl(".summary\\.csv$", ysd_list) | grepl("_R2\\.csv$", ysd_list))]
ysd_names <- paste("Model", 1:10)
ysd_names <- factor(ysd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

osdplot <- plot_forest_models(osd_list, osd_names)
osdplot
osdplot <- osdplot +
  labs(title = "Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("EMT", "EXT", "MWMT", "RH", "DOY", "TD", "Year"))+
  theme(legend.position = "none")
osdplot

osfplot <- plot_forest_models(osf_list, osf_names)
osfplot
osfplot <- osfplot +
  labs(title = "Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("AHM", "EXT", "MAP", "MWMT", "DOY", "TD","Year"))+
  theme(legend.position = "none")
osfplot

#young leaves
ysdplot <- plot_forest_models(ysd_list, ysd_names)
ysdplot
ysdplot <- ysdplot +
  labs(title = "Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ysdplot

ysfplot <- plot_forest_models(ysf_list, ysf_names)
ysfplot
ysfplot <- ysfplot +
  labs(title = "Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ysfplot

specplot <- plot_grid(ysfplot,osfplot, 
                      ysdplot,osdplot, 
                      nrow = 2,
                      ncol = 2)
specplot <- plot_grid(specplot, modellegend,
                      nrow = 1, ncol = 2,
                      rel_widths = c(1,.10))
specplot

#80 Specialized plots----
osdplot80 <- plot80_forest_models(osd_list, osd_names)
osdplot80
osdplot80 <- osdplot80 +
  labs(title = "Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("EMT", "MWMT", "TD", "Year"))+
  theme(legend.position = "none")
osdplot80

osfplot80 <- plot80_forest_models(osf_list, osf_names)
osfplot80
osfplot80 <- osfplot80 +
  labs(title = "Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("TD","Year"))+
  theme(legend.position = "none")
osfplot80

#young leaves
ysdplot80 <- plot80_forest_models(ysd_list, ysd_names)
ysdplot80
ysdplot80 <- ysdplot80 +
  labs(title = "Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("AHM", "MAP", "Year"))+
  theme(legend.position = "none")
ysdplot80

ysfplot80 <- plot80_forest_models(ysf_list, ysf_names)
ysfplot80
ysfplot80 <- ysfplot80 +
  labs(title = "Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("AHM"))+
  theme(legend.position = "none")
ysfplot80

specplot80 <- plot_grid(ysfplot80,osfplot80, 
                      ysdplot80,osdplot80, 
                      nrow = 2,
                      ncol = 2)
specplot80 <- plot_grid(specplot80, modellegend,
                      nrow = 1, ncol = 2,
                      rel_widths = c(1,.10))
specplot80

#Mine plot----
opam_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                          pattern = "^leafperc_area_mine", 
                          full.names = TRUE)
)
opam_list <- opam_list[!(grepl(".summary\\.csv$", opam_list) | grepl("_R2\\.csv$", opam_list))]
opam_names <- paste("Model", 1:10)
opam_names <- factor(opam_names, 
                     levels = paste("Model", 1:10),
                     ordered = TRUE)

ypam_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                          pattern = "^leafperc_area_mine", 
                          full.names = TRUE)
)
ypam_list <- ypam_list[!(grepl(".summary\\.csv$", ypam_list) | grepl("_R2\\.csv$", ypam_list))]
ypam_names <- paste("Model", 1:9)
ypam_names <- factor(ypam_names, 
                    levels = paste("Model", 1:9),
                    ordered = TRUE)

omf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafmine.freq", 
                         full.names = TRUE)
)
omf_list <- omf_list[!(grepl(".summary\\.csv$", omf_list) | grepl("_R2\\.csv$", omf_list))]
omf_names <- paste("Model", 1:10)
omf_names <- factor(omf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ymf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafmine.freq", 
                         full.names = TRUE)
)
ymf_list <- ymf_list[!(grepl(".summary\\.csv$", ymf_list) | grepl("_R2\\.csv$", ymf_list))]
ymf_names <- paste("Model", 1:10)
ymf_names <- factor(ymf_names, 
                     levels = paste("Model", 1:10),
                     ordered = TRUE)

omd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafmine.rich", 
                         full.names = TRUE)
)
omd_list <- omd_list[!(grepl(".summary\\.csv$", omd_list) | grepl("_R2\\.csv$", omd_list))]
omd_names <- paste("Model", 1:10)
omd_names <- factor(omd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ymd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafmine.rich", 
                         full.names = TRUE)
)
ymd_list <- ymd_list[!(grepl(".summary\\.csv$", ymd_list) | grepl("_R2\\.csv$", ymd_list))]
ymd_names <- paste("Model", 1:10)
ymd_names <- factor(ymd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

# pamplot <- plot_forest_models(pam_list, pam_names)
# pamplot
# pamplot <- pamplot +
#   labs(title = "Percent Leaf Area \nMined") +
#   scale_y_discrete(labels = c("CMT", "Leaf Age", "DD18", "MWMT", "RH", "TD", "Year"))+
#   theme(legend.position = "none")
# pamplot

opamplot <- plot_forest_models(opam_list, opam_names)
opamplot
opamplot <- opamplot +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("CMD", "eFFP", "EMT", "MWMT", "DOY","TD", "Year"))+
  theme(legend.position = "none")
opamplot

ypamplot <- plot_forest_models(ypam_list, ypam_names)
ypamplot
ypamplot <- ypamplot +
  labs(title = "Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+ #MAR likely didn't converge
  theme(legend.position = "none")
ypamplot

omdplot <- plot_forest_models(omd_list, omd_names)
omdplot
omdplot <- omdplot +
  labs(title = "Richness of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DD_18", "DD5", "MAT", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
omdplot

omfplot <- plot_forest_models(omf_list, omf_names)
omfplot
omfplot <- omfplot +
  labs(title = "Frequency of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("bFFP", "DD_18", "DD1040", "DD5", "MAT", "DOY","Year"))+
  theme(legend.position = "none")
omfplot

#young leaves
ymdplot <- plot_forest_models(ymd_list, ymd_names)
ymdplot
ymdplot <- ymdplot +
  labs(title = "Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ymdplot

ymfplot <- plot_forest_models(ymf_list, ymf_names)
ymfplot
ymfplot <- ymfplot +
  labs(title = "Frequency of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ymfplot

mineplot <- plot_grid(ypamplot,opamplot, 
                      ymfplot,omfplot, 
                      ymdplot,omdplot, 
                      nrow = 3,
                      ncol = 2)
mineplot
mineplot <- plot_grid(mineplot,
                      modellegend,
                      nrow = 1,
                      ncol = 2, 
                      rel_widths = c(1,.10))
mineplot

#80 Mine plot----
opamplot80 <- plot80_forest_models(opam_list, opam_names)
opamplot80
opamplot80 <- opamplot80 +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("CMD", "MWMT", "DOY","TD"))+
  theme(legend.position = "none")
opamplot80

ypamplot80 <- plot80_forest_models(ypam_list, ypam_names)
ypamplot80
ypamplot80 <- ypamplot80 +
  labs(title = "Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("AHM", "MAP", "MSP", "SHM"))+ 
  theme(legend.position = "none")
ypamplot80

omdplot80 <- plot80_forest_models(omd_list, omd_names)
omdplot80
omdplot80 <- omdplot80 +
  labs(title = "Richness of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DD_18", "DD5", "MAT", "SHM"))+
  theme(legend.position = "none")
omdplot80

omfplot80 <- plot80_forest_models(omf_list, omf_names)
omfplot80
omfplot80 <- omfplot80 +
  labs(title = "Frequency of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("bFFP", "DD_18", "DD1040", "DD5", "MAT"))+
  theme(legend.position = "none")
omfplot80

#young leaves
ymdplot80 <- plot80_forest_models(ymd_list, ymd_names)
ymdplot80
#all cross 0

ymfplot80 <- plot80_forest_models(ymf_list, ymf_names)
ymfplot80
#all cross zero

mineplot80 <- plot_grid(ypamplot80,opamplot80, 
                       ymfplot80,omfplot80,
                       ymdplot80,omdplot80,
                      nrow = 3,
                      ncol = 2)
mineplot80
mineplot80 <- plot_grid(mineplot80,
                      modellegend,
                      nrow = 1,
                      ncol = 2, 
                      rel_widths = c(1,.10))
mineplot80


# Gall plot----
opag_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                          pattern = "^leafperc_area_gall", 
                          full.names = TRUE)
)
opag_list <- opag_list[!(grepl(".summary\\.csv$", opag_list) | grepl("_R2\\.csv$", opag_list))]
opag_names <- paste("Model", 1:10)
opag_names <- factor(opag_names, 
                     levels = paste("Model", 1:10),
                     ordered = TRUE)

ypag_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE),
               list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                          pattern = "^leafperc_area_gall", 
                          full.names = TRUE)
)
ypag_list <- ypag_list[!(grepl(".summary\\.csv$", ypag_list) | grepl("_R2\\.csv$", ypag_list))]
ypag_names <- paste("Model", 1:8)
ypag_names <- factor(ypag_names, 
                    levels = paste("Model", 1:8),
                    ordered = TRUE)

ogf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafgall.freq", 
                         full.names = TRUE)
)
ogf_list <- ogf_list[!(grepl(".summary\\.csv$", ogf_list) | grepl("_R2\\.csv$", ogf_list))]
ogf_names <- paste("Model", 1:10)
ogf_names <- factor(ogf_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ygf_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafgall.freq", 
                         full.names = TRUE)
)
ygf_list <- ygf_list[!(grepl(".summary\\.csv$", ygf_list) | grepl("_R2\\.csv$", ygf_list))]
ygf_names <- paste("Model", 1:10)
ygf_names <- factor(ygf_names, 
                     levels = paste("Model", 1:10),
                     ordered = TRUE)

ogd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/old/year/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/old/doy/", 
                         pattern = "^leafgall.rich", 
                         full.names = TRUE)
)
ogd_list <- ogd_list[!(grepl(".summary\\.csv$", ogd_list) | grepl("_R2\\.csv$", ogd_list))]
ogd_names <- paste("Model", 1:10)
ogd_names <- factor(ogd_names, 
                    levels = paste("Model", 1:10),
                    ordered = TRUE)

ygd_list <- c(list.files("Script_output/Brms/climate_models/l.lvl/young/year/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE),
              list.files("Script_output/Brms/climate_models/l.lvl/young/doy/", 
                         pattern = "^leafgall.rich", 
                         full.names = TRUE)
)
ygd_list <- ygd_list[!(grepl(".summary\\.csv$", ygd_list) | grepl("_R2\\.csv$", ygd_list))]
ygd_names <- paste("Model", 1:6)
ygd_names <- factor(ygd_names, 
                    levels = paste("Model", 1:6),
                    ordered = TRUE)

#plot--
# pagplot <- plot_forest_models(pag_list, pag_names)
# pagplot
# pagplot <- pagplot +
#   labs(title = "Percent Area Damaged by Galls")+
#   scale_y_discrete(labels = c("Leaf Age", "AHM", "CMI", "MAP", "MSP", "SHM", "Year"))+
#   theme(legend.position = "none")
# pagplot

opagplot <- plot_forest_models(opag_list, opag_names)
opagplot
opagplot <- opagplot +
  labs(title = "Percent Area \nDamaged by Galls (late-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MSP", "SHM", "DOY", "TD","Year"))+
  theme(legend.position = "none")
opagplot

#young
ypagplot <- plot_forest_models(ypag_list, ypag_names)
ypagplot
ypagplot <- ypagplot +
  labs(title = "Percent Area Damaged by \nGalls (early-season)")+
  scale_y_discrete(labels = c("AHM", "CMD", "DD18", "Eref", "MAR", "DOY","Year"))+ #AHM likely didn't converge
  theme(legend.position = "none")
ypagplot

ogdplot <- plot_forest_models(ogd_list, ogd_names)
ogdplot
ogdplot <- ogdplot +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c("EMT", "EXT", "MAP", "MWMT", "DOY", "TD","Year"))+
  theme(legend.position = "none")
ogdplot

ogfplot <- plot_forest_models(ogf_list, ogf_names)
ogfplot
ogfplot <- ogfplot +
  labs(title = "Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("EXT", "MAP", "MSP", "MWMT", "DOY", "TD","Year"))+
  theme(legend.position = "none")
ogfplot

#young
ygdplot <- plot_forest_models(ygd_list, ygd_names)
ygdplot
ygdplot <- ygdplot +
  labs(title = "Richness of \nGalls (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ygdplot

ygfplot <- plot_forest_models(ygf_list, ygf_names)
ygfplot
ygfplot <- ygfplot +
  labs(title = "Frequency \nof Galls (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "MAR", "MSP", "SHM", "DOY","Year"))+
  theme(legend.position = "none")
ygfplot


gallplot <- plot_grid(ypagplot,opagplot, 
                      ygfplot,ogfplot,
                      ygdplot,ogdplot,
                      nrow = 3,
                      ncol = 2)
gallplot <- plot_grid(gallplot,modellegend,
                      nrow = 1,
                      ncol = 2,
                      rel_widths = c(1,.10))
gallplot

#80 Gall plot----
opagplot80 <- plot80_forest_models(opag_list, opag_names)
opagplot80
#all overlap 0

#young
ypagplot80 <- plot80_forest_models(ypag_list, ypag_names)
ypagplot80
ypagplot80 <- ypagplot80 +
  labs(title = "Percent Area Damaged by \nGalls (early-season)")+
  scale_y_discrete(labels = c("AHM", "CMD", "Eref", "MAR", "DOY","Year"))+ #AHM likely didn't converge
  theme(legend.position = "none")
ypagplot80

ogdplot80 <- plot80_forest_models(ogd_list, ogd_names)
ogdplot80
ogdplot80 <- ogdplot80 +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c( "EXT", "MAP", "MWMT", "DOY", "TD","Year"))+
  theme(legend.position = "none")
ogdplot80

ogfplot80 <- plot80_forest_models(ogf_list, ogf_names)
ogfplot80
ogfplot80 <- ogfplot80 +
  labs(title = "Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("Year"))+
  theme(legend.position = "none")
ogfplot80

#young
ygdplot80 <- plot80_forest_models(ygd_list, ygd_names)
ygdplot80
ygdplot80 <- ygdplot80 +
  labs(title = "Richness of \nGalls (early-season)")+
  scale_y_discrete(labels = c("AHM", "MAP", "Year"))+
  theme(legend.position = "none")
ygdplot80

ygfplot80 <- plot80_forest_models(ygf_list, ygf_names)
ygfplot80
#all overlap 0


gallplot80 <- plot_grid(ypagplot80,NULL, 
                      NULL,ogfplot80,
                      ygdplot80,ogdplot80,
                      nrow = 3,
                      ncol = 2)
gallplot80 <- plot_grid(gallplot80,modellegend,
                      nrow = 1,
                      ncol = 2,
                      rel_widths = c(1,.10))
gallplot80


#saving manuscript (80% CI) plots----
ggsave("figures/l.lvl/totalbrms80_plot.pdf", totalplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/chewbrms80_plot.pdf", ypacplot80)
ggsave("figures/l.lvl/specbrms80_plot.pdf", specplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/minebrms80_plot.pdf", mineplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/gallbrms80_plot.pdf", gallplot80, height = 15, width = 13, units = "in")

#saving SI plots----
ggsave("figures/l.lvl/SI/totalbrms_plot.pdf", totalplot, height = 15, width = 13, units = "in") #convert back to .png when the are "final" and ready for a manuscript
ggsave("figures/l.lvl/SI/chewbrms_plot.pdf", chewplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/specbrms_plot.pdf", specplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/minebrms_plot.pdf", mineplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/gallbrms_plot.pdf", gallplot, height = 15, width = 13, units = "in")



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



