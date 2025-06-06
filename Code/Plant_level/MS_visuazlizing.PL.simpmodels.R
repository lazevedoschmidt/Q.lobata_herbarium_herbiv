

#load libraries
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

#SI Function----
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
  
  # Define model labels
  model_labels <- c("Model 1" = "+ DOY",
                    "Model 2" = "+ YearxDOY",
                    "Model 3" = "+ Year"
  )
  
  # Color definitions for parameters (climate variables)
  parameter_color_def <- c(
    "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scalestartDayOfYear" = "#C1AE7C" 
  )
  
  # Model line type definitions
  model_linetype_def <- c(
    "Model 1" = "dashed", 
    "Model 2" = "longdash", 
    "Model 3" = "solid" 
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
    scale_linetype_manual(values = model_linetype_def, 
                          labels = model_labels[names(model_linetype_def)],
                          name = "Model")+
    #scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
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
  
  # Determine model type (year, doy, yearxdoy) from csv_paths
  model_types <- rep(NA, length(csv_paths))
  
  for (i in seq_along(csv_paths)) {
    if (grepl("_year_model\\.csv$", csv_paths[i]) || grepl("/year/", csv_paths[i])) {
      model_types[i] <- "year"
    } else if (grepl("_doy_model\\.csv$", csv_paths[i]) || grepl("/doy/", csv_paths[i])) {
      model_types[i] <- "doy"
    } else if (grepl("_yearxdoy_model\\.csv$", csv_paths[i]) || grepl("/yearxdoy/", csv_paths[i])) {
      model_types[i] <- "yearxdoy"
    } else {
      # If we can't determine from path, use existing mechanism with model names
      if (model_names[i] == "Model 1") {
        model_types[i] <- "doy"
      } else if (model_names[i] == "Model 2") {
        model_types[i] <- "yearxdoy"
      } else if (model_names[i] == "Model 3") {
        model_types[i] <- "year"
      } else {
        # Default fallback
        model_types[i] <- "unknown"
      }
    }
  }
  
  # Create model labels based on detected model types
  model_labels <- setNames(
    ifelse(model_types == "year", "+ Year",
           ifelse(model_types == "doy", "+ DOY", 
                  ifelse(model_types == "yearxdoy", "+ YearxDOY", "Unknown"))),
    model_names
  )
  
  # Create model line types based on detected model types
  model_linetype_def <- setNames(
    ifelse(model_types == "year", "solid",
           ifelse(model_types == "doy", "dashed", 
                  ifelse(model_types == "yearxdoy", "longdash", "dotted"))),
    model_names
  )
  
  # Color definitions for parameters (climate variables)
  parameter_color_def <- c(
    "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scalestartDayOfYear" = "#C1AE7C" 
  )
  
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
    scale_linetype_manual(values = model_linetype_def, 
                          labels = model_labels,
                          name = "Model") +
    # Improved layout with legend at the right
    theme(
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  return(forest_plot)
}

#Total plot----
otf_list <- 
  list.files("Script_output/Brms/simp_models/p.lvl/old/", 
             pattern = "^perc.dam", 
             full.names = TRUE)
otf_list <- otf_list[grepl("_model\\.csv$", otf_list)]
otf_names <- paste("Model", 1:3)
otf_names <- factor(otf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ytf_list <- 
  list.files("Script_output/Brms/simp_models/p.lvl/young/", 
             pattern = "^perc.dam", 
             full.names = TRUE)
ytf_list <- ytf_list[grepl("_model\\.csv$", ytf_list)]
ytf_names <- paste("Model", 1:3)
ytf_names <- factor(ytf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

otd_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^total.div", 
                       full.names = TRUE)
otd_list <- otd_list[grepl("_model\\.csv$", otd_list)]
otd_names <- paste("Model", 1:3)
otd_names <- factor(otd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ytd_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^total.div", 
                       full.names = TRUE)
ytd_list <- ytd_list[grepl("_model\\.csv$", ytd_list)]
ytd_names <- paste("Model", 1:3)
ytd_names <- factor(ytd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

otdplot <- plot_forest_models(otd_list, otd_names)
otdplot
otdplot <- otdplot +
  labs(title = "Plant Level Richness of \nTotal Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
otdplot

otfplot <- plot_forest_models(otf_list, otf_names)
otfplot
otfplot <- otfplot +
  labs(title = "Plant Level Frequency of \nTotal Damage (late-season)")+
  scale_y_discrete (labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
otfplot

#young leaves
ytdplot <- plot_forest_models(ytd_list, ytd_names)
ytdplot
#pulling legend
modellegend <- get_legend(ytdplot)
ytdplot <- ytdplot +
  labs(title = "Plant Level Richness of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ytdplot

ytfplot <- plot_forest_models(ytf_list, ytf_names)
ytfplot
ytfplot <- ytfplot +
  labs(title = "Plant Level Frequency of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year","YearxDOY"))+
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

#80 Total plot----
otdplot80 <- plot80_forest_models(otd_list, otd_names)
otdplot80
otdplot80 <- otdplot80 +
  labs(title = "Plant Level Richness of \nTotal Damage (late-season)")+
  scale_y_discrete(labels = c("Year", "YearxDOY"))+
  theme(legend.position = "none")
otdplot80

otfplot80 <- plot80_forest_models(otf_list, otf_names)
otfplot80
otfplot80 <- otfplot80 +
  labs(title = "Plant Level Frequency of \nTotal Damage (late-season)")+
  scale_y_discrete (labels = c("YearxDOY"))+
  theme(legend.position = "none")
otfplot80

#young leaves
ytdplot80 <- plot80_forest_models(ytd_list, ytd_names)
ytdplot80
#all cross 0

ytfplot80 <- plot80_forest_models(ytf_list, ytf_names)
ytfplot80
#all cross 0

totalplot80 <- plot_grid(NULL,otfplot80,
                         NULL,otdplot80,
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

opac_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
opac_list <- opac_list[grepl("_model\\.csv$", opac_list)]
opac_names <- paste("Model", 1:3)
opac_names <- factor(opac_names, 
                     levels = paste("Model", 1:3),
                     ordered = TRUE)

ypac_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
ypac_list <- ypac_list[grepl("_model\\.csv$", ypac_list)]
ypac_names <- paste("Model", 1:3)
ypac_names <- factor(ypac_names, 
                     levels = paste("Model", 1:3),
                     ordered = TRUE)

#chewing herbs 
ochew_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ochew_list <- ochew_list[grepl("_model\\.csv$", ochew_list)]
ochew_names <- paste("Model", 1:3)
ochew_names <- factor(ochew_names, 
                      levels = paste("Model", 1:3),
                      ordered = TRUE)

ychew_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ychew_list <- ychew_list[grepl("_model\\.csv$", ychew_list)]
ychew_names <- paste("Model", 1:3)
ychew_names <- factor(ychew_names, 
                      levels = paste("Model", 1:3),
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
  labs(title = "Plant Level Percent Leaf \nArea Chewed (late-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
opacplot

ypacplot <- plot_forest_models(ypac_list, ypac_names)
ypacplot
ypacplot <- ypacplot +
  labs(title = "Plant Level Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot

#chewing
ochewplot <- plot_forest_models(ochew_list, ochew_names)
ochewplot
ochewplot <- ochewplot +
  labs(title = "Plant Level Chewing \nDamage (proportion; late-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot

ychewplot <- plot_forest_models(ychew_list, ychew_names)
ychewplot
ychewplot <- ychewplot +
  labs(title = "Plant Level Chewing Damage \n(proportion; early-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ychewplot

#chew plot
chewplot <- plot_grid(ypacplot,opacplot,
                      ychewplot,ochewplot,
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
#all cross 0

ypacplot80 <- plot80_forest_models(ypac_list, ypac_names)
ypacplot80
ypacplot80 <- ypacplot80 +
  labs(title = "Plant Level Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("DOY","YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot80

#chewing
ochewplot80 <- plot80_forest_models(ochew_list, ochew_names)
ochewplot80
ochewplot80 <- ochewplot80 +
  labs(title = "Plant Level Chewing \nDamage (proportion; late-season)") +
  scale_y_discrete(labels=c("DOY", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot80

ychewplot80 <- plot80_forest_models(ychew_list, ychew_names)
ychewplot80
#all overlap 0

#chew plot
chewplot80 <- plot_grid(ypacplot80,NULL,
                        NULL,ochewplot80,
                        nrow = 2,
                        ncol = 2)
chewplot80
chewplot80 <- plot_grid(chewplot80,
                        modellegend,
                        nrow = 1,
                        ncol = 2, 
                        rel_widths = c(1, .10))
chewplot80

#Specialized plot----
osf_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^perc.spec", 
                       full.names = TRUE)
osf_list <- osf_list[grepl("_model\\.csv$", osf_list)]
osf_names <- paste("Model", 1:3)
osf_names <- factor(osf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ysf_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^perc.spec", 
                       full.names = TRUE)
ysf_list <- ysf_list[grepl("_model\\.csv$", ysf_list)]
ysf_names <- paste("Model", 1:3)
ysf_names <- factor(ysf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

osd_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^spec.div", 
                       full.names = TRUE)
osd_list <- osd_list[grepl("_model\\.csv$", osd_list)]
osd_names <- paste("Model", 1:3)
osd_names <- factor(osd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ysd_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^spec.div", 
                       full.names = TRUE)
ysd_list <- ysd_list[grepl("_model\\.csv$", ysd_list)]
ysd_names <- paste("Model", 1:3)
ysd_names <- factor(ysd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

osdplot <- plot_forest_models(osd_list, osd_names)
osdplot
osdplot <- osdplot +
  labs(title = "Plant Level Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
osdplot

osfplot <- plot_forest_models(osf_list, osf_names)
osfplot
osfplot <- osfplot +
  labs(title = "Plant Level Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
osfplot

#young leaves
ysdplot <- plot_forest_models(ysd_list, ysd_names)
ysdplot
ysdplot <- ysdplot +
  labs(title = "Plant Level Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ysdplot

ysfplot <- plot_forest_models(ysf_list, ysf_names)
ysfplot
ysfplot <- ysfplot +
  labs(title = "Plant Level Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("DOY", "Year", "YearxDOY"))+
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

#80 Specialized plot----
osdplot80 <- plot80_forest_models(osd_list, osd_names)
osdplot80
osdplot80 <- osdplot80 +
  labs(title = "Plant Level Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("YearxDOY"))+
  theme(legend.position = "none")
osdplot80

osfplot80 <- plot80_forest_models(osf_list, osf_names)
osfplot80
#All cross 0

#young leaves
ysdplot80 <- plot80_forest_models(ysd_list, ysd_names)
ysdplot80
ysdplot80 <- ysdplot80 +
  labs(title = "Plant Level Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("Year"))+
  theme(legend.position = "none")
ysdplot80

ysfplot80 <- plot80_forest_models(ysf_list, ysf_names)
ysfplot80
#all cross 0

specplot80 <- plot_grid(ysdplot80,osdplot80, 
                        nrow = 1,
                        ncol = 2)
specplot80 <- plot_grid(specplot80, modellegend,
                        nrow = 1, ncol = 2,
                        rel_widths = c(1,.10))
specplot80

#Mine plot----
opam_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
opam_list <- opam_list[grepl("_model\\.csv$", opam_list)]
opam_names <- paste("Model", 1:3)
opam_names <- factor(opam_names, 
                     levels = paste("Model", 1:3),
                     ordered = TRUE)

ypam_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
ypam_list <- ypam_list[grepl("_model\\.csv$", ypam_list)]
ypam_names <- paste("Model", 1:2)
ypam_names <- factor(ypam_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)

omf_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^perc.mine", 
                       full.names = TRUE)
omf_list <- omf_list[grepl("_model\\.csv$", omf_list)]
omf_names <- paste("Model", 1:3)
omf_names <- factor(omf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ymf_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^perc.mine", 
                       full.names = TRUE)
ymf_list <- ymf_list[grepl("_model\\.csv$", ymf_list)]
ymf_names <- paste("Model", 1:3)
ymf_names <- factor(ymf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

omd_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^mine.div", 
                       full.names = TRUE)
omd_list <- omd_list[grepl("_model\\.csv$", omd_list)]
omd_names <- paste("Model", 1:3)
omd_names <- factor(omd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ymd_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^mine.div", 
                       full.names = TRUE)
ymd_list <- ymd_list[grepl("_model\\.csv$", ymd_list)]
ymd_names <- paste("Model", 1:3)
ymd_names <- factor(ymd_names, 
                    levels = paste("Model", 1:3),
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
  labs(title = "Plant Level Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
opamplot

ypamplot <- plot_forest_models(ypam_list, ypam_names)
ypamplot
ypamplot <- ypamplot +
  labs(title = "Plant Level Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ypamplot

omdplot <- plot_forest_models(omd_list, omd_names)
omdplot
omdplot <- omdplot +
  labs(title = "Plant Level Richness of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
omdplot

omfplot <- plot_forest_models(omf_list, omf_names)
omfplot
omfplot <- omfplot +
  labs(title = "Plant Level Frequency of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
omfplot

#young leaves
ymdplot <- plot_forest_models(ymd_list, ymd_names)
ymdplot
ymdplot <- ymdplot +
  labs(title = "Plant Level Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ymdplot

ymfplot <- plot_forest_models(ymf_list, ymf_names)
ymfplot
ymfplot <- ymfplot +
  labs(title = "Plant Level Frequency of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ymfplot

mineplot <- plot_grid( ypamplot,opamplot,
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
#all cross 0

ypamplot80 <- plot80_forest_models(ypam_list, ypam_names)
ypamplot80
ypamplot80 <- ypamplot80 +
  labs(title = "Plant Level Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ypamplot80

omdplot80 <- plot80_forest_models(omd_list, omd_names)
omdplot80
#all cross 0

omfplot80 <- plot80_forest_models(omf_list, omf_names)
omfplot80
#all cross 0

#young leaves
ymdplot80 <- plot80_forest_models(ymd_list, ymd_names)
ymdplot80
ymdplot80 <- ymdplot80 +
  labs(title = "Plant Level Mine \nRichness (early-season)") +
  scale_y_discrete(labels = c("DOY"))+
  theme(legend.position = "none")
ymdplot80

ymfplot80 <- plot80_forest_models(ymf_list, ymf_names)
ymfplot80
#all cross 0

mineplot80 <- plot_grid(ypamplot80, ymdplot80,
                        nrow = 1,
                        ncol = 2)
mineplot80
mineplot80 <- plot_grid(mineplot80,
                        modellegend,
                        nrow = 1,
                        ncol = 2, 
                        rel_widths = c(1,.10))
mineplot80


# Gall plot----
opag_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
opag_list <- opag_list[grepl("_model\\.csv$", opag_list)]
opag_names <- paste("Model", 1:2)
opag_names <- factor(opag_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)

ypag_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
ypag_list <- ypag_list[grepl("_model\\.csv$", ypag_list)]
ypag_names <- paste("Model", 1:2)
ypag_names <- factor(ypag_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)

ogf_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^perc.gall", 
                       full.names = TRUE)
ogf_list <- ogf_list[grepl("_model\\.csv$", ogf_list)]
ogf_names <- paste("Model", 1:3)
ogf_names <- factor(ogf_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ygf_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^perc.gall", 
                       full.names = TRUE)
ygf_list <- ygf_list[grepl("_model\\.csv$", ygf_list)]
ygf_names <- paste("Model", 1:2)
ygf_names <- factor(ygf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)

ogd_list <- list.files("Script_output/Brms/simp_models/p.lvl/old/", 
                       pattern = "^gall.div", 
                       full.names = TRUE)
ogd_list <- ogd_list[grepl("_model\\.csv$", ogd_list)]
ogd_names <- paste("Model", 1:3)
ogd_names <- factor(ogd_names, 
                    levels = paste("Model", 1:3),
                    ordered = TRUE)

ygd_list <- list.files("Script_output/Brms/simp_models/p.lvl/young/", 
                       pattern = "^gall.div", 
                       full.names = TRUE)
ygd_list <- ygd_list[grepl("_model\\.csv$", ygd_list)]
ygd_names <- paste("Model", 1:3)
ygd_names <- factor(ygd_names, 
                    levels = paste("Model", 1:3),
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
  labs(title = "Plant Level Percent Area \nDamaged by Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
opagplot

#young
ypagplot <- plot_forest_models(ypag_list, ypag_names)
ypagplot
ypagplot <- ypagplot +
  labs(title = "Plant Level Percent Area Damaged by \nGalls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year"))+
  theme(legend.position = "none")
ypagplot

ogdplot <- plot_forest_models(ogd_list, ogd_names)
ogdplot
ogdplot <- ogdplot +
  labs(title = "Plant Level Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogdplot

ogfplot <- plot_forest_models(ogf_list, ogf_names)
ogfplot
ogfplot <- ogfplot +
  labs(title = "Plant Level Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogfplot

#young
ygdplot <- plot_forest_models(ygd_list, ygd_names)
ygdplot
ygdplot <- ygdplot +
  labs(title = "Plant Level Richness of \nGalls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygdplot

ygfplot <- plot_forest_models(ygf_list, ygf_names)
ygfplot
ygfplot <- ygfplot +
  labs(title = "Plant level  Frequency \nof Galls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygfplot


gallplot <- plot_grid( ypagplot,opagplot,
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
#all cross 0

#young
ypagplot80 <- plot80_forest_models(ypag_list, ypag_names)
ypagplot80
#all cross 0

ogdplot80 <- plot80_forest_models(ogd_list, ogd_names)
ogdplot80
#all cross 0

ogfplot80 <- plot80_forest_models(ogf_list, ogf_names)
ogfplot80
#all cross 0

#young
ygdplot80 <- plot80_forest_models(ygd_list, ygd_names)
ygdplot80
#all cross 0

ygfplot80 <- plot80_forest_models(ygf_list, ygf_names)
ygfplot80
#all cross 0


# gallplot80 <- plot_grid(ogfplot80,NULL,
#                         ogdplot80,ygdplot80,
#                         nrow = 2,
#                         ncol = 2)
# gallplot80 <- plot_grid(gallplot80,modellegend,
#                         nrow = 1,
#                         ncol = 2,
#                         rel_widths = c(1,.10))
# gallplot80

#saving manuscript plots----
ggsave("figures/p.lvl/simp.totalbrms80_plot.pdf", totalplot80, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/simp.chewbrms80_plot.pdf", chewplot80, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/simp.specbrms80_plot.pdf", specplot80, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/simp.minebrms80_plot.pdf", mineplot80, height = 15, width = 13, units = "in")
#ggsave("figures/p.lvl/simp.gallbrms80_plot.pdf", gallplot80, height = 15, width = 13, units = "in")

#saving SI plots----
ggsave("figures/p.lvl/SI/year_doy.totalbrms_plot.pdf", totalplot, height = 15, width = 13, units = "in") 
ggsave("figures/p.lvl/SI/year_doy.chewbrms_plot.pdf", chewplot, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/SI/year_doy.specbrms_plot.pdf", specplot, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/SI/year_doy.minebrms_plot.pdf", mineplot, height = 15, width = 13, units = "in")
ggsave("figures/p.lvl/SI/year_doy.gallbrms_plot.pdf", gallplot, height = 15, width = 13, units = "in")


