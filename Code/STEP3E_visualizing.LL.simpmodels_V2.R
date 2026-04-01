#Purpose: Visualize Leaf Level model results for "simple" models
#Date started: 03/21/2025
#Author: LAS
#R version: 4.3.3. 

#load packages
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

# Simplified forest plot function for temporal models only (year.x + doy.clean + interaction)

simpplot_forest_models <- function(csv_paths, model_names = NULL, show_all = TRUE) {
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
  
  # Simple parameter ordering for temporal variables only
  # Order: interaction first, then year, then doy (so they appear top to bottom in that order)
  all_params <- unique(all_summaries$parameter)
  param_order <- c("scaleyear.x.scaledoy.clean", "scaleyear.x", "scaledoy.clean")
  param_order <- param_order[param_order %in% all_params]
  
  # Add any unexpected parameters at the end
  remaining_params <- setdiff(all_params, param_order)
  param_order <- c(param_order, remaining_params)
  
  # Convert to ordered factor (reverse order so top items appear at top of plot)
  all_summaries$parameter <- factor(all_summaries$parameter, 
                                    levels = rev(param_order), 
                                    ordered = TRUE)
  
  # Simple color scheme - all temporal variables in one color
  parameter_color_def <- c(
    "scaleyear.x.scaledoy.clean" = "#C1AE7C",
    "scaleyear.x" = "#C1AE7C", 
    "scaledoy.clean" = "#C1AE7C"
  )
  
  # Create line type mapping for models
  model_linetype_def <- setNames(
    rep(c("solid", "dashed", "dotted", "dotdash", "longdash"), 
        length.out = length(unique(all_summaries$model))),
    sort(unique(all_summaries$model))
  )
  
  # Filter data based on show_all parameter
  plot_data <- if (show_all) {
    all_summaries
  } else {
    all_summaries %>% filter(significant_80)
  }
  
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

# # Convenience function for showing only significant estimates (80% CI) ----
# simpplot_forest_significant <- function(csv_paths, model_names = NULL) {
#   return(simpplot_forest_models(csv_paths, model_names, show_all = FALSE))
# }
# 
# # Alternative function with separate plots for comparison ----
# simpplot_forest_comparison <- function(csv_paths, model_names = NULL) {
#   # Create both plots
#   full_plot <- simpplot_forest_models(csv_paths, model_names, show_all = TRUE)
#   sig_plot <- simpplot_forest_models(csv_paths, model_names, show_all = FALSE)
#   
#   # Return a list of both plots
#   return(list(
#     full = full_plot + ggtitle("All Parameter Estimates"),
#     significant = sig_plot + ggtitle("Significant Parameter Estimates (80% CI)")
#   ))
# }

#Total plot----
otf_list <- 
  list.files("Script_output/Brms/simp_models/l.lvl/old/", 
             pattern = "^leaftotal.freq", 
             full.names = TRUE)
otf_list <- otf_list[grepl("_model\\.csv$", otf_list)]
otf_names <- paste("Model", 1:length(otf_list))

ytf_list <- 
  list.files("Script_output/Brms/simp_models/l.lvl/young/", 
             pattern = "^leaftotal.freq", 
             full.names = TRUE)
ytf_list <- ytf_list[grepl("_model\\.csv$", ytf_list)]
ytf_names <- paste("Model", 1:length(ytf_list))

otd_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leaftotal.rich", 
                       full.names = TRUE)
otd_list <- otd_list[grepl("_model\\.csv$", otd_list)]
otd_names <- paste("Model", 1:length(otd_list))

ytd_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leaftotal.rich", 
                       full.names = TRUE)
ytd_list <- ytd_list[grepl("_model\\.csv$", ytd_list)]
ytd_names <- paste("Model", 1:length(ytd_list))

otdplot <- simpplot_forest_models(otd_list, otd_names,show_all = TRUE)
otdplot
otdplot <- otdplot +
  labs(title = "Richness of \nAll Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
otdplot

otfplot <- simpplot_forest_models(otf_list, otf_names,show_all = TRUE)
otfplot
otfplot <- otfplot +
  labs(title = "Frequency of \nAll Damage (late-season)")+
  scale_y_discrete (labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
otfplot

#young leaves
ytdplot <- simpplot_forest_models(ytd_list, ytd_names, show_all = TRUE)
ytdplot
#pulling legend
#modellegend <- get_legend(ytdplot)
ytdplot <- ytdplot +
  labs(title = "Richness of All \nDamage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ytdplot

ytfplot <- simpplot_forest_models(ytf_list, ytf_names, show_all = TRUE)
ytfplot
ytfplot <- ytfplot +
  labs(title = "Frequency of All \nDamage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year","YearxDOY"))+
  theme(legend.position = "none")
ytfplot

totalplot <- plot_grid(ytfplot,otfplot,
                       ytdplot,otdplot,
                       nrow = 2,
                       ncol = 2)
totalplot

#80 Total plot----
otdplot80 <- simpplot_forest_models(otd_list, otd_names, show_all = FALSE)
otdplot80
otdplot80 <- otdplot80 +
  labs(title = "Richness of \nAll Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "YearxDOY"))+
  theme(legend.position = "none")
otdplot80

otfplot80 <- simpplot_forest_models(otf_list, otf_names, show_all = FALSE)
otfplot80
otfplot80 <- otfplot80 +
  labs(title = "Frequency of \nAll Damage (late-season)")+
  scale_y_discrete (labels = c("YearxDOY"))+
  theme(legend.position = "none")
otfplot80

#young leaves
ytdplot80 <- simpplot_forest_models(ytd_list, ytd_names, show_all = FALSE)
ytdplot80
#all cross 0

ytfplot80 <- simpplot_forest_models(ytf_list, ytf_names,show_all = FALSE)
ytfplot80
#all cross 0

totalplot80 <- plot_grid(otfplot80,otdplot80,
                         nrow = 1,
                         ncol = 2)

totalplot80

#Chewing plot----
opac_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
opac_list <- opac_list[grepl("_model\\.csv$", opac_list)]
opac_names <- paste("Model", 1:length(opac_list))

ypac_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                        pattern = "^leafperc_area_chew", 
                        full.names = TRUE)
ypac_list <- ypac_list[grepl("_model\\.csv$", ypac_list)]
ypac_names <- paste("Model", 1:length(ypac_list))

#chewing herbs 
ochew_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ochew_list <- ochew_list[grepl("_model\\.csv$", ochew_list)]
ochew_names <- paste("Model", 1:length(ochew_list))

ychew_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                         pattern = "^chewvalue", 
                         full.names = TRUE)
ychew_list <- ychew_list[grepl("_model\\.csv$", ychew_list)]
ychew_names <- paste("Model", 1:length(ychew_list))

opacplot <- simpplot_forest_models(opac_list, opac_names, show_all = TRUE)
opacplot
opacplot <- opacplot +
  labs(title = "Percent Leaf \nArea Chewed (late-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
opacplot

ypacplot <- simpplot_forest_models(ypac_list, ypac_names,show_all = TRUE)
ypacplot
ypacplot <- ypacplot +
  labs(title = "Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot

#chewing
ochewplot <- simpplot_forest_models(ochew_list, ochew_names,show_all = TRUE)
ochewplot
ochewplot <- ochewplot +
  labs(title = "Chewing \nDamage (proportion; late-season)") +
  scale_y_discrete(labels=c("DOY", "Year", "YearxDOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot

ychewplot <- simpplot_forest_models(ychew_list, ychew_names,show_all = TRUE)
ychewplot
ychewplot <- ychewplot +
  labs(title = "Chewing Damage \n(proportion; early-season)") +
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
opacplot80 <- simpplot_forest_models(opac_list, opac_names,show_all = FALSE)
opacplot80
#all cross 0

ypacplot80 <- simpplot_forest_models(ypac_list, ypac_names,show_all = FALSE)
ypacplot80
#all cross 0

#chewing
ochewplot80 <- simpplot_forest_models(ochew_list, ochew_names,show_all = FALSE)
ochewplot80
# all cross 0

ychewplot80 <- simpplot_forest_models(ychew_list, ychew_names,show_all = FALSE)
ychewplot80
#all overlap 0

#Generalized plot----
ogenf_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafgen.freq", 
                       full.names = TRUE)
ogenf_list <- ogenf_list[grepl("_model\\.csv$", ogenf_list)]
ogenf_names <- paste("Model", 1:length(ogenf_list))

ygenf_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafgen.freq", 
                       full.names = TRUE)
ygenf_list <- ygenf_list[grepl("_model\\.csv$", ygenf_list)]
ygenf_names <- paste("Model", 1:length(ygenf_list))

ogend_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafgen.rich", 
                       full.names = TRUE)
ogend_list <- ogend_list[grepl("_model\\.csv$", ogend_list)]
ogend_names <- paste("Model", 1:length(ogend_list))

ygend_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafgen.rich", 
                       full.names = TRUE)
ygend_list <- ygend_list[grepl("_model\\.csv$", ygend_list)]
ygend_names <- paste("Model", 1:length(ygend_list))

ogendplot <- simpplot_forest_models(ogend_list, ogend_names,show_all = TRUE)
ogendplot
ogendplot <- ogendplot +
  labs(title = "Richness of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogendplot

ogenfplot <- simpplot_forest_models(ogenf_list, ogenf_names,show_all = TRUE)
ogenfplot
ogenfplot <- ogenfplot +
  labs(title = "Frequency of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogenfplot

#young leaves
ygendplot <- simpplot_forest_models(ygend_list, ygend_names,show_all = TRUE)
ygendplot
ygendplot <- ygendplot +
  labs(title = "Richness of Generalized \nDamage (early-season)")+
  scale_y_discrete(labels =c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygendplot

ygenfplot <- simpplot_forest_models(ygenf_list, ygenf_names,show_all = TRUE)
ygenfplot
ygenfplot <- ygenfplot +
  labs(title = "Frequency of Generalized \nDamage (early-season)")+
  scale_y_discrete(labels= c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygenfplot

genplot <- plot_grid(ygenfplot,ogenfplot, 
                      ygendplot,ogendplot, 
                      nrow = 2,
                      ncol = 2)
genplot <- plot_grid(genplot, modellegend,
                      nrow = 1, ncol = 2,
                      rel_widths = c(1,.10))
genplot

#80 Generalized plot----
ogendplot80 <- simpplot_forest_models(ogend_list, ogend_names,show_all = FALSE)
ogendplot80
ogendplot80 <- ogendplot80 +
  labs(title = "Richness of \nGeneralized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY"))+
  theme(legend.position = "none")
ogendplot80

ogenfplot80 <- simpplot_forest_models(ogenf_list, ogenf_names,show_all = FALSE)
ogenfplot80
#all cross zero

#young leaves
ygendplot80 <- simpplot_forest_models(ygend_list, ygend_names,show_all = FALSE)
ygendplot80
#all cross zero

ygenfplot80 <- simpplot_forest_models(ygenf_list, ygenf_names,show_all = FALSE)
ygenfplot80
#all cross 0

#Only ogenplot80 has results

#Specialized plot----
osf_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE)
osf_list <- osf_list[grepl("_model\\.csv$", osf_list)]
osf_names <- paste("Model", 1:length(osf_list))

ysf_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafspec.freq", 
                       full.names = TRUE)
ysf_list <- ysf_list[grepl("_model\\.csv$", ysf_list)]
ysf_names <- paste("Model", 1:length(ysf_list))

osd_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE)
osd_list <- osd_list[grepl("_model\\.csv$", osd_list)]
osd_names <- paste("Model", 1:length(osd_list))

ysd_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafspec.rich", 
                       full.names = TRUE)
ysd_list <- ysd_list[grepl("_model\\.csv$", ysd_list)]
ysd_names <- paste("Model", 1:length(ysd_list))

osdplot <- simpplot_forest_models(osd_list, osd_names,show_all = TRUE)
osdplot
osdplot <- osdplot +
  labs(title = "Richness of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
osdplot

osfplot <- simpplot_forest_models(osf_list, osf_names,show_all = TRUE)
osfplot
osfplot <- osfplot +
  labs(title = "Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
osfplot

#young leaves
ysdplot <- simpplot_forest_models(ysd_list, ysd_names,show_all = TRUE)
ysdplot
ysdplot <- ysdplot +
  labs(title = "Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ysdplot

ysfplot <- simpplot_forest_models(ysf_list, ysf_names,show_all = TRUE)
ysfplot
ysfplot <- ysfplot +
  labs(title = "Frequency of Specialized \nDamage (early-season)")+
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
osdplot80 <- simpplot_forest_models(osd_list, osd_names,show_all = FALSE)
osdplot80
#all cross zero

osfplot80 <- simpplot_forest_models(osf_list, osf_names,show_all = FALSE)
osfplot80
osfplot80 <- osfplot80 +
  labs(title = "Frequency of \nSpecialized Damage (late-season)")+
  scale_y_discrete(labels = c("YearxDOY"))+
  theme(legend.position = "none")
osfplot80

#young leaves
ysdplot80 <- simpplot_forest_models(ysd_list, ysd_names,show_all = FALSE)
ysdplot80
#all cross zero

ysfplot80 <- simpplot_forest_models(ysf_list, ysf_names)
ysfplot80
#all cross 0

#only osfplot80 is significant

#Mine plot----
opam_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
opam_list <- opam_list[grepl("_model\\.csv$", opam_list)]
opam_names <- paste("Model", 1:length(opam_list))

ypam_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                        pattern = "^leafperc_area_mine", 
                        full.names = TRUE)
ypam_list <- ypam_list[grepl("_model\\.csv$", ypam_list)]
ypam_names <- paste("Model", 1:length(ypam_list))

omf_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE)
omf_list <- omf_list[grepl("_model\\.csv$", omf_list)]
omf_names <- paste("Model", 1:length(omf_list))

ymf_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafmine.freq", 
                       full.names = TRUE)
ymf_list <- ymf_list[grepl("_model\\.csv$", ymf_list)]
ymf_names <- paste("Model", 1:length(ymf_list))

omd_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE)
omd_list <- omd_list[grepl("_model\\.csv$", omd_list)]
omd_names <- paste("Model", 1:length(omd_list))

ymd_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafmine.rich", 
                       full.names = TRUE)
ymd_list <- ymd_list[grepl("_model\\.csv$", ymd_list)]
ymd_names <- paste("Model", 1:length(ymd_list))

opamplot <- simpplot_forest_models(opam_list, opam_names,show_all = TRUE)
opamplot
opamplot <- opamplot +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
opamplot

ypamplot <- simpplot_forest_models(ypam_list, ypam_names,show_all = TRUE)
ypamplot
ypamplot <- ypamplot +
  labs(title = "Percent Leaf Area \nMined (early-season)") +
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ypamplot

omdplot <- simpplot_forest_models(omd_list, omd_names,show_all = TRUE)
omdplot
omdplot <- omdplot +
  labs(title = "Richness of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
omdplot

omfplot <- simpplot_forest_models(omf_list, omf_names,show_all = TRUE)
omfplot
omfplot <- omfplot +
  labs(title = "Frequency of \nMine Damage (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
omfplot

#young leaves
ymdplot <- simpplot_forest_models(ymd_list, ymd_names,show_all = TRUE)
ymdplot
ymdplot <- ymdplot +
  labs(title = "Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ymdplot

ymfplot <- simpplot_forest_models(ymf_list, ymf_names,show_all = TRUE)
ymfplot
ymfplot <- ymfplot +
  labs(title = "Frequency of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ymfplot

mineplot <- plot_grid(ypamplot,opamplot, 
                      ymfplot,omfplot,
                      ymdplot,omdplot,
                      nrow = 3,
                      ncol = 2)
mineplot

#80 Mine plot----
opamplot80 <- simpplot_forest_models(opam_list, opam_names,show_all = FALSE)
opamplot80
opamplot80 <- opamplot80 +
  labs(title = "Percent Leaf \nArea Mined (late-season)") +
  scale_y_discrete(labels = c("Year"))+
  theme(legend.position = "none")
opamplot80

ypamplot80 <- simpplot_forest_models(ypam_list, ypam_names,show_all = FALSE)
ypamplot80
#all cross zero

omdplot80 <- simpplot_forest_models(omd_list, omd_names,show_all = FALSE)
omdplot80
#all cross 0

omfplot80 <- simpplot_forest_models(omf_list, omf_names,show_all = FALSE)
omfplot80
#all cross 0

#young leaves
ymdplot80 <- simpplot_forest_models(ymd_list, ymd_names,show_all = FALSE)
ymdplot80
#all cross 0

ymfplot80 <- simpplot_forest_models(ymf_list, ymf_names,show_all = FALSE)
ymfplot80
#all cross 0


# Gall plot----
opag_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
opag_list <- opag_list[grepl("_model\\.csv$", opag_list)]
opag_names <- paste("Model", 1:length(opag_list))

ypag_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                        pattern = "^leafperc_area_gall", 
                        full.names = TRUE)
ypag_list <- ypag_list[grepl("_model\\.csv$", ypag_list)]
ypag_names <- paste("Model", 1:length(ypag_list))

ogf_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE)
ogf_list <- ogf_list[grepl("_model\\.csv$", ogf_list)]
ogf_names <- paste("Model", 1:length(ogf_list))

ygf_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafgall.freq", 
                       full.names = TRUE)
ygf_list <- ygf_list[grepl("_model\\.csv$", ygf_list)]
ygf_names <- paste("Model", 1:length(ygf_list))

ogd_list <- list.files("Script_output/Brms/simp_models/l.lvl/old/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE)
ogd_list <- ogd_list[grepl("_model\\.csv$", ogd_list)]
ogd_names <- paste("Model", 1:length(ogd_list))

ygd_list <- list.files("Script_output/Brms/simp_models/l.lvl/young/", 
                       pattern = "^leafgall.rich", 
                       full.names = TRUE)
ygd_list <- ygd_list[grepl("_model\\.csv$", ygd_list)]
ygd_names <- paste("Model", 1:length(ygd_list))

#plot--
opagplot <- simpplot_forest_models(opag_list, opag_names,show_all = TRUE)
opagplot
opagplot <- opagplot +
  labs(title = "Percent Area \nDamaged by Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
opagplot

#young
ypagplot <- simpplot_forest_models(ypag_list, ypag_names,show_all = TRUE)
ypagplot
ypagplot <- ypagplot +
  labs(title = "Percent Area Damaged by \nGalls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ypagplot

ogdplot <- simpplot_forest_models(ogd_list, ogd_names,show_all = TRUE)
ogdplot
ogdplot <- ogdplot +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogdplot

ogfplot <- simpplot_forest_models(ogf_list, ogf_names,show_all = TRUE)
ogfplot
ogfplot <- ogfplot +
  labs(title = "Frequency \nof Galls (late-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ogfplot

#young
ygdplot <- simpplot_forest_models(ygd_list, ygd_names,show_all = TRUE)
ygdplot
ygdplot <- ygdplot +
  labs(title = "Richness of \nGalls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygdplot

ygfplot <- simpplot_forest_models(ygf_list, ygf_names,show_all = TRUE)
ygfplot
ygfplot <- ygfplot +
  labs(title = "Frequency \nof Galls (early-season)")+
  scale_y_discrete(labels = c("DOY", "Year", "YearxDOY"))+
  theme(legend.position = "none")
ygfplot


gallplot <- plot_grid(ypagplot,opagplot, 
                      ygfplot,ogfplot,
                      ygdplot,ogdplot,
                      nrow = 3,
                      ncol = 2)
gallplot

#80 Gall plot----
opagplot80 <- simpplot_forest_models(opag_list, opag_names,show_all = FALSE)
opagplot80
#all cross 0

#young
ypagplot80 <- simpplot_forest_models(ypag_list, ypag_names,show_all = FALSE)
ypagplot80
#all cross 0

ogdplot80 <- simpplot_forest_models(ogd_list, ogd_names,show_all = FALSE)
ogdplot80
ogdplot80 <- ogdplot80 +
  labs(title = "Richness \nof Galls (late-season)")+
  scale_y_discrete(labels = c("DOY"))+
  theme(legend.position = "none")
ogdplot80

ogfplot80 <- simpplot_forest_models(ogf_list, ogf_names,show_all = FALSE)
ogfplot80
#all cross 0

#young
ygdplot80 <- simpplot_forest_models(ygd_list, ygd_names,show_all = FALSE)
ygdplot80
#all cross 0

ygfplot80 <- simpplot_forest_models(ygf_list, ygf_names,show_all = FALSE)
ygfplot80
#all cross 0

#only ogdplot80 is significant 

#saving manuscript plots----
ggsave("figures/l.lvl/simp.totalbrms80_plot.pdf", totalplot80, height = 15, width = 13, units = "in")
#ggsave("figures/l.lvl/simp.chewbrms80_plot.pdf", chewplot80, height = 15, width = 13, units = "in")
#no significant chewing damage
ggsave("figures/l.lvl/simp.specbrms80_plot.pdf", osfplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/simp.minebrms80_plot.pdf", opamplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/simp.gallbrms80_plot.pdf", ogdplot80, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/simp.genbrms80_plot.pdf", ogendplot80, height = 15, width = 13, units = "in")


#saving SI plots----
ggsave("figures/l.lvl/SI/year_doy.totalbrms_plot.pdf", totalplot, height = 15, width = 13, units = "in") 
ggsave("figures/l.lvl/SI/year_doy.chewbrms_plot.pdf", chewplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/year_doy.specbrms_plot.pdf", specplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/year_doy.minebrms_plot.pdf", mineplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/year_doy.gallbrms_plot.pdf", gallplot, height = 15, width = 13, units = "in")
ggsave("figures/l.lvl/SI/year_doy.genbrms_plot.pdf", genplot, height = 15, width = 13, units = "in")


#function boneyard----
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
#   # Define model labels
#   model_labels <- c("Model 1" = "+ DOY",
#                     "Model 2" = "+ YearxDOY",
#                     "Model 3" = "+ Year"
#   )
#   
#   # Color definitions for parameters (climate variables)
#   parameter_color_def <- c(
#     "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
#     "scaleyear.x" = "#C1AE7C", 
#     "scalestartDayOfYear" = "#C1AE7C" 
#   )
#   
#   # Model line type definitions
#   model_linetype_def <- c(
#     "Model 1" = "dashed", 
#     "Model 2" = "longdash", 
#     "Model 3" = "solid" 
#   )
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
#     scale_linetype_manual(values = model_linetype_def, 
#                           labels = model_labels[names(model_linetype_def)],
#                           name = "Model")+
#     #scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
#     # Improved layout with legend at the right
#     theme(
#       legend.position = "right",
#       legend.box = "vertical",
#       panel.grid.minor = element_blank()
#     )
#   
#   return(forest_plot)
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
#   # Define model labels
#   model_labels <- c("Model 1" = "+ DOY",
#                     "Model 2" = "+ YearxDOY",
#                     "Model 3" = "+ Year"
#   )
#   
#   # Color definitions for parameters (climate variables)
#   parameter_color_def <- c(
#     "scaleyear.x.scalestartDayOfYear" = "#C1AE7C",
#     "scaleyear.x" = "#C1AE7C", 
#     "scalestartDayOfYear" = "#C1AE7C" 
#   )
#   
#   # Model line type definitions
#   model_linetype_def <- c(
#     "Model 1" = "dashed", 
#     "Model 2" = "longdash", 
#     "Model 3" = "solid" 
#   )
#   
#   # Filter data to only keep estimates that don't cross 0 (80% CI)
#   significant_summaries <- all_summaries %>%
#     filter(significant_80)
#   
#   # Create the forest plot with only significant estimates
#   forest_plot <- significant_summaries %>%
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
#     geom_vline(xintercept = 0, linetype = "dashed") +
#     theme_minimal() +
#     labs(x = "Estimate", y = "Parameter") +
#     # Colors by parameter, line types by model
#     scale_color_manual(values = parameter_color_def, name = "Parameter") +
#     scale_linetype_manual(values = model_linetype_def, 
#                           labels = model_labels[names(model_linetype_def)],
#                           name = "Model")+
#     #scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +
#     # Improved layout with legend at the right
#     theme(
#       legend.position = "right",
#       legend.box = "vertical",
#       panel.grid.minor = element_blank()
#     )
#   
#   return(forest_plot)
# }
