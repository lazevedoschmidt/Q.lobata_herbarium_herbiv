#Purpose: Visualize PLANT level timespace correlation models
#Author: LAS
#Date started: 03.20.2025
#R version: 4.3.3.

#Load packages
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

#path to data folder: Script_output/Brms/time_spac_models/p.lvl/young/

#Function----
tspa_forest_models <- function(csv_paths, model_names = NULL) {  
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
        lower_80 = quantile(value, 0.1),  # 80% CI lower bound  
        upper_80 = quantile(value, 0.9),  # 80% CI upper bound  
        positive = sum(value > 0) / length(value),  
        negative = sum(value < 0) / length(value)  
      ) %>%  
      mutate(  
        parameter = gsub("b_", "", parameter),  
        model = model_name,  
        significant_80 = (lower_80 > 0) | (upper_80 < 0) # Flag estimates that don’t cross zero (80% CI)  
      )  
    
    return(summary_stats)  
  }  
  
  # Process all models  
  all_summaries <- map2_df(models_data, model_names, process_model)  
  
  # Define model labels  
  model_labels <- c("Model 1" = "+ Latitude", "Model 2" = "+ Longitude")  
  
  # Model color definitions  
  model_color_def <- c("Model 1" = "#5e548e", "Model 2" = "#f08080")  
  
  # Model line type definitions  
  model_linetype_def <- c("Model 1" = "solid", "Model 2" = "dashed")  
  
  # Filter data to only keep estimates that don’t cross 0 (80% CI)  
  significant_summaries <- all_summaries %>% filter(significant_80)  
  
  # Create the forest plot with only significant estimates  
  forest_plot <- significant_summaries %>%  
    ggplot(aes(y = parameter, x = mean, color = model, linetype = model)) +  
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
    # Add model labels next to lines  
    # geom_text(aes(x = upper_95 + abs(upper_95) * 0.1, label = model_labels[model]),  
    #           position = position_dodge(width = 0.5),  
    #           hjust = 0,   
    #           size = 3,   
    #           show.legend = FALSE) +  
    theme_minimal() +  
    labs(x = "Estimate", y = "Parameter") +  
    # Update legend labels  
    scale_color_manual(values = model_color_def, name = "Model", labels = model_labels) +  
    scale_linetype_manual(values = model_linetype_def, name = "Model", labels = model_labels) +  
    # Improved layout with legend at the top  
    theme(
      legend.position = "right",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
  
  return(forest_plot)  
}

#Total data----
#frequency
otf_list <- c(
  list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
             pattern = "^perc.dam", 
             full.names = TRUE)
)
otf_list <- otf_list[!(grepl(".summary\\.csv$", otf_list) | grepl("_R2\\.csv$", otf_list))]
otf_names <- paste("Model", 1:2)
otf_names <- factor(otf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)

#OTF plot----
otfplot80 <- tspa_forest_models(otf_list, otf_names)
otfplot80
#All cross 0

ytf_list <- c(
  list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
             pattern = "^perc.dam", 
             full.names = TRUE)
)
ytf_list <- ytf_list[!(grepl(".summary\\.csv$", ytf_list) | grepl("_R2\\.csv$", ytf_list))]
ytf_names <- paste("Model", 1:2)
ytf_names <- factor(ytf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)

#YTF plot----
ytfplot80 <- tspa_forest_models(ytf_list, ytf_names)
ytfplot80
tspalegend <- get_legend(ytfplot80)
ytfplot80 <- ytfplot80 +
  labs(title = "Plant Level Frequency of \nTotal Damage (early-season)")+
  scale_y_discrete (labels = c("Latitude","Longitude"))+
  theme(legend.position = "none")
ytfplot80

#richness
otd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^total.div", 
                         full.names = TRUE)
)
otd_list <- otd_list[!(grepl(".summary\\.csv$", otd_list) | grepl("_R2\\.csv$", otd_list))]
otd_names <- paste("Model", 1:2)
otd_names <- factor(otd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OTD plot----
otdplot80 <- tspa_forest_models(otd_list, otd_names)
otdplot80
otdplot80 <- otdplot80 +
  labs(title = "Plant Level Richness of \nTotal Damage (late-season)")+
  scale_y_discrete(labels = c("Year"))+
  theme(legend.position = "none")
otdplot80

ytd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^total.div", 
                         full.names = TRUE)
)
ytd_list <- ytd_list[!(grepl(".summary\\.csv$", ytd_list) | grepl("_R2\\.csv$", ytd_list))]
ytd_names <- paste("Model", 1:2)
ytd_names <- factor(ytd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YTD plot----
ytdplot80 <- tspa_forest_models(ytd_list, ytd_names)
ytdplot80
ytdplot80 <- ytdplot80 +
  labs(title = "Plant Level Richness of Total \nDamage (early-season)")+
  scale_y_discrete(labels = c("Latitude", "Longitude"))+
  theme(legend.position = "none")
ytdplot80

#combined total plots----
totalplot80 <- plot_grid(ytfplot80,tspalegend, #using blank space for legend
                         ytdplot80,otdplot80,
                         nrow = 2,
                         ncol = 2)
totalplot80

#Chewing data----
opac_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                          pattern = "^leafperc_area_chew", 
                          full.names = TRUE)
)
opac_list <- opac_list[!(grepl(".summary\\.csv$", opac_list) | grepl("_R2\\.csv$", opac_list))]
opac_names <- paste("Model", 1:2)
opac_names <- factor(opac_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)

#OPAC plot----
opacplot <- tspa_forest_models(opac_list, opac_names)
opacplot
#all cross 0

ypac_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                          pattern = "^leafperc_area_chew", 
                          full.names = TRUE)
)
ypac_list <- ypac_list[!(grepl(".summary\\.csv$", ypac_list) | grepl("_R2\\.csv$", ypac_list))]
ypac_names <- paste("Model", 1:2)
ypac_names <- factor(ypac_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)

#YPAC plot----
ypacplot <- tspa_forest_models(ypac_list, ypac_names)
ypacplot
ypacplot <- ypacplot +
  labs(title = "Plant Level Percent Leaf Area \nChewed (early-season)") +
  scale_y_discrete(labels=c("Latitude", "Longitude", "DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ypacplot

#pres/abs chewing 
ochew_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                           pattern = "^chewvalue", 
                           full.names = TRUE)
)
ochew_list <- ochew_list[!(grepl(".summary\\.csv$", ochew_list) | grepl("_R2\\.csv$", ochew_list))]
ochew_names <- paste("Model", 1:2)
ochew_names <- factor(ochew_names, 
                      levels = paste("Model", 1:2),
                      ordered = TRUE)

#OCHEW plot----
ochewplot <- tspa_forest_models(ochew_list, ochew_names)
ochewplot
ochewplot <- ochewplot +
  labs(title = "Plant Level Chewing \nDamage (proportion; late-season)") +
  scale_y_discrete(labels=c("DOY"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ochewplot


ychew_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                           pattern = "^chewvalue", 
                           full.names = TRUE)
)
ychew_list <- ychew_list[!(grepl(".summary\\.csv$", ychew_list) | grepl("_R2\\.csv$", ychew_list))]
ychew_names <- paste("Model", 1:2)
ychew_names <- factor(ychew_names, 
                      levels = paste("Model", 1:2),
                      ordered = TRUE)
#YCHEW plot----
ychewplot <- tspa_forest_models(ychew_list, ychew_names)
ychewplot 
ychewplot <- ychewplot +
  labs(title = "Plant Level Chewing Damage \n(proportion; early-season)") +
  scale_y_discrete(labels=c("Longitude"))+ #we are replacing the y-axis (cord flip in function)
  theme(legend.position = "none") 
ychewplot

#combined chew plot----
chewplot <- plot_grid(ypacplot,tspalegend, 
                      ychewplot, ochewplot,
                      nrow = 2,
                      ncol = 2)
chewplot

#Specialized data----
osf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^perc.spec", 
                         full.names = TRUE)
)
osf_list <- osf_list[!(grepl(".summary\\.csv$", osf_list) | grepl("_R2\\.csv$", osf_list))]
osf_names <- paste("Model", 1:2)
osf_names <- factor(osf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OSF plot----
osfplot <- tspa_forest_models(osf_list, osf_names)
osfplot
#all cross zero

ysf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^perc.spec", 
                         full.names = TRUE)
)
ysf_list <- ysf_list[!(grepl(".summary\\.csv$", ysf_list) | grepl("_R2\\.csv$", ysf_list))]
ysf_names <- paste("Model", 1:2)
ysf_names <- factor(ysf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YSF plot----
ysfplot <- tspa_forest_models(ysf_list, ysf_names)
ysfplot
ysfplot <- ysfplot +
  labs(title = "Plant Level Frequency of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels= c("Latitude", "Longitude"))+
  theme(legend.position = "none")
ysfplot

osd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^spec.div", 
                         full.names = TRUE)
)
osd_list <- osd_list[!(grepl(".summary\\.csv$", osd_list) | grepl("_R2\\.csv$", osd_list))]
osd_names <- paste("Model", 1:2)
osd_names <- factor(osd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OSD plot----
osdplot <- tspa_forest_models(osd_list, osd_names)
osdplot
osdplot <- osdplot +
  labs(title = "Plant Level Richness of Specialized \nDamage (late-season)")+
  scale_y_discrete(labels =c("Year"))+
  theme(legend.position = "none")
osdplot

ysd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^spec.div", 
                         full.names = TRUE)
)
ysd_list <- ysd_list[!(grepl(".summary\\.csv$", ysd_list) | grepl("_R2\\.csv$", ysd_list))]
ysd_names <- paste("Model", 1:2)
ysd_names <- factor(ysd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YSD plot----
ysdplot <- tspa_forest_models(ysd_list, ysd_names)
ysdplot
ysdplot <- ysdplot +
  labs(title = "Plant Level Richness of Specialized \nDamage (early-season)")+
  scale_y_discrete(labels =c("Latitude", "Longitude", "Year"))+
  theme(legend.position = "none")
ysdplot

#combinded Specialized plot----
specplot <- plot_grid(ysfplot,tspalegend, 
                      ysdplot,osdplot, 
                      nrow = 2,
                      ncol = 2)
specplot

#Mining data----
opam_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                          pattern = "^leafperc_area_mine", 
                          full.names = TRUE)
)
opam_list <- opam_list[!(grepl(".summary\\.csv$", opam_list) | grepl("_R2\\.csv$", opam_list))]
opam_names <- paste("Model", 1:2)
opam_names <- factor(opam_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)
#OPAM plot----
opamplot <- tspa_forest_models(opam_list, opam_names)
opamplot
#all cross zero

ypam_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                          pattern = "^leafperc_area_mine", 
                          full.names = TRUE)
)
ypam_list <- ypam_list[!(grepl(".summary\\.csv$", ypam_list) | grepl("_R2\\.csv$", ypam_list))]
ypam_names <- paste("Model", 1:2)
ypam_names <- factor(ypam_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)
#YPAM plot----
ypamplot <- tspa_forest_models(ypam_list, ypam_names)
ypamplot 
ypamplot <- ypamplot +
  labs(title = "Plant Level Percent Leaf Area \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("Longitude", "DOY"))+
  theme(legend.position = "none")
ypamplot


omf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^perc.mine", 
                         full.names = TRUE)
)
omf_list <- omf_list[!(grepl(".summary\\.csv$", omf_list) | grepl("_R2\\.csv$", omf_list))]
omf_names <- paste("Model", 1:2)
omf_names <- factor(omf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OMF plot----
omfplot <- tspa_forest_models(omf_list, omf_names)
omfplot
#all cross zero

ymf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^perc.mine", 
                         full.names = TRUE)
)
ymf_list <- ymf_list[!(grepl(".summary\\.csv$", ymf_list) | grepl("_R2\\.csv$", ymf_list))]
ymf_names <- paste("Model", 1:2)
ymf_names <- factor(ymf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YMF plot----
ymfplot <- tspa_forest_models(ymf_list, ymf_names)
ymfplot
ymfplot <- ymfplot +
  labs(title = "Plant Level Frequency of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("Latitude", "Longitude"))+
  theme(legend.position = "none")
ymfplot

omd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^mine.div", 
                         full.names = TRUE)
)
omd_list <- omd_list[!(grepl(".summary\\.csv$", omd_list) | grepl("_R2\\.csv$", omd_list))]
omd_names <- paste("Model", 1:2)
omd_names <- factor(omd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OMD plot----
omdplot <- tspa_forest_models(omd_list, omd_names)
omdplot
#all cross zero

ymd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^mine.div", 
                         full.names = TRUE)
)
ymd_list <- ymd_list[!(grepl(".summary\\.csv$", ymd_list) | grepl("_R2\\.csv$", ymd_list))]
ymd_names <- paste("Model", 1:2)
ymd_names <- factor(ymd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YMD plot----
ymdplot <- tspa_forest_models(ymd_list, ymd_names)
ymdplot
ymdplot <- ymdplot +
  labs(title = "Plant Level Richness of \nMine Damage (early-season)")+
  scale_y_discrete(labels = c("DOY"))+
  theme(legend.position = "none")
ymdplot

#combined Mine plot----
mineplot <- plot_grid(ypamplot, 
                      ymfplot, 
                      ymdplot, 
                      nrow = 3,
                      ncol = 1)
mineplot <- plot_grid(mineplot,
                      tspalegend,
                      ncol = 2, 
                      nrow = 1,
                      rel_widths = c(1, 0.15))
mineplot

#Gall data----
opag_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                          pattern = "^leafperc_area_gall", 
                          full.names = TRUE)
)
opag_list <- opag_list[!(grepl(".summary\\.csv$", opag_list) | grepl("_R2\\.csv$", opag_list))]
opag_names <- paste("Model", 1:2)
opag_names <- factor(opag_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)
#OPAG plot----
opagplot <- tspa_forest_models(opag_list, opag_names)
opagplot
#all cross zero

ypag_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                          pattern = "^leafperc_area_gall", 
                          full.names = TRUE)
)
ypag_list <- ypag_list[!(grepl(".summary\\.csv$", ypag_list) | grepl("_R2\\.csv$", ypag_list))]
ypag_names <- paste("Model", 1:2)
ypag_names <- factor(ypag_names, 
                     levels = paste("Model", 1:2),
                     ordered = TRUE)
#YPAG plot----
#No early season percent area galled data
# ypagplot <- tspa_forest_models(ypag_list, ypag_names)
# ypagplot
#all cross zero

ogf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^perc.gall", 
                         full.names = TRUE)
)
ogf_list <- ogf_list[!(grepl(".summary\\.csv$", ogf_list) | grepl("_R2\\.csv$", ogf_list))]
ogf_names <- paste("Model", 1:2)
ogf_names <- factor(ogf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OGF plot----
ogfplot <- tspa_forest_models(ogf_list, ogf_names)
ogfplot
#all cross zero

ygf_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^perc.gall", 
                         full.names = TRUE)
)
ygf_list <- ygf_list[!(grepl(".summary\\.csv$", ygf_list) | grepl("_R2\\.csv$", ygf_list))]
ygf_names <- paste("Model", 1:2)
ygf_names <- factor(ygf_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YGF plot----
ygfplot <- tspa_forest_models(ygf_list, ygf_names)
ygfplot
#all cross zero

ogd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/old/", 
                         pattern = "^gall.div", 
                         full.names = TRUE)
)
ogd_list <- ogd_list[!(grepl(".summary\\.csv$", ogd_list) | grepl("_R2\\.csv$", ogd_list))]
ogd_names <- paste("Model", 1:2)
ogd_names <- factor(ogd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#OGD plot----
ogdplot <- tspa_forest_models(ogd_list, ogd_names)
ogdplot
#all cross zero

ygd_list <- c(list.files("Script_output/Brms/time_spac_models/p.lvl/young/", 
                         pattern = "^gall.div", 
                         full.names = TRUE)
)
ygd_list <- ygd_list[!(grepl(".summary\\.csv$", ygd_list) | grepl("_R2\\.csv$", ygd_list))]
ygd_names <- paste("Model", 1:2)
ygd_names <- factor(ygd_names, 
                    levels = paste("Model", 1:2),
                    ordered = TRUE)
#YGD plot----
ygdplot <- tspa_forest_models(ygd_list, ygd_names)
ygdplot
#all cross zero

#combined Gall plot----
# gallplot <- plot_grid(ypagplot,opagplot, 
#                       ygfplot,ogfplot,
#                       ygdplot,ogdplot,
#                       nrow = 3,
#                       ncol = 2)
# gallplot

#saving plots----
ggsave("figures/p.lvl/PL.tspa.total_plot.pdf", totalplot80, height = 11, width = 8.5, units = "in")
ggsave("figures/p.lvl/PL.tspa.chew_plot.pdf", chewplot, height = 11, width = 8.5, units = "in")
ggsave("figures/p.lvl/PL.tspa.spec_plot.pdf", specplot, height = 11, width = 8.5, units = "in")
ggsave("figures/p.lvl/PL.tspa.mine_plot.pdf", mineplot, height = 11, width = 8.5, units = "in")
#ggsave("figures/p.lvl/PL.tspa.gall_plot.pdf", gallplot, height = 11, width = 8.5, units = "in") #nothing is significant




