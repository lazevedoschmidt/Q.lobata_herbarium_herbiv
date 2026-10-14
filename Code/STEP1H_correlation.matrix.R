#Purpose: Checking correlation of temporal and spatial variables
#Date started: 11.20.2024
#Date updated:03/10/2026
  #Adding addition correlation tests for climate variables and year
#Date updated: 03/31/2025
  #Ran correlation matrices for time/space/climate variables
#R version: 4.3.3.
#author: LAS

#Load packages
require(tidyverse)
require(lubridate) #used for splitting date information into separate columns 
require(ggplot2)
require(pwr) #needed for pwr.t.test function
require(corrplot)
require(arsenal) #need for comparef() 
library(knitr) #need these two for making the SI tables at the end for climate and year
library(kableExtra)

oldleaves <- read_csv("Cleaned_data/reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/reclass.young.leaves.csv")

#Correlation of occurrence data:----
oldQL_corr_matrix <- oldleaves %>% 
  select(c(catalogNumber, year.x, month.x, day.x, doy.clean,
           decimalLatitude, decimalLongitude)) %>%  # selecting the columns I'm interested in
  group_by(catalogNumber) %>%
  summarise(across(everything(), first)) %>%  # Use summarise() instead of mutate() with summarise_all()
  remove_rownames %>%  # removing rownames
  column_to_rownames(var = "catalogNumber")

M_o <-   cor(oldQL_corr_matrix) #creating correlation matrix to then plot 
old_QL_cor.plot <- corrplot(M_o, method = 'number') # colorful number
cor.test(oldQL_corr_matrix$year.x, oldQL_corr_matrix$doy.clean)
cor.test(oldQL_corr_matrix$decimalLatitude, oldQL_corr_matrix$decimalLongitude)
#save output

youngQL_corr_matrix <- youngleaves %>% 
  select(c(catalogNumber, year.x, month.x, day.x, doy.clean,
           decimalLatitude, decimalLongitude)) %>%  # selecting the columns I'm interested in
  group_by(catalogNumber) %>%
  summarise(across(everything(), first)) %>%  # Use summarise() instead of mutate() with summarise_all()
  remove_rownames %>%  # removing rownames
  column_to_rownames(var = "catalogNumber")

M_y <-   cor(youngQL_corr_matrix) #creating correlation matrix to then plot 
young_QL_cor.plot <- corrplot(M_y, method = 'number') # colorful number
cor.test(youngQL_corr_matrix$year.x, youngQL_corr_matrix$doy.clean)
cor.test(youngQL_corr_matrix$decimalLatitude, youngQL_corr_matrix$decimalLongitude)
#save output

#time/space/climate corr----
oldtspcorr <- oldleaves %>% 
  select(c(catalogNumber:DD1040,year.x)) %>% 
  select(!c(el, DD_0, PAS, EMT, FFP, MCMT)) %>% 
  group_by(catalogNumber) %>% 
  summarise(across(everything(), first)) %>%  # Use summarise() instead of mutate() with summarise_all()
  remove_rownames %>%  # removing rownames
  column_to_rownames(var = "catalogNumber")

M_otsp <- cor(oldtspcorr)
otsp.plot <- corrplot(M_otsp, method = 'number', type = 'lower', diag = FALSE)
otsp.plot + title("Spatiotemporal Correlation Old Leaves")

youngtspcorr <- youngleaves %>% 
  select(c(catalogNumber:DD1040,year.x)) %>% 
  select(!c(el, DD_0, PAS, EMT, FFP, MCMT)) %>% 
  group_by(catalogNumber) %>% 
  summarise(across(everything(), first)) %>%  # Use summarise() instead of mutate() with summarise_all()
  remove_rownames %>%  # removing rownames
  column_to_rownames(var = "catalogNumber")

M_ytsp <- cor(youngtspcorr)
ytsp.plot <- corrplot(M_ytsp, method = 'number', type = 'lower', diag = FALSE)
ytsp.plot + title("Spatiotemporal Correlation Young Leaves")

####Testing correlation between year and climate variables###
# Year/Climate correlation with p-values ----
plot_year_cors <- function(df, label) {
  
  yr_clim <- df %>%
    select(catalogNumber, MAT:CMI_12, year.x) %>%
    group_by(catalogNumber) %>%
    summarise(across(everything(), first)) %>%
    ungroup() %>%
    distinct(catalogNumber, .keep_all = TRUE) %>%
    remove_rownames %>%
    column_to_rownames(var = "catalogNumber")
  
  message(paste("Running correlations on", nrow(yr_clim), "unique specimens for", label))
  
  climate_vars <- colnames(yr_clim)[colnames(yr_clim) != "year.x"]
  
  results <- map_dfr(climate_vars, function(var) {
    x <- yr_clim[[var]]
    y <- yr_clim[["year.x"]]
    
    complete <- complete.cases(x, y)
    if (sum(complete) < 3) return(NULL)
    
    # Skip variables with near-zero standard deviation
    if (sd(x[complete], na.rm = TRUE) < 1e-10) {
      message(paste("Skipping", var, "- no variation"))
      return(NULL)
    }
    
    test <- cor.test(x[complete], y[complete], method = "pearson")
    
    data.frame(
      variable = var,
      r        = round(test$estimate, 3),
      p_value  = round(test$p.value, 4),
      sig      = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ ".",
        TRUE                 ~ "ns"
      )
    )
  })
  
  results_sig <- results %>% 
    filter(p_value < 0.05) %>% 
    arrange(r)
  
  par(mar = c(4, 8, 4, 2))
  barplot(results_sig$r,
          names.arg = paste0(results_sig$variable, " ", results_sig$sig),
          horiz = TRUE,
          las = 1,
          col = ifelse(results_sig$r > 0, "steelblue", "tomato"),
          main = paste("Variables Significantly Correlated with Year (p < 0.05)\n", label),
          xlab = "Pearson r",
          xlim = c(-0.6, 0.6),
          cex.names = 0.75)
  abline(v = 0, lty = 2, col = "grey40")
  
  return(results)
}

## Run for both datasets
old_yr_results   <- plot_year_cors(oldleaves,   "Old Leaves")
young_yr_results <- plot_year_cors(youngleaves, "Young Leaves")

## Print full tables
print(old_yr_results)
print(young_yr_results)

##Making tables that can be saved----
save_year_cor_table <- function(old_results, young_results, filename = "Table_S1_year_climate_correlations.html") {
  
  # Format old leaves table
  old_table <- old_results %>%
    arrange(p_value) %>%
    mutate(
      Dataset = "Old Leaves",
      r = sprintf("%.3f", r),
      p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))
    )
  
  # Format young leaves table
  young_table <- young_results %>%
    arrange(p_value) %>%
    mutate(
      Dataset = "Young Leaves",
      r = sprintf("%.3f", r),
      p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))
    )
  
  # Combine and select columns
  combined <- bind_rows(old_table, young_table) %>%
    select(Dataset, variable, r, p_value, sig) %>%
    rename(
      "Variable"     = variable,
      "Pearson r"    = r,
      "P-value"      = p_value,
      "Significance" = sig
    )
  
  # Row indices for each dataset (for grouping)
  n_old   <- nrow(old_table)
  n_young <- nrow(young_table)
  
  # Build table
  tbl <- combined %>%
    kbl(
      caption  = "Table S1. Pearson correlations between collection year and climate variables. Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1, ns = not significant.",
      align    = c("l", "l", "c", "c", "c"),
      booktabs = TRUE
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "bordered"),
      full_width        = FALSE,
      font_size         = 12
    ) %>%
    # Group rows by dataset
    pack_rows("Old Leaves",   1,       n_old) %>%
    pack_rows("Young Leaves", n_old+1, n_old + n_young) %>%
    # Bold significant rows
    row_spec(which(combined$Significance %in% c("*", "**", "***")),
             bold = TRUE) %>%
    # Style header
    row_spec(0, bold = TRUE, background = "#2E75B6", color = "white")
  
  # Save as html
  save_kable(tbl, file = filename, self_contained = TRUE)
  message(paste("Table saved to:", filename))
  
  return(tbl)
}

# Run it - saves as HTML by default
save_year_cor_table(old_yr_results, young_yr_results)

# Or save as PDF (requires webshot2 package)
# save_year_cor_table(old_yr_results, young_yr_results,
#                     filename = "Table_S1_year_climate_correlations.pdf")

##### Only variables for the models we ran----
# Define variables for each dataset
young_climate_vars <- c("AHM", "MAR", "RH", "MSP", "SHM", "EXT")
old_climate_vars   <- c("AHM", "MAR", "RH", "MSP", "SHM", "EXT", "MAP", "TD")

# Correlation + p-values function for specific variables ----
get_specific_cor <- function(df, climate_vars, label) {
  
  yr_clim <- df %>%
    select(catalogNumber, all_of(climate_vars), year.x) %>%
    group_by(catalogNumber) %>%
    summarise(across(everything(), first)) %>%
    ungroup() %>%
    distinct(catalogNumber, .keep_all = TRUE) %>%
    remove_rownames %>%
    column_to_rownames(var = "catalogNumber")
  
  message(paste("Running correlations on", nrow(yr_clim), "unique specimens for", label))
  
  results <- map_dfr(climate_vars, function(var) {
    x <- yr_clim[[var]]
    y <- yr_clim[["year.x"]]
    
    complete <- complete.cases(x, y)
    if (sum(complete) < 3) return(NULL)
    
    if (sd(x[complete], na.rm = TRUE) < 1e-10) {
      message(paste("Skipping", var, "- no variation"))
      return(NULL)
    }
    
    test <- cor.test(x[complete], y[complete], method = "pearson")
    
    data.frame(
      variable = var,
      r        = round(test$estimate, 3),
      p_value  = round(test$p.value, 4),
      sig      = case_when(
        test$p.value < 0.001 ~ "***",
        test$p.value < 0.01  ~ "**",
        test$p.value < 0.05  ~ "*",
        test$p.value < 0.1   ~ ".",
        TRUE                 ~ "ns"
      )
    )
  })
  
  return(results)
}

# Run for both datasets
old_key_results   <- get_specific_cor(oldleaves,   old_climate_vars,   "Old Leaves")
young_key_results <- get_specific_cor(youngleaves, young_climate_vars, "Young Leaves")

# Print results
print(old_key_results)
print(young_key_results)

# Bar plots ----
plot_specific_cors <- function(results, label) {
  results_sorted <- results %>% arrange(r)
  
  par(mar = c(4, 8, 4, 2))
  barplot(results_sorted$r,
          names.arg = paste0(results_sorted$variable, " ", results_sorted$sig),
          horiz     = TRUE,
          las       = 1,
          col       = ifelse(results_sorted$r > 0, "steelblue", "tomato"),
          main      = paste("Year-Climate Correlation\n", label),
          xlab      = "Pearson r",
          xlim      = c(-0.6, 0.6),
          cex.names = 0.9)
  abline(v = 0, lty = 2, col = "grey40")
}

plot_specific_cors(old_key_results,   "Old Leaves")
plot_specific_cors(young_key_results, "Young Leaves")

# Supplementary table ----
save_key_cor_table <- function(old_results, young_results, 
                               filename = "Table_S_key_climate_year_correlations.html") {
  
  old_table <- old_results %>%
    arrange(p_value) %>%
    mutate(
      Dataset = "Late-Season Leaves",
      r       = sprintf("%.3f", r),
      p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))
    )
  
  young_table <- young_results %>%
    arrange(p_value) %>%
    mutate(
      Dataset = "Early-Season Leaves",
      r       = sprintf("%.3f", r),
      p_value = ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))
    )
  
  combined <- bind_rows(old_table, young_table) %>%
    select(Dataset, variable, r, p_value, sig) %>%
    rename(
      "Climate Variable \n(correlated with year)"     = variable,
      "Pearson r"    = r,
      "P-value"      = p_value,
      "Significance" = sig
    )
  
  n_old   <- nrow(old_table)
  n_young <- nrow(young_table)
  
  tbl <- combined %>%
    kbl(
      #caption  = "Table S5. Pearson correlations between collection year and key climate predictors of insect herbivory. Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1, ns = not significant.",
      align    = c("l", "l", "c", "c", "c"),
      booktabs = TRUE
    ) %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "bordered"),
      full_width        = FALSE,
      font_size         = 12
    ) %>%
    pack_rows("Late-Season Leaves",   1,       n_old) %>%
    pack_rows("Early-Season Leaves", n_old+1, n_old + n_young) %>%
    row_spec(which(combined$Significance %in% c("*", "**", "***")),
             bold = TRUE) %>%
    row_spec(0, bold = TRUE, background = "#2E75B6", color = "white")
  
  save_kable(tbl, file = filename, self_contained = TRUE)
  message(paste("Table saved to:", filename))
  
  return(tbl)
}

save_key_cor_table(old_key_results, young_key_results)


