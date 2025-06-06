#Purpose: Plant Level correlation models to understand time (year) and space (lat/long)
#Create simple correlation models for 1) herbivory ~ year + lat
# 2) herbivory ~ year + long. Add in DOY into the models if they will allow it. 
#Author: LAS
#Date started: 03.19.2025
#R version: 4.3.3.

#Load packages
require(rstan) #need rstan for brms package
require(brms) #bay. models
require(tidyverse)
require(bayesplot)
require(shinystan)

tspcorr_brm_models <- function(data, iter = 4000, cores = 2, init = 0, 
                            warmup = 2000, chains = 4, 
                            control = list(adapt_delta = 0.95)) {
  
  # Set seed for reproducibility
  set.seed(888)
  
  # Create a list to store the models
  models <- list()
  
  # Get response variable names from the data
  freq_vars <- grep("^perc.", names(data), value = TRUE) 
  rich_vars <- grep("\\.div$", names(data), value = TRUE)
  perc_vars <- grep("^leafperc_", names(data), value = TRUE)
  chew_vars <- grep("chewvalue", names(data), value = TRUE)  
  
  # Combine all response variables
  response_vars <- c(freq_vars, rich_vars, perc_vars, chew_vars)
  
  # Function to determine the family based on the response variable name
  get_family <- function(response) {
    if (grepl("^perc.", response) || grepl("^leafperc_", response) || grepl("chewvalue", response)) {
      return("zero_inflated_beta")
    } else if (grepl("\\.div$", response)) {
      return("gaussian") 
    }
  }
  
  # Loop through each response variable
  for (response in response_vars) {
    # Create formulas for different model combinations
    formulas <- list(
      year.lat = paste(response, "~ scale(year.x) + scale(startDayOfYear) + scale(lat)"), #might have to remove this if DOY breaks the model 
      #doy = paste(response, "~ scale(startDayOfYear)"),
      #lat = paste(response, "~scale(lat)"),
      year.long = paste(response, "~ scale(year.x) + scale(startDayOfYear) + scale(lon)")
      #interaction = paste(response, "~ scale(year.x) * scale(startDayOfYear)")
    )
    
    # Loop through each formula
    for (model_type in names(formulas)) {
      formula <- as.formula(formulas[[model_type]])
      
      # Determine the family
      family <- get_family(response)
      
      # Fit the model
      model <- brm(formula, data = data, family = family,
                   iter = iter, cores = cores, init = init, 
                   warmup = warmup, chains = chains, control = control)
      
      # Calculate Bayes R2
      r2 <- bayes_R2(model)
      
      # Store the model
      model_name <- paste(response, model_type, sep = "_")
      models[[model_name]] <- model
      
      # Get the model summary
      model_summary <- summary(model)$fixed
      
      # Write the model summary to a file
      write.csv(model_summary,
                file = paste0("Script_output/Brms/time_spac_models/p.lvl/young/", model_name, "_model_summary.csv"))
      
      # Write full model output 
      write.csv(model,
                file = paste0("Script_output/Brms/time_spac_models/p.lvl/young/", model_name, "_model.csv"))
      
      # Write the R2 information in the same directory but with the summary in the filename
      r2_summary <- data.frame(
        Mean_R2 = mean(r2),
        Median_R2 = median(r2),
        Lower_CI = quantile(r2, probs = 0.025),
        Upper_CI = quantile(r2, probs = 0.975),
        row.names = "Bayes_R2"  # Set an explicit row name
      )
      
      write.csv(r2_summary,
                file = paste0("Script_output/Brms/time_spac_models/p.lvl/young/", model_name, "_model_summary_R2.csv"))
    }
  }
  
  return(models)
}

#Load data----
oldleaves <- read_csv("Cleaned_data/reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/reclass.young.leaves.csv")

#Below is the old not reclassified leaf_rel_ages
#QL_leafdata <- read_csv("Cleaned_data/QL_leaflevel.df.csv") #LEAF LEVEL DATA ONLY
#clmdata <- read_csv("Cleaned_data/climate.data.csv") 

#Clean data: Plant Level----
oldleaves2 <- oldleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "perc.dam", "perc.spec", "perc.mine", "perc.gall"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12", "lat", "lon","startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "perc.dam", "perc.spec", "perc.gall", "perc.mine", "total.div", "spec.div", 
           "mine.div", "gall.div")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) #converting percent leaf area measurements to proportions (0:1)


# #converting all 0 = 0.0001 and 1 = 0.9999
oldleaves2$perc.dam <- pmax(0.0001, pmin(0.9999, oldleaves2$perc.dam))
oldleaves2$perc.spec <- pmax(0.0001, pmin(0.9999, oldleaves2$perc.spec))
oldleaves2$perc.mine <- pmax(0.0001, pmin(0.9999, oldleaves2$perc.mine))
oldleaves2$perc.gall <- pmax(0.0001, pmin(0.9999, oldleaves2$perc.gall))

youngleaves2 <- youngleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "perc.dam", "perc.spec", "perc.mine", "perc.gall"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12","lat", "lon","startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "perc.dam", "perc.spec", "perc.gall", "perc.mine", "total.div", "spec.div", 
           "mine.div", "gall.div")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) #converting percent leaf area measurements to proportions (0:1)

# #converting all 0 = 0.0001 and 1 = 0.9999
youngleaves2$perc.dam <- pmax(0.0001, pmin(0.9999, youngleaves2$perc.dam))
youngleaves2$perc.spec <- pmax(0.0001, pmin(0.9999, youngleaves2$perc.spec))
youngleaves2$perc.mine <- pmax(0.0001, pmin(0.9999, youngleaves2$perc.mine))
youngleaves2$perc.gall <- pmax(0.0001, pmin(0.9999, youngleaves2$perc.gall))

#creating a new column of data that is 0:1 for chewing damage to then run models on
oldleaves3 <- oldleaves2
oldleaves3$chewbinary <- oldleaves3$leafperc_area_chew 
oldleaves3 <- oldleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% 
  select(!chewbinary)
oldleaves3$chewvalue <- pmax(0.0001, pmin(0.9999, oldleaves3$chewvalue))
oldleaves3$leafperc_area_chew <- round(oldleaves3$leafperc_area_chew, 4) #rounding percent area damaged columns to only have 4 decimal places
oldleaves3$leafperc_area_gall <- round(oldleaves3$leafperc_area_gall, 4) #this will make all columns standarized
oldleaves3$leafperc_area_mine <- round(oldleaves3$leafperc_area_mine, 4)

#removing columns with incorrect climate data (i.e., -9999.0)
oldleaves3 <- oldleaves3%>% 
  select(where(~if_else(any(. == -9999.0, na.rm = TRUE), FALSE, TRUE)))

#Summarized data (grouped by catalogNumber) for models
oldleaves.sum <- oldleaves3 %>% 
  select(!"leaf_rel_age") %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>% 
  mutate(across(-catalogNumber, ~round(., 4))) %>% 
  column_to_rownames(var = "catalogNumber") %>% 
  mutate(across(c(total.div:gall.div), ~replace_na(., 0)))

youngleaves3 <- youngleaves2
youngleaves3$chewbinary <- youngleaves3$leafperc_area_chew 
youngleaves3 <- youngleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% 
  select(!chewbinary)
youngleaves3$chewvalue <- pmax(0.0001, pmin(0.9999, youngleaves3$chewvalue))
youngleaves3$leafperc_area_chew <- round(youngleaves3$leafperc_area_chew, 4)
youngleaves3$leafperc_area_gall <- round(youngleaves3$leafperc_area_gall, 4)
youngleaves3$leafperc_area_mine <- round(youngleaves3$leafperc_area_mine, 4)

youngleaves3<- youngleaves3%>% 
  select(where(~if_else(any(. == -9999.0, na.rm = TRUE), FALSE, TRUE)))

#Summarized data (grouped by catalogNumber) for models
youngleaves.sum <- youngleaves3 %>%
  select(!"leaf_rel_age") %>% 
  group_by(catalogNumber) %>% 
  summarise_all(mean) %>% 
  mutate(across(-catalogNumber, ~round(., 4))) %>% 
  column_to_rownames(var = "catalogNumber") %>% 
  mutate(across(c(total.div:gall.div), ~replace_na(., 0)))

#recombine old and young datasets to make a full dataset that is the same number of leaves as the seperated datasets
leafdata2 <- youngleaves.sum %>% 
  rbind(oldleaves.sum)

oldleaves4 <- oldleaves.sum %>% 
  select(c(starts_with("perc."), ends_with(".div"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear",
           "lon", "lat"))

youngleaves4 <- youngleaves.sum %>% 
  select(c(starts_with("perc."), ends_with(".div"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear",
           "lat", "lon"))

# Using function----
tspcorr_oldmodels <- tspcorr_brm_models(oldleaves4)
tspcorr_youngmodels <- tspcorr_brm_models(youngleaves4)






