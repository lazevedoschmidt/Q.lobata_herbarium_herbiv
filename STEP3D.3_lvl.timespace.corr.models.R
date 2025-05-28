#Purpose: Leaf Level correlation models to understand time (year) and space (lat/long)
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

#Function----
tspcorr_brm_models <- function(data, iter = 4000, cores = 2, init = 0, 
                               warmup = 2000, chains = 4, 
                               control = list(adapt_delta = 0.95)) {
  
  # Set seed for reproducibility
  set.seed(888)
  
  # Create a list to store the models
  models <- list()
  
  # Get response variable names from the data
  freq_vars <- grep("\\.freq$", names(data), value = TRUE)
  rich_vars <- grep("\\.rich$", names(data), value = TRUE)
  perc_vars <- grep("^leafperc_", names(data), value = TRUE)
  chew_vars <- grep("chewvalue", names(data), value = TRUE)  # Added chew value
  
  # Combine all response variables
  response_vars <- c(freq_vars, rich_vars, perc_vars, chew_vars)
  
  # Function to determine the family based on the response variable name
  get_family <- function(response) {
    if (grepl("\\.freq$", response) || grepl("^leafperc_", response) || grepl("chewvalue", response)) {
      return("zero_inflated_beta")
    } else if (grepl("\\.rich$", response)) {
      return("zero_inflated_poisson")
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
                file = paste0("Script_output/Brms/time_spac_models/l.lvl/young/", model_name, "_model_summary.csv"))
      
      # Write full model output 
      write.csv(model,
                file = paste0("Script_output/Brms/time_spac_models/l.lvl/young/", model_name, "_model.csv"))
      
      # Write the R2 information in the same directory but with the summary in the filename
      r2_summary <- data.frame(
        Mean_R2 = mean(r2),
        Median_R2 = median(r2),
        Lower_CI = quantile(r2, probs = 0.025),
        Upper_CI = quantile(r2, probs = 0.975),
        row.names = "Bayes_R2"  # Set an explicit row name
      )
      
      write.csv(r2_summary,
                file = paste0("Script_output/Brms/time_spac_models/l.lvl/young/", model_name, "_model_summary_R2.csv"))
    }
  }
  
  return(models)
}

#Load leaf level data----
oldleaves <- read_csv("Cleaned_data/reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/reclass.young.leaves.csv")

#Below is the old not reclassified leaf_rel_ages
#QL_leafdata <- read_csv("Cleaned_data/QL_leaflevel.df.csv") #LEAF LEVEL DATA ONLY
#clmdata <- read_csv("Cleaned_data/climate.data.csv") 

#Clean data----
oldleaves2 <- oldleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12", "lat", "lon","startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) #converting percent leaf area measurements to proportions (0:1)

# #converting all 0 = 0.0001 and 1 = 0.9999
oldleaves2$leaftotal.freq <- pmax(0.0001, pmin(0.9999, oldleaves2$leaftotal.freq))
oldleaves2$leafspec.freq <- pmax(0.0001, pmin(0.9999, oldleaves2$leafspec.freq))
oldleaves2$leafmine.freq <- pmax(0.0001, pmin(0.9999, oldleaves2$leafmine.freq))
oldleaves2$leafgall.freq <- pmax(0.0001, pmin(0.9999, oldleaves2$leafgall.freq))

youngleaves2 <- youngleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12","lat", "lon","startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) %>% 
  mutate(across(starts_with("leafperc_area_"), ~ . / 100)) #converting percent leaf area measurements to proportions (0:1)
# #converting all 0 = 0.0001 and 1 = 0.9999
youngleaves2$leaftotal.freq <- pmax(0.0001, pmin(0.9999, youngleaves2$leaftotal.freq))
youngleaves2$leafspec.freq <- pmax(0.0001, pmin(0.9999, youngleaves2$leafspec.freq))
youngleaves2$leafmine.freq <- pmax(0.0001, pmin(0.9999, youngleaves2$leafmine.freq))
youngleaves2$leafgall.freq <- pmax(0.0001, pmin(0.9999, youngleaves2$leafgall.freq))

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


#recombine old and young datasets to make a full dataset that is the same number of leaves as the seperated datasets
leafdata2 <- youngleaves3 %>% 
  rbind(oldleaves3)

#Subset data for simple correlation models
oldleaves4 <- oldleaves3 %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear", 
           "lat", "lon"))

youngleaves4 <- youngleaves3 %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear",
           "lat", "lon"))

# Using function----
tspcorr_oldmodels <- tspcorr_brm_models(oldleaves4)
tspcorr_youngmodels <- tspcorr_brm_models(youngleaves4)


