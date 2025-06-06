#Purpose: LEAF level models with LEAF level PCA importance
  #This uses a function to import the PCA importance files and run through each brm()
  #for each unique response and predictor variable. Then, it writes a .csv file for the summary
  #output and the full model output. 
#R version: 4.3.3.
#Aurhors: LAS 

#load packages
require(rstan) #need rstan for brms package
require(brms) #bay. models
require(tidyverse)
require(bayesplot)
require(shinystan)


#Leaf Level brms() function
#functions----
#year and doy models only. No climate variables
simp_brm_models <- function(data, iter = 2000, cores = 2, init = 0, 
                              warmup = 200, chains = 2, 
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
      year = paste(response, "~ scale(year.x)"),
      doy = paste(response, "~ scale(startDayOfYear)"),
      interaction = paste(response, "~ scale(year.x) * scale(startDayOfYear)")
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
      r2 <- bayes_R2(model) #Added on 02/13/2025 to give us the R2 for each model
      
      # Store the model
      model_name <- paste(response, model_type, sep = "_")
      models[[model_name]] <- model
      
      
      # Write full model output
      write.csv(model,
                file = paste0("Script_output/Brms/simp_models/old/", model_name, "_model.csv"))

      # Write the model summary
      write.csv(summary(model)$fixed,
                file = paste0("Script_output/Brms/simp_models/old/", model_name, "_model_summary.csv"))

      #Bays R2 files
      write.csv(as.data.frame(r2),
                paste0("Script_output/Brms/simp_models/old/", model_name, "_bayes_R2.csv")) #here too too (added 02/13/2025)
    }
  }
  
  return(models)
}

##########################
#climate variable models 
create_brm_models <- function(data_file, response_data, additional_predictors = NULL, 
                              iter = 2000, cores = 2, init = 0, 
                              warmup = 200, chains = 2, 
                              control = list(adapt_delta = 0.95)) {
  
  #set seed
  set.seed(888)
  
  # Read the CSV file
  data <- read.csv(data_file)
  
  # Create a list to store the models
  models <- list()
  
  # Function to determine the family based on the response variable name
  get_family <- function(response) {
    if (grepl("^leaf.*_area_", response) || grepl("\\.freq$", response) ||  grepl("chewvalue", response)) {
      return("zero_inflated_beta")
    } else if (grepl("\\.rich$", response)) {
      return("zero_inflated_poisson")
    } else {
      stop(paste("Unknown family for response variable:", response))
    }
  }
  
  # Get unique response variables
  response_vars <- unique(data$Response)
  
  # Define model types
  model_types <- list(
    year = function(response) paste(response, "~ scale(year.x)"),
    doy = function(response) paste(response, "~ scale(startDayOfYear)") #add comma back if the interaction is wanted
    #interaction = function(response) paste(response, "~ scale(year.x) * scale(startDayOfYear)")
  )
  
  # Loop through each response variable
  for (response in response_vars) {
    # Get predictor variables for the current response
    predictors <- data %>%
      filter(Response == response) %>%
      pull(Variable)
    
    # Combine with additional predictors if provided
    all_predictors <- unique(c(predictors, additional_predictors))
    
    # Loop through each predictor for the current response
    for (predictor in all_predictors) {
      # Loop through each model type
      for (model_type_name in names(model_types)) {
        # Create the formula using the appropriate model type function
        base_formula <- model_types[[model_type_name]](response)
        
        # Add the predictor to the formula if it's not already part of the model type
        if (!grepl(predictor, base_formula, fixed = TRUE)) {
          formula_str <- paste(base_formula, "+", paste0("scale(", predictor, ")"))
        } else {
          formula_str <- base_formula
        }
        
        formula <- as.formula(formula_str)
        
        # Determine the family
        family <- get_family(response)
        
        # Fit the model
        model <- brm(formula, data = response_data, family = family,  
                     iter = iter, cores = cores, init = init, warmup = warmup,
                     chains = chains, control = control)
        
        # Calculate Bayes R2
        r2 <- bayes_R2(model)
        
        # Store the model
        model_name <- paste(response, predictor, model_type_name, sep = "_")
        models[[model_name]] <- model
        
        # Create output directory if it doesn't exist
        output_dir <- paste0("Script_output/Brms/climate_models/l.lvl/young/", model_type_name, "/") #switch directory between old/young runs
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Write full model output
        write.csv(model, paste0(output_dir, model_name, "_model.csv"))
        
        # Write the model summary to a CSV file
        write.csv(summary(model)$fixed, 
                  paste0(output_dir, model_name, "_model.summary.csv"))
        
        # Write the Bayes R2 to a CSV file
        write.csv(as.data.frame(r2), 
                  paste0(output_dir, model_name, "_bayes_R2.csv"))
      }
    }
  }
  
  return(models)
}

#load data and subset----
oldleaves <- read_csv("Cleaned_data/reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/reclass.young.leaves.csv")

oldleaves2 <- oldleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12", "startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
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
  select(c("MAT":"CMI_12","startDayOfYear","catalogNumber", "leaf_rel_age","year.x","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
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

#removing columns with incorrect climate data (i.e., -9999.0)
oldleaves3 <- oldleaves3%>% 
  select(where(~if_else(any(. == -9999.0, na.rm = TRUE), FALSE, TRUE)))

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


#recombine old and young datasets to make a full dataset that is the same number of leaves as the seperated datasets
leafdata2 <- youngleaves3 %>% 
  rbind(oldleaves3)

#Subset data for simple models
oldleaves4 <- oldleaves3 %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear"))

youngleaves4 <- youngleaves3 %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "startDayOfYear"))

# Usage of function----
#simple year and doy models below: 
simp_oldmodels <- simp_brm_models(oldleaves4)
simp_youngmodels <- simp_brm_models(youngleaves4)

 #Climate variable models below: 
#additional_predictors <- c("year.x") #can update this to run with any additional predictor

#prior to running these models, edit the function to not include "leaf_rel_age" in the models because they aren't appropriate
  #for these (already age specific) models. 
#ALSO make sure you update the data frame call and how the summary and model outputs are being saved (i.e., add the age qualifier)
oldYmodels <- create_brm_models("Script_output/PCA/l.lvl/old/Ltop5oldY.csv", oldleaves3)

youngYmodels <- create_brm_models("Script_output/PCA/l.lvl/young/Ltop5youngY.csv", youngleaves3)

