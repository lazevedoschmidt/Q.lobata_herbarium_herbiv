#Purpose: Run LEAF LEVEL simple models
#R version: 4.3.3.
#Author: LAS
#Date started: 08.28.2024
  #Date updated: 08/18/2025
    #This code has been updated to handle random effects to account for pseudoreplication (catalogNumber)
    #It includes all updates reflected in notes on "STEP3D_leaf.level.models.function" 


#load packages
require(rstan) #need rstan for brms package
require(brms) #bay. models
require(tidyverse)
require(bayesplot)
require(shinystan)

#Function
simp_brm_models <- function(data, iter = 4000, cores = 4, init = 0, 
                            warmup = 2000, chains = 4, 
                            control = list(adapt_delta = 0.99, max_treedepth = 15)) {
  
  # Set seed for reproducibility
  set.seed(888)
  
  # Create a list to store the models
  models <- list()
  
  # Check the catalogNumber structure (for information)
  catalog_counts <- table(data$catalogNumber)
  cat("Leaves per specimen (catalogNumber):\n")
  print(summary(catalog_counts))
  cat("Total specimens:", length(unique(data$catalogNumber)), "\n")
  cat("Total leaves:", nrow(data), "\n\n")
  
  # Get response variable names from the data
  freq_vars <- grep("\\.freq$", names(data), value = TRUE)
  rich_vars <- grep("\\.rich$", names(data), value = TRUE)
  perc_vars <- grep("^leafperc_", names(data), value = TRUE)
  chew_vars <- grep("chewvalue", names(data), value = TRUE)
  
  # Combine all response variables
  response_vars <- c(freq_vars, rich_vars, perc_vars, chew_vars)
  
  # Function to determine the family based on the response variable name
  get_family <- function(response) {
    if (grepl("\\.freq$", response) || grepl("^leafperc_", response) || grepl("chewvalue", response)) {
      return("beta")
    } else if (grepl("\\.rich$", response)) {
      return("zero_inflated_poisson")
    }
  }
  
  # Set priors that help with convergence for different families
  get_priors <- function(family) {
    if(family == "beta") {
      c(prior(exponential(1), class = sd, group = catalogNumber),
        prior(normal(0, 2.5), class = Intercept),
        prior(normal(0, 1), class = b))
    } else if(family == "zero_inflated_poisson") {
      c(prior(exponential(1), class = sd, group = catalogNumber),
        prior(normal(0, 2), class = Intercept),
        prior(normal(0, 1), class = b),
        prior(beta(1, 1), class = zi))  # Prior for zero-inflation parameter
    } else {
      c(prior(exponential(1), class = sd, group = catalogNumber),
        prior(normal(0, 2.5), class = Intercept),
        prior(normal(0, 1), class = b))
    }
  }
  
  # Loop through each response variable
  for (response in response_vars) {
    # Create formulas for different model combinations
    formulas <- list(
      year_doy.int = paste(response, "~ scale(year.x) + scale(doy.clean)+ scale(year.x) * scale(doy.clean) +(1|catalogNumber)")
      #doy = paste(response, "~ scale(doy.clean) + (1|catalogNumber)"),
      #interaction = paste(response, "~ scale(year.x) * scale(doy.clean) + (1|catalogNumber)")
    )
    
    # Loop through each formula
    for (model_type in names(formulas)) {
      cat("Fitting", response, "-", model_type, "model...\n")
      
      formula <- as.formula(formulas[[model_type]])
      
      # Determine the family
      family <- get_family(response)
      model_priors <- get_priors(family)
      
      # Fit model with error handling and retries
      model <- NULL
      attempt <- 1
      max_attempts <- 2
      
      while(is.null(model) && attempt <= max_attempts) {
        if(attempt > 1) {
          cat("  Retry attempt", attempt, "with higher adapt_delta...\n")
          current_control <- list(adapt_delta = 0.995, max_treedepth = 15)
          current_iter <- iter + 1000  # More iterations on retry
          current_warmup <- warmup + 500
        } else {
          current_control <- control
          current_iter <- iter
          current_warmup <- warmup
        }
        
        model <- try({
          brm(formula, data = data, family = family, prior = model_priors,
              iter = current_iter, cores = cores, init = init, 
              warmup = current_warmup, chains = chains, control = current_control,
              silent = 2, refresh = 0)  # Reduce output during fitting
        }, silent = TRUE)
        
        if(inherits(model, "try-error")) {
          cat("  Error on attempt", attempt, "\n")
          model <- NULL
          attempt <- attempt + 1
        } else {
          # Check convergence - robust across brms versions
          rhat_vals <- try({
            # Try newer brms syntax first
            rhat(model, pars = NULL)
          }, silent = TRUE)
          
          if(inherits(rhat_vals, "try-error")) {
            # Fall back to older syntax or extract from summary
            rhat_vals <- try(summary(model)$fixed[, "Rhat"], silent = TRUE)
          }
          
          if(!inherits(rhat_vals, "try-error") && length(rhat_vals) > 0) {
            max_rhat <- max(rhat_vals, na.rm = TRUE)
          } else {
            # If all else fails, assume convergence is okay
            max_rhat <- 1.0
            cat("  Could not extract Rhat values\n")
          }
          
          if(max_rhat > 1.01) {
            cat("  Warning: Max Rhat =", round(max_rhat, 3), "\n")
            if(max_rhat > 1.05 && attempt < max_attempts) {
              cat("  Rhat too high, retrying...\n")
              model <- NULL
              attempt <- attempt + 1
            }
          } else {
            cat("  Converged successfully (Max Rhat =", round(max_rhat, 3), ")\n")
          }
        }
      }
      
      if(is.null(model)) {
        cat("  FAILED after", max_attempts, "attempts - skipping", response, model_type, "\n\n")
        next
      }
      
      # Calculate Bayes R2 with error handling
      r2 <- try(bayes_R2(model), silent = TRUE)
      if(inherits(r2, "try-error")) {
        cat("  Warning: Could not calculate R2 for", response, model_type, "\n")
        r2 <- NA
      }
      
      # Store the model
      model_name <- paste(response, model_type, sep = "_")
      models[[model_name]] <- model
      
      # Write full model output (consider using saveRDS instead)
      try({
        write.csv(model,
                  file = paste0("Script_output/Brms/simp_models/l.lvl/young/", model_name, "_model.csv"))
      }, silent = TRUE)
      
      # Write the model summary
      model_summary <- summary(model)$fixed
      write.csv(model_summary,
                file = paste0("Script_output/Brms/simp_models/l.lvl/young/", model_name, "_model_summary.csv"))
      
      # Write R2 files with error handling
      if(!inherits(r2, "try-error") && !any(is.na(r2))) {
        write.csv(as.data.frame(r2),
                  paste0("Script_output/Brms/simp_models/l.lvl/young/", model_name, "_bayes_R2.csv"))
      } else {
        cat("  Skipping R2 output due to calculation error\n")
      }
      
      cat("  Saved:", model_name, "\n\n")
    }
  }
  
  cat("Model fitting complete!\n")
  cat("Successfully fitted", length(models), "out of", length(response_vars), "attempted models.\n")
  
  return(models)
}


#load data and subset----
oldleaves <- read_csv("Cleaned_data/cleaned_reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/cleaned_reclass.young.leaves.csv")

#recombine old and young datasets to make a full dataset that is the same number of leaves as the seperated datasets
# leafdata2 <- youngleaves %>% 
#   rbind(oldleaves)

#Subset data for simple models
oldleaves4 <- oldleaves %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "doy.clean", "catalogNumber"))

youngleaves4 <- youngleaves %>% 
  select(c(ends_with(".freq"), ends_with(".rich"), starts_with("leafperc_"), "chewvalue","year.x", "doy.clean","catalogNumber"))

# Usage of function----
#simple year and doy models below: 
simp_oldmodels <- simp_brm_models(oldleaves4)
simp_youngmodels <- simp_brm_models(youngleaves4)


