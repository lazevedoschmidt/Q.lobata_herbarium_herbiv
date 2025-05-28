#Purpose: Variable selection for future models
#Using PCA instead of random forests which have now been abandoned (RIP 07/22/2024)
#How output is used: 
  #PCA importance at the leaf level are used to justify/explain which climate variables are used in the 
  #bayesian models in future steps. The output here is then fed into step3D
#Author: LAS
#R version: 4.2.2.
#Date started: 07/23/2024
  #Date last revised: 10/22/2024

library(stats)
library(caret)
library(tidyverse)
library(reshape2)

#Functions
#PCA function----
create_multiple_pca_models <- function(data, n_components = NULL, na.action = "omit") {
  # Identify response variables (columns that begin with "perc.", "perc_area", or end with ".div")
  #response_vars <- grep("^ffg", names(data), value = TRUE)
  response_vars <- grep("^leaf|^chew|^ffg", names(data), value = TRUE)
  # Identify predictor variables (all other columns except catalogNumber)
  predictor_vars <- setdiff(names(data), c(response_vars, "catalogNumber"))
  
  # Create a list to store the PCA results
  pca_results <- list()
  
  # Handle NA values based on the specified method
  if (na.action == "omit") {
    data <- na.omit(data)
  } else if (na.action == "impute") {
    # Simple imputation: replace NA with column mean for numeric, mode for categorical
    for (col in names(data)) {
      if (is.numeric(data[[col]])) {
        data[[col]][is.na(data[[col]])] <- mean(data[[col]], na.rm = TRUE)
      } else {
        mode_val <- names(sort(table(data[[col]]), decreasing = TRUE))[1]
        data[[col]][is.na(data[[col]])] <- mode_val
      }
    }
  } #Simple imputation that takes the mean or if it's categorical, take the most common value. Can't just drop NAs because it would remove the entire row so this is our work around
  
  # Perform PCA on predictor variables
  #set seed
  set.seed(888)
  X <- data[, predictor_vars, drop = FALSE]
  
  # Store original variable names before removing constant columns
  original_var_names <- colnames(X)
  
  # Remove constant columns
  non_constant_cols <- apply(X, 2, function(col) length(unique(col)) > 1)
  X <- X[, non_constant_cols, drop = FALSE]
  #removing columns where every value is the same (double check this)
  
  # Store the names of variables actually used in PCA
  used_var_names <- colnames(X)
  
  # Check if there are any columns left after removing constant ones
  if (ncol(X) == 0) {
    stop("All predictor variables have zero variance.")
  }
  
  X_scaled <- scale(X)
  
  # Perform PCA
  #Set seed
  set.seed(888)
  pca_model <- prcomp(X_scaled, center = TRUE, scale. = TRUE)
  
  # Determine number of components to keep
  max_components <- min(nrow(X) - 1, ncol(X))
  if (is.null(n_components)) {
    n_components <- max_components
  } else {
    n_components <- min(n_components, max_components)
  }
  
  # Extract PCA components
  pca_components <- pca_model$x[, 1:n_components, drop = FALSE]
  
  # Calculate explained variance ratio
  explained_variance_ratio <- pca_model$sdev^2 / sum(pca_model$sdev^2)
  explained_variance_ratio <- explained_variance_ratio[1:n_components]
  
  # Calculate cumulative explained variance ratio
  cumulative_variance_ratio <- cumsum(explained_variance_ratio)
  
  # Store general PCA results
  pca_results$general <- list(
    pca_components = pca_components,
    explained_variance_ratio = explained_variance_ratio,
    cumulative_variance_ratio = cumulative_variance_ratio,
    feature_importance = pca_model$rotation[, 1:n_components, drop = FALSE],
    removed_constant_columns = setdiff(original_var_names, used_var_names),
    n_components_used = n_components,
    variable_names = used_var_names  # Store the names of variables actually used in PCA
  )
  
  # Analyze each response variable
  for (response_var in response_vars) {
    tryCatch({
      # Check if response variable has zero variance
      if (length(unique(data[[response_var]])) == 1) {
        warning(paste("Response variable", response_var, "has zero variance. Skipping correlation calculation."))
        correlations <- rep(NA, n_components)
      } else {
        # Calculate correlations between PCA components and the response variable
        correlations <- cor(pca_components, data[[response_var]], use = "complete.obs")
      }
      
      # Store results for this response variable
      pca_results[[response_var]] <- list(
        correlations_with_response = correlations
      )
      
      cat("PCA analysis completed for:", response_var, "\n")
    }, error = function(e) {
      warning(paste("Error analyzing", response_var, ":", e$message))
    })
  }
  
  pca_results$predictor_variables <- predictor_vars
  pca_results$response_variables <- response_vars
  
  return(pca_results)
}

#creating .csv function----
generate_pca_csv_output <- function(pca_results) {
  # Initialize an empty list to store all rows
  all_rows <- list()
  
  # Extract necessary information
  loadings <- pca_results$general$feature_importance
  var_names <- pca_results$general$variable_names
  n_components <- ncol(loadings)
  
  # Add loadings
  for (i in 1:nrow(loadings)) {
    row <- c(Variable = var_names[i],
             Type = "Predictor",
             setNames(as.list(loadings[i,]), paste0("PC", 1:n_components)),
             Value = NA)
    all_rows[[length(all_rows) + 1]] <- row
  }
  
  # Add explained variance
  for (i in 1:n_components) {
    row <- c(Variable = paste0("Explained_Variance_PC", i),
             Type = "Metadata",
             setNames(as.list(rep(NA, n_components)), paste0("PC", 1:n_components)),
             Value = pca_results$general$explained_variance_ratio[i])
    all_rows[[length(all_rows) + 1]] <- row
  }
  
  # Add cumulative explained variance
  for (i in 1:n_components) {
    row <- c(Variable = paste0("Cumulative_Explained_Variance_PC", i),
             Type = "Metadata",
             setNames(as.list(rep(NA, n_components)), paste0("PC", 1:n_components)),
             Value = pca_results$general$cumulative_variance_ratio[i])
    all_rows[[length(all_rows) + 1]] <- row
  }
  
  # Add correlations with response variables
  for (response_var in pca_results$response_variables) {
    if (!is.null(pca_results[[response_var]]$correlations_with_response)) {
      correlations <- pca_results[[response_var]]$correlations_with_response
      for (i in 1:length(correlations)) {
        row <- c(Variable = paste0(response_var, "_Correlation_PC", i),
                 Type = "Response",
                 setNames(as.list(rep(NA, n_components)), paste0("PC", 1:n_components)),
                 Value = correlations[i])
        all_rows[[length(all_rows) + 1]] <- row
      }
    }
  }
  
  # Convert list of rows to data frame
  final_df <- do.call(rbind, lapply(all_rows, data.frame))
  
  return(final_df)
}

#load files----
oldleaves <- read_csv("Cleaned_data/reclass.old.leaves.csv")
youngleaves <- read_csv("Cleaned_data/reclass.young.leaves.csv")

#Below is the old not reclassified leaf_rel_ages
#QL_leafdata <- read_csv("Cleaned_data/QL_leaflevel.df.csv") #LEAF LEVEL DATA ONLY
#clmdata <- read_csv("Cleaned_data/climate.data.csv") 

oldleaves2 <- oldleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12", "startDayOfYear","catalogNumber", "leaf_rel_age","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) 

youngleaves2 <- youngleaves%>%
  select(!c("perc_area_PS", "perc_area_ovi")) %>% 
  mutate_at(c("leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
              "leaftotal.freq", "leafspec.freq", "leafgall.freq",
              "leafmine.freq"), as.numeric) %>% #R is importing these as characters which isn't correct so updating them here
  mutate(leaf_rel_age = as.factor(leaf_rel_age)) %>% 
  select(c("MAT":"CMI_12","startDayOfYear","catalogNumber", "leaf_rel_age","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall",
           "leaftotal.freq", "leafspec.freq", "leafgall.freq",
           "leafmine.freq", "leaftotal.rich", "leafspec.rich", "leafgall.rich", "leafmine.rich")) 


#creating a new column of data that is 0:1 for chewing damage to then run models on
oldleaves3 <- oldleaves2
oldleaves3$chewbinary <- oldleaves3$leafperc_area_chew 
oldleaves3 <- oldleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% #keep column chewvalue in PCA, we dont need the binary anymore
  select(!chewbinary)
#Rounding percent leaf area numbers
oldleaves3$leafperc_area_chew <- round(oldleaves3$leafperc_area_chew, 4)
oldleaves3$leafperc_area_gall <- round(oldleaves3$leafperc_area_gall, 4)
oldleaves3$leafperc_area_mine <- round(oldleaves3$leafperc_area_mine, 4)

youngleaves3 <- youngleaves2
youngleaves3$chewbinary <- youngleaves3$leafperc_area_chew 
youngleaves3 <- youngleaves3 %>% 
  mutate(value = if_else(chewbinary < 0.0001, 0, 1)) %>% 
  rename(chewvalue = value) %>% 
  select(!chewbinary)
#rounding
youngleaves3$leafperc_area_chew <- round(youngleaves3$leafperc_area_chew, 4)
youngleaves3$leafperc_area_gall <- round(youngleaves3$leafperc_area_gall, 4)
youngleaves3$leafperc_area_mine <- round(youngleaves3$leafperc_area_mine, 4)

#recombine old and young datasets to make a full dataset that is the same number of leaves as the seperated datasets
leafdata2 <- youngleaves3 %>% 
  rbind(oldleaves3)


#subsetting data into annual data
Y_data <- leafdata2 %>% 
  select(c("catalogNumber", "MAT":"DD1040", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#seasonal data
S_data <- leafdata2 %>% 
  select(c("catalogNumber", "Tmax_wt":"CMI_at", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#monthly data
M_data <- leafdata2 %>% 
  select(c("catalogNumber", "Tmax_01":"CMI_12","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#old leaves
Yold_data <- oldleaves3 %>% 
  select(c("catalogNumber", "MAT":"DD1040", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#seasonal data
Sold_data <- oldleaves3 %>% 
  select(c("catalogNumber", "Tmax_wt":"CMI_at", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#monthly data
Mold_data <- oldleaves3 %>% 
  select(c("catalogNumber", "Tmax_01":"CMI_12","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#young leaves
Yyoung_data <- youngleaves3 %>% 
  select(c("catalogNumber", "MAT":"DD1040", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#seasonal data
Syoung_data <- youngleaves3 %>% 
  select(c("catalogNumber", "Tmax_wt":"CMI_at", "leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue"))%>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 

#monthly data
Myoung_data <- youngleaves3 %>% 
  select(c("catalogNumber", "Tmax_01":"CMI_12","leafperc_area_chew", "leafperc_area_mine", "leafperc_area_gall", 
           "leaftotal.freq":"chewvalue")) %>% 
  select(where(~!all(is.nan(.)))) %>%  #removes the columns (seasonal or monthly) with NaN values due to no data
  select(where(~!all(is.na(.)))) 


# running models on all datasets
Y_pca_results <- create_multiple_pca_models(Y_data, n_components = 10) #annual
Y_pca_results <- generate_pca_csv_output(Y_pca_results)
S_pca_results <- create_multiple_pca_models(S_data, n_components = 10) #seasonal
S_pca_results <- generate_pca_csv_output(S_pca_results)
M_pca_results <- create_multiple_pca_models(M_data, n_components = 10) #monthly
M_pca_results <- generate_pca_csv_output(M_pca_results)

#old leaves
oldY_pca_results <- create_multiple_pca_models(Yold_data, n_components = 10) #annual
oldY_pca_results <- generate_pca_csv_output(oldY_pca_results)
oldS_pca_results <- create_multiple_pca_models(Sold_data, n_components = 10) #seasonal
oldS_pca_results <- generate_pca_csv_output(oldS_pca_results)
oldM_pca_results <- create_multiple_pca_models(Mold_data, n_components = 10) #monthly
oldM_pca_results <- generate_pca_csv_output(oldM_pca_results)

#young leaves
youngY_pca_results <- create_multiple_pca_models(Yyoung_data, n_components = 10) #annual
youngY_pca_results <- generate_pca_csv_output(youngY_pca_results)
youngS_pca_results <- create_multiple_pca_models(Syoung_data, n_components = 10) #seasonal
youngS_pca_results <- generate_pca_csv_output(youngS_pca_results)
youngM_pca_results <- create_multiple_pca_models(Myoung_data, n_components = 10) #monthly
youngM_pca_results <- generate_pca_csv_output(youngM_pca_results)

#write/save outputs
# write.csv(Y_pca_results, "Script_output/PCA/LY_pca_results.csv", row.names = FALSE)
# write.csv(S_pca_results, "Script_output/PCA/LS_pca_results.csv", row.names = FALSE)
# write.csv(M_pca_results, "Script_output/PCA/LM_pca_results.csv", row.names = FALSE)
#old leaves
write.csv(oldY_pca_results, "Script_output/PCA/l.lvl/old/LoldY_pca_results.csv", row.names = FALSE)
write.csv(oldS_pca_results, "Script_output/PCA/l.lvl/old/LoldS_pca_results.csv", row.names = FALSE)
write.csv(oldM_pca_results, "Script_output/PCA/l.lvl/old/LoldM_pca_results.csv", row.names = FALSE)
#young leaves
write.csv(youngY_pca_results, "Script_output/PCA/l.lvl/young/LyoungY_pca_results.csv", row.names = FALSE)
write.csv(youngS_pca_results, "Script_output/PCA/l.lvl/young/LyoungS_pca_results.csv", row.names = FALSE)
write.csv(youngM_pca_results, "Script_output/PCA/l.lvl/young/LyoungM_pca_results.csv", row.names = FALSE)


# functions for calculating the variable contributions and the top 5 important climate variables
#Must save the above output and then read it back in for the contributions functions below
#read back in each saved pca result dataframe
# Y_results <- read.csv("Script_output/PCA/LY_pca_results.csv", stringsAsFactors = FALSE)
# S_results <- read.csv("Script_output/PCA/LS_pca_results.csv", stringsAsFactors = FALSE)
# M_results <- read.csv("Script_output/PCA/LM_pca_results.csv", stringsAsFactors = FALSE)

#old leaves
oldY_results <- read.csv("Script_output/PCA/l.lvl/old/LoldY_pca_results.csv", stringsAsFactors = FALSE)
oldS_results <- read.csv("Script_output/PCA/l.lvl/old/LoldS_pca_results.csv", stringsAsFactors = FALSE)
oldM_results <- read.csv("Script_output/PCA/l.lvl/old/LoldM_pca_results.csv", stringsAsFactors = FALSE)

#young leaves
youngY_results <- read.csv("Script_output/PCA/l.lvl/young/LyoungY_pca_results.csv", stringsAsFactors = FALSE)
youngS_results <- read.csv("Script_output/PCA/l.lvl/young/LyoungS_pca_results.csv", stringsAsFactors = FALSE)
youngM_results <- read.csv("Script_output/PCA/l.lvl/young/LyoungM_pca_results.csv", stringsAsFactors = FALSE)


# Function to calculate variable contributions----
calculate_contributions <- function(data) {
  predictors <- data[data$Type == "Predictor", ]
  pc_cols <- c("PC1", "PC2", "PC3", "PC4")
  
  # Convert PC columns to numeric
  predictors[, pc_cols] <- lapply(predictors[, pc_cols], as.numeric)
  
  # Square the loadings
  squared_loadings <- predictors[, pc_cols]^2
  
  # Calculate the sum of squared loadings for each PC
  sum_squared_loadings <- colSums(squared_loadings)
  
  # Calculate the contributions
  contributions <- sweep(squared_loadings, 2, sum_squared_loadings, "/") * 100
  
  # Add the variable names
  contributions <- data.frame(Variable = predictors$Variable, contributions)
  
  return(contributions)
}

#calculating contribtion for less than PC4 
calculate_contributions2 <- function(data) {
  predictors <- data[data$Type == "Predictor", ]
  pc_cols <- c("PC1", "PC2", "PC3")
  
  # Convert PC columns to numeric
  predictors[, pc_cols] <- lapply(predictors[, pc_cols], as.numeric)
  
  # Square the loadings
  squared_loadings <- predictors[, pc_cols]^2
  
  # Calculate the sum of squared loadings for each PC
  sum_squared_loadings <- colSums(squared_loadings)
  
  # Calculate the contributions
  contributions <- sweep(squared_loadings, 2, sum_squared_loadings, "/") * 100
  
  # Add the variable names
  contributions <- data.frame(Variable = predictors$Variable, contributions)
  
  return(contributions)
}

# Function to calculate importance for a single response variable----
calculate_importance <- function(response_name, contributions, data) {
  responses <- data[data$Type == "Response", ]
  response_corr <- responses[grep(paste0(response_name, "_Correlation_PC"), responses$Variable), "Value"]
  response_corr <- as.numeric(response_corr)
  
  importance <- contributions
  pc_cols <- c("PC1", "PC2", "PC3", "PC4")
  for (i in 1:4) {
    importance[, pc_cols[i]] <- importance[, pc_cols[i]] * abs(response_corr[i]) #i.e., how far away is this response from zero. 
  }
  importance$Total <- rowSums(importance[, pc_cols])
  importance <- importance[order(-importance$Total), ]
  importance$Response <- response_name
  return(importance)
}

#Function to calculate importance when there are less than 4 PC columns
calculate_importance2 <- function(response_name, contributions, data) {
  responses <- data[data$Type == "Response", ]
  response_corr <- responses[grep(paste0(response_name, "_Correlation_PC"), responses$Variable), "Value"]
  response_corr <- as.numeric(response_corr)
  
  importance <- contributions[contributions$Type == "Predictor", ]
  pc_cols <- c("PC1", "PC2", "PC3")
  
  for (i in seq_along(pc_cols)) {
    importance[, pc_cols[i]] <- importance[, pc_cols[i]] * abs(response_corr[i])
  }
  
  importance$Total <- rowSums(importance[, pc_cols]) #total is the sum of the importance (percentage)
  importance <- importance[order(-importance$Total), ]
  importance$Response <- response_name
  
  # Select only relevant columns for the output
  importance <- importance[, c("Variable", "Type", pc_cols, "Total", "Response")]
  
  return(importance)
}


# Calculate contributions for each result dataframe that was read back into R
#full data set
# Ycon <- calculate_contributions(Y_results)
# Scon <- calculate_contributions(S_results)
# Mcon <- calculate_contributions(M_results)

#old leaves
oldYcon <- calculate_contributions2(oldY_results) 
oldScon <- calculate_contributions(oldS_results)
oldMcon <- calculate_contributions(oldM_results)

#young leaves
youngYcon <- calculate_contributions(youngY_results)
youngScon <- calculate_contributions(youngS_results)
youngMcon <- calculate_contributions(youngM_results)


#top important function----
#run this LAST function to create the necessary outputs
calculate_top_important_variables <- function(data, contributions, output_file = "top_important_variables.csv") {
  
  # Helper function to calculate importance for a single response variable
  calculate_importance <- function(response_name, contributions, data) {
    responses <- data[data$Type == "Response", ]
    response_corr <- responses[grep(paste0(response_name, "_Correlation_PC"), responses$Variable), "Value"]
    response_corr <- as.numeric(response_corr)
    
    importance <- contributions
    pc_cols <- c("PC1", "PC2", "PC3", "PC4")
    for (i in 1:4) {
      importance[, pc_cols[i]] <- importance[, pc_cols[i]] * abs(response_corr[i])
    }
    importance$Total <- rowSums(importance[, pc_cols])
    importance <- importance[order(-importance$Total), ]
    importance$Response <- response_name
    return(importance)
  }
  
  # Get list of response variables
  response_vars <- unique(gsub("_Correlation_PC[1-4]", "", data$Variable[data$Type == "Response"]))
  
  # Calculate importance for all response variables
  all_importance <- do.call(rbind, lapply(response_vars, function(x) calculate_importance(x, contributions, data)))
  
  # Create an empty list to store results
  results_list <- list()
  
  # Loop through each unique response
  for (response in unique(all_importance$Response)) {
    # Get top 5 important variables for this response
    top_5 <- head(all_importance[all_importance$Response == response, c("Variable", "Total")], 5)
    
    # Add response column
    top_5$Response <- response
    
    # Add to results list
    results_list[[response]] <- top_5
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, results_list)
  
  # Reorder columns
  final_results <- final_results[, c("Response", "Variable", "Total")]
  
  # Write to CSV
  write.csv(final_results, output_file, row.names = FALSE)
  
  return(final_results)
}

#Secondary function when less than PC4 
calculate_top_important_variables2 <- function(data, contributions, output_file = "top_important_variables.csv") {
  
  # Helper function to calculate importance for a single response variable
  calculate_importance <- function(response_name, contributions, data) {
    responses <- data[data$Type == "Response", ]
    response_corr <- responses[grep(paste0(response_name, "_Correlation_PC"), responses$Variable), "Value"]
    response_corr <- as.numeric(response_corr)
    
    importance <- contributions
    pc_cols <- c("PC1", "PC2", "PC3")
    for (i in 1:3) {
      importance[, pc_cols[i]] <- importance[, pc_cols[i]] * abs(response_corr[i])
    }
    importance$Total <- rowSums(importance[, pc_cols])
    importance <- importance[order(-importance$Total), ]
    importance$Response <- response_name
    return(importance)
  }
  
  # Get list of response variables
  response_vars <- unique(gsub("_Correlation_PC[1-3]", "", data$Variable[data$Type == "Response"]))
  
  # Calculate importance for all response variables
  all_importance <- do.call(rbind, lapply(response_vars, function(x) calculate_importance(x, contributions, data)))
  
  # Create an empty list to store results
  results_list <- list()
  
  # Loop through each unique response
  for (response in unique(all_importance$Response)) {
    # Get top 5 important variables for this response
    top_5 <- head(all_importance[all_importance$Response == response, c("Variable", "Total")], 5)
    
    # Add response column
    top_5$Response <- response
    
    # Add to results list
    results_list[[response]] <- top_5
  }
  
  # Combine all results into a single data frame
  final_results <- do.call(rbind, results_list)
  
  # Reorder columns
  final_results <- final_results[, c("Response", "Variable", "Total")]
  
  # Write to CSV
  write.csv(final_results, output_file, row.names = FALSE)
  
  return(final_results)
}

#Run contributions function on each result input dataframe with the corresponding contributions calculated above
# topY <- calculate_top_important_variables(Y_results, Ycon) %>% 
#   drop_na()
# topS <- calculate_top_important_variables(S_results, Scon) %>% 
#   drop_na()
# topM <- calculate_top_important_variables(M_results, Mcon) %>% 
#   drop_na()

topoldY <- calculate_top_important_variables2(oldY_results, oldYcon) %>% #same as above
  drop_na()
topoldS <-calculate_top_important_variables(oldS_results, oldScon) %>% 
  drop_na()
topoldM <- calculate_top_important_variables(oldM_results, oldMcon) %>% 
  drop_na()

topyoungY <- calculate_top_important_variables(youngY_results, youngYcon) %>% 
  drop_na()
topyoungS <- calculate_top_important_variables(youngS_results, youngScon) %>% 
  drop_na()
topyoungM <- calculate_top_important_variables(youngM_results, youngMcon) %>% 
  drop_na()

#Notes: 
#remember that the "total" that is being calculate is already a percentage and is the sum of 
#the contributions per PC by resopnse variable 

#write the final .csv files
# write.csv(topY, "Script_Output/PCA/Ltop5Y.csv", row.names = FALSE)
# write.csv(topS, "Script_Output/PCA/Ltop5S.csv", row.names = FALSE)
# write.csv(topM, "Script_Output/PCA/Ltop5M.csv", row.names = FALSE)

#old leaves
write.csv(topoldY, "Script_Output/PCA/l.lvl/old/Ltop5oldY.csv", row.names = FALSE)
write.csv(topoldS, "Script_Output/PCA/l.lvl/old/Ltop5oldS.csv", row.names = FALSE)
write.csv(topoldM, "Script_Output/PCA/l.lvl/old/Ltop5oldM.csv", row.names = FALSE)

#young leaves
write.csv(topyoungY, "Script_Output/PCA/l.lvl/young/Ltop5youngY.csv", row.names = FALSE)
write.csv(topyoungS, "Script_Output/PCA/l.lvl/young/Ltop5youngS.csv", row.names = FALSE)
write.csv(topyoungM, "Script_Output/PCA/l.lvl/young/Ltop5youngM.csv", row.names = FALSE)

#NOTE: All non-reclassified PCA outputs were moved into the "archive" folder within the PCA folder.
  #All PCAs are now based on the reclassified ages of the leaves as of 10/24/2024





