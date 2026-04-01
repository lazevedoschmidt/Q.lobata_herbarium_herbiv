library(tidyverse)
library(stringr)
library(patchwork)
library(cowplot)

# FUNCTION 1: Consolidate time-space model data - FIXED for your filename pattern
consolidate_timespace_data <- function(directory_path, pattern = "_model\\.csv$") {
  
  # Get all model CSV files
  model_files <- list.files(directory_path, 
                            pattern = pattern, 
                            full.names = TRUE)
  
  if (length(model_files) == 0) {
    stop("No model files found in directory: ", directory_path)
  }
  
  message("Found ", length(model_files), " model files")
  
  # Function to extract herbivory category and spatial variable from filename
  extract_file_info <- function(filepath) {
    filename <- basename(filepath)
    
    # Remove "_model.csv" from end
    name_core <- str_remove(filename, "_model\\.csv$")
    
    # Pattern: herbivory_category_year.spatial (e.g., "chewvalue_year.lat" or "leafgall.freq_year.long")
    # Split by "_year." to separate herbivory category from spatial variable
    if (grepl("_year\\.", name_core)) {
      parts <- str_split(name_core, "_year\\.")[[1]]
      
      if (length(parts) == 2) {
        herbivory_cat <- parts[1]
        spatial_var <- parts[2]
        
        # Convert "lat" to "latitude" and "long" to "longitude" for clarity
        if (spatial_var == "lat") {
          spatial_var <- "latitude"
        } else if (spatial_var == "long") {
          spatial_var <- "longitude"
        }
        
        return(list(
          herbivory_category = herbivory_cat,
          spatial_variable = spatial_var
        ))
      }
    }
    
    warning(paste("Could not parse filename:", filename))
    return(list(
      herbivory_category = NA,
      spatial_variable = NA
    ))
  }
  
  # Process each file
  all_data <- map_df(model_files, function(filepath) {
    # Extract metadata from filename
    file_info <- extract_file_info(filepath)
    
    # Skip if parsing failed
    if (is.na(file_info$herbivory_category)) {
      return(NULL)
    }
    
    # Read the CSV
    model_data <- read_csv(filepath, show_col_types = FALSE)
    
    # We want the year slope estimate - try different possible column names
    possible_year_cols <- c("b_scaleyear.x", "b_scaleyear", "b_year")
    year_col <- possible_year_cols[possible_year_cols %in% colnames(model_data)][1]
    
    if (is.na(year_col)) {
      warning(paste("No year column found in", basename(filepath)))
      return(NULL)
    }
    
    # Create a simplified dataframe with metadata
    model_data <- model_data %>%
      select(all_of(year_col)) %>%
      mutate(
        herbivory_category = file_info$herbivory_category,
        spatial_variable = file_info$spatial_variable,
        year_estimate = .data[[year_col]],
        .before = 1
      ) %>%
      select(-all_of(year_col))
    
    return(model_data)
  })
  
  # Remove any NULL entries from failed file processing
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No data was successfully processed.")
  }
  
  all_data <- all_data %>% filter(!is.na(herbivory_category))
  
  message("Successfully processed data for:")
  message("  Herbivory categories: ", paste(unique(all_data$herbivory_category), collapse = ", "))
  message("  Spatial variables: ", paste(unique(all_data$spatial_variable), collapse = ", "))
  
  return(all_data)
}


# FUNCTION 2: Create forest plot for year slopes by spatial variable
plot_year_slopes_by_spatial <- function(consolidated_data, 
                                        show_all = TRUE) {
  
  # Calculate summary statistics for the year estimates
  summary_data <- consolidated_data %>%
    group_by(herbivory_category, spatial_variable) %>%
    summarize(
      mean = mean(year_estimate, na.rm = TRUE),
      lower_95 = quantile(year_estimate, 0.025, na.rm = TRUE),
      upper_95 = quantile(year_estimate, 0.975, na.rm = TRUE),
      lower_80 = quantile(year_estimate, 0.1, na.rm = TRUE),
      upper_80 = quantile(year_estimate, 0.9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Flag significant estimates (80% CI doesn't cross zero)
      significant_80 = (lower_80 > 0) | (upper_80 < 0)
    )
  
  # Filter based on show_all parameter
  if (!show_all) {
    plot_data <- summary_data %>% filter(significant_80)
  } else {
    plot_data <- summary_data
  }
  
  # Calculate x-axis limits
  x_limits <- c(
    min(plot_data$lower_95, na.rm = TRUE),
    max(plot_data$upper_95, na.rm = TRUE)
  )
  
  # Get unique herbivory categories for ordering
  herbivory_categories <- unique(plot_data$herbivory_category)
  
  # Create ordered factor for herbivory categories (reverse for top-to-bottom display)
  plot_data$herbivory_category <- factor(
    plot_data$herbivory_category,
    levels = rev(sort(herbivory_categories)),
    ordered = TRUE
  )
  
  category_colors <- c(
    "chewvalue" = "#E16036",
    "leafgall.freq" = "#5EB1BF",
    "leafgall.rich" = "#6D4F92",
    "leafgen.freq" = "#C1AE7C",
    "leafgen.rich" = "#2A9D8F",
    "leafmine.freq" = "#E76F51",
    "leafmine.rich" = "#F4A261",
    "leafperc_area_chew" = "#E63946",
    "leafperc_area_gall" = "#457B9D",
    "leafperc_area_mine" = "#A8DADC",
    "leafspec.freq" = "#8338EC",
    "leafspec.rich" = "#FF006E",
    "leaftotal.freq" = "#FB5607",
    "leaftotal.rich" = "#FFBE0B"
  )
  
  # Custom labels for herbivory categories
  category_labels <- c(
    "chewvalue" = "Chewing Frequency",
    "leafgall.freq" = "Gall Frequency",
    "leafgall.rich" = "Gall Richness",
    "leafgen.freq" = "Generalized Frequency",
    "leafgen.rich" = "Generalized Richness",
    "leafmine.freq" = "Mine Frequency",
    "leafmine.rich" = "Mine Richness",
    "leafperc_area_chew" = "% Area Chewed",
    "leafperc_area_gall" = "% Area Galled",
    "leafperc_area_mine" = "% Area Mined",
    "leafspec.freq" = "Specialized Frequency",
    "leafspec.rich" = "Specialized Richness",
    "leaftotal.freq" = "Total Frequency",
    "leaftotal.rich" = "Total Richness"
  )
  
  # Apply labels to the plot data
  plot_data <- plot_data %>%
    mutate(herbivory_label = factor(
      herbivory_category,
      levels = rev(sort(herbivory_categories)),
      labels = category_labels[rev(sort(herbivory_categories))]
    ))
  
  # Create a list to store plots
  plot_list <- list()
  
  # Create a plot for each spatial variable (latitude/longitude)
  spatial_vars <- unique(plot_data$spatial_variable)
  
  for (spatial_var in spatial_vars) {
    var_data <- plot_data %>%
      filter(spatial_variable == spatial_var)
    
    # Skip if no data for this variable
    if (nrow(var_data) == 0) {
      message(paste("No data found for spatial variable:", spatial_var))
      next
    }
    
    p <- var_data %>%
      ggplot(aes(y = herbivory_label, x = mean, color = herbivory_category)) +
      # 95% CI with thin line
      geom_linerange(
        aes(xmin = lower_95, xmax = upper_95),
        linewidth = 0.5
      ) +
      # 80% CI with thick line
      geom_linerange(
        aes(xmin = lower_80, xmax = upper_80),
        linewidth = 2.5
      ) +
      # Point estimate
      geom_point(size = 3) +
      # Zero reference line
      geom_vline(xintercept = 0, linetype = "dotted", color = "gray50") +
      theme_minimal() +
      labs(
        title = paste("Year Slope with", tools::toTitleCase(spatial_var)),
        x = "Year Slope Estimate",
        y = NULL
      ) +
      scale_x_continuous(limits = x_limits) +
      scale_color_manual(values = category_colors, labels = category_labels, name = "Herbivory Category") +
      theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      )
    
    plot_list[[spatial_var]] <- p
  }
  
  return(plot_list)
}


# USING FUNCTIONS ----
# OLD LEAVES
old_timespace_data <- consolidate_timespace_data(
  directory_path = "Script_output/Brms/time_spac_models/l.lvl/old/"
)

old_timespace_plots <- plot_year_slopes_by_spatial(
  consolidated_data = old_timespace_data,
  show_all = FALSE  # Set to TRUE to see all effects
)
old_timespace_plots$latitude
old_timespace_plots$longitude

# YOUNG LEAVES
young_timespace_data <- consolidate_timespace_data(
  directory_path = "Script_output/Brms/time_spac_models/l.lvl/young/"
)

young_timespace_plots <- plot_year_slopes_by_spatial(
  consolidated_data = young_timespace_data,
  show_all = TRUE #FALSE = no plots for young leaves here
)


#BELOW IF FOR NON-SIGNIFICANT BECAUSE THERE ARE NO SIGNIFICANT RELATIONSHIPS WITH YOUNG LEAVES

# Combine young and old plots
combined_timespace_plots <- list()

# Get all spatial variables that exist in either dataset
all_spatial_vars <- unique(c(names(young_timespace_plots), names(old_timespace_plots)))

for (spatial_var in all_spatial_vars) {
  has_young <- spatial_var %in% names(young_timespace_plots)
  has_old <- spatial_var %in% names(old_timespace_plots)
  
  if (has_young && has_old) {
    combined_timespace_plots[[spatial_var]] <- 
      young_timespace_plots[[spatial_var]] + 
      labs(title = paste("Year Slope with", tools::toTitleCase(spatial_var), "- Young Leaves")) +
      old_timespace_plots[[spatial_var]] + 
      labs(title = paste("Year Slope with", tools::toTitleCase(spatial_var), "- Old Leaves")) +
      plot_annotation(
        title = paste("Temporal Trends with", tools::toTitleCase(spatial_var), "in Model"),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
  } else if (has_young) {
    combined_timespace_plots[[spatial_var]] <- 
      young_timespace_plots[[spatial_var]] + 
      labs(title = paste("Year Slope with", tools::toTitleCase(spatial_var), "- Young Leaves Only"))
  } else if (has_old) {
    combined_timespace_plots[[spatial_var]] <- 
      old_timespace_plots[[spatial_var]] + 
      labs(title = paste("Year Slope with", tools::toTitleCase(spatial_var), "- Old Leaves Only"))
  }
}

# View plots
if ("latitude" %in% names(combined_timespace_plots)) {
  plot_latitude <- combined_timespace_plots$latitude
  print(plot_latitude)
}

if ("longitude" %in% names(combined_timespace_plots)) {
  plot_longitude <- combined_timespace_plots$longitude
  print(plot_longitude)
}

# Create legend (same as before)
category_colors <- c(
  "chewvalue" = "#E16036",
  "leafgall.freq" = "#5EB1BF",
  "leafgall.rich" = "#6D4F92",
  "leafgen.freq" = "#C1AE7C",
  "leafgen.rich" = "#2A9D8F",
  "leafmine.freq" = "#E76F51",
  "leafmine.rich" = "#F4A261",
  "leafperc_area_chew" = "#E63946",
  "leafperc_area_gall" = "#457B9D",
  "leafperc_area_mine" = "#A8DADC",
  "leafspec.freq" = "#8338EC",
  "leafspec.rich" = "#FF006E",
  "leaftotal.freq" = "#FB5607",
  "leaftotal.rich" = "#FFBE0B"
)

category_labels <- c(
  "chewvalue" = "Chew Value",
  "leafgall.freq" = "Gall Frequency",
  "leafgall.rich" = "Gall Richness",
  "leafgen.freq" = "Generalized Frequency",
  "leafgen.rich" = "Generalized Richness",
  "leafmine.freq" = "Mine Frequency",
  "leafmine.rich" = "Mine Richness",
  "leafperc_area_chew" = "% Area Chewed",
  "leafperc_area_gall" = "% Area Galled",
  "leafperc_area_mine" = "% Area Mined",
  "leafspec.freq" = "Specialized Frequency",
  "leafspec.rich" = "Specialized Richness",
  "leaftotal.freq" = "Total Frequency",
  "leaftotal.rich" = "Total Richness"
)

all_categories <- names(category_colors)
dummy_data <- data.frame(
  herbivory_category = all_categories,
  x = 1:length(all_categories),
  y = 1:length(all_categories)
)

legend_plot <- ggplot(dummy_data, aes(x = x, y = y, color = herbivory_category)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = category_colors,
    labels = category_labels,
    name = "Herbivory Category"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +
  guides(colour = guide_legend(nrow = 2))

full_legend <- get_legend(legend_plot)

# Combine all plots with legend
if (length(combined_timespace_plots) > 0) {
  all_spatial_combined <- plot_grid(
    plotlist = combined_timespace_plots,
    nrow = length(combined_timespace_plots),
    ncol = 1,
    labels = LETTERS[1:length(combined_timespace_plots)]
  )
  
  final_plot <- plot_grid(
    all_spatial_combined,
    full_legend,
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 0.15)
  )
  
  print(final_plot)
  
  # Save
  ggsave("figures/l.lvl/timespace_year_nonsig.pdf", final_plot, #FIGURE FOR NON-"SIGNIFICANT" OUTPUT
         width = 14, height = 10, dpi = 300)
}

#Plotting significant relationships for old leaves only
oldsig.lat <- old_timespace_plots$latitude
oldsig.long <- old_timespace_plots$longitude

oldsig.tsp <- plot_grid(oldsig.lat,
                        oldsig.long,
                        nrow = 1,
                        ncol = 2
                        )
oldsig.tsp2 <- plot_grid(oldsig.tsp,
                        full_legend,
                        nrow = 2,
                        ncol = 1,
                        rel_heights = c(1,.25))
oldsig.tsp2

ggsave("figures/l.lvl/timespace_year_sig.pdf", oldsig.tsp2, width = 14, height = 10, dpi = 300)

