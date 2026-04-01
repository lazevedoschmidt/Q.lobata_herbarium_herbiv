library(tidyverse)
library(stringr)
library(patchwork)
library(cowplot)

# FUNCTION 1: Consolidate all model CSV files into a single dataframe
consolidate_model_data <- function(directory_path, pattern = "_model\\.csv$") {
  
  # Get all model CSV files
  model_files <- list.files(directory_path, 
                            pattern = pattern, 
                            full.names = TRUE)
  
  # Function to extract herbivory category and climate variable from filename
  extract_file_info <- function(filepath) {
    filename <- basename(filepath)
    # Remove "_year_doy_model.csv" from end
    name_core <- str_remove(filename, "_year_doy_model\\.csv$")
    
    # Split by underscore to get herbivory_category and climate variable
    # Pattern: herbivory_category_CLIMATE (e.g., "chewvalue_EXT" or "leafgall.freq_MAP")
    parts <- str_split(name_core, "_")[[1]]
    
    # The climate variable is the last part
    # Everything before it is the herbivory category
    if (length(parts) >= 2) {
      climate_var <- parts[length(parts)]
      herbivory_cat <- paste(parts[1:(length(parts)-1)], collapse = "_")
      
      return(list(
        herbivory_category = herbivory_cat,
        climate_variable = climate_var
      ))
    } else {
      warning(paste("Could not parse filename:", filename))
      return(list(
        herbivory_category = NA,
        climate_variable = NA
      ))
    }
  }
  
  # Process each file
  all_data <- map_df(model_files, function(filepath) {
    # Extract metadata from filename
    file_info <- extract_file_info(filepath)
    
    # Read the CSV
    model_data <- read_csv(filepath, show_col_types = FALSE)
    
    # Find the climate variable column (should match the climate variable from filename)
    climate_col <- paste0("b_scale", file_info$climate_variable)
    
    # Check if the expected column exists
    if (!climate_col %in% colnames(model_data)) {
      warning(paste("Expected column", climate_col, "not found in", basename(filepath)))
      return(NULL)
    }
    
    # Select only the relevant columns: intercept, year, doy, and the climate variable
    # We'll keep b_scaleyear.x and b_scaledoy.clean as they're likely in all models
    keep_cols <- c("b_Intercept", "b_scaleyear.x", "b_scaledoy.clean", climate_col)
    keep_cols <- keep_cols[keep_cols %in% colnames(model_data)]
    
    # Create a simplified dataframe with metadata
    model_data <- model_data %>%
      select(all_of(keep_cols)) %>%
      mutate(
        herbivory_category = file_info$herbivory_category,
        climate_variable = file_info$climate_variable,
        # Rename the climate column to a standard name
        climate_estimate = .data[[climate_col]],
        .before = 1
      ) %>%
      select(-all_of(climate_col))
    
    return(model_data)
  })
  
  # Remove any NULL entries from failed file processing
  all_data <- all_data %>% filter(!is.na(herbivory_category))
  
  return(all_data)
}


# FUNCTION 2: Create forest plots organized by climate variable
plot_forest_by_climate <- function(consolidated_data, 
                                   climate_vars = c("AHM", "MAR", "RH", "MSP", "MAP", "SHM", "EXT", "TD"),
                                   show_all = TRUE) {
  
  # Calculate summary statistics for the climate estimates
  summary_data <- consolidated_data %>%
    group_by(herbivory_category, climate_variable) %>%
    summarize(
      mean = mean(climate_estimate, na.rm = TRUE),
      lower_95 = quantile(climate_estimate, 0.025, na.rm = TRUE),
      upper_95 = quantile(climate_estimate, 0.975, na.rm = TRUE),
      lower_80 = quantile(climate_estimate, 0.1, na.rm = TRUE),
      upper_80 = quantile(climate_estimate, 0.9, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Flag significant estimates (80% CI doesn't cross zero)
      significant_80 = (lower_80 > 0) | (upper_80 < 0)
    )
  
  # Check what climate variables are actually in the data
  actual_climate_vars <- unique(summary_data$climate_variable)
  message("Climate variables found in data: ", paste(actual_climate_vars, collapse = ", "))
  
  # Standardize climate variable names (case-insensitive match)
  summary_data <- summary_data %>%
    mutate(climate_variable = toupper(climate_variable))
  
  # Filter to only include specified climate variables
  plot_data <- summary_data %>%
    filter(climate_variable %in% climate_vars)
  
  # Filter based on show_all parameter
  if (!show_all) {
    plot_data <- plot_data %>% filter(significant_80)
  }
  
  # Debug: show what we have after filtering
  message("After filtering, climate variables: ", paste(unique(plot_data$climate_variable), collapse = ", "))
  message("Number of rows in plot_data: ", nrow(plot_data))
  
  # Calculate x-axis limits BEFORE the loop
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
  
  # Create a plot for each climate variable
  for (climate_var in climate_vars) {
    var_data <- plot_data %>%
      filter(climate_variable == climate_var)
    
    # Skip if no data for this variable
    if (nrow(var_data) == 0) {
      message(paste("No data found for climate variable:", climate_var))
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
        title = climate_var,
        x = "Slope Estimate",
        y = "Herbivory Category"
      ) +
      scale_x_continuous(limits = x_limits) +
      scale_color_manual(values = category_colors, labels = category_labels ,name = "Herbivory Category") +
      theme(legend.position = "none")
      # theme(
      #   legend.position = "right",
      #   panel.grid.minor = element_blank(),
      #   plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
      #)
    
    plot_list[[climate_var]] <- p
  }
  
  return(plot_list)
}


# USING FUNCTIONS----
# OLD LEAVES
old_model_data <- consolidate_model_data(
  directory_path = "Script_output/Brms/climate_models/l.lvl/old/year_doy/"
)

old_climate_plots <- plot_forest_by_climate(
  consolidated_data = old_model_data,
  climate_vars = c("AHM", "MAR", "RH", "MSP", "MAP", "SHM", "EXT", "TD"),
  show_all = TRUE  # or FALSE to see only significant effects
)

# YOUNG LEAVES
young_model_data <- consolidate_model_data(
  directory_path = "Script_output/Brms/climate_models/l.lvl/young/year_doy/"
)

young_climate_plots <- plot_forest_by_climate(
  consolidated_data = young_model_data,
  climate_vars = c("AHM", "MAR", "RH", "MSP", "MAP", "SHM", "EXT", "TD"),
  show_all = TRUE  # or FALSE
)

# Create a list to store combined plots
combined_plots <- list()

# Get all climate variables that exist in either dataset
all_climate_vars <- unique(c(names(young_climate_plots), names(old_climate_plots)))

# Loop through each climate variable and combine young + old
for (climate_var in all_climate_vars) {
  
  # Check if the variable exists in both datasets
  has_young <- climate_var %in% names(young_climate_plots)
  has_old <- climate_var %in% names(old_climate_plots)
  
  if (has_young && has_old) {
    # Both exist - combine side by side
    combined_plots[[climate_var]] <- 
      young_climate_plots[[climate_var]] + 
      labs(title = paste(climate_var, "- Young Leaves")) +
      old_climate_plots[[climate_var]] + 
      labs(title = paste(climate_var, "- Old Leaves")) +
      plot_annotation(title = paste("Climate Variable:", climate_var),
                      theme = theme(plot.title = element_text(size = 16, face = "bold")))
    
  } else if (has_young) {
    # Only young exists
    combined_plots[[climate_var]] <- 
      young_climate_plots[[climate_var]] + 
      labs(title = paste(climate_var, "- Young Leaves Only"))
    
  } else if (has_old) {
    # Only old exists
    combined_plots[[climate_var]] <- 
      old_climate_plots[[climate_var]] + 
      labs(title = paste(climate_var, "- Old Leaves Only"))
  }
}

# View individual combined plots
combined_plots
plotEXT <- combined_plots$EXT
plotMAR <- combined_plots$MAR
plotAHM <- combined_plots$AHM
plotRH <- combined_plots$RH
plotMSP <- combined_plots$MSP
plotMAP <- combined_plots$MAP 
plotSHM <- combined_plots$SHM
plotTD <- combined_plots$TD

# Create a dummy plot with all herbivory categories
# Define your color palette (copy from the function)
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

# Define your labels (copy from the function)
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

# Now create the dummy plot
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
  theme(legend.position = "right")+
  guides(colour = guide_legend(nrow = 2)) #making the legend horizonal (2 rows because there are a lot of names)

# Extract the legend
full_legend <- get_legend(legend_plot)
# Display it
plot(full_legend)

youngcombined <- plot_grid(plotEXT,
                         plotMAR,
                        nrow = 2, 
                        ncol = 1,
                        labels = c("A", "B"))
youngcombined
oldcombined <- plot_grid(plotRH,
                         plotMSP,
                         plotMAP,
                         plotSHM,
                         plotTD,
                         nrow = 5,
                         ncol = 1,
                         labels = c("D","E","F", "G", "H"))
oldcombined
allcombined1 <- plot_grid(youngcombined,
                          oldcombined,
                          nrow = 1,
                          ncol = 2,
                          rel_widths = c(1,1.15))
allcombined1
allcombined2 <- plot_grid(allcombined1,
                          plotAHM,
                          full_legend,
                          nrow = 3,
                          ncol = 1,
                          rel_heights = c(1,.65,.15),
                          rel_widths = c(1,1.15,1))
allcombined2

#nonsig plot for SI
allnonsig <- plot_grid(plotAHM,
                       plotMAR,
                       plotRH,
                       plotMSP,
                       plotSHM,
                       plotEXT,
                       nrow = 6,
                       ncol = 1)
allnonsig
allnonsig2 <- plot_grid(NULL,plotMAP,
                        NULL,plotTD,
                        nrow = 2,
                        ncol = 2)
allnonsig2
allnonsig3 <- plot_grid(allnonsig,
                        allnonsig2,
                        nrow = 2,
                        ncol = 1,
                        rel_heights = c(1,.25))
allnonsig3

# Save individual plots
ggsave("figures/l.lvl/EXT_plot.png", plotEXT, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/MAR_plot.png", plotMAR, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/RH_plot.png", plotRH, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/MSP_plot.png", plotMSP, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/MAP_plot.png", plotMAP, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/SHM_plot.png", plotSHM, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/TD_plot.png", plotTD, width = 7, height = 7, dpi = 300)
ggsave("figures/l.lvl/AHM_plot.png", plotAHM, width = 14, height = 7, dpi = 300)
ggsave("figures/l.lvl/allclimate.combined_plot.pdf", allcombined2, width = 14, height = 16, dpi = 300)
ggsave("figures/l.lvl/allclimate.combinednon.sig_plot.pdf", allnonsig3, width = 14, height = 30, dpi = 300)






