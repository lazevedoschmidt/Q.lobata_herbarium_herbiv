#Reviewer response checks----
#Added 10/14/2026 in response to reviewer comments.
#These checks do NOT change the PCAs, importance scores, or top-5 variables above.
#They use objects created earlier in this script, so run the full script first.
#  1) Temporal resolution comparison: R2 of each herbivory category on the 10 retained PCs
#  2) Comparison of PCA-selected variables with a ranking by direct correlation (|r|)
#All outputs are written to Script_output/PCA/l.lvl/reviewer_checks/

out_dir <- "Script_output/PCA/l.lvl/reviewer_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#datasets used for the PCAs above
clim_data <- list(oldY = Yold_data,   oldS = Sold_data,   oldM = Mold_data,
                  youngY = Yyoung_data, youngS = Syoung_data, youngM = Myoung_data)
pca_out   <- list(oldY = oldY_results,   oldS = oldS_results,   oldM = oldM_results,
                  youngY = youngY_results, youngS = youngS_results, youngM = youngM_results)

#number of leaves used in each PCA (rows remaining after na.omit)
n_rows <- map_int(clim_data, ~ nrow(na.omit(.x)))
print(n_rows)


#1) Temporal resolution comparison----
#PC scores are uncorrelated, so the R2 of a regression of a herbivory category on the
#10 retained PCs equals the sum of its squared correlations with those PCs.
#Adjusted R2 accounts for the 10 predictors: 1 - (1 - R2)(n - 1)/(n - 11)
r2 <- imap_dfr(pca_out, function(res, name) {
  res %>%
    filter(Type == "Response") %>%
    mutate(Response = sub("_Correlation_PC[0-9]+$", "", Variable),
           r        = as.numeric(Value)) %>%
    group_by(Response) %>%
    summarise(R2 = sum(r^2), .groups = "drop") %>%
    mutate(n          = n_rows[[name]],
           adjR2      = 1 - (1 - R2) * (n - 1) / (n - 11),
           leaf_age   = sub("[YSM]$", "", name),
           resolution = recode(sub("^(old|young)", "", name),
                               Y = "annual", S = "seasonal", M = "monthly"))
})

#one row per herbivory category, resolutions side by side
r2_wide <- r2 %>%
  select(Response, leaf_age, resolution, R2, adjR2) %>%
  pivot_wider(names_from = resolution, values_from = c(R2, adjR2)) %>%
  mutate(best = c("annual", "seasonal", "monthly")[
    max.col(cbind(R2_annual, R2_seasonal, R2_monthly), ties.method = "first")])

r2_wide %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% print(n = 28)

r2_summary <- r2_wide %>%
  group_by(leaf_age) %>%
  summarise(median_R2_annual      = median(R2_annual,      na.rm = TRUE),
            median_R2_seasonal    = median(R2_seasonal,    na.rm = TRUE),
            median_R2_monthly     = median(R2_monthly,     na.rm = TRUE),
            median_adjR2_annual   = median(adjR2_annual,   na.rm = TRUE),
            median_adjR2_seasonal = median(adjR2_seasonal, na.rm = TRUE),
            median_adjR2_monthly  = median(adjR2_monthly,  na.rm = TRUE),
            n_best_annual         = sum(best == "annual"),
            n_best_seasonal       = sum(best == "seasonal"),
            n_best_monthly        = sum(best == "monthly"),
            .groups = "drop")
print(r2_summary)

write.csv(r2_wide %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
          file.path(out_dir, "resolution_R2_by_response.csv"), row.names = FALSE)
write.csv(r2_summary, file.path(out_dir, "resolution_R2_summary.csv"), row.names = FALSE)


#2) PCA selection vs. direct correlation----
#Uses the same rows (na.omit) and constant-column removal as create_multiple_pca_models()

#prepare data exactly as the PCA function did
prep_data <- function(data) {
  response_vars  <- grep("^leaf|^chew|^ffg", names(data), value = TRUE)
  predictor_vars <- setdiff(names(data), c(response_vars, "catalogNumber"))
  data <- na.omit(data)
  X <- as.matrix(data[, predictor_vars])
  X <- X[, apply(X, 2, function(col) length(unique(col)) > 1), drop = FALSE]
  list(data = data, X = X, response_vars = response_vars)
}

#rank every climate variable by |r| with each herbivory category (Pearson, as in the PCA)
direct_ranks <- function(prepped) {
  map_dfr(prepped$response_vars, function(rv) {
    y <- prepped$data[[rv]]
    if (length(unique(y)) == 1) return(NULL)              #zero-variance response, skipped as in the PCA
    r <- cor(prepped$X, y)[, 1]
    tibble(Response = rv, Variable = names(r), r = r) %>%
      arrange(desc(abs(r))) %>%
      mutate(direct_rank = row_number())
  })
}

#largest |r| between any two variables in a selected set (how redundant the set is)
max_pairwise_r <- function(vars, clim_cor) {
  if (length(vars) < 2) return(NA_real_)
  m <- abs(clim_cor[vars, vars])
  max(m[upper.tri(m)])
}

compare_one <- function(name, data, top5) {
  prepped  <- prep_data(data)
  clim_cor <- cor(prepped$X)
  direct   <- direct_ranks(prepped)
  pca      <- top5 %>%
    group_by(Response) %>%
    mutate(pca_rank = row_number()) %>%                   #top-5 tables are sorted by Total
    ungroup()
  
  details <- pca %>%
    left_join(direct, by = c("Response", "Variable")) %>%
    mutate(dataset = name, in_direct_top5 = direct_rank <= 5) %>%
    select(dataset, Response, Variable, pca_rank, Total, direct_rank, r, in_direct_top5)
  
  by_response <- map_dfr(unique(pca$Response), function(rv) {
    pca_vars    <- pca$Variable[pca$Response == rv]
    direct_vars <- direct$Variable[direct$Response == rv & direct$direct_rank <= 5]
    rv_direct   <- direct[direct$Response == rv, ]
    tibble(
      dataset            = name,
      Response           = rv,
      n_shared           = length(intersect(pca_vars, direct_vars)),
      pca_top5           = paste(pca_vars, collapse = ", "),
      direct_top5        = paste(direct_vars, collapse = ", "),
      maxr_within_pca    = max_pairwise_r(pca_vars, clim_cor),
      maxr_within_direct = max_pairwise_r(direct_vars, clim_cor),
      top_absr_direct    = max(abs(rv_direct$r)),                                  #strongest single correlate
      top_absr_pca       = max(abs(rv_direct$r[rv_direct$Variable %in% pca_vars])),
      mean_absr_direct   = mean(abs(rv_direct$r[rv_direct$Variable %in% direct_vars])),
      mean_absr_pca      = mean(abs(rv_direct$r[rv_direct$Variable %in% pca_vars])),
      coverage_r         = max(abs(clim_cor[direct_vars[1], pca_vars]))           #stand-in for the strongest correlate?
    )
  })
  
  list(details = details, by_response = by_response)
}

top5_list <- list(oldY = topoldY,   oldS = topoldS,   oldM = topoldM,
                  youngY = topyoungY, youngS = topyoungS, youngM = topyoungM)

direct_results <- imap(clim_data, ~ compare_one(.y, .x, top5_list[[.y]]))
direct_details     <- map_dfr(direct_results, "details")
direct_by_response <- map_dfr(direct_results, "by_response")

direct_summary <- direct_by_response %>%
  group_by(dataset) %>%
  summarise(
    n_categories            = n(),
    mean_shared_of_5        = mean(n_shared),
    n_distinct_vars_pca     = n_distinct(unlist(strsplit(pca_top5, ", "))),
    n_distinct_vars_direct  = n_distinct(unlist(strsplit(direct_top5, ", "))),
    median_maxr_pca         = median(maxr_within_pca, na.rm = TRUE),
    median_maxr_direct      = median(maxr_within_direct, na.rm = TRUE),
    median_top_absr_direct  = median(top_absr_direct),
    median_top_absr_pca     = median(top_absr_pca),
    median_mean_absr_direct = median(mean_absr_direct),
    median_mean_absr_pca    = median(mean_absr_pca),
    median_coverage_r       = median(coverage_r),
    n_coverage_above_0.8    = sum(coverage_r > 0.8),
    .groups = "drop"
  )
print(direct_summary)

write.csv(direct_by_response %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
          file.path(out_dir, "direct_check_by_response.csv"), row.names = FALSE)
write.csv(direct_details %>% mutate(across(c(Total, r), ~ round(.x, 3))),
          file.path(out_dir, "direct_check_details.csv"), row.names = FALSE)
write.csv(direct_summary, file.path(out_dir, "direct_check_summary.csv"), row.names = FALSE)

#Notes on the outputs:
#n_shared: how many of the PCA top 5 are also in the direct-correlation top 5 (0-5)
#n_distinct_vars_*: number of different variables selected across all herbivory categories
#maxr_within_*: largest |r| between any two variables in a top-5 set (near 1 = near-duplicates)
#top_absr_* / mean_absr_*: strongest and average |r| with herbivory within each top-5 set
#coverage_r: largest |r| between the #1 direct-correlation variable and any PCA-selected variable