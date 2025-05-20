# KG_Intervals: Korn and Graubard Confidence Intervals + NCHS Standards
# Adapted from NCHS SAS macro by Crescent Martin

library(dplyr)
library(readr)

#' Calculate Korn and Graubard Confidence Intervals
calculate_KG_CI <- function(df,
                            prop_var = "proportion",
                            se_var = "se",
                            n_var = "counts",
                            atlev1 = NULL,
                            atlev2 = NULL,
                            proportionTotal = 100,
                            roundingIncrement = 0.1,
                            df_adjust = TRUE) {
  
  # Automatically disable df_adjust if PSU/Strata not available
  if (df_adjust) {
    if (is.null(atlev1) || is.null(atlev2) ||
        !(atlev1 %in% names(df)) || !(atlev2 %in% names(df))) {
      message("Note: atlev1 or atlev2 not found. df_adjust set to FALSE.")
      df_adjust <- FALSE
    }
  }
  
  df <- df %>%
    mutate(
      p = .data[[prop_var]] / proportionTotal,
      sep = .data[[se_var]] / proportionTotal,
      percent = round(p * 100, roundingIncrement),
      sepercent = round(sep * 100, roundingIncrement),
      n_eff = ifelse(p %in% c(0, 1), .data[[n_var]], p * (1 - p) / (sep^2)),
      deffmean = ifelse(p %in% c(0, 1), NA, (sep^2) / (p * (1 - p) / .data[[n_var]])),
      df = if (df_adjust) .data[[atlev2]] - .data[[atlev1]] else NA,
      t_adj = ifelse(df_adjust & !is.na(df) & df > 0,
                     (qt(.975, .data[[n_var]] - 1) / qt(.975, df))^2,
                     1),
      n_eff_df = pmin(.data[[n_var]], n_eff * t_adj),
      x = p * n_eff_df
    ) %>%
    rowwise() %>%
    mutate(
      lowerCL = case_when(
        p == 0 ~ 0,
        p == 1 ~ (1 + ((n_eff_df - x + 1) / (x * qf(0.025, 2 * x, 2 * (n_eff_df - x + 1)))))^-1,
        TRUE ~ (1 + ((n_eff_df - x + 1) / (x * qf(0.025, 2 * x, 2 * (n_eff_df - x + 1)))))^-1
      ),
      upperCL = case_when(
        p == 0 ~ (1 + ((n_eff_df - x) / ((x + 1) * qf(0.975, 2 * (x + 1), 2 * (n_eff_df - x)))))^-1,
        p == 1 ~ 1,
        TRUE ~ (1 + ((n_eff_df - x) / ((x + 1) * qf(0.975, 2 * (x + 1), 2 * (n_eff_df - x)))))^-1
      ),
      CIW = round((upperCL - lowerCL) * 100, roundingIncrement),
      RCIW = round(((upperCL - lowerCL) / p) * 100, roundingIncrement),
      RCIW_complement = round(((upperCL - lowerCL) / (1 - p)) * 100, roundingIncrement),
      suppress1 = pmin(n_eff, .data[[n_var]]) < 30,
      suppress2 = CIW >= 30,
      suppress3 = CIW > 5 & RCIW > 130,
      statReview1 = p %in% c(0, 1),
      statReview2 = if (df_adjust) df < 8 else FALSE,
      suppress = suppress1 | suppress2 | suppress3,
      statReview = !suppress & (statReview1 | statReview2),
      present = !suppress & !statReview,
      fn_complement = present & CIW > 5 & RCIW_complement > 130,
      pct_se = sprintf("%.1f (%.1f)", percent, sepercent),
      pct_ci = sprintf("%.1f (%.1f, %.1f)", percent, lowerCL * 100, upperCL * 100)
    ) %>%
    ungroup()
  
  return(df)
}

# ---- Example usage ----
example <- data.frame(
  Sex = c("Men", "Women"),
  counts = c(1654, 1674),
  proportion = c(5.51, 10.05),
  se = c(0.646, 0.804)
)

# Run the function
result <- calculate_KG_CI(example, df_adjust = FALSE)

# Show key output
print(result %>% select(Sex, pct_se, pct_ci, suppress, statReview, present))

# Export to CSV
write_csv(result, "KG_CI_output.csv")
