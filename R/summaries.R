# Method performance within one scenario. Every method is compared with the
# generating within-site treatment effect.

# Scenario columns repeated in every summary row.
summary_scenario_columns <- c(
  "scenario_id", "scenario_set", "scenario_family", "scenario_family_label",
  "task_id", "trial_size_sites", "total_n", "mean_site_size",
  "median_site_size", "min_site_size", "max_site_size", "site_size_cv",
  "shifted_site_fraction", "shifted_site_count", "site_size_pattern",
  "size_allocation_association", "prognosis_strength", "outcome_frequency",
  "scenario_list_hash", "outcome", "effect_scenario", "residual_variance"
)

# One row per method. `results` holds one row per replicate and method.
summarize_task <- function(results, scenario) {
  critical_value <- qnorm(0.975)
  true_effect <- scenario$true_effect

  performance <- results |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      true_effect_model = true_effect,
      n_attempted = dplyr::n(),
      # n_converged is the original name for fits with status "ok".
      n_converged = sum(status == "ok"),
      n_warning = sum(status == "warning"),
      n_nonconvergence = sum(status == "nonconvergence"),
      n_error = sum(status == "error"),
      n_successful = sum(is.finite(estimate) & is.finite(se)),
      modal_status_message = most_frequent_message(message),
      mean_estimate_model = mean(estimate, na.rm = TRUE),
      empirical_se = stats::sd(estimate, na.rm = TRUE),
      # The mean of the reported SEs, including sandwich SEs where used.
      model_based_se = mean(se, na.rm = TRUE),
      coverage = mean(
        estimate - critical_value * se <= true_effect &
          estimate + critical_value * se >= true_effect,
        na.rm = TRUE
      ),
      rejection_rate = mean(abs(estimate / se) >= critical_value, na.rm = TRUE)
    ) |>
    dplyr::mutate(
      # Monte Carlo SEs of the coverage and rejection proportions.
      mc_se_coverage = sqrt(coverage * (1 - coverage) / n_successful),
      mc_se_rejection = sqrt(
        rejection_rate * (1 - rejection_rate) / n_successful
      ),
      # HRs and ORs are exponentiated mean log estimates.
      mean_estimate_display = to_display_scale(mean_estimate_model, scenario),
      true_effect_display = to_display_scale(true_effect, scenario)
    )
  dplyr::bind_cols(
    as.data.frame(scenario[summary_scenario_columns]), performance
  )
}

to_display_scale <- function(value, scenario) {
  if (scenario$outcome == "Continuous") value else exp(value)
}

# The most frequent nonempty message, shortened to 120 characters.
most_frequent_message <- function(messages) {
  messages <- messages[!is.na(messages) & nzchar(messages)]
  if (length(messages) == 0L) {
    return(NA_character_)
  }
  message_counts <- sort(table(messages), decreasing = TRUE)
  substr(names(message_counts)[1], 1, 120)
}
