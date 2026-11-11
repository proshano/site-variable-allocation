# Interval checks for IPTW with cluster-robust SEs (Table S41). In the six
# scenarios of Tables 3 and 4, this regenerates the simulated trials and
# refits IPTW with cluster-robust SEs.
#
#   Part 1 compares normal and t(29) critical values, in all six scenarios.
#   Part 2 uses the CR2 covariance with Satterthwaite degrees of freedom, in
#          the four binary and continuous scenarios. CR2 was not applied to
#          the Cox models.
#
# Run it with: Rscript run_interval_checks.R
#
# The default uses 100 replicates per scenario and takes under a minute. The
# paper used 10,000, which takes about half an hour on one core:
#
#   Rscript run_interval_checks.R --replicates=10000 --cores=4
#
# Results are written to output/interval_checks/. The paper's results are in
# data/saved_results/.
source("R/load.R")
if (!requireNamespace("clubSandwich", quietly = TRUE)) {
  stop("Install the clubSandwich package first. See README.md.")
}

settings <- parse_options(list(
  replicates = "100", cores = "1", output = "output/interval_checks"
))
n_replicates <- parse_count(settings$replicates, 10000L)
cores <- parse_count(settings$cores, 128L)
output <- output_directory(settings$output)

scenarios <- read_scenarios()
saved_summary <- read_saved_results("simulation_summary.csv")

# Part 1: normal and t(29) intervals ------------------------------------

# Coverage of both intervals in one scenario, with the saved simulation
# results for the same method alongside.
check_normal_and_t <- function(task) {
  scenario <- get_scenario(scenarios, task)
  fit_iptw_cluster <-
    analysis_methods[[scenario$outcome]][["IPTW (cluster-robust SE)"]]
  estimates <- ses <- rep(NA_real_, n_replicates)
  for (i in seq_len(n_replicates)) {
    trial <- simulate_trial(scenario, replicate_seed(scenario, i))
    fit <- fit_iptw_cluster(trial)
    estimates[i] <- fit$estimate
    ses[i] <- fit$se
  }
  stopifnot(all(is.finite(estimates) & is.finite(ses) & ses > 0))

  error <- abs(estimates - scenario$true_effect)
  z_coverage <- mean(error <= qnorm(0.975) * ses)
  t_coverage <- mean(
    error <= qt(0.975, df = scenario$trial_size_sites - 1L) * ses
  )
  saved <- saved_summary[
    saved_summary$task_id == task &
      saved_summary$method == "IPTW (cluster-robust SE)",
  ]
  summary_difference <- max(abs(
    c(mean(estimates), sd(estimates), mean(ses), z_coverage) -
      c(saved$mean_estimate_model, saved$empirical_se, saved$model_based_se,
        saved$coverage)
  ))
  data.frame(
    task_id = task, outcome = scenario$outcome,
    effect = scenario$effect_scenario,
    replicates = n_replicates, sites = scenario$trial_size_sites,
    estimate = mean(estimates), empirical_sd = sd(estimates),
    mean_se = mean(ses),
    normal_coverage = z_coverage, t_coverage = t_coverage,
    normal_mcse = sqrt(z_coverage * (1 - z_coverage) / n_replicates),
    t_mcse = sqrt(t_coverage * (1 - t_coverage) / n_replicates),
    submitted_estimate = saved$mean_estimate_model,
    submitted_normal_coverage = saved$coverage,
    max_abs_summary_difference = summary_difference
  )
}

normal_and_t <- parallel_tasks(
  c(82L, 83L, 88L, 89L, 94L, 95L), check_normal_and_t, cores,
  exports = c("scenarios", "n_replicates", "saved_summary")
)
normal_and_t <- dplyr::bind_rows(normal_and_t)
readr::write_csv(normal_and_t, file.path(output, "interval_checks.csv"))
print(normal_and_t[, c(
  "task_id", "outcome", "effect", "replicates", "normal_coverage", "t_coverage"
)], row.names = FALSE, digits = 3)


# Part 2: CR2 covariance with Satterthwaite degrees of freedom ----------

# The same weighted models as the IPTW methods in R/models.R.
fit_weighted_model <- function(trial, outcome) {
  if (outcome == "Binary") {
    return(glm(
      binary_outcome ~ treatment, data = trial,
      family = quasibinomial(), weights = iptw_weight
    ))
  }
  lm(continuous_outcome ~ treatment, data = trial, weights = iptw_weight)
}

# CR2 on all participants is slow. Within each site and arm, predictors and
# weights are identical, so this fits the cell means instead, with weights
# multiplied by the cell size and working variances divided by it. That keeps
# the estimating equations, CR2 covariance, and degrees of freedom unchanged.
# check_one_replicate() confirms this against the full data for 20 replicates.
fit_cell_model <- function(trial, outcome, starting_values) {
  trial$y <- if (outcome == "Binary") {
    trial$binary_outcome
  } else {
    trial$continuous_outcome
  }
  cells <- trial |>
    dplyr::group_by(site_id, treatment, iptw_weight) |>
    dplyr::summarise(n = dplyr::n(), y = mean(y), .groups = "drop")
  if (outcome == "Binary") {
    fit <- glm(
      y ~ treatment, data = cells,
      family = quasibinomial(), weights = n * iptw_weight,
      start = starting_values, control = glm.control(epsilon = 1e-12)
    )
    target <- fit$family$variance(fitted(fit)) / cells$n
  } else {
    fit <- lm(y ~ treatment, data = cells, weights = n * iptw_weight)
    target <- 1 / cells$n
  }
  list(fit = fit, target = target, cluster = cells$site_id)
}

check_one_replicate <- function(scenario, replicate_index) {
  outcome <- scenario$outcome
  trial <- simulate_trial(scenario, replicate_seed(scenario, replicate_index))
  # `fit` is the original analysis: it gives the estimate and original SE.
  # `cells$fit` gives the same coefficients and is used for CR2.
  fit <- fit_weighted_model(trial, outcome)
  original_se <- sandwich_se(fit, trial$site_id)
  cells <- fit_cell_model(trial, outcome, coef(fit))
  coefficient_difference <- max(abs(coef(fit) - coef(cells$fit)))
  stopifnot(coefficient_difference < 1e-7)
  # These are probability weights. Do not use the glm default working
  # variance divided by prior weights, which treats them as precision weights.
  cr2 <- clubSandwich::vcovCR(
    cells$fit, cluster = cells$cluster, type = "CR2",
    inverse_var = FALSE, target = cells$target
  )
  inference <- clubSandwich::coef_test(
    cells$fit, vcov = cr2, coefs = "treatment", test = "Satterthwaite"
  )

  # For the first 20 replicates, confirm the cell aggregation against the
  # full individual-level data.
  full_data <- c(cr0 = NA_real_, cr2 = NA_real_, df = NA_real_)
  if (replicate_index <= 20L) {
    full_data <- full_data_differences(trial, outcome, fit, cr2, inference)
  }

  values <- c(
    estimate = unname(coef(fit)["treatment"]),
    original_se = original_se, cr2_se = inference$SE,
    cr2_df = inference$df_Satt
  )
  stopifnot(all(is.finite(values)), all(values[-1] > 0))
  data.frame(
    task_id = scenario$task_id, replicate = replicate_index,
    estimate = values["estimate"], original_se = values["original_se"],
    cr2_se = values["cr2_se"], cr2_df = values["cr2_df"],
    coefficient_difference = coefficient_difference,
    cr0_relative_difference = full_data[["cr0"]],
    cr2_relative_difference = full_data[["cr2"]],
    df_difference = full_data[["df"]]
  )
}

# Repeat the calculation on the individual-level data and return how far the
# cell-mean results are from it: the CR2 covariance (relative difference) and
# the degrees of freedom. The first value compares clubSandwich's CR0
# covariance with sandwich's cluster covariance on the same data. All three
# must agree.
full_data_differences <- function(trial, outcome, fit, cr2, inference) {
  if (outcome == "Binary") {
    # Restarting at the fitted coefficients refreshes glm's saved IRLS
    # weights, which otherwise describe the penultimate iteration.
    full_fit <- glm(
      binary_outcome ~ treatment, data = trial, family = quasibinomial(),
      weights = iptw_weight, start = coef(fit),
      control = glm.control(epsilon = 1e-12)
    )
    target <- full_fit$family$variance(fitted(full_fit))
  } else {
    full_fit <- fit
    target <- rep(1, nrow(trial))
  }
  cr0 <- clubSandwich::vcovCR(
    full_fit, cluster = trial$site_id, type = "CR0",
    inverse_var = FALSE, target = target
  )
  hc0 <- sandwich::vcovCL(
    full_fit, cluster = trial$site_id, type = "HC0", cadjust = FALSE
  )
  full_cr2 <- clubSandwich::vcovCR(
    full_fit, cluster = trial$site_id, type = "CR2",
    inverse_var = FALSE, target = target
  )
  full_inference <- clubSandwich::coef_test(
    full_fit, vcov = full_cr2, coefs = "treatment", test = "Satterthwaite"
  )
  differences <- c(
    cr0 = max(abs(cr0 - hc0)) / max(abs(hc0)),
    cr2 = max(abs(full_cr2 - cr2)) / max(abs(full_cr2)),
    df = abs(full_inference$df_Satt - inference$df_Satt)
  )
  stopifnot(all(differences < 1e-7))
  differences
}

summarize_cr2_task <- function(records, scenario) {
  saved <- saved_summary[
    saved_summary$task_id == scenario$task_id &
      saved_summary$method == "IPTW (cluster-robust SE)",
  ]
  error <- abs(records$estimate - scenario$true_effect)
  z_coverage <- mean(error <= qnorm(0.975) * records$original_se)
  t_coverage <- mean(
    error <= qt(0.975, scenario$trial_size_sites - 1L) * records$original_se
  )
  cr2_coverage <- mean(error <= qt(0.975, records$cr2_df) * records$cr2_se)
  difference <- max(abs(
    c(mean(records$estimate), sd(records$estimate),
      mean(records$original_se), z_coverage) -
      c(saved$mean_estimate_model, saved$empirical_se,
        saved$model_based_se, saved$coverage)
  ))
  # A full run must reproduce the saved summaries.
  if (n_replicates == 10000L) stopifnot(difference < 1e-8)
  data.frame(
    task_id = scenario$task_id, outcome = scenario$outcome,
    effect = scenario$effect_scenario, replicates = n_replicates,
    sites = scenario$trial_size_sites, estimate = mean(records$estimate),
    empirical_sd = sd(records$estimate), mean_se = mean(records$original_se),
    normal_coverage = z_coverage, t_coverage = t_coverage,
    cr2_mean_se = mean(records$cr2_se), cr2_coverage = cr2_coverage,
    cr2_mcse = sqrt(cr2_coverage * (1 - cr2_coverage) / n_replicates),
    cr2_df_mean = mean(records$cr2_df), cr2_df_min = min(records$cr2_df),
    cr2_df_max = max(records$cr2_df), max_abs_summary_difference = difference
  )
}

check_cr2 <- function(task) {
  scenario <- get_scenario(scenarios, task)
  records <- vector("list", n_replicates)
  for (i in seq_len(n_replicates)) {
    records[[i]] <- check_one_replicate(scenario, i)
  }
  records <- dplyr::bind_rows(records)
  list(summary = summarize_cr2_task(records, scenario), records = records)
}

cr2_checks <- parallel_tasks(
  c(88L, 89L, 94L, 95L), check_cr2, cores,
  exports = c(
    "scenarios", "n_replicates", "saved_summary", "fit_weighted_model",
    "fit_cell_model", "check_one_replicate", "full_data_differences",
    "summarize_cr2_task"
  )
)
cr2_summary <- dplyr::bind_rows(lapply(cr2_checks, `[[`, "summary"))
readr::write_csv(cr2_summary, file.path(output, "interval_checks_cr2.csv"))
print(cr2_summary[, c(
  "task_id", "outcome", "effect", "replicates", "normal_coverage",
  "cr2_coverage", "cr2_df_mean"
)], row.names = FALSE, digits = 3)

# How far the cell-mean calculation was from the individual-level data in the
# replicates where both were done.
records <- dplyr::bind_rows(lapply(cr2_checks, `[[`, "records"))
message(
  "Largest difference from the individual-level calculation: ",
  signif(max(records$cr2_relative_difference, na.rm = TRUE), 2),
  " (CR2 covariance, relative) and ",
  signif(max(records$df_difference, na.rm = TRUE), 2), " (degrees of freedom)"
)

