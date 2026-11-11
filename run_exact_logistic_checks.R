# Exact conditional logistic regression checks (Supplementary Appendix,
# section 6.4). In the simulation,
# exact conditional logistic regression stratified by allocation ratio
# returned no estimate in the scenarios with large strata. This script
# regenerates replicates 1 to 3 of five binary scenarios and fits each trial
# three ways:
#
#   1. the exact method as in the simulation, survival::clogit();
#   2. the same conditional likelihood, evaluated in log space; and
#   3. logistic regression with allocation ratio as a covariate.
#
# Run it with: Rscript run_exact_logistic_checks.R
#
# It takes a few seconds and writes two files to output/exact_logistic_checks/.
# The paper's results are in data/saved_results/.
source("R/load.R")

settings <- parse_options(list(output = "output/exact_logistic_checks"))
output <- output_directory(settings$output)
scenarios <- read_scenarios()

# Participants (n) and events (x) by arm (1 or 0) within each allocation ratio.
binary_counts <- function(trial) {
  ratio_groups <- split(trial, trial$allocation_ratio)
  dplyr::bind_rows(lapply(ratio_groups, function(group) {
    data.frame(
      ratio = as.character(group$allocation_ratio[1]),
      n1 = sum(group$treatment == 1L), n0 = sum(group$treatment == 0L),
      x1 = sum(group$binary_outcome[group$treatment == 1L]),
      x0 = sum(group$binary_outcome[group$treatment == 0L])
    )
  }))
}

# The same conditional likelihood as clogit(method = "exact"), evaluated in
# log space. With one binary predictor, each stratum's likelihood depends only
# on how many of its events are in the intervention arm, so the calculation
# lists those counts rather than subsets of participants. Subtracting the
# largest log term before exponentiating avoids forming the enormous
# denominator.
#
# Returns the negative log likelihood at `beta` and the information, which is
# the variance of the intervention-arm event count, summed over strata.
conditional_likelihood <- function(beta, counts) {
  log_denominator <- numeric(nrow(counts))
  information <- numeric(nrow(counts))
  for (i in seq_len(nrow(counts))) {
    n1 <- counts$n1[i]
    n0 <- counts$n0[i]
    events <- counts$x1[i] + counts$x0[i]
    # Every possible number of events in the intervention arm.
    k <- seq.int(max(0L, events - n0), min(events, n1))
    log_term <- lchoose(n1, k) + lchoose(n0, events - k) + beta * k
    largest <- max(log_term)
    term <- exp(log_term - largest)
    probability <- term / sum(term)
    mean_k <- sum(probability * k)
    log_denominator[i] <- largest + log(sum(term))
    information[i] <- sum(probability * (k - mean_k)^2)
  }
  list(
    negative_loglik = sum(log_denominator) - beta * sum(counts$x1),
    information = sum(information)
  )
}

# The estimate minimizes the negative log likelihood.
stable_conditional_fit <- function(counts) {
  fit <- optimize(
    function(beta) conditional_likelihood(beta, counts)$negative_loglik,
    c(-15, 15), tol = 1e-10
  )
  information <- conditional_likelihood(fit$minimum, counts)$information
  stopifnot(abs(fit$minimum) < 14.9, information > 0)
  c(estimate = fit$minimum, se = sqrt(1 / information), loglik = -fit$objective)
}

fit_rows <- list()
stratum_rows <- list()
for (task in c(16L, 18L, 52L, 88L, 89L)) {
  scenario <- get_scenario(scenarios, task)
  for (replicate_index in 1:3) {
    seed <- replicate_seed(scenario, replicate_index)
    trial <- simulate_trial(scenario, seed)
    counts <- binary_counts(trial)
    stable <- stable_conditional_fit(counts)
    native <- capture_fit_with_status(clogit(
      binary_outcome ~ treatment + strata(allocation_ratio),
      data = trial, method = "exact"
    ))
    native_estimate <- unname(coef(native$model)[1])
    native_se <- sqrt(vcov(native$model)[1, 1])
    native_loglik <- tail(native$model$loglik, 1)
    # Where the native fit succeeded, the log-space fit must agree with it.
    if (is.finite(native_estimate) && is.finite(native_se) && native_se > 0) {
      stopifnot(abs(stable["estimate"] - native_estimate) < 1e-5)
      stopifnot(abs(stable["se"] - native_se) < 1e-5)
      stopifnot(abs(stable["loglik"] - native_loglik) < 1e-5)
    }
    companion <- glm(
      binary_outcome ~ treatment + allocation_ratio,
      data = trial, family = binomial()
    )

    counts$task_id <- task
    counts$replicate <- replicate_index
    counts$log_denominator_at_zero <- lchoose(
      counts$n1 + counts$n0, counts$x1 + counts$x0
    )
    stratum_rows[[length(stratum_rows) + 1L]] <- counts
    fit_rows[[length(fit_rows) + 1L]] <- data.frame(
      task_id = task, replicate = replicate_index,
      sites = scenario$trial_size_sites,
      native_status = native$status, native_message = native$message,
      native_estimate = native_estimate, native_se = native_se,
      native_loglik = native_loglik,
      native_iterations = native$model$iter,
      stable_estimate = stable["estimate"], stable_se = stable["se"],
      stable_loglik = stable["loglik"],
      glm_estimate = unname(coef(companion)["treatment"]),
      glm_se = summary(companion)$coef["treatment", "Std. Error"],
      all_four_cells_positive = all(
        counts$x1 > 0 & counts$x0 > 0 &
          counts$x1 < counts$n1 & counts$x0 < counts$n0
      )
    )
  }
}
fits <- dplyr::bind_rows(fit_rows)
readr::write_csv(fits, file.path(output, "exact_logistic_fits.csv"))
readr::write_csv(
  dplyr::bind_rows(stratum_rows),
  file.path(output, "exact_logistic_strata.csv")
)

# The odds ratio from each of the three fits. Tables S5 and S9 are tasks 16
# and 18, Table S17 is task 52, and Tables 3 and 4 are tasks 88 and 89.
odds_ratios <- data.frame(
  task_id = fits$task_id,
  replicate = fits$replicate,
  exact_clogit = ifelse(
    is.finite(fits$native_estimate),
    sprintf("%.4f", exp(fits$native_estimate)), "Not estimated"
  ),
  exact_log_scale = sprintf("%.4f", exp(fits$stable_estimate)),
  ratio_covariate = sprintf("%.4f", exp(fits$glm_estimate))
)
print(odds_ratios, row.names = FALSE)
message("Results: ", output)
