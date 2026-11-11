# Generate one simulated trial. These are the data-generating equations in the
# Supplementary Appendix. For patient i at site s:
#
#   time to event:  hazard = h0 * exp(delta_s + b_s + beta * A_is)
#   binary:         logit(risk) = alpha + delta_s + b_s + beta * A_is
#   continuous:     Y_is = mu + delta_s + b_s + beta * A_is + residual error
#
# A_is is treatment (1 = intervention), beta the true treatment effect, b_s the
# site's random intercept, and delta_s its prognosis shift. Only 1:3 sites have
# a shift, and only under the Strong prognosis setting. Event times are
# exponential and censored uniformly over follow-up.
#
# Each site's intervention and control totals are fixed by its allocation
# ratio. Randomness enters through the site effects, the order of assignments
# within each site, and the outcomes.

simulate_trial <- function(scenario, seed) {
  set.seed(seed)
  sites <- scenario$sites

  # Keep these random draws in this order, and keep the terms of the sums in
  # simulate_outcome() in their order. Otherwise a replicate no longer
  # reproduces the same trial as the original study.
  sites$prognosis_shift <- draw_prognosis_shifts(scenario)
  sites$random_intercept <- rnorm(nrow(sites), 0, scenario$site_intercept_sd)

  site_trials <- lapply(seq_len(nrow(sites)), function(i) {
    site <- sites[i, ]
    # The site's fixed arm totals, in random order.
    treatment <- sample(c(rep(1L, site$n_treatment), rep(0L, site$n_control)))
    data.frame(
      site_id = site$site_id,
      allocation_ratio = site$allocation_ratio,
      treatment_share = site$treatment_share,
      treatment = treatment,
      simulate_outcome(site, treatment, scenario)
    )
  })

  trial <- dplyr::bind_rows(site_trials)
  trial$site_id <- factor(trial$site_id)
  trial$allocation_ratio <- factor(trial$allocation_ratio)
  # Assignment probabilities are known by design, so the inverse probability
  # weights use each site's treatment share, not an estimated propensity score.
  trial$iptw_weight <- ifelse(
    trial$treatment == 1L,
    1 / trial$treatment_share,
    1 / (1 - trial$treatment_share)
  )
  trial
}

draw_prognosis_shifts <- function(scenario) {
  sites <- scenario$sites
  shifts <- numeric(nrow(sites))
  if (scenario$prognosis_strength == "Strong") {
    at_1_3 <- sites$allocation_ratio == "1:3"
    shifts[at_1_3] <- rnorm(sum(at_1_3), scenario$shift$mean, scenario$shift$sd)
  }
  shifts
}

# The outcome column or columns for the patients at one site.
simulate_outcome <- function(site, treatment, scenario) {
  n <- length(treatment)
  effect <- treatment * scenario$true_effect

  if (scenario$outcome == "Time-to-Event") {
    hazard <- scenario$baseline *
      exp(site$prognosis_shift + site$random_intercept + effect)
    event_time <- rexp(n, rate = hazard)
    censoring_time <- runif(n, 0, max_follow_up)
    return(data.frame(
      time = pmin(event_time, censoring_time),
      status = as.integer(event_time <= censoring_time)
    ))
  }
  if (scenario$outcome == "Binary") {
    log_odds <- qlogis(scenario$baseline) +
      site$prognosis_shift + site$random_intercept + effect
    return(data.frame(binary_outcome = rbinom(n, 1, plogis(log_odds))))
  }
  mean_outcome <- scenario$baseline +
    site$prognosis_shift + site$random_intercept + effect
  residual_sd <- scenario$residual_sd[[site$allocation_ratio]]
  data.frame(continuous_outcome = rnorm(n, mean_outcome, residual_sd))
}
