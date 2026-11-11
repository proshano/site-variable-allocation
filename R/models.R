# The ten analysis methods for each outcome. Each method is a function of one
# simulated trial that returns the treatment estimate, its SE, and the fit
# status (see R/fit_helpers.R). The list names are the method names stored in
# the results; R/tables.R gives the labels printed in the paper.
#
# The two IPTW methods weight participants by the inverse of their known
# assignment probability and differ only in their variance estimator.

# The two random-intercept models use gamma frailty, the default of frailty(),
# although site effects were generated as normal on the log-hazard scale. They
# also keep coxph's default limit of 10 outer iterations, which it applies
# without a warning; see the notes in README.md.
time_to_event_methods <- list(
  "Unadjusted" = function(data) {
    fit_method(coxph(Surv(time, status) ~ treatment, data = data))
  },
  "Unadjusted (cluster-robust SE)" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment,
      data = data, robust = TRUE, cluster = site_id
    ))
  },
  "IPTW (robust SE)" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment,
      data = data, weights = iptw_weight, robust = TRUE
    ))
  },
  "IPTW (cluster-robust SE)" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment,
      data = data, weights = iptw_weight, robust = TRUE, cluster = site_id
    ))
  },
  "Allocation ratio covariate" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment + allocation_ratio,
      data = data
    ))
  },
  "Stratified by allocation ratio" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment + strata(allocation_ratio),
      data = data
    ))
  },
  "Allocation ratio covariate + site random intercept" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment + allocation_ratio + frailty(site_id),
      data = data
    ))
  },
  "Site covariate" = function(data) {
    fit_method(coxph(Surv(time, status) ~ treatment + site_id, data = data))
  },
  "Stratified by site" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment + strata(site_id),
      data = data
    ))
  },
  "Site random intercept" = function(data) {
    fit_method(coxph(
      Surv(time, status) ~ treatment + frailty(site_id),
      data = data
    ))
  }
)

# The two stratified binary analyses use exact conditional logistic regression.
# Both random-intercept models use the bobyqa optimizer with a raised limit on
# function evaluations.
binary_methods <- list(
  "Unadjusted" = function(data) {
    fit_method(glm(
      binary_outcome ~ treatment,
      data = data, family = binomial()
    ))
  },
  "Unadjusted (cluster-robust SE)" = function(data) {
    fit_with_sandwich_se(
      glm(binary_outcome ~ treatment, data = data, family = binomial()),
      clusters = data$site_id
    )
  },
  "IPTW (robust SE)" = function(data) {
    fit_with_sandwich_se(glm(
      binary_outcome ~ treatment,
      data = data, family = quasibinomial(), weights = iptw_weight
    ))
  },
  "IPTW (cluster-robust SE)" = function(data) {
    fit_with_sandwich_se(
      glm(
        binary_outcome ~ treatment,
        data = data, family = quasibinomial(), weights = iptw_weight
      ),
      clusters = data$site_id
    )
  },
  "Allocation ratio covariate" = function(data) {
    fit_method(glm(
      binary_outcome ~ treatment + allocation_ratio,
      data = data, family = binomial()
    ))
  },
  "Stratified by allocation ratio" = function(data) {
    fit_method(clogit(
      binary_outcome ~ treatment + strata(allocation_ratio),
      data = data, method = "exact"
    ))
  },
  "Allocation ratio covariate + site random intercept" = function(data) {
    fit_method(lme4::glmer(
      binary_outcome ~ treatment + allocation_ratio + (1 | site_id),
      data = data, family = binomial(), control = glmer_control()
    ))
  },
  "Site covariate" = function(data) {
    fit_method(glm(
      binary_outcome ~ treatment + site_id,
      data = data, family = binomial()
    ))
  },
  "Stratified by site" = function(data) {
    fit_method(clogit(
      binary_outcome ~ treatment + strata(site_id),
      data = data, method = "exact"
    ))
  },
  "Site random intercept" = function(data) {
    fit_method(lme4::glmer(
      binary_outcome ~ treatment + (1 | site_id),
      data = data, family = binomial(), control = glmer_control()
    ))
  }
)

glmer_control <- function() {
  lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
}

# The two "(Pooled)" methods combine within-stratum mean differences rather
# than fitting one stratified model; see fit_pooled_within_strata().
continuous_methods <- list(
  "Unadjusted" = function(data) {
    fit_method(lm(continuous_outcome ~ treatment, data = data))
  },
  "Unadjusted (cluster-robust SE)" = function(data) {
    fit_with_sandwich_se(
      lm(continuous_outcome ~ treatment, data = data),
      clusters = data$site_id
    )
  },
  "IPTW (robust SE)" = function(data) {
    fit_with_sandwich_se(lm(
      continuous_outcome ~ treatment,
      data = data, weights = iptw_weight
    ))
  },
  "IPTW (cluster-robust SE)" = function(data) {
    fit_with_sandwich_se(
      lm(continuous_outcome ~ treatment, data = data, weights = iptw_weight),
      clusters = data$site_id
    )
  },
  "Allocation ratio covariate" = function(data) {
    fit_method(lm(
      continuous_outcome ~ treatment + allocation_ratio,
      data = data
    ))
  },
  "Stratified by allocation ratio (Pooled)" = function(data) {
    fit_pooled_within_strata(data, data$allocation_ratio)
  },
  "Allocation ratio covariate + site random intercept" = function(data) {
    fit_method(lme4::lmer(
      continuous_outcome ~ treatment + allocation_ratio + (1 | site_id),
      data = data, REML = FALSE
    ))
  },
  "Site covariate" = function(data) {
    fit_method(lm(continuous_outcome ~ treatment + site_id, data = data))
  },
  "Stratified by site (Pooled)" = function(data) {
    fit_pooled_within_strata(data, data$site_id)
  },
  "Site random intercept" = function(data) {
    fit_method(lme4::lmer(
      continuous_outcome ~ treatment + (1 | site_id),
      data = data, REML = FALSE
    ))
  }
)

analysis_methods <- list(
  `Time-to-Event` = time_to_event_methods,
  Binary = binary_methods,
  Continuous = continuous_methods
)

# Fit every method for the trial's outcome. One row per method.
analyze_trial <- function(trial, outcome) {
  fits <- lapply(analysis_methods[[outcome]], function(method) method(trial))
  dplyr::bind_rows(fits, .id = "method")
}

# Simulate one replicate of a scenario and fit every method.
run_replicate <- function(scenario, replicate_index) {
  trial <- simulate_trial(scenario, replicate_seed(scenario, replicate_index))
  cbind(replicate = replicate_index, analyze_trial(trial, scenario$outcome))
}
