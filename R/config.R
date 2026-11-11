# Fixed parameters of the simulated trials. Each row of data/scenarios.csv
# selects among them, and R/scenarios.R combines the two.

# The true treatment effect on the model scale: log hazard ratio, log odds
# ratio, or mean difference. "Realistic" is the stored name of the paper's
# moderate effect. Scenario names enter the seeds, so they cannot be changed.
treatment_effects <- list(
  `Time-to-Event` = c(Null = log(1), Realistic = log(0.8), Large = log(0.5)),
  Binary = c(Null = log(1), Realistic = log(0.8), Large = log(0.5)),
  Continuous = c(Null = 0, Realistic = -20, Large = -40)
)

# Outcome level at 1:1 sites, by event frequency.
baseline_hazard <- c(Common = 0.10, Sparse = 0.03)
baseline_risk <- c(Common = 0.20, Sparse = 0.05)
baseline_mean <- 100

# Every site receives a random intercept with this SD.
site_intercept_sd <- c(`Time-to-Event` = 0.35, Binary = 0.50, Continuous = 10)

# Under the Strong prognosis setting, each 1:3 site also receives its own
# adverse shift, centered on these values with these SDs.
strong_hazard_ratio <- 2.5
strong_hazard_shift_sd <- 0.40
strong_risk <- 0.40
strong_risk_shift_sd <- 0.60
strong_mean <- 120
strong_mean_shift_sd <- 15

# Residual SD of continuous outcomes, by allocation ratio.
residual_sds <- list(
  Equal = c(`1:1` = 20, `1:3` = 20),
  Differential = c(`1:1` = 20, `1:3` = 30)
)

# Censoring times are uniform between 0 and this time.
max_follow_up <- 5
