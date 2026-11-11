# The 99 simulation scenarios. Each row of data/scenarios.csv is one scenario
# with one outcome. A scenario's sites are fixed: only the replicate seed
# changes from one simulated trial to the next.

# The saved list defines the study, and its hash enters every replicate's seed
# (see replicate_seed()), so stop if the list has changed. The hash is
# recomputed exactly as the original study serialized the list, with the
# field names it used.
read_scenarios <- function(path = "data/scenarios.csv") {
  expected_hash <-
    "24a4b165f97decf6dbe9335d09ba505f203aec2c95983ec28a4714d88d9267c1"
  scenarios <- as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  scenario_json <- jsonlite::toJSON(
    list(
      manifest_version = "expanded_two_ratio_v2_targeted",
      tasks = scenarios[, names(scenarios) != "scenario_list_hash"]
    ),
    auto_unbox = TRUE, dataframe = "rows", null = "null", digits = NA
  )
  actual_hash <- digest::digest(
    scenario_json, algo = "sha256", serialize = FALSE
  )
  if (actual_hash != expected_hash) {
    stop("The saved scenario list has changed.")
  }
  scenarios
}

# One scenario as a list: its row of data/scenarios.csv, plus its sites and
# the parameter values for its outcome.
get_scenario <- function(scenarios, task_id) {
  row <- scenarios[scenarios$task_id == task_id, ]
  if (nrow(row) != 1L) {
    stop("Unknown task ID: ", task_id)
  }
  outcome <- row$outcome

  scenario <- as.list(row)
  scenario$sites <- site_layout(row)
  scenario$true_effect <- treatment_effects[[outcome]][[row$effect_scenario]]
  scenario$site_intercept_sd <- site_intercept_sd[[outcome]]
  scenario$baseline <- switch(outcome,
    `Time-to-Event` = baseline_hazard[[row$outcome_frequency]],
    Binary = baseline_risk[[row$outcome_frequency]],
    Continuous = baseline_mean
  )
  # On the model scale: log hazard, log odds, or outcome units.
  scenario$shift <- switch(outcome,
    `Time-to-Event` = list(
      mean = log(strong_hazard_ratio), sd = strong_hazard_shift_sd
    ),
    Binary = list(
      mean = qlogis(strong_risk) - qlogis(scenario$baseline),
      sd = strong_risk_shift_sd
    ),
    Continuous = list(
      mean = strong_mean - baseline_mean, sd = strong_mean_shift_sd
    )
  )
  if (outcome == "Continuous") {
    scenario$residual_sd <- residual_sds[[row$residual_variance]]
  }
  scenario
}

# One row per site: its size, allocation ratio, and fixed arm totals.
site_layout <- function(row) {
  sizes <- build_site_sizes(row$trial_size_sites, row$site_size_pattern)
  uses_1_3 <- choose_1_3_sites(
    sizes, row$size_allocation_association, row$shifted_site_fraction
  )
  treatment_share <- ifelse(uses_1_3, 0.25, 0.5)
  layout <- data.frame(
    site_id = seq_along(sizes),
    site_size = sizes,
    allocation_ratio = ifelse(uses_1_3, "1:3", "1:1"),
    treatment_share = treatment_share,
    n_treatment = as.integer(sizes * treatment_share),
    n_control = as.integer(sizes * (1 - treatment_share))
  )
  stopifnot(
    sum(layout$site_size) == row$total_n,
    sum(uses_1_3) == row$shifted_site_count,
    all(layout$n_treatment + layout$n_control == layout$site_size)
  )
  layout
}

# Variable site sizes are evenly spaced quantiles of a log-normal distribution
# with mean 40 and coefficient of variation 1.2. They are rounded to multiples
# of 4 so that 1:1 and 1:3 sites both split into whole numbers of intervention
# and control participants. Evenly spaced quantiles miss the far right tail, so
# the total falls short. Four participants at a time then go to the smallest
# site until the mean is exactly 40. The realized coefficient of variation is
# 0.98 to 1.09, depending on the number of sites.
build_site_sizes <- function(num_sites, site_size_pattern) {
  if (site_size_pattern == "Equal") {
    return(rep(40L, num_sites))
  }
  stopifnot(site_size_pattern == "HighCV")

  sdlog <- sqrt(log(1 + 1.2^2))
  meanlog <- log(40) - sdlog^2 / 2
  quantiles <- (seq_len(num_sites) - 0.5) / num_sites
  raw_sizes <- qlnorm(quantiles, meanlog = meanlog, sdlog = sdlog)
  sizes <- pmax(4L, as.integer(round(raw_sizes / 4) * 4L))

  while (sum(sizes) < 40L * num_sites) {
    smallest <- which.min(sizes)
    sizes[smallest] <- sizes[smallest] + 4L
  }
  sizes
}

# Choose which sites use 1:3 allocation. When sites differ in size, these are
# the largest or the smallest sites, with ties kept in site order. When sites
# are of equal size, they are the last sites. The stored names of this fraction
# use "shifted" because 1:3 sites are the ones given a prognosis shift under
# the Strong setting.
choose_1_3_sites <- function(sizes, size_allocation_association, fraction) {
  num_sites <- length(sizes)
  count <- as.integer(round(num_sites * fraction))
  size_order <- order(sizes, seq_along(sizes))
  chosen <- switch(size_allocation_association,
    Not_Applicable = seq.int(num_sites - count + 1L, num_sites),
    Ratio13Larger = utils::tail(size_order, count),
    Ratio13Smaller = utils::head(size_order, count),
    stop("Unknown association: ", size_allocation_association)
  )
  seq_len(num_sites) %in% chosen
}

# A replicate's seed depends only on the scenario and the replicate number, so
# results do not depend on the number of workers or on how a run is split.
# Seven hexadecimal digits of the hash give a 28-bit integer, as in the
# original study.
replicate_seed <- function(scenario, replicate_index) {
  seed_key <- paste(
    scenario$scenario_list_hash, scenario$scenario_id, scenario$task_id,
    replicate_index, sep = "::"
  )
  seed_hash <- digest::digest(seed_key, algo = "xxhash64", serialize = FALSE)
  strtoi(substr(seed_hash, 1L, 7L), base = 16L)
}
