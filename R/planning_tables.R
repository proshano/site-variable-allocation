# Planning calculations for Tables 1, 2, and S1.
#
# A site can recruit at most `full_capacity` participants and can deliver the
# intervention to at most `intervention_limit` of them. With intervention
# probability (treatment share) p, it recruits min(full_capacity,
# intervention_limit / p) participants. Information is counted in
# participants under 1:1 allocation, 4 n p (1 - p), which equals n under 1:1
# allocation. This benchmark is exact for continuous outcomes with a common
# residual variance. It holds near the null for binary and time-to-event
# outcomes when the event risk or event fraction is common across sites. It is
# not an endpoint-specific power calculation.

information_target <- 10000

site_recruitment <- function(full_capacity, intervention_limit,
                             treatment_share) {
  pmin(full_capacity, intervention_limit / treatment_share)
}

site_information <- function(full_capacity, intervention_limit,
                             treatment_share) {
  recruitment <- site_recruitment(
    full_capacity, intervention_limit, treatment_share
  )
  4 * recruitment * treatment_share * (1 - treatment_share)
}

# The most balanced treatment share that still fills a site. This is 0.5
# (1:1) when the intervention limit allows it and intervention_limit /
# full_capacity otherwise.
capacity_filling_share <- function(full_capacity, intervention_limit) {
  pmin(0.5, intervention_limit / full_capacity)
}

# Sites needed to reach the information target, rounded up. The offset keeps
# floating-point error from rounding an exact whole number up.
sites_required <- function(information_per_site) {
  ceiling(information_target / information_per_site - 1e-10)
}

# The whole-number ratio with the smallest total whose intervention share is
# closest to `treatment_share`; for example, "3:5" for 0.375.
ratio_label <- function(treatment_share, max_total = 40) {
  totals <- 2:max_total
  intervention <- round(treatment_share * totals)
  valid <- intervention > 0 & intervention < totals
  error <- abs(intervention / totals - treatment_share)
  best <- which(valid)[which.min(error[valid])]
  paste(intervention[best], totals[best] - intervention[best], sep = ":")
}

# Expected intervention:control recruitment, such as "10.0:30.0".
format_recruitment <- function(recruitment, treatment_share) {
  paste0(
    format_fixed(recruitment * treatment_share, 1), ":",
    format_fixed(recruitment * (1 - treatment_share), 1)
  )
}

# Table 1 ---------------------------------------------------------------

# One site that can recruit 40 participants but treat only 10, under ratios
# 1:1 to 1:5.
planning_table_1 <- function() {
  p <- 1 / (2:6)
  recruitment <- site_recruitment(40, 10, p)
  information <- site_information(40, 10, p)
  data.frame(
    Ratio = paste0("1:", 1:5),
    `Intervention proportion` = format_fixed(p, 3),
    `Expected recruitment (Intervention:control)` =
      format_recruitment(recruitment, p),
    `Total recruitment/site` = format_fixed(recruitment, 1),
    `Required total n` = format_count(information_target / (4 * p * (1 - p))),
    `Benchmark information/site` = format_fixed(information, 1),
    `Sites required` = sites_required(information),
    check.names = FALSE
  )
}

# Table 2 ---------------------------------------------------------------

# A constrained site paired with an unconstrained site, for three pairs of
# site capacities.
planning_table_2 <- function() {
  rbind(
    table_2_scenario(constrained_capacity = 40, unconstrained_capacity = 40),
    table_2_scenario(constrained_capacity = 20, unconstrained_capacity = 60),
    table_2_scenario(constrained_capacity = 60, unconstrained_capacity = 20)
  )
}

# Allocation strategies for one pair of sites. The constrained site can
# treat one quarter of its capacity, so 1:3 allocation fills it. Each vector
# below holds the constrained site, then the unconstrained site.
table_2_scenario <- function(constrained_capacity, unconstrained_capacity) {
  capacity <- c(constrained_capacity, unconstrained_capacity)
  intervention_limit <- c(constrained_capacity / 4, unconstrained_capacity)
  filling_share <- capacity_filling_share(capacity[1], intervention_limit[1])
  filling_ratio <- ratio_label(filling_share)

  # The uniform compromise is the overall intervention proportion of the
  # site-specific strategy.
  specific_shares <- c(filling_share, 0.5)
  compromise_share <- weighted.mean(
    specific_shares,
    site_recruitment(capacity, intervention_limit, specific_shares)
  )

  rows <- rbind(
    table_2_row("Uniform 1:1", c(0.5, 0.5), capacity, intervention_limit),
    table_2_row(
      paste("Uniform", filling_ratio), rep(filling_share, 2),
      capacity, intervention_limit
    ),
    table_2_row(
      paste("Uniform compromise", ratio_label(compromise_share)),
      rep(compromise_share, 2), capacity, intervention_limit
    ),
    table_2_row(
      paste("Site-specific 1:1 or", filling_ratio), specific_shares,
      capacity, intervention_limit
    )
  )
  if (constrained_capacity == unconstrained_capacity) {
    rows <- rbind(
      rows,
      table_2_row(
        "Mixed 1:3 and 3:1 allocation (overall 1:1)", c(0.25, 0.75),
        capacity, intervention_limit
      )
    )
  }
  rows
}

# One strategy for one pair of sites: recruitment and information averaged
# over the two sites, and the participants and sites needed to reach the
# information target.
table_2_row <- function(strategy, treatment_shares, capacity,
                        intervention_limit) {
  recruitment <- site_recruitment(
    capacity, intervention_limit, treatment_shares
  )
  information <- site_information(
    capacity, intervention_limit, treatment_shares
  )
  data.frame(
    `Max capacity constrained/unconstrained` = paste(capacity, collapse = "/"),
    Strategy = strategy,
    Constrained = format_recruitment(recruitment[1], treatment_shares[1]),
    Unconstrained = format_recruitment(recruitment[2], treatment_shares[2]),
    n = format_count(
      information_target * mean(recruitment) / mean(information)
    ),
    `Capacity per site` = format_fixed(mean(recruitment), 1),
    `Information per site` = format_fixed(mean(information), 2),
    `Sites required` = sites_required(mean(information)),
    check.names = FALSE
  )
}

# Table S1 --------------------------------------------------------------

# Site-specific allocation compared with one compromise ratio for all sites.
planning_table_s1 <- function() {
  rows <- list()
  for (intervention_limit in c(5, 8, 10, 12, 15)) {
    for (constrained_site_share in c(0.25, 0.50, 0.75)) {
      rows[[length(rows) + 1L]] <- tradeoff_row(
        intervention_limit, constrained_site_share
      )
    }
  }
  do.call(rbind, rows)
}

# Constrained sites use their capacity-filling ratio and other sites use 1:1;
# the compromise is the single ratio with the same overall intervention
# proportion. The percentages compare site-specific allocation with the
# compromise. Each vector below holds the constrained sites, then the others.
tradeoff_row <- function(intervention_limit, constrained_site_share,
                         full_capacity = 40) {
  site_weights <- c(constrained_site_share, 1 - constrained_site_share)
  limits <- c(intervention_limit, full_capacity)
  filling_share <- capacity_filling_share(full_capacity, intervention_limit)
  specific_shares <- c(filling_share, 0.5)

  specific_recruitment <- site_recruitment(
    full_capacity, limits, specific_shares
  )
  participant_weights <- site_weights * specific_recruitment
  compromise_share <- weighted.mean(specific_shares, participant_weights)
  compromise_recruitment <- site_recruitment(
    full_capacity, limits, compromise_share
  )

  # Information per participant when every site uses the compromise, and when
  # each site keeps its own share.
  information_per_participant_compromise <-
    4 * compromise_share * (1 - compromise_share)
  information_per_participant_specific <- weighted.mean(
    4 * specific_shares * (1 - specific_shares), participant_weights
  )

  # Site-specific allocation needs more participants but recruits more per
  # site. Sites required change by the ratio of these two quantities.
  participant_ratio <- information_per_participant_compromise /
    information_per_participant_specific
  recruitment_ratio <- sum(site_weights * specific_recruitment) /
    sum(site_weights * compromise_recruitment)
  site_ratio <- participant_ratio / recruitment_ratio

  data.frame(
    `Constrained site intervention limit` = paste0(
      format_fixed(100 * intervention_limit / full_capacity, 1), "%"
    ),
    `Tailored constrained ratio` = ratio_label(filling_share, full_capacity),
    `Constrained sites` = paste0(
      format_fixed(100 * constrained_site_share, 0), "%"
    ),
    `Average compromise p` = format_fixed(compromise_share, 3),
    `Participants required` = format_signed_percent(participant_ratio),
    `Recruitment per site` = format_signed_percent(recruitment_ratio),
    `Number of sites` = format_signed_percent(site_ratio),
    check.names = FALSE
  )
}

# A ratio of 1.071 prints as "+7.1%".
format_signed_percent <- function(ratio) {
  percent <- 100 * (ratio - 1)
  paste0(if (percent > 0) "+", format_fixed(percent, 1), "%")
}
