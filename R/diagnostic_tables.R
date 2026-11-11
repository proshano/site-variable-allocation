# Supplementary Tables S39-S41. S39 states what each method estimates, S40
# compares empirical SDs with estimated SEs, and S41 compares interval methods
# for IPTW with cluster-robust SEs. S40 uses the simulation summary and S41
# the saved results of run_interval_checks.R.

# Table S39: what each method estimates.
method_interpretation_table <- function() {
  rows <- data.frame(
    Method = c(
      "No adjustment, with either SE",
      "IPTW, with either SE",
      "Site covariate or site stratification",
      "Allocation ratio covariate or stratification",
      "Continuous pooling within strata",
      "Site random intercept",
      "Allocation ratio covariate + site random intercept"
    ),
    Interpretation = c(
      paste(
        "The difference between all intervention and all control",
        "participants, with sites pooled. Mean differences agree with the",
        "within-site effect when prognosis is unrelated to allocation ratio.",
        "Hazard and odds ratios then agree only for a null effect, because",
        "pooling sites that differ in prognosis moves a non-null ratio toward",
        "1. A cluster-robust SE changes the interval, not the estimate."
      ),
      paste(
        "The effect averaged over all participants (a marginal effect), with",
        "sites balanced between the arms by weighting. It agrees with the",
        "within-site effect for mean differences and for a null effect. A",
        "non-null marginal hazard or odds ratio is generally closer to 1, and",
        "a marginal hazard ratio can change over follow-up."
      ),
      paste(
        "The effect within site, which is the target of all comparisons here.",
        "Stratification allows each site its own baseline risk or hazard. A",
        "site covariate in a Cox model assumes proportional site hazards, and",
        "estimating a parameter for every site moved hazard and odds ratios",
        "up to 3% further from 1."
      ),
      paste(
        "The effect within each ratio group, averaged over the sites in that",
        "group. It agrees with the within-site effect for mean differences",
        "and for a null effect. A non-null hazard or odds ratio is generally",
        "closer to 1."
      ),
      paste(
        "The mean difference within sites, or within ratio groups, combined",
        "by inverse-variance weighting. It agrees with the within-site",
        "effect. The SE treats each estimated variance as known."
      ),
      paste(
        "The effect within site, provided site prognosis is unrelated to",
        "allocation ratio. Otherwise part of the confounding between sites",
        "stays in the estimate."
      ),
      paste(
        "The effect within site, after allowing the ratio groups to differ in",
        "average prognosis. It agreed with the within-site effect in every",
        "scenario. In all of them, prognosis also varied between sites within",
        "each ratio group."
      )
    ),
    check.names = FALSE
  )
  note <- paste(
    "Every comparison in this paper is with the true within-site effect, not",
    "with each method's own target. IPTW = inverse probability of treatment",
    "weighting; SE = standard error."
  )
  list(table_id = "S39", rows = rows, note = note)
}

# Table S40 has one block of rows for each of these tables of the paper. Each
# block starts with a row that names the table.
sd_and_se_blocks <- c(
  `3` = "Table 3 (30 sites, null effect)",
  `4` = "Table 4 (30 sites, moderate effect)",
  S5 = "Table S5 (120 sites, null effect)",
  S9 = "Table S9 (120 sites, large effect)",
  S17 = "Table S17 (30 sites, null effect)"
)

# Table S40: empirical SD and mean estimated SE of each method. A block is
# left out when `summary` lacks one of its scenarios, so a short run of
# run_simulation.R shows the blocks it covers.
sd_and_se_table <- function(summary, table_map) {
  display <- NULL
  for (source_table in names(sd_and_se_blocks)) {
    task_ids <- table_tasks(table_map, source_table)
    if (!all(task_ids %in% summary$task_id)) {
      next
    }
    rows <- summary[summary$task_id %in% task_ids, ]
    heading <- data.frame(
      Method = sd_and_se_blocks[[source_table]],
      `Log HR` = "", `Log OR` = "", `Mean difference` = "",
      check.names = FALSE
    )
    block <- data.frame(
      Method = table_rows,
      `Log HR` = sd_and_se_cells(rows, "Time-to-Event"),
      `Log OR` = sd_and_se_cells(rows, "Binary"),
      `Mean difference` = sd_and_se_cells(rows, "Continuous"),
      check.names = FALSE
    )
    display <- rbind(display, heading, block)
  }
  if (is.null(display)) {
    return(NULL)
  }

  note <- paste(c(
    paste(
      "Cells show the empirical SD of the estimates / the mean estimated SE,",
      "on the coefficient scale (log HR, log OR, or mean difference)."
    ),
    paste(
      "Each block repeats the scenario of the table named in its first row.",
      "That table gives the design, the number of replicates with finite",
      "estimates, and coverage."
    ),
    if (any(display == "Not estimated")) {
      paste(
        "Not estimated indicates that the method produced no finite estimate",
        "for that outcome."
      )
    },
    paste(
      "HR = hazard ratio; IPTW = inverse probability of treatment weighting;",
      "OR = odds ratio; SD = standard deviation; SE = standard error."
    )
  ), collapse = " ")
  list(table_id = "S40", rows = display, note = note)
}

# The cells of one outcome: empirical SD / mean estimated SE.
sd_and_se_cells <- function(rows, outcome) {
  in_order <- rows_in_table_order(rows, outcome)
  ifelse(
    is.finite(in_order$empirical_se),
    sprintf("%.3f / %.3f", in_order$empirical_se, in_order$model_based_se),
    "Not estimated"
  )
}

# Table S41: normal, t(29), and CR2/Satterthwaite intervals for IPTW with
# cluster-robust SEs in the six main-table settings. CR2 was not applied to
# Cox models.
interval_table <- function(intervals, cr2) {
  stopifnot(
    nrow(intervals) == 6L, all(intervals$replicates == 10000L),
    all(intervals$sites == 30L),
    all(intervals$max_abs_summary_difference < 1e-8)
  )
  stopifnot(
    nrow(cr2) == 4L, !anyDuplicated(cr2$task_id),
    setequal(cr2$task_id, c(88L, 89L, 94L, 95L)),
    all(cr2$replicates == 10000L), all(cr2$sites == 30L),
    all(cr2$max_abs_summary_difference < 1e-8)
  )
  # CR2 results exist for the binary and continuous settings only.
  matched <- cr2[match(intervals$task_id, cr2$task_id), ]
  has_cr2 <- !is.na(matched$task_id)
  cr2_mean_se <- ifelse(
    has_cr2, sprintf("%.3f", matched$cr2_mean_se), "Not applied"
  )
  cr2_coverage <- ifelse(
    has_cr2, sprintf("%.1f%%", 100 * matched$cr2_coverage), "Not applied"
  )
  cr2_df <- ifelse(
    has_cr2, sprintf("%.1f", matched$cr2_df_mean), "Not applied"
  )

  rows <- data.frame(
    Outcome = ifelse(
      intervals$outcome == "Time-to-Event", "Time-to-event", intervals$outcome
    ),
    Effect = ifelse(intervals$effect == "Realistic", "Moderate", "Null"),
    `Empirical\nSD` = sprintf("%.3f", intervals$empirical_sd),
    `Mean SE\noriginal` = sprintf("%.3f", intervals$mean_se),
    `Mean SE\nCR2` = cr2_mean_se,
    `Coverage\nnormal` = sprintf("%.1f%%", 100 * intervals$normal_coverage),
    `Coverage\nt(29)` = sprintf("%.1f%%", 100 * intervals$t_coverage),
    `Coverage\nCR2` = cr2_coverage,
    `Mean df\nCR2` = cr2_df,
    check.names = FALSE
  )
  mcse <- 100 * c(intervals$normal_mcse, intervals$t_mcse, cr2$cr2_mcse)
  note <- paste(
    "Each row summarizes 10,000 replicates of a scenario with 30 sites.",
    "Normal is the original interval, which used the cluster-robust SE and",
    "the normal quantile 1.960.",
    "t(29) used the same SE and the t quantile on 29 degrees of freedom,",
    "2.045.",
    "CR2 used the bias-reduced cluster-robust SE and a t quantile on",
    "Satterthwaite degrees of freedom, and was not applied to Cox models.",
    "Coverage is of the true within-site effect. SDs and SEs are on the",
    "coefficient scale (log HR, log OR, or mean difference).",
    sprintf(
      "Monte Carlo SEs of coverage were %.2f to %.2f percentage points.",
      min(mcse), max(mcse)
    ),
    "CR2 = bias-reduced cluster-robust variance estimator; df = degrees of",
    "freedom; HR = hazard ratio; IPTW = inverse probability of treatment",
    "weighting; OR = odds ratio; SD = standard deviation; SE = standard",
    "error."
  )
  list(table_id = "S41", rows = rows, note = note)
}
