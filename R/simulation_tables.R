# Simulation result tables: main Tables 3 and 4 and Supplementary Tables
# S2-S38. Each table shows one scenario, with the analysis methods as rows and
# the outcomes as columns. data/table_tasks.csv lists the simulation tasks in
# each table. This file also writes the tables to files, and
# R/diagnostic_tables.R reuses the footnote helpers below.

outcome_order <- c("Time-to-Event", "Binary", "Continuous")
outcome_labels <- c(
  `Time-to-Event` = "Time-to-event", Binary = "Binary",
  Continuous = "Continuous"
)
effect_symbols <- c(`Time-to-Event` = "HR", Binary = "OR", Continuous = "MD")

# The rows of the simulation tables, in the paper's order.
table_rows <- c(
  "No adjustment",
  "No adjustment + cluster-robust SE",
  "IPTW (robust SE)",
  "IPTW (cluster-robust SE)",
  "Site covariate",
  "Allocation ratio covariate",
  "Stratified by site",
  "Stratified by allocation ratio",
  "Site random intercept",
  "Allocation ratio covariate + site random intercept"
)

# The table row of a stored method name. The paper prints "No adjustment" for
# "Unadjusted" and shows the pooled continuous analyses as "Stratified".
method_label <- function(method) {
  label <- sub(" (Pooled)", "", method, fixed = TRUE)
  label[label == "Unadjusted"] <- "No adjustment"
  label[label == "Unadjusted (cluster-robust SE)"] <-
    "No adjustment + cluster-robust SE"
  label
}

# Number formatting ----------------------------------------------------

format_fixed <- function(value, digits) {
  formatC(value, format = "f", digits = digits)
}

format_count <- function(value) {
  format(round(value), big.mark = ",", scientific = FALSE, trim = TRUE)
}

# Monte Carlo SEs too small to show at this precision print as, for example,
# "<0.001".
format_mcse <- function(value, digits) {
  if (round(value, digits) == 0) {
    return(paste0("<", format_fixed(10^(-digits), digits)))
  }
  format_fixed(value, digits)
}

# Tables ---------------------------------------------------------------

# The paper's results from the full simulation and the two sets of checks.
read_saved_results <- function(file_name) {
  as.data.frame(readr::read_csv(
    file.path("data/saved_results", file_name), show_col_types = FALSE
  ))
}

# Write every table of the paper, with its footnote, to `output`. Tables 1,
# 2, and S1 are recalculated from their formulas. The others use the paper's
# saved results.
export_paper_tables <- function(output = "output/tables") {
  summary <- read_saved_results("simulation_summary.csv")
  tables <- c(
    list(
      list(table_id = "1", rows = planning_table_1()),
      list(table_id = "2", rows = planning_table_2()),
      list(table_id = "S1", rows = planning_table_s1())
    ),
    simulation_tables(summary),
    list(
      method_interpretation_table(),
      interval_table(
        read_saved_results("interval_checks.csv"),
        read_saved_results("interval_checks_cr2.csv")
      )
    )
  )
  write_tables(tables, output)
}

# Write the tables that a new simulation run covers.
export_run_tables <- function(summary_file, output) {
  summary <- as.data.frame(
    readr::read_csv(summary_file, show_col_types = FALSE)
  )
  write_tables(simulation_tables(summary), output)
}

# The tables built from a simulation summary. Tables 3, 4, and S2-S38 show
# estimates and coverage, and Table S40 shows SDs and SEs. A table is left out
# when the summary lacks one of its scenarios.
simulation_tables <- function(summary) {
  table_map <- read_table_map()
  tables <- list()
  for (i in seq_len(nrow(table_map))) {
    if (table_map$source[i] != "Simulation") {
      next
    }
    table_id <- table_map$table[i]
    task_ids <- table_tasks(table_map, table_id)
    if (!all(task_ids %in% summary$task_id)) {
      next
    }
    rows <- summary[summary$task_id %in% task_ids, ]
    tables[[length(tables) + 1L]] <- build_result_table(rows, table_id)
  }
  sd_and_se <- sd_and_se_table(summary, table_map)
  if (!is.null(sd_and_se)) {
    tables[[length(tables) + 1L]] <- sd_and_se
  }
  tables
}

# One CSV file per table and one text file per footnote.
write_tables <- function(tables, output) {
  output <- output_directory(output)
  for (table in tables) {
    file_stem <- file.path(output, paste0("Table_", table$table_id))
    readr::write_csv(table$rows, paste0(file_stem, ".csv"))
    if (!is.null(table$note)) {
      writeLines(table$note, paste0(file_stem, "_note.txt"))
    }
  }
  message("Wrote ", length(tables), " tables to ", output)
}

# The title, simulation tasks, and source of every table.
read_table_map <- function(path = "data/table_tasks.csv") {
  as.data.frame(readr::read_csv(path, col_types = "cccc"))
}

# The simulation tasks in one table, such as 82, 88, and 94 for Table 3.
table_tasks <- function(table_map, table_id) {
  task_ids <- table_map$task_ids[table_map$table == table_id]
  if (is.na(task_ids)) {
    return(integer(0))
  }
  as.integer(strsplit(task_ids, ",", fixed = TRUE)[[1]])
}

# Tables 3, 4, and S2-S38. `rows` holds the summary rows of the table's tasks.
build_result_table <- function(rows, table_id) {
  effect_scenario <- rows$effect_scenario[1]
  display <- result_cells(rows)
  names(display) <- c(
    "Analysis method", outcome_column_headers(effect_scenario)
  )

  # Main Tables 3 and 4 share one scenario, described in this fixed sentence,
  # and their titles state the effect. The supplementary tables describe their
  # scenario and effect from the summary.
  if (table_id %in% c("3", "4")) {
    scenario_note <- paste(
      "The scenario has 30 sites of variable size; 50% use 1:3 allocation and",
      "are larger, with strongly and heterogeneously worse prognosis.",
      "Time-to-event and binary outcomes are sparse, and continuous residual",
      "variance is higher in 1:3 sites."
    )
    other_tables_note <- paste(
      "Supplementary Appendix Tables S2-S38 summarize results for the",
      "additional scenarios."
    )
  } else {
    scenario_note <- c(
      scenario_design_note(rows), true_effect_note(effect_scenario, rows)
    )
    other_tables_note <- NULL
  }

  note <- paste(c(
    paste(
      "Cells show the exponentiated mean log HR or log OR, or arithmetic mean",
      "MD, followed by 95% confidence interval coverage of the true",
      "within-site effect."
    ),
    scenario_note,
    "HR = hazard ratio. OR = odds ratio. MD = mean difference.",
    "IPTW = inverse probability of treatment weighting. SE = standard error.",
    cell_label_notes(display),
    other_tables_note,
    replicate_count_note(rows),
    finite_result_note(rows),
    fit_status_note(rows),
    estimate_mcse_note(rows),
    coverage_mcse_note(rows)
  ), collapse = " ")
  list(table_id = table_id, rows = display, note = note)
}

# One row per method and one column per outcome.
result_cells <- function(rows) {
  data.frame(
    table_rows,
    estimate_coverage_cells(rows, "Time-to-Event"),
    estimate_coverage_cells(rows, "Binary"),
    estimate_coverage_cells(rows, "Continuous")
  )
}

# The cells of one outcome. A cell shows the displayed mean estimate and
# coverage, "Not estimated", or "Not simulated" when the scenario has no such
# outcome.
estimate_coverage_cells <- function(rows, outcome) {
  if (!any(rows$outcome == outcome)) {
    return(rep("Not simulated", length(table_rows)))
  }
  in_order <- rows_in_table_order(rows, outcome)
  digits <- if (outcome == "Continuous") 2 else 3
  ifelse(
    is.finite(in_order$mean_estimate_display) & is.finite(in_order$coverage),
    paste0(
      format_fixed(in_order$mean_estimate_display, digits), " / ",
      format_fixed(100 * in_order$coverage, 1), "%"
    ),
    "Not estimated"
  )
}

# The summary rows of one outcome, one per table row and in that order.
rows_in_table_order <- function(rows, outcome) {
  outcome_rows <- rows[rows$outcome == outcome, ]
  outcome_rows[match(table_rows, method_label(outcome_rows$method)), ]
}

# The outcome above its true effect, such as "Binary\nOR = 1.0".
outcome_column_headers <- function(effect_scenario) {
  effects <- true_effect_values(effect_scenario)
  c(
    paste0("Time-to-event\nHR = ", format_fixed(effects$hr, 1)),
    paste0("Binary\nOR = ", format_fixed(effects$or, 1)),
    paste0("Continuous\nMD = ", format_fixed(effects$md, 1))
  )
}

true_effect_values <- function(effect_scenario) {
  list(
    hr = exp(treatment_effects$`Time-to-Event`[[effect_scenario]]),
    or = exp(treatment_effects$Binary[[effect_scenario]]),
    md = treatment_effects$Continuous[[effect_scenario]]
  )
}

# Footnote sentences ---------------------------------------------------

scenario_design_note <- function(rows) {
  first <- rows[1, ]
  site_sizes <- if (first$site_size_pattern == "Equal") {
    "sites of equal size"
  } else {
    paste0(
      "sites of variable size generated with mean ",
      format_fixed(first$mean_site_size, 0), " patients/site, median ",
      format_fixed(first$median_site_size, 0), ", and range ",
      format_fixed(first$min_site_size, 0), "-",
      format_fixed(first$max_site_size, 0), " (coefficient of variation ",
      format_fixed(first$site_size_cv, 2), ")"
    )
  }
  # The site-size phrase already describes equal-size sites.
  size_ratio <- switch(first$size_allocation_association,
    Ratio13Larger = "1:3 sites larger",
    Ratio13Smaller = "1:3 sites smaller"
  )
  prognosis <- switch(first$prognosis_strength,
    Independent = "site prognosis independent of allocation ratio",
    Strong = "strong heterogeneous worse prognosis in 1:3 sites"
  )
  frequency <- switch(first_applicable(rows$outcome_frequency),
    Common = "common time-to-event and binary outcome frequencies",
    Sparse = "sparse time-to-event and binary outcome frequencies"
  )
  variance <- switch(first_applicable(rows$residual_variance),
    Equal = "equal continuous residual variance",
    Differential = "higher continuous residual variance in 1:3 sites"
  )
  pieces <- c(
    paste(first$trial_size_sites, "sites"),
    site_sizes,
    paste0(
      format_fixed(100 * first$shifted_site_fraction, 0),
      "% of sites use 1:3 allocation"
    ),
    size_ratio, prognosis, frequency, variance,
    outcome_design_phrase(rows)
  )
  paste0("Scenario design: ", paste(pieces, collapse = "; "), ".")
}

# The first value other than "Not_Applicable", if any.
first_applicable <- function(values) {
  values <- values[values != "Not_Applicable"]
  if (length(values) == 0L) "Not_Applicable" else values[[1]]
}

outcome_design_phrase <- function(rows) {
  if (all(rows$outcome == "Continuous")) {
    return("continuous outcome only")
  }
  if (!any(rows$outcome == "Continuous")) {
    return("time-to-event and binary outcomes only")
  }
  "time-to-event, binary, and continuous outcomes"
}

# The mean difference is also given in units of the 1:1-site residual SD (20).
true_effect_note <- function(effect_scenario, rows) {
  effects <- true_effect_values(effect_scenario)
  parts <- character()
  if ("Time-to-Event" %in% rows$outcome) {
    parts <- c(parts, paste0("HR = ", format_fixed(effects$hr, 1)))
  }
  if ("Binary" %in% rows$outcome) {
    parts <- c(parts, paste0("OR = ", format_fixed(effects$or, 1)))
  }
  if ("Continuous" %in% rows$outcome) {
    parts <- c(parts, paste0(
      "MD = ", format_fixed(effects$md, 1),
      " (", format_fixed(abs(effects$md) / 20, 1), " residual SD)"
    ))
  }
  paste0("True treatment effect: ", paste(parts, collapse = "; "), ".")
}

# Define the "Not estimated" and "Not simulated" labels when a table uses them.
cell_label_notes <- function(display) {
  notes <- character()
  if (any(display == "Not estimated")) {
    notes <- c(notes, paste(
      "Not estimated indicates that the fitted method produced no",
      "finite estimate for that outcome."
    ))
  }
  if (any(display == "Not simulated")) {
    notes <- c(notes, paste(
      "Not simulated indicates that the outcome type was outside",
      "this scenario design."
    ))
  }
  notes
}

replicate_count_note <- function(rows) {
  paste("The simulation used", format_count(rows$n_attempted[1]), "replicates.")
}

# Which methods returned fewer finite estimates and SEs than attempted.
finite_result_note <- function(rows) {
  attempted <- format_count(rows$n_attempted[1])
  some <- rows[rows$n_successful > 0 & rows$n_successful < rows$n_attempted, ]
  none <- rows[rows$n_successful == 0, ]
  if (nrow(some) == 0L && nrow(none) == 0L) {
    return(paste(
      "All", attempted,
      "replicates returned finite estimates and SEs for each method."
    ))
  }

  if (nrow(some) == 0L) {
    note <- paste(
      "All", attempted,
      "replicates returned finite estimates and SEs for each estimable cell."
    )
  } else {
    some_counts <- paste0(
      outcome_method_labels(some), " had ", format_count(some$n_successful),
      " of ", attempted, " replicates with finite estimates and SEs"
    )
    note <- c(
      paste0(paste(some_counts, collapse = "; "), "."),
      paste(
        "All other estimable cells had", attempted, "of", attempted,
        "replicates with finite estimates and SEs."
      )
    )
  }
  if (nrow(none) > 0L) {
    none_labels <- paste0(
      outcome_method_labels(none),
      " produced no finite estimates in any of ", attempted, " replicates"
    )
    note <- c(note, paste0(paste(none_labels, collapse = "; "), "."))
  }
  note
}

# Counts of fits with a warning or a nonconvergence status, by method. See
# detect_nonconvergence() for what each status covers.
fit_status_note <- function(rows) {
  affected <- rows[rows$n_warning > 0 | rows$n_nonconvergence > 0, ]
  if (nrow(affected) == 0L) {
    return("No warning, convergence, or singularity flag was recorded.")
  }

  note <- character()
  if (any(affected$n_warning > 0)) {
    note <- c(note, paste0(
      "Fits with a warning: ",
      method_count_list(affected, affected$n_warning), "."
    ))
  }
  if (any(affected$n_nonconvergence > 0)) {
    note <- c(
      note,
      paste0(
        "Fits with convergence or singularity flags: ",
        method_count_list(affected, affected$n_nonconvergence), "."
      ),
      paste(
        "The flag counts combine convergence problems and singular fits,",
        "whose estimated site variance is at or near zero."
      )
    )
  }
  note <- c(note, paste(
    "Counts of finite results include these fits when the treatment",
    "estimate and SE were finite."
  ))

  site_cox <- affected[
    affected$outcome == "Time-to-Event" &
      affected$method == "Site covariate" & affected$n_warning > 0,
  ]
  if (nrow(site_cox) > 0L &&
      all(concerns_site_terms_only(site_cox$modal_status_message))) {
    note <- c(note, paste(
      "The most frequent warning for Time-to-event / Site covariate",
      "concerned possibly infinite coefficients for site indicators, not the",
      "treatment coefficient. This can occur when sites have few or no",
      "events."
    ))
  }
  note
}

# The methods with a positive count, such as
# "Binary / Site random intercept 30; Continuous / Site random intercept 2".
method_count_list <- function(rows, counts) {
  keep <- counts > 0
  paste0(
    outcome_method_labels(rows[keep, ]), " ",
    format(counts[keep], big.mark = ",", trim = TRUE),
    collapse = "; "
  )
}

# coxph warns, for example, "Loglik converged before variable 3,5 ;
# coefficient may be infinite." It lists variables in increasing order. In the
# site-covariate model, variable 1 is treatment and the others are sites.
concerns_site_terms_only <- function(messages) {
  grepl("Loglik converged before variable", messages) &
    !grepl("before variable +1[ ,]", messages)
}

# Labels such as "Binary / Stratified by allocation ratio", named as in the
# table rows.
outcome_method_labels <- function(rows) {
  paste0(outcome_labels[rows$outcome], " / ", method_label(rows$method))
}

# The Monte Carlo SE of a mean estimate is its empirical SD divided by the
# square root of the number of finite replicates. For exponentiated HRs and
# ORs, the delta method multiplies this by the displayed estimate.
estimate_mcse_note <- function(rows) {
  estimable <- rows[rows$n_successful > 0 & is.finite(rows$empirical_se), ]
  mcse <- estimable$empirical_se / sqrt(estimable$n_successful)
  ratio_scale <- estimable$outcome != "Continuous"
  mcse[ratio_scale] <- mcse[ratio_scale] *
    estimable$mean_estimate_display[ratio_scale]

  labels <- character()
  for (outcome in intersect(outcome_order, estimable$outcome)) {
    largest <- max(mcse[estimable$outcome == outcome])
    digits <- if (outcome == "Continuous") 2 else 3
    labels <- c(labels, paste0(
      format_mcse(largest, digits), " for ", effect_symbols[[outcome]]
    ))
  }
  paste0(
    "The largest Monte Carlo standard errors for the displayed mean ",
    "treatment effect estimates were ", paste(labels, collapse = ", "), "."
  )
}

# The Monte Carlo SE of a coverage proportion is largest at 50%.
coverage_mcse_note <- function(rows) {
  fewest <- min(rows$n_successful[rows$n_successful > 0])
  max_mcse <- 100 * sqrt(0.25 / fewest)
  near_95_mcse <- 100 * sqrt(0.95 * 0.05 / fewest)
  paste(
    "For coverage percentages, the maximum Monte Carlo standard error was",
    paste0(format_fixed(max_mcse, 1), " percentage points,"),
    "and near 95% coverage it was",
    paste0(format_fixed(near_95_mcse, 2), " percentage points.")
  )
}
