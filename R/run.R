# Command-line options, parallel workers, and running and saving the
# simulation.

# Parse --name=value options. Every option must be named in `defaults`.
parse_options <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  option_names <- sub("=.*$", "", sub("^--", "", args))
  if (any(!grepl("^--[^=]+=.+$", args)) || anyDuplicated(option_names) ||
      any(!option_names %in% names(defaults))) {
    stop(
      "Use each option once as --name=value. Options: ",
      paste(names(defaults), collapse = ", ")
    )
  }
  defaults[option_names] <- sub("^--[^=]+=", "", args)
  defaults
}

# Integers and ascending ranges, such as "82,88,94" or "1-99".
integer_sequence <- function(text, maximum) {
  if (!grepl("^[0-9]+(-[0-9]+)?(,[0-9]+(-[0-9]+)?)*$", text)) {
    stop("Expected integers or ascending ranges, for example 82,88,94 or 1-99.")
  }
  values <- integer()
  for (part in strsplit(text, ",", fixed = TRUE)[[1]]) {
    # A single number is a range that starts and ends at that number.
    limits <- as.numeric(strsplit(part, "-", fixed = TRUE)[[1]])
    first <- limits[1]
    last <- limits[length(limits)]
    if (first < 1 || last > maximum || first > last) {
      stop("Invalid range: ", part, ". Use values from 1 to ", maximum, ".")
    }
    values <- c(values, seq.int(first, last))
  }
  if (anyDuplicated(values)) {
    stop("Ranges must not overlap.")
  }
  as.integer(values)
}

# One whole number, such as a replicate or worker count.
parse_count <- function(text, maximum) {
  value <- integer_sequence(text, maximum)
  if (length(value) != 1L) {
    stop("Expected one whole number from 1 to ", maximum, ".")
  }
  value
}

# Create an output folder. It must be outside data/, which holds the paper's
# results.
output_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path <- normalizePath(path, mustWork = TRUE)
  data_folder <- normalizePath("data", mustWork = TRUE)
  if (path == data_folder || startsWith(path, paste0(data_folder, "/"))) {
    stop("Write new results outside data/, which holds the paper's results.")
  }
  path
}

# Apply `worker` to each task, in parallel when cores > 1. Each worker process
# starts in the release folder and loads the same code. `...` passes further
# arguments to `worker`. `exports` names objects that a script defines and its
# workers also need. Tasks are handed out one at a time because their running
# times differ greatly.
parallel_tasks <- function(tasks, worker, cores, ..., exports = character()) {
  workers <- min(cores, length(tasks))
  if (workers <= 1L) {
    return(lapply(tasks, worker, ...))
  }
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, setwd, getwd())
  parallel::clusterEvalQ(cluster, source("R/load.R"))
  parallel::clusterExport(cluster, exports, envir = globalenv())
  parallel::parLapplyLB(cluster, tasks, worker, ..., chunk.size = 1)
}

# The file that holds one scenario's results.
task_results_path <- function(output, task_id, n_replicates) {
  file.path(output, "tasks", sprintf(
    "task_%03d_%d_replicates.rds", task_id, n_replicates
  ))
}

# Run replicates 1 to `n_replicates` of one scenario and save the results of
# every method, one row per replicate and method. The file is written under a
# temporary name and then renamed, so that an interrupted run never leaves a
# partial file.
run_scenario <- function(scenario, n_replicates, output) {
  records <- lapply(seq_len(n_replicates), function(replicate_index) {
    run_replicate(scenario, replicate_index)
  })
  path <- task_results_path(output, scenario$task_id, n_replicates)
  temporary <- tempfile(tmpdir = dirname(path))
  saveRDS(dplyr::bind_rows(records), temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not save task ", scenario$task_id)
  }
  invisible(path)
}

# Summarize the saved results of the given tasks, one row per task and
# method, and write them to summary.csv.
summarize_run <- function(output, scenarios, tasks, n_replicates) {
  task_summaries <- lapply(tasks, function(task) {
    results <- readRDS(task_results_path(output, task, n_replicates))
    summarize_task(results, get_scenario(scenarios, task))
  })
  run_summary <- dplyr::bind_rows(task_summaries)
  readr::write_csv(run_summary, file.path(output, "summary.csv"))
  invisible(run_summary)
}
