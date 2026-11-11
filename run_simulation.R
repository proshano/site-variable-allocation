# Run the simulation. For each scenario this simulates trials, analyzes every
# trial with the ten methods, and summarizes the estimates, SEs, and coverage.
#
# Run it with: Rscript run_simulation.R
#
# The default is 100 replicates of the six scenarios in Tables 3 and 4, which
# takes 1 to 2 minutes. The paper used 10,000 replicates of all 99 scenarios,
# which takes about 230 core-hours:
#
#   Rscript run_simulation.R --tasks=1-99 --replicates=10000 --cores=8 \
#     --output=output/full_simulation
#
# Tasks are rows of data/scenarios.csv, and data/table_tasks.csv lists the
# tasks behind every table. Each scenario is saved when it finishes. If a run
# stops, repeat the same command and the finished scenarios are skipped.
source("R/load.R")

settings <- parse_options(list(
  tasks = "82,83,88,89,94,95", replicates = "100", cores = "1",
  output = "output/simulation"
))
tasks <- integer_sequence(settings$tasks, 99L)
n_replicates <- parse_count(settings$replicates, 100000L)
cores <- parse_count(settings$cores, 128L)
output <- output_directory(settings$output)
dir.create(file.path(output, "tasks"), showWarnings = FALSE)

# Run the scenarios that are not yet saved in the output folder.
scenarios <- read_scenarios()
is_saved <- file.exists(task_results_path(output, tasks, n_replicates))
if (any(is_saved)) {
  message("Skipping ", sum(is_saved), " scenarios already saved in ", output)
}
scenarios_to_run <- lapply(tasks[!is_saved], function(task) {
  get_scenario(scenarios, task)
})
invisible(parallel_tasks(
  scenarios_to_run, run_scenario, cores,
  n_replicates = n_replicates, output = output
))

# One row per scenario and method, then the paper's tables for the scenarios
# in this run.
summarize_run(output, scenarios, tasks, n_replicates)
message("Summary: ", file.path(output, "summary.csv"))
export_run_tables(
  file.path(output, "summary.csv"), file.path(output, "tables")
)
