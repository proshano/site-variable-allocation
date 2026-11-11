# Site-variable allocation ratios in multicentre randomized trials

Code and saved results for the paper **Tailoring allocation ratios to recruitment conditions at each site in multicentre randomized trials**, by Pavel S. Roshanov.

The paper compares recruitment and statistical information under different allocation ratios. It then compares ten analysis methods when allocation ratios and prognosis vary between sites. All patient data are simulated.

## What you need

R and eight packages. The paper used R 4.6.1 and the versions below.

```r
install.packages(c(
  "survival", "lme4", "sandwich", "clubSandwich",
  "dplyr", "readr", "jsonlite", "digest"
))
```

| Package | Version | Use |
| --- | --- | --- |
| survival | 3.8-6 | Cox models, frailty models, conditional logistic regression |
| lme4 | 2.0-1 | Logistic and linear mixed models |
| sandwich | 3.1-1 | Robust and cluster-robust standard errors |
| clubSandwich | 0.7.0 | CR2 standard errors in `run_interval_checks.R` |
| dplyr | 1.2.1 | Summaries |
| readr | 2.2.0 | Reading and writing CSV files |
| jsonlite, digest | 2.0.0, 0.6.39 | Seeds for the simulated trials |

Other versions should work, but results can differ slightly. To install exactly the versions above, use the lockfile. `renv.lock` records the version of R and of every package used for the paper, including the packages that these eight depend on.

```r
install.packages("renv")
renv::activate()   # gives this folder its own package library
renv::restore()    # in a new R session: installs the recorded versions
```

Your other R projects are not affected. Use a folder whose path has no apostrophe, which can stop packages from compiling.

## The four scripts

Run them from this folder, either with `Rscript` in a terminal or with `source()` in R after setting the working directory to this folder. Results go to `output/`.

| Script | What it does | Default | Time |
| --- | --- | --- | --- |
| `make_tables.R` | Recreates every table in the paper, with its footnote, from the saved results | All 45 tables | Seconds |
| `run_simulation.R` | Simulates trials, analyzes each with the ten methods, and summarizes estimates, standard errors, and coverage | 100 replicates of the six scenarios in Tables 3 and 4 | 1 to 2 minutes |
| `run_interval_checks.R` | Table S41: normal, t, and CR2 intervals for IPTW with cluster-robust standard errors | 100 replicates of each scenario | Under a minute |
| `run_exact_logistic_checks.R` | Supplementary Appendix section 6.4: exact conditional logistic regression in three replicates of five scenarios | The 15 fits in the paper | Seconds |

The defaults are small so that everything runs quickly. Results from 100 replicates are close to the paper's but not the same. The paper used 10,000 replicates, and `data/saved_results/` holds those results.

## Running more

Options are written as `--name=value`.

```sh
Rscript run_simulation.R --replicates=1000 --cores=4
Rscript run_simulation.R --tasks=10,16,22 --replicates=1000 --output=output/table_s5
Rscript run_interval_checks.R --replicates=10000 --cores=4
```

`--tasks` selects scenarios by their row in `data/scenarios.csv`, as single numbers or ranges such as `1-99`. `data/table_tasks.csv` lists the tasks behind every table. `run_simulation.R` also writes the paper's tables for the scenarios it ran, in the `tables` folder of its output.

The complete simulation is 10,000 replicates of all 99 scenarios:

```sh
Rscript run_simulation.R --tasks=1-99 --replicates=10000 --cores=8 --output=output/full_simulation
```

This takes about 230 core-hours. Scenarios with 30 sites take 1 to 5 minutes per 1,000 replicates, and binary scenarios with 120 sites take about 45 minutes. Each scenario is saved when it finishes. If a run stops, repeat the same command and the finished scenarios are skipped.

A simulated trial depends only on its scenario and replicate number, so the number of cores does not change the results.

## Which code makes which table

| Tables | Code | Based on |
| --- | --- | --- |
| 1, 2, S1 | `R/planning_tables.R` | Planning formulas |
| 3, 4, S2–S38 | `R/simulation_tables.R` | `run_simulation.R` |
| S39 | `R/diagnostic_tables.R` | Text only |
| S40 | `R/diagnostic_tables.R` | `run_simulation.R` |
| S41 | `R/diagnostic_tables.R` | `run_interval_checks.R` |

## What is in the folders

`R/` holds the functions that the scripts use.

| File | Contents |
| --- | --- |
| `config.R` | Fixed parameters of the simulated trials |
| `scenarios.R` | Turns a row of `data/scenarios.csv` into its sites and parameter values, and sets each replicate's seed |
| `simulate.R` | Generates one trial, following the equations in the Supplementary Appendix |
| `models.R` | The ten analysis methods for each outcome, with every model formula written out |
| `fit_helpers.R` | Fits one model and keeps its treatment estimate, standard error, and any warning |
| `summaries.R` | Mean estimate, empirical SD, mean standard error, and coverage for each method |
| `run.R` | Command-line options, parallel workers, and saving results |
| `planning_tables.R`, `simulation_tables.R`, `diagnostic_tables.R` | The tables and their footnotes |
| `load.R` | Loads the packages and the files above |

`data/` holds the inputs and the paper's results.

| File | Contents |
| --- | --- |
| `scenarios.csv` | The 99 simulation scenarios, one outcome each |
| `table_tasks.csv` | The title of every table and the scenarios it shows |
| `saved_results/simulation_summary.csv` | Results of the complete simulation: one row per scenario and method |
| `saved_results/interval_checks.csv`, `interval_checks_cr2.csv` | Results of `run_interval_checks.R` at 10,000 replicates |
| `saved_results/exact_logistic_fits.csv`, `exact_logistic_strata.csv` | Results of `run_exact_logistic_checks.R` |

To try the functions directly:

```r
source("R/load.R")
scenarios <- read_scenarios()
scenario <- get_scenario(scenarios, 88)           # binary outcome, Table 3
trial <- simulate_trial(scenario, replicate_seed(scenario, 1))
analyze_trial(trial, scenario$outcome)            # all ten methods
```

## Notes

- Coverage always refers to the true within-site treatment effect. Weighted analyses and adjustment for allocation ratio alone estimate effects averaged across sites, which differ from it for non-null hazard and odds ratios.
- In the summary, `mean_estimate_display` is the exponentiated mean log hazard or odds ratio, or the mean difference. `empirical_se` is the SD of the estimates across replicates, and `model_based_se` is the mean of the reported standard errors. `n_successful` counts fits with a finite estimate and standard error. `n_warning` and `n_nonconvergence` count fits that the software flagged. For mixed models, the flags include singular fits.
- The scenario labels `Realistic` (the moderate effect) and `HighCV` (variable site sizes) enter the seeds, so they cannot be renamed.
- A complete rerun reproduces the saved results. Mixed models can differ very slightly between computers. In a rerun of 15 scenarios on a second computer, 1 printed cell changed, by 0.1 percentage points of coverage.
- The Cox models with a site random intercept use the gamma frailty that `frailty()` fits by default, with its default limit of 10 outer iterations.
- The planning benchmark in `R/planning_tables.R` counts participants. It is exact for continuous outcomes with a common residual variance and is not a power calculation for a specific endpoint.

## License and citation

The code is available under the [MIT license](LICENSE). If you use it, please cite the paper.
