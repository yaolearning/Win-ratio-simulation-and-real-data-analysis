\name{run_win_simulation}
\alias{run_win_simulation}
\alias{write_win_simulation_results}
\title{Run and Export MaxWin Simulation Studies}
\description{
Run repeated scenario-based simulation studies and summarize empirical rejection proportions, selection behavior, tie summaries, and log-rank comparator results.
}
\usage{
run_win_simulation(scenarios = NULL,
  scenario_grid = default_win_scenarios(),
  nsim = 1000L, B = 500L, seed = 2026L, alpha = 0.05,
  output_dir = NULL, checkpoint_every = 25L, resume = FALSE,
  verbose = TRUE, save_plots = TRUE,
  save_scenario_tables = TRUE)

write_win_simulation_results(x, output_dir, prefix = "GLOBAL",
  save_plots = TRUE, save_scenario_tables = TRUE)
}
\arguments{
\item{scenarios}{Selected scenario indices, identifiers, or a scenario data frame. Use \code{NULL} for all scenarios.}
\item{scenario_grid}{Scenario table.}
\item{nsim}{Number of simulated trials per scenario.}
\item{B}{Number of treatment-label permutations per simulated trial.}
\item{seed}{Master random seed.}
\item{alpha}{Significance level used for empirical rejection summaries.}
\item{output_dir}{Optional output directory.}
\item{checkpoint_every}{Interval for simulation-level checkpoint files.}
\item{resume}{Whether an existing simulation-level checkpoint is resumed.}
\item{verbose}{Whether progress is printed.}
\item{save_plots}{Whether simulation plots are written.}
\item{save_scenario_tables}{Whether per-scenario output tables are written.}
\item{x}{A \code{win_simulation} object.}
\item{prefix}{File prefix for exported simulation results.}
}
\value{
\code{run_win_simulation()} returns an object of class \code{win_simulation}. \code{write_win_simulation_results()} invisibly returns the output directory.
}
