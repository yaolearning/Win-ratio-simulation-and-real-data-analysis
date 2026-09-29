\name{MaxWin-simulation-data}
\alias{generate_win_dataset}
\alias{default_win_scenarios}
\alias{generate_win_scenario_dataset}
\title{Generate Simulation Data and Scenarios}
\description{
Generate simulated randomized-trial datasets with survival and recurrent-event outcomes, inspect the default scenario grid, or generate one dataset from a selected scenario.
}
\usage{
generate_win_dataset(n_control = 50L, n_treatment = 50L,
  mort_rate_control = -log(0.6), hazard_ratio = 1,
  hosp_shape = 5, hosp_scale_control = 1,
  hosp_scale_treatment = 1, followup = 1,
  censor_rate = 0, seed = 2026L)

default_win_scenarios()

generate_win_scenario_dataset(scenario = 1L,
  scenario_grid = default_win_scenarios(),
  sim_index = 1L, seed = 2026L)
}
\arguments{
\item{n_control}{Number of control subjects.}
\item{n_treatment}{Number of treatment subjects.}
\item{mort_rate_control}{Control-group exponential mortality rate.}
\item{hazard_ratio}{Treatment-to-control mortality hazard ratio.}
\item{hosp_shape}{Gamma shape parameter for subject-specific recurrent-event rates.}
\item{hosp_scale_control}{Control-group gamma scale parameter.}
\item{hosp_scale_treatment}{Treatment-group gamma scale parameter.}
\item{followup}{Administrative follow-up duration.}
\item{censor_rate}{Independent exponential censoring rate.}
\item{seed}{Random seed.}
\item{scenario}{Scenario index, scenario identifier, or one-row scenario data frame.}
\item{scenario_grid}{Scenario table.}
\item{sim_index}{Simulation replicate index.}
}
\value{
\code{generate_win_dataset()} and \code{generate_win_scenario_dataset()} return simulated dataset objects. \code{default_win_scenarios()} returns the default 15-scenario data frame.
}
