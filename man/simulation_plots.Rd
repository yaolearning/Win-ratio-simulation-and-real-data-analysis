\name{MaxWin-simulation-plots}
\alias{plot_win_simulation_power}
\alias{plot_win_simulation_gain}
\alias{plot_win_simulation_scenario_gain}
\alias{save_win_simulation_plots}
\title{Plot MaxWin Simulation Results}
\description{
Plot empirical rejection proportions and gains relative to a baseline method from a \code{win_simulation} object.
}
\usage{
plot_win_simulation_power(x,
  measure = c("WR", "WO", "Log-rank"),
  side = c("one", "two"), main = NULL)

plot_win_simulation_gain(x, measure = c("WR", "WO"),
  side = c("one", "two"), baseline_method = "original",
  main = NULL)

plot_win_simulation_scenario_gain(x, scenario,
  measure = c("WR", "WO"), side = c("one", "two"),
  baseline_method = "original", main = NULL)

save_win_simulation_plots(x, output_dir, prefix = "SIM")
}
\arguments{
\item{x}{A \code{win_simulation} object.}
\item{measure}{Statistic or comparator to display.}
\item{side}{One-sided or two-sided result.}
\item{main}{Optional plot title.}
\item{baseline_method}{Method identifier used as the gain reference.}
\item{scenario}{Scenario index or identifier.}
\item{output_dir}{Directory in which figures are saved.}
\item{prefix}{Output file prefix.}
}
\value{
The plotting functions invisibly return the plotted summary object. \code{save_win_simulation_plots()} invisibly returns the output directory.
}
