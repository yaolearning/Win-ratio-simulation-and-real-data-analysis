\name{MaxWin-analysis-plots}
\alias{plot_win_pvalues}
\alias{plot_win_km}
\alias{plot_win_ties}
\alias{save_win_plots}
\title{Plot MaxWin Real-Data Results}
\description{
Create permutation p-value plots, Kaplan-Meier plots, and tie summaries from a \code{win_analysis} object.
}
\usage{
plot_win_pvalues(x, measure = c("WR", "WO"),
  side = c("one", "two"), alpha = 0.05, main = NULL)

plot_win_km(x, endpoint = 1L, main = NULL)

plot_win_ties(x, measure = c("WR", "WO"),
  side = c("one", "two"), main = NULL)

save_win_plots(x, output_dir, prefix = "WIN", alpha = 0.05)
}
\arguments{
\item{x}{A \code{win_analysis} object.}
\item{measure}{Statistic to plot, WR or WO.}
\item{side}{One-sided or two-sided result.}
\item{alpha}{Reference significance level.}
\item{main}{Optional plot title.}
\item{endpoint}{Time-to-event endpoint number for the Kaplan-Meier plot.}
\item{output_dir}{Directory in which figures are saved.}
\item{prefix}{Output file prefix.}
}
\value{
The plotting functions invisibly return the data or fitted object used for plotting. \code{save_win_plots()} invisibly returns the output directory.
}
