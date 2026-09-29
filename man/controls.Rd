\name{MaxWin-controls}
\alias{threshold_control}
\alias{weight_control}
\alias{permutation_control}
\alias{method_control}
\alias{engine_control}
\title{Control MaxWin Analyses}
\description{
Construct control objects for threshold search, endpoint weights, permutation inference, method selection, and computational backend.
}
\usage{
threshold_control(time_endpoints = "auto", candidates = NULL,
  unit = "months", data_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
  max_event_times = 800L, max_thresholds = 10L, seed = 999L)

weight_control(mode = c("vertices", "grid"), step = 0.05,
  constraint = c("ordered", "simplex"))

permutation_control(enabled = TRUE, B = 500L, seed = 2026L)

method_control(run_wr = TRUE, run_wo = TRUE,
  report_all_fixed_orders = TRUE)

engine_control(backend = c("auto", "general", "two_endpoint_recurrent"),
  eps = 1e-8)
}
\arguments{
\item{time_endpoints}{Time endpoints eligible for thresholding, or \code{"auto"}.}
\item{candidates}{User-supplied threshold candidates.}
\item{unit}{Unit for threshold candidates.}
\item{data_probs}{Quantiles used for data-driven threshold generation.}
\item{max_event_times}{Maximum number of event times used to construct data-driven candidates.}
\item{max_thresholds}{Maximum number of threshold candidates.}
\item{seed}{Random seed.}
\item{mode}{Weight-grid mode.}
\item{step}{Grid spacing when \code{mode = "grid"}.}
\item{constraint}{Weight constraint.}
\item{enabled}{Whether permutation inference is run.}
\item{B}{Number of treatment-label permutations.}
\item{run_wr}{Whether Win Ratio methods are evaluated.}
\item{run_wo}{Whether Win Odds methods are evaluated.}
\item{report_all_fixed_orders}{Whether all fixed endpoint orders are reported.}
\item{backend}{Computational backend.}
\item{eps}{Small positive constant used by the specialized two-endpoint engine for Win Ratio stabilization.}
}
\value{
A control object used by other MaxWin functions.
}
