\name{win_analysis}
\alias{win_analysis}
\alias{preview_win_plan}
\alias{write_win_results}
\title{Run Win Ratio and Win Odds Analyses}
\description{
Run standard and adaptive Win Ratio and Win Odds analyses for two- or three-endpoint hierarchical composite outcomes.
}
\usage{
preview_win_plan(data, endpoints, format = c("wide", "event_long"),
  id = "SUBJID", treatment = "ARM", event_time = NULL,
  event_status = NULL, event_time_unit = "years",
  threshold = threshold_control(), weight = weight_control(),
  methods = method_control(), engine = engine_control())

win_analysis(data, endpoints, format = c("wide", "event_long"),
  id = "SUBJID", treatment = "ARM", event_time = NULL,
  event_status = NULL, event_time_unit = "years",
  threshold = threshold_control(), weight = weight_control(),
  methods = method_control(), permutation = permutation_control(),
  engine = engine_control(), run_logrank = TRUE,
  run_composite_logrank_if_possible = TRUE, output_dir = NULL,
  output_prefix = "WIN")

write_win_results(x, output_dir, prefix = "WIN", save_plots = TRUE)
}
\arguments{
\item{data}{A \code{win_data} object or supported input data object.}
\item{endpoints}{List of endpoint specifications created with \code{endpoint_time()}, \code{endpoint_count()}, \code{endpoint_binary()}, or \code{endpoint_continuous()}.}
\item{format}{Input format, either \code{"wide"} or \code{"event_long"}.}
\item{id}{Subject identifier column.}
\item{treatment}{Treatment-assignment column coded 1 for treatment and 0 for control.}
\item{event_time}{Event-time column for event-long data.}
\item{event_status}{Event-status column for event-long data.}
\item{event_time_unit}{Time unit for event-long data.}
\item{threshold}{Threshold control object.}
\item{weight}{Weight control object.}
\item{methods}{Method control object.}
\item{permutation}{Permutation control object.}
\item{engine}{Engine control object.}
\item{run_logrank}{Whether log-rank comparator analyses are run for time-to-event endpoints.}
\item{run_composite_logrank_if_possible}{Whether a time-to-first-event composite log-rank test is added when all endpoints are time-to-event.}
\item{output_dir}{Optional output directory.}
\item{output_prefix}{Prefix for automatically written output files.}
\item{x}{A \code{win_analysis} object.}
\item{prefix}{Prefix for files written by \code{write_win_results()}.}
\item{save_plots}{Whether figures are saved with the tabular output.}
}
\value{
\code{preview_win_plan()} returns the planned candidate grids. \code{win_analysis()} returns an object of class \code{win_analysis}. \code{write_win_results()} invisibly returns the output directory.
}
\details{
One-sided adaptive selection uses the raw WR or WO statistic. Two-sided selection uses the absolute log statistic. For adaptive methods, endpoint order, weights, and thresholds are re-selected within each treatment-label permutation.
}
