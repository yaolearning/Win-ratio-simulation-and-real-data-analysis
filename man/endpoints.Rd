\name{MaxWin-endpoints}
\alias{endpoint_time}
\alias{endpoint_count}
\alias{endpoint_binary}
\alias{endpoint_continuous}
\title{Define Hierarchical Endpoints}
\description{
Create endpoint specifications for use in \code{win_analysis()}.
}
\usage{
endpoint_time(name, time = NULL, event = NULL, unit = "years",
  status_code = NULL)

endpoint_count(name, count = NULL,
  comparison = c("observed_total", "pairwise_common_followup"),
  recurrent_time = NULL, recurrent_id = NULL, followup = NULL,
  unit = "years", followup_unit = NULL,
  recurrent_time_type = c("absolute", "gap"), status_code = NULL)

endpoint_binary(name, value, adverse = TRUE)

endpoint_continuous(name, value, higher_better = TRUE)
}
\arguments{
\item{name}{Endpoint display name.}
\item{time}{Column containing time-to-event or follow-up time.}
\item{event}{Column containing the event indicator, coded 1 for event and 0 for censored.}
\item{unit}{Time unit for the endpoint.}
\item{status_code}{Event code used for event-long data.}
\item{count}{Column containing an observed count.}
\item{comparison}{Count comparison rule. Use \code{"pairwise_common_followup"} to compare recurrent-event counts within the common follow-up time of each treatment-control pair.}
\item{recurrent_time}{Column containing recurrent-event times.}
\item{recurrent_id}{Subject identifier column in a recurrent-event table.}
\item{followup}{Subject-level follow-up column used for pairwise common follow-up.}
\item{followup_unit}{Time unit for the follow-up column.}
\item{recurrent_time_type}{Whether recurrent-event times are absolute times or gap times.}
\item{value}{Column containing a binary or continuous endpoint.}
\item{adverse}{For binary endpoints, whether larger values represent worse outcomes.}
\item{higher_better}{For continuous endpoints, whether larger values represent better outcomes.}
}
\value{
An object of class \code{win_endpoint}.
}
