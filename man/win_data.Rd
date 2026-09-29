\name{win_data}
\alias{win_data}
\title{Create a MaxWin Data Input}
\description{
Create a data source object from an in-memory object, external file, or installed R package, optionally with a separate recurrent-event table.
}
\usage{
win_data(data = NULL, file = NULL, package = NULL, dataset = NULL,
  object_name = NULL, recurrent_data = NULL, recurrent_file = NULL,
  recurrent_package = NULL, recurrent_dataset = NULL,
  recurrent_object_name = NULL)
}
\arguments{
\item{data}{Primary in-memory data object.}
\item{file}{Path to a primary CSV, TXT, RDS, RData, or RDA file.}
\item{package}{Installed R package containing the primary dataset.}
\item{dataset}{Dataset name in \code{package}.}
\item{object_name}{Object name to extract from an RData/RDA file.}
\item{recurrent_data}{Optional in-memory recurrent-event table.}
\item{recurrent_file}{Optional recurrent-event data file.}
\item{recurrent_package}{Optional package containing the recurrent-event dataset.}
\item{recurrent_dataset}{Dataset name for recurrent-event data.}
\item{recurrent_object_name}{Object name for recurrent-event data in an RData/RDA file.}
}
\value{
An object of class \code{win_data}.
}
