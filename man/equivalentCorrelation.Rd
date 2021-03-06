% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tests.R
\name{equivalentCorrelation}
\alias{equivalentCorrelation}
\title{Run equivalence tests on a bivariate correlation.}
\usage{
equivalentCorrelation(x, y, lleq, uleq)
}
\arguments{
\item{x}{A numeric vector with no missing data.}

\item{y}{A numeric vector with no missing data; same length as x.}

\item{lleq}{A value between -1 and 1 indicating the lower limit of the equivalence interval.}

\item{uleq}{A value between -1 and 1 indicating the upper limit of the equivalence interval.}
}
\value{
NULL.
}
\description{
This procedure works with the raw data for the two variables. No provision
has been made for missing data and the two vectors of data must be the
same length.
}
\examples{
equivalentCorrelation(x=iris$Sepal.Length, y=iris$Sepal.Width, lleq=-0.25, uleq=0.25)
}
