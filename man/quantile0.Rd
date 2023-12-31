% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{quantile0}
\alias{quantile0}
\title{Produce a quantile number using Rcpp.}
\usage{
quantile0(x, q)
}
\arguments{
\item{x}{a numeric vector}

\item{q}{the order of the quantile}
}
\value{
the quantile of order q of vector x
}
\description{
Produce the quantile of order q of a numeric vector.
}
\examples{
\dontrun{
x<-rcauchy(5,0,1)
quantile(x,0.5)
}
}
