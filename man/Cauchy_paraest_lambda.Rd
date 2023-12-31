% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{Cauchy_paraest_lambda}
\alias{Cauchy_paraest_lambda}
\title{Predict the scale parameter of a Cauchy distribution using Rcpp.}
\usage{
Cauchy_paraest_lambda(x, b)
}
\arguments{
\item{x}{the sorted i.i.d Cauchy variables (NumericVector)}

\item{b}{the method (1, 2 or 3) (int)}
}
\value{
the point estimation of the scale parameter
}
\description{
Predict the scale parameter of a Cauchy distribution using Rcpp.
}
\examples{
\dontrun{
x<-rcauchy(50,0,1)
Cauchy_paraest_lambda(x,1)
}
}
