% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estlambda.R
\name{estlambda}
\alias{estlambda}
\title{Estimate the inflation factor for a distribution of P-values}
\usage{
estlambda(
  data,
  plot = FALSE,
  proportion = 1,
  method = "regression",
  filter = TRUE,
  df = 1,
  ...
)
}
\arguments{
\item{data}{A vector of reals. If all are <=1, it is assumed that this is a
vector of P-values, else it is treated as a vector of chi-squares}

\item{plot}{Whether the plot should be shown or not (default).}

\item{proportion}{The proportion of lowest P (or
\eqn{\chi^2}{chi^2}) values to be used when estimating the inflation
factor \eqn{\lambda}{lambda}.}

\item{method}{"regression" (default), "median", or "KS": method to
be used for \eqn{\lambda}{lambda} estimation.}

\item{filter}{if the test statistics with 0-value of
\eqn{\chi^2}{chi^2} should be excluded prior to estimation of
\eqn{\lambda}{lambda}.}

\item{df}{Number of degrees of freedom.}

\item{...}{arguments passed to the \code{\link{plot}} function.}
}
\value{
A list with elements 
  \item{estimate}{Estimate of \eqn{\lambda}{lambda}}
  \item{se}{Standard error of the estimate}
}
\description{
NOTE: This is copied (with slight modifications) from the now-defunct GenABEL 
package. The original 'KS' method has been removed as it had many internal
GenABEL dependencies. All credit goes to the function author Yurii Aulchenko.
}
\details{
Estimate the inflation factor for a distribution of P-values or 1df
chi-square test. The major use of this procedure is the Genomic Control, but
can also be used to visualise the distribution of P-values coming from other
tests. Methods implemented include 'median' (median(chi2)/0.455...),
regression (of observed onto expected) and 'KS' (optimizing the
chi2.1df distribution fit by use of Kolmogorov-Smirnov test)
}
\author{
Yurii Aulchenko
}
