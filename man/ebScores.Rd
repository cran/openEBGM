% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_createObject.R
\name{ebScores}
\alias{ebScores}
\title{Construct an openEBGM object}
\usage{
ebScores(processed, hyper_estimate, quantiles = c(5, 95), digits = 2)
}
\arguments{
\item{processed}{A data frame resulting from running \code{\link{processRaw}}.}

\item{hyper_estimate}{A list resulting from running \code{\link{autoHyper}}.}

\item{quantiles}{Either a numeric vector of desired quantiles to be
calculated from the posterior distribution or NULL for no calculation of
quantiles.}

\item{digits}{A whole number scalar specifying how many decimal places to
use for rounding \emph{EBGM} and the quantiles scores.}
}
\value{
An openEBGM object (class S3) containing:
  \itemize{
    \item{\emph{data}: }{A data frame containing the results (scores, etc.).}
    \item{\emph{hyper_parameters}: }{A list containing the hyperparameter
                                     estimation results.}
    \item{\emph{quantiles}: }{The chosen percentiles.}
  }
}
\description{
\code{ebScores} calculates EBGM scores as well as the quantiles from the
posterior distribution and returns an object of class openEBGM.
}
\details{
This function takes the processed data as well as the hyperparameter
  estimates and instantiates an object of class openEBGM. This object then
  contains additional calculations, such as the EBGM score, and the quantiles
  that are supplied by the quantiles parameter at the time of calling the
  function.

The function allows for the choice of an arbitrary amount of
  quantiles or no quantiles at all to be calculated. This may be helpful for
  large datasets, or when the EBGM score is the only metric of interest.
}
\examples{
data.table::setDTthreads(2)  #only needed for CRAN checks
theta_init <- data.frame(
  alpha1 = c(0.5, 1),
  beta1  = c(0.5, 1),
  alpha2 = c(2,   3),
  beta2  = c(2,   3),
  p      = c(0.1, 0.2)
)
data(caers)
proc <- processRaw(caers)
squashed <- squashData(proc, bin_size = 300, keep_pts = 10)
squashed <- squashData(squashed, count = 2, bin_size = 13, keep_pts = 10)
suppressWarnings(
  hypers <- autoHyper(data = squashed, theta_init = theta_init)
)
obj <- ebScores(processed = proc, hyper_estimate = hypers, quantiles = 5)
print(obj)

}
\keyword{openEBGM}
