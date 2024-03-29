% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_transformInput.R
\name{processRaw}
\alias{processRaw}
\title{Process raw data}
\usage{
processRaw(
  data,
  stratify = FALSE,
  zeroes = FALSE,
  digits = 2,
  max_cats = 10,
  list_ids = FALSE
)
}
\arguments{
\item{data}{A data frame containing columns named: \emph{id}, \emph{var1},
and \emph{var2}. Possibly includes columns for stratification variables
identified by the substring \emph{'strat'} (e.g. \emph{strat_age}). Other
columns will be ignored.}

\item{stratify}{A logical scalar (TRUE or FALSE) specifying if stratification
is to be used.}

\item{zeroes}{A logical scalar specifying if zero counts should be included.
Using zero counts is only recommended for small data sets because it will
dramatically increase run time.}

\item{digits}{A whole number scalar specifying how many decimal places to
use for rounding \emph{RR} and \emph{PRR}.}

\item{max_cats}{A whole number scalar specifying the maximum number of
categories to allow in any given stratification variable. Used to help
prevent a situation where the user forgets to categorize a continuous
variable, such as age.}

\item{list_ids}{A logical scalar specifying if a column for pipe-concatenated
IDs should be returned.}
}
\value{
A data frame with actual counts (\emph{N}), expected counts
  (\emph{E}), relative reporting ratio (\emph{RR}), and proportional
  reporting ratio (\emph{PRR}) for \emph{var1-var2} pairs. Also includes a
  column for IDs (\emph{ids}) if \code{list_ids = TRUE}.
}
\description{
\code{processRaw} finds the actual and expected counts using the methodology
described by DuMouchel (1999); however, an adjustment is applied to expected
counts to prevent double-counting (i.e., using unique marginal ID counts
instead of contingency table marginal totals). Also calculates the relative
reporting ratio (\emph{RR}) and the proportional reporting ratio
(\emph{PRR}).
}
\details{
An \emph{id} column must be included in \code{data}. If your data
  set does not include IDs, make a column of unique IDs using \code{df$id <-
  1:nrow(df)}. However, unique IDs should only be constructed if the cells in
  the contingency table are mutually exclusive. For instance, unique IDs for
  each row in \code{data} are not appropriate with CAERS data since a given
  report can include multiple products and/or adverse events.

Stratification variables are identified by searching for the
  substring \emph{'strat'}. Only variables containing \emph{'strat'} (case
  sensitive) will be used as stratification variables. \emph{PRR}
  calculations ignore stratification, but \emph{E} and \emph{RR} calculations
  are affected. A warning will be displayed if any stratum contains less than
  50 unique IDs.

If a \emph{PRR} calculation results in division by zero, \code{Inf}
  is returned.
}
\section{Warnings}{

  Use of the \code{zeroes = TRUE} option will result in a considerable
  increase in runtime. Using zero counts is not recommended if the contingency
  table is moderate or large in size (~500K cells or larger). However, using
  zeroes could be useful if the optimization algorithm fails to converge
  when estimating hyperparameters.


Any columns in \code{data} containing the substring \emph{'strat'} in the
  column name will be considered stratification variables, so verify that you
  do not have any extraneous columns with that substring.
}

\examples{
data.table::setDTthreads(2)  #only needed for CRAN checks
var1 <- c("product_A", rep("product_B", 3), "product_C",
          rep("product_A", 2), rep("product_B", 2), "product_C")
var2 <- c("event_1", rep("event_2", 2), rep("event_3", 2),
          "event_2", rep("event_3", 3), "event_1")
strat1 <- c(rep("Male", 5), rep("Female", 3), rep("Male", 2))
strat2 <- c(rep("age_cat1", 5), rep("age_cat1", 3), rep("age_cat2", 2))
dat <- data.frame(
  var1 = var1, var2 = var2, strat1 = strat1, strat2 = strat2,
  stringsAsFactors = FALSE
)
dat$id <- 1:nrow(dat)
processRaw(dat)
suppressWarnings(
  processRaw(dat, stratify = TRUE)
)
processRaw(dat, zeroes = TRUE)
suppressWarnings(
  processRaw(dat, stratify = TRUE, zeroes = TRUE)
)
processRaw(dat, list_ids = TRUE)

}
\references{
DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
  Tables, With an Application to the FDA Spontaneous Reporting System."
  \emph{The American Statistician}, 53(3), 177-190.
}
\keyword{openEBGM}
