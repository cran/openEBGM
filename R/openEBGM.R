#' \pkg{openEBGM}: EBGM Scores for Mining Large Contingency Tables
#'
#' \pkg{openEBGM} is a Bayesian data mining package for calculating Empirical
#' Bayes scores based on the \emph{Gamma-Poisson Shrinker} (\emph{GPS}) model
#' for large, sparse contingency (frequency) tables. \pkg{openEBGM} includes
#' several important functions implementing DuMouchel's (1999, 2001) methods for
#' calculating the EBGM (Empirical Bayes Geometric Mean) score and the quantile
#' scores used to create credibility intervals. Some simple disproportionality
#' scores (relative report rate and proportional reporting ratio) are also
#' included. Adverse event report data are used as an example application. Much
#' of \pkg{openEBGM}'s code is derived from the \pkg{PhViD} and \pkg{mederrRank}
#' packages.
#'
#' @section Data preparation & squashing functions:
#' The data preparation function, \code{\link{processRaw}}, converts raw data
#'   into actual and expected counts for product/event pairs.
#'   \code{\link{processRaw}} also adds the relative reporting ratio (RR) and
#'   proportional reporting ratio (PRR). The data squashing function,
#'   \code{\link{squashData}}, implements the simple version of data squashing
#'   described in DuMouchel et al. (2001). Data squashing can be used to reduce
#'   computational burden.
#'
#' @section Negative log-likelihood functions:
#' The negative log-likelihood functions (\code{\link{negLL}},
#'   \code{\link{negLLsquash}}, \code{\link{negLLzero}}, and
#'   \code{\link{negLLzeroSquash}}) provide the means of calculating the
#'   negative log-likelihoods as mentioned in the DuMouchel papers. DuMouchel
#'   uses the likelihood function, based on the marginal distributions of the
#'   counts, to estimate the hyperparameters of the prior distribution.
#'
#' @section Hyperparameter estimation functions:
#' The hyperparameter estimation functions (\code{\link{exploreHypers}} and
#'   \code{\link{autoHyper}}) use gradient-based approaches to estimate the
#'   hyperparameters, \eqn{\theta}, of the prior distribution (gamma mixture)
#'   using the negative log-likelihood functions from the marginal distributions
#'   of the counts (negative binomial). \eqn{\theta} is a vector containing five
#'   parameters (\eqn{\alpha_1}, \eqn{\beta_1}, \eqn{\alpha_2}, \eqn{\beta_2},
#'   and \eqn{P}). \code{\link{hyperEM}} estimates \eqn{\theta} using a version
#'   of the EM algorithm.
#'
#' @section Posterior distribution functions:
#' The posterior distribution functions calculate the mixture fraction
#'   (\code{\link{Qn}}), geometric mean (\code{\link{ebgm}}), and quantiles
#'   (\code{\link{quantBisect}}) of the posterior distribution. Alternatively,
#'   \code{\link{ebScores}} can be used to create an object of class openEBGM
#'   that contains the EBGM and quantiles scores. Appropriate methods exist for
#'   the generic functions \code{\link[base]{print}},
#'   \code{\link[base]{summary}}, and \code{\link[graphics]{plot}} for openEBGM
#'   objects.
#'
#' @references Ahmed I, Poncet A (2016). \pkg{PhViD}: an R package for
#'   PharmacoVigilance signal Detection. \emph{R package version 1.0.8}.
#'
#' @references Venturini S, Myers J (2015). \pkg{mederrRank}: Bayesian Methods
#'   for Identifying the Most Harmful Medication Errors. \emph{R package version
#'   0.0.8}.
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#'
#' @references Evans SJW, Waller P, Davis S (2001). "Use of Proportional
#'   Reporting Ratios (PRRs) for Signal Generation from Spontaneous Adverse Drug
#'   Reaction Reports." \emph{Pharmacoepidemiology and Drug Safety}, 10(6),
#'   483-486.
#'
#' @references FDA (2017). "CFSAN Adverse Event Reporting System (CAERS)."
#'   URL \url{https://www.fda.gov/Food/ComplianceEnforcement/ucm494015.htm}.
#'
#' @docType package
#' @name openEBGM
NULL
