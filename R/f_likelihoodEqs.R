#' Likelihood with zero counts
#'
#' \code{negLLzero} computes the negative log-likelihood based on the
#' unconditional marginal distribution of \emph{N} (equation 12 in DuMouchel
#' 1999, except taking negative natural log). This function is minimized to
#' estimate the hyperparameters of the prior distribution. Use this function if
#' including zero counts but not squashing data. Generally this function is not
#' recommended (\code{\link{negLLsquash}} is typically more efficient).
#'
#' @param theta A numeric vector of hyperparameters ordered as:
#'   \eqn{\alpha_1, \beta_1, \alpha_2, \beta_2, P}.
#' @param N A whole number vector of actual counts from
#'   \code{\link{processRaw}}.
#' @param E A numeric vector of expected counts from \code{\link{processRaw}}.
#' @return A scalar negative log-likelihood value.
#'
#' @details The marginal distribution of the counts, \emph{N}, is a mixture of
#'   two negative binomial distributions. The hyperparameters for the prior
#'   distribution (mixture of gammas) are estimated by optimizing the likelihood
#'   equation from this marginal distribution.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameters of the second component}
#'     \item{\eqn{P}: }{Mixture fraction}
#' }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N} actually contains zeroes before using this function. You
#'   should have used the \code{zeroes = TRUE} option when calling the
#'   \code{\link{processRaw}} function.
#' @section Warnings:
#'   Make sure the data were not squashed before using this function.
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family negative log-likelihood functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}, and
#'   \code{\link[stats]{optim}} for optimization
#' @importFrom stats dnbinom
#' @export
negLLzero <- function(theta, N, E) {

  .checkInputs_negLLzero(theta, N, E)

  size_f1 <- theta[1]  #alpha1
  prob_f1 <- theta[2] / (theta[2] + E)  #beta1 / (beta1 + E)
  size_f2 <- theta[3]  #alpha2
  prob_f2 <- theta[4] / (theta[4] + E)  #beta2 / (beta2 + E)

  f1 <- dnbinom(N, size = size_f1, prob = prob_f1)
  f2 <- dnbinom(N, size = size_f2, prob = prob_f2)

  P <- theta[5]
  L <- (P * f1) + ((1 - P) * f2)
  sum(-log(L))
}


#' Likelihood with data squashing & zero counts
#'
#' \code{negLLzeroSquash} computes the negative log-likelihood based on the
#' unconditional marginal distribution of \emph{N} (DuMouchel et al. 2001).
#' This function is minimized to estimate the hyperparameters of the prior
#' distribution. Use this function if including zero counts and using data
#' squashing. Generally this function is not recommended unless convergence
#' issues occur without zero counts (\code{\link{negLLsquash}} is typically more
#' efficient).
#'
#' @inheritParams negLLzero
#' @param ni A whole number vector of squashed actual counts from
#'   \code{\link{squashData}}.
#' @param ei A numeric vector of squashed expected counts from
#'   \code{\link{squashData}}.
#' @param wi A whole number vector of bin weights from \code{\link{squashData}}.
#' @return A scalar negative log-likelihood value.
#'
#' @details The marginal distribution of the counts, \emph{N}, is a mixture of
#'   two negative binomial distributions. The hyperparameters for the prior
#'   distribution (mixture of gammas) are estimated by optimizing the likelihood
#'   equation from this marginal distribution.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameters of the second component}
#'     \item{\eqn{P}: }{Mixture fraction}
#' }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{ni} actually contains zeroes before using this function.
#'   You should have used the \code{zeroes = TRUE} option when calling the
#'   \code{processRaw} function.
#' @section Warnings:
#'   Make sure the data were actually squashed (see \code{\link{squashData}})
#'   before using this function.
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#' @family negative log-likelihood functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}, and
#'   \code{\link[stats]{optim}} for optimization and \code{\link{squashData}}
#'   for data squashing
#' @importFrom stats dnbinom
#' @export
negLLzeroSquash <- function(theta, ni, ei, wi) {

  .checkInputs_negLLzeroSquash(theta, ni, ei, wi)

  size_f1 <- theta[1]  #alpha1
  prob_f1 <- theta[2] / (theta[2] + ei)  #beta1 / (beta1 + E)
  size_f2 <- theta[3]  #alpha2
  prob_f2 <- theta[4] / (theta[4] + ei)  #beta2 / (beta2 + E)

  f1 <- dnbinom(ni, size = size_f1, prob = prob_f1)
  f2 <- dnbinom(ni, size = size_f2, prob = prob_f2)

  P <- theta[5]
  logL <- wi * (log((P * f1) + ((1 - P) * f2)))
  sum(-logL)
}


#' Likelihood without zero counts
#'
#' \code{negLL} computes the negative log-likelihood based on the conditional
#' marginal distribution of the counts, \emph{N}, given that \emph{N >= N*},
#' where \emph{N*} is the smallest count used for estimating the hyperparameters
#' (DuMouchel et al. 2001). This function is minimized to estimate the
#' hyperparameters of the prior distribution. Use this function when neither
#' zero counts nor data squashing are being used. Generally this function is not
#' recommended unless using a small data set since data squashing (see
#' \code{\link{squashData}} and \code{\link{negLLsquash}}) can increase
#' efficiency (DuMouchel et al. 2001).
#'
#' @inheritParams negLLzero
#' @param N_star A scalar whole number for the minimum count size used.
#' @return A scalar negative log-likelihood value
#'
#' @details The conditional marginal distribution for the counts, \emph{N},
#'   given that \emph{N >= N*}, is based on a mixture of two negative binomial
#'   distributions. The hyperparameters for the prior distribution (mixture of
#'   gammas) are estimated by optimizing the likelihood equation from this
#'   conditional marginal distribution. It is recommended to use \code{N_star =
#'   1} when practical.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameters of the second component}
#'     \item{\eqn{P}: }{Mixture fraction}
#' }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N_star} matches the smallest actual count in \emph{N}
#'   before using this function. Filter \emph{N} and \emph{E} if needed.
#' @section Warnings:
#'   Make sure the data were not squashed before using this function.
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#' @family negative log-likelihood functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}, and
#'   \code{\link[stats]{optim}} for optimization
#' @importFrom stats pnbinom
#' @export
negLL <- function(theta, N, E, N_star = 1) {

  .checkInputs_negLL(theta, N, E, N_star)

  size_f1 <- theta[1]  #alpha1
  prob_f1 <- theta[2] / (theta[2] + E)  #beta1 / (beta1 + E)
  size_f2 <- theta[3]  #alpha2
  prob_f2 <- theta[4] / (theta[4] + E)  #beta2 / (beta2 + E)

  f1 <- dnbinom(N, size = size_f1, prob = prob_f1)
  f2 <- dnbinom(N, size = size_f2, prob = prob_f2)

  sum_limit <- N_star - 1
  f1_cumul  <- pnbinom(sum_limit, size = size_f1, prob = prob_f1)
  f2_cumul  <- pnbinom(sum_limit, size = size_f2, prob = prob_f2)
  f_star1   <- f1 / (1 - f1_cumul)
  f_star2   <- f2 / (1 - f2_cumul)

  P <- theta[5]
  L <- (P * f_star1) + ((1 - P) * f_star2)
  sum(-log(L))
}


#' Likelihood with data squashing and no zero counts
#'
#' \code{negLLsquash} computes the negative log-likelihood based on the
#' conditional marginal distribution of the counts, \emph{N}, given that
#' \emph{N >= N*}, where \emph{N*} is the smallest count used for estimating the
#' hyperparameters. This function is minimized to estimate the hyperparameters
#' of the prior distribution. Use this function when zero counts are not used
#' and data squashing is used as described by DuMouchel et al. (2001). This
#' function is the likelihood function that should usually be chosen.
#'
#' @inheritParams negLLzeroSquash
#' @inheritParams negLL
#' @return A scalar negative log-likelihood value
#'
#' @details The conditional marginal distribution for the counts, \emph{N},
#'   given that \emph{N >= N*}, is based on a mixture of two negative binomial
#'   distributions. The hyperparameters for the prior distribution (mixture of
#'   gammas) are estimated by optimizing the likelihood equation from this
#'   conditional marginal distribution. It is recommended to use \code{N_star =
#'   1} when practical.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameters of the second component}
#'     \item{\eqn{P}: }{Mixture fraction}
#' }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N_star} matches the smallest actual count in \emph{ni}
#'   before using this function. Filter \emph{ni}, \emph{ei}, and \emph{wi} if
#'   needed.
#' @section Warnings:
#'   Make sure the data were actually squashed (see \code{\link{squashData}})
#'   before using this function.
#'
#' @examples
#' theta_init <- c(0.2, 0.1, 2, 4, 1/3)  #initial guess
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, count = 1, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' negLLsquash(theta = theta_init, ni = squashed$N, ei = squashed$E,
#'             wi = squashed$weight)
#' #For hyperparameter estimation...
#' stats::nlminb(start = theta_init, objective = negLLsquash, ni = squashed$N,
#'               ei = squashed$E, wi = squashed$weight)
#'
#' @references DuMouchel W, Pregibon D (2001). "Empirical Bayes Screening for
#'   Multi-item Associations." In \emph{Proceedings of the Seventh ACM SIGKDD
#'   International Conference on Knowledge Discovery and Data Mining}, KDD '01,
#'   pp. 67-76. ACM, New York, NY, USA. ISBN 1-58113-391-X.
#' @family negative log-likelihood functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}, and
#'   \code{\link[stats]{optim}} for optimization and \code{\link{squashData}}
#'   for data squashing
#' @importFrom stats dnbinom
#' @importFrom stats pnbinom
#' @export
negLLsquash <- function(theta, ni, ei, wi, N_star = 1) {

  .checkInputs_negLLsquash(theta, ni, ei, wi, N_star)

  size_f1 <- theta[1]  #alpha1
  prob_f1 <- theta[2] / (theta[2] + ei)  #beta1 / (beta1 + E)
  size_f2 <- theta[3]  #alpha2
  prob_f2 <- theta[4] / (theta[4] + ei)  #beta2 / (beta2 + E)

  f1 <- dnbinom(ni, size = size_f1, prob = prob_f1)
  f2 <- dnbinom(ni, size = size_f2, prob = prob_f2)

  sum_limit <- N_star - 1
  f1_cumul  <- pnbinom(sum_limit, size = size_f1, prob = prob_f1)
  f2_cumul  <- pnbinom(sum_limit, size = size_f2, prob = prob_f2)
  f_star1   <- f1 / (1 - f1_cumul)
  f_star2   <- f2 / (1 - f2_cumul)

  P <- theta[5]
  logL <- wi * log((P * f_star1) + ((1 - P) * f_star2))
  sum(-logL)
}
