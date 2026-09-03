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
#' @param transformed A logical scalar indicating whether or not \eqn{theta}
#'   is being optimized on the transformed scale (natural log for \eqn{\alpha})s
#'   and \eqn{\beta}s; logit scale for the mixture fraction). The transformed
#'   scale can be especially useful if using unconstrained optimization.
#' @returns A scalar negative log-likelihood value.
#'
#' @details The marginal distribution of the counts, \emph{N}, is a mixture of
#'   two negative binomial distributions. The hyperparameters for the prior
#'   distribution (mixture of gammas) are estimated by optimizing the likelihood
#'   equation from this marginal distribution.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item \eqn{\alpha_1, \beta_1}: Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)
#'     \item \eqn{\alpha_2, \beta_2}: Parameters of the second component
#'     \item \eqn{P}: Mixture fraction
#'   }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N} actually contains zeroes before using this function. You
#'   should have used the \code{zeroes = TRUE} option when calling the
#'   \code{\link{processRaw}} function.
#' @section Warnings:
#'   Make sure the data were not squashed before using this function.
#' @section Warnings:
#'   Make sure there are no missing or infinite counts (actual or expected)
#'   before using this function.
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family negative log-likelihood functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlm}}, \code{\link[stats]{nlminb}}, and
#'   \code{\link[stats]{optim}} for optimization
#' @importFrom stats dnbinom
#' @importFrom stats plogis
#' @export
negLLzero <- function(theta, N, E, transformed = FALSE) {
  .checkInputs_negLLzero(theta, N, E)
  if (transformed) {
    theta[1:4] <- exp(theta[1:4])
    theta[5] <- plogis(theta[5])  #inverse logit
  } else {
    if (any(theta == 0) || theta[5] == 1) {
      return(NA)
    }
  }
  L <- .negLLzero_Lvec(theta, N, E)
  L_is_infinite <- is.infinite(L)
  if (all(L_is_infinite)) {
    neg_LL <- 1e+10
    neg_LL <- neg_LL + 100 * stats::rbinom(1, size = 1, prob = 0.5)  #prevent false convergence
    return(neg_LL)
  }
  L[L < 1e-100] <- 1e-100
  L[L_is_infinite & L > 0] <- 1e+100
  return(sum(-log(L)))
}

.negLLzero_Lvec <- function(theta, N, E) {
  # Returns a vector for L calculations
  log_f1 <- .negLLzero_logfstar(theta[1:2], N, E)
  log_f2 <- .negLLzero_logfstar(theta[3:4], N, E)
  log_P <- log(theta[5])
  log_1_minus_P <- log(1 - theta[5])
  # (P * f1) + ((1 - P) * f2)
  exp(log_P + log_f1) + exp(log_1_minus_P + log_f2)
}

.negLLzero_logfstar <- function(alpha_beta, N, E) {
  # Returns log(f) for one of the mixture components
  # alpha_beta: A vector of length 2
  size_f <- alpha_beta[1]  #alpha
  prob_f <- exp(log(alpha_beta[2]) - log(alpha_beta[2] + E))  #beta / (beta + E)
  prob_f[prob_f < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f[prob_f > max_prob] <- max_prob
  dnbinom(N, size = size_f, prob = prob_f, log = TRUE)
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
#' @returns A scalar negative log-likelihood value.
#'
#' @details The marginal distribution of the counts, \emph{N}, is a mixture of
#'   two negative binomial distributions. The hyperparameters for the prior
#'   distribution (mixture of gammas) are estimated by optimizing the likelihood
#'   equation from this marginal distribution.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item \eqn{\alpha_1, \beta_1}: Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)
#'     \item \eqn{\alpha_2, \beta_2}: Parameters of the second component
#'     \item \eqn{P}: Mixture fraction
#'   }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{ni} actually contains zeroes before using this function.
#'   You should have used the \code{zeroes = TRUE} option when calling the
#'   \code{processRaw} function.
#' @section Warnings:
#'   Make sure the data were actually squashed (see \code{\link{squashData}})
#'   before using this function.
#' @section Warnings:
#'   Make sure there are no missing or infinite counts (actual or expected)
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
#' @importFrom stats plogis
#' @export
negLLzeroSquash <- function(theta, ni, ei, wi, transformed = FALSE) {
  .checkInputs_negLLzeroSquash(theta, ni, ei, wi)
  if (transformed) {
    theta[1:4] <- exp(theta[1:4])
    theta[5] <- plogis(theta[5])  #inverse logit
  } else {
    if (any(theta == 0) || theta[5] == 1) {
      return(NA)
    }
  }
  logL <- .negLLzeroSquash_logLvec(theta, ni, ei, wi) 
  logL_is_infinite <- is.infinite(logL)
  if (all(logL_is_infinite)) {
    neg_LL <- 1e+10
    neg_LL <- neg_LL + 100 * stats::rbinom(1, size = 1, prob = 0.5)  #prevent false convergence
    return(neg_LL)
  }
  logL[logL_is_infinite & logL < 0] <- -230  #log(1e-100)
  logL[logL_is_infinite & logL > 0] <- 230  #log(1e+100)
  return(sum(-logL))
}

.negLLzeroSquash_logLvec <- function(theta, ni, ei, wi) {
  # Returns a vector for logL calculations
  log_f1 <- .negLLzeroSquash_logfstar(theta[1:2], ni, ei)
  log_f2 <- .negLLzeroSquash_logfstar(theta[3:4], ni, ei)
  log_P <- log(theta[5])
  log_1_minus_P <- log(1 - theta[5])
  # wi * log[(P * f1) + ((1 - P) * f2)]
  wi * log(exp(log_P + log_f1) + exp(log_1_minus_P + log_f2))
}

.negLLzeroSquash_logfstar <- function(alpha_beta, ni, ei) {
  # Returns log(f) for one of the mixture components
  # alpha_beta: A vector of length 2
  size_f <- alpha_beta[1]  #alpha
  prob_f <- exp(log(alpha_beta[2]) - log(alpha_beta[2] + ei))  #beta / (beta + E)
  prob_f[prob_f < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f[prob_f > max_prob] <- max_prob
  dnbinom(ni, size = size_f, prob = prob_f, log = TRUE)
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
#' @returns A scalar negative log-likelihood value
#'
#' @details The conditional marginal distribution for the counts, \emph{N},
#'   given that \emph{N >= N*}, is based on a mixture of two negative binomial
#'   distributions. The hyperparameters for the prior distribution (mixture of
#'   gammas) are estimated by optimizing the likelihood equation from this
#'   conditional marginal distribution. It is recommended to use \code{N_star =
#'   1} when practical.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item \eqn{\alpha_1, \beta_1}: Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)
#'     \item \eqn{\alpha_2, \beta_2}: Parameters of the second component
#'     \item \eqn{P}: Mixture fraction
#'   }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N_star} matches the smallest actual count in \emph{N}
#'   before using this function. Filter \emph{N} and \emph{E} if needed.
#' @section Warnings:
#'   Make sure the data were not squashed before using this function.
#' @section Warnings:
#'   Make sure there are no missing or infinite counts (actual or expected)
#'   before using this function.
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
#' @importFrom stats plogis
#' @export
negLL <- function(theta, N, E, N_star = 1, transformed = FALSE) {
  .checkInputs_negLL(theta, N, E, N_star)
  if (transformed) {
    theta[1:4] <- exp(theta[1:4])
    theta[5] <- plogis(theta[5])  #inverse logit
  } else {
    if (any(theta == 0) || theta[5] == 1) {
      return(NA)
    }
  }
  L <- .negLL_Lvec(theta, N, E, N_star)
  L_is_infinite <- is.infinite(L)
  if (all(L_is_infinite)) {
    neg_LL <- 1e+10
    neg_LL <- neg_LL + 100 * stats::rbinom(1, size = 1, prob = 0.5)  #prevent false convergence
    return(neg_LL)
  }
  L[L < 1e-100] <- 1e-100
  L[L_is_infinite & L > 0] <- 1e+100
  return(sum(-log(L)))
}

.negLL_Lvec <- function(theta, N, E, N_star) {
  # Returns a vector for L calculations
  sum_limit <- N_star - 1
  log_f_star1 <- .negLL_logfstar(theta[1:2], N, E, sum_limit)
  log_f_star2 <- .negLL_logfstar(theta[3:4], N, E, sum_limit)
  log_P <- log(theta[5])
  log_1_minus_P <- log(1 - theta[5])
  # (P * f_star1) + ((1 - P) * f_star2)
  exp(log_P + log_f_star1) + exp(log_1_minus_P + log_f_star2)
}

.negLL_logfstar <- function(alpha_beta, N, E, sum_limit) {
  # Returns log(f_star) for one of the mixture components
  # alpha_beta: A vector of length 2
  size_f <- alpha_beta[1]  #alpha
  prob_f <- exp(log(alpha_beta[2]) - log(alpha_beta[2] + E))  #beta / (beta + E)
  prob_f[prob_f < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f[prob_f > max_prob] <- max_prob
  log_f <- dnbinom(N, size = size_f, prob = prob_f, log = TRUE)
  log_comp_f_cumul <- pnbinom(sum_limit, size = size_f, prob = prob_f,
    lower.tail = FALSE, log.p = TRUE
  )
  log_f - log_comp_f_cumul  #log(f / (1-f_cumul))
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
#' @returns A scalar negative log-likelihood value
#'
#' @details The conditional marginal distribution for the counts, \emph{N},
#'   given that \emph{N >= N*}, is based on a mixture of two negative binomial
#'   distributions. The hyperparameters for the prior distribution (mixture of
#'   gammas) are estimated by optimizing the likelihood equation from this
#'   conditional marginal distribution. It is recommended to use \code{N_star =
#'   1} when practical.
#' @details The hyperparameters are:
#'   \itemize{
#'     \item \eqn{\alpha_1, \beta_1}: Parameters of the first component of the
#'       marginal distribution of the counts (also the prior distribution)
#'     \item \eqn{\alpha_2, \beta_2}: Parameters of the second component
#'     \item \eqn{P}: Mixture fraction
#'   }
#' @details This function will not need to be called directly if using
#'   \code{\link{exploreHypers}} or \code{\link{autoHyper}}.
#' @section Warnings:
#'   Make sure \emph{N_star} matches the smallest actual count in \emph{ni}
#'   before using this function. Filter \emph{ni}, \emph{ei}, and \emph{wi} if
#'   needed.
#' @section Warnings:
#'   Make sure the data were actually squashed (see \code{\link{squashData}})
#'   before using this function.
#' @section Warnings:
#'   Make sure there are no missing or infinite counts (actual or expected)
#'   before using this function.
#'
#' @examples
#' data.table::setDTthreads(2)  #only needed for CRAN checks
#' theta_init <- c(1, 1, 3, 3, .2)  #initial guess
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 10)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 10)
#' negLLsquash(theta = theta_init, ni = squashed$N, ei = squashed$E,
#'   wi = squashed$weight
#' )
#' # For hyperparameter estimation...
#' # Log & logit transformations can help with unconstrained optimization techniques
#' theta_init_transf <- theta_init
#' theta_init_transf[1:4] <- log(theta_init_transf[1:4])
#' theta_init_transf[5] <- log(theta_init_transf[5] / (1 - theta_init_transf[5]))
#' theta_transf <- stats::nlminb(start = theta_init_transf, objective = negLLsquash,
#'   ni = squashed$N, ei = squashed$E, wi = squashed$weight,
#'   transformed = TRUE
#' )
#' theta_transf
#' theta <- theta_transf$par
#' # Backtransform
#' theta[1:4] <- exp(theta[1:4])
#' theta[5] <- stats::plogis(theta[5])
#' theta
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
#' @importFrom stats plogis
#' @export
negLLsquash <- function(theta, ni, ei, wi, N_star = 1, transformed = FALSE) {
  .checkInputs_negLLsquash(theta, ni, ei, wi, N_star)
  if (transformed) {
    theta[1:4] <- exp(theta[1:4])
    theta[5] <- plogis(theta[5])  #inverse logit
  } else {
    if (any(theta == 0) || theta[5] == 1) {
      return(NA)
    }
  }
  logL <- .negLLsquash_logLvec(theta, ni, ei, wi, N_star)
  logL_is_infinite <- is.infinite(logL)
  if (all(logL_is_infinite)) {
    neg_LL <- 1e+10
    neg_LL <- neg_LL + 100 * stats::rbinom(1, size = 1, prob = 0.5)  #prevent false convergence
    return(neg_LL)
  }
  logL[logL_is_infinite & logL < 0] <- -230  #log(1e-100)
  logL[logL_is_infinite & logL > 0] <- 230  #log(1e+100)
  return(sum(-logL))
}

.negLLsquash_logLvec <- function(theta, ni, ei, wi, N_star) {
  # Return a vector for logL calculations
  sum_limit <- N_star - 1
  log_f_star1 <- .negLLsquash_logfstar(theta[1:2], ni, ei, sum_limit)
  log_f_star2 <- .negLLsquash_logfstar(theta[3:4], ni, ei, sum_limit)
  log_P <- log(theta[5])
  log_1_minus_P <- log(1 - theta[5])
  # wi * log((P * f_star1) + ((1 - P) * f_star2))
  wi * log(exp(log_P + log_f_star1) + exp(log_1_minus_P + log_f_star2))
}

.negLLsquash_logfstar <- function(alpha_beta, ni, ei, sum_limit) {
  # Returns log(f_star) for one of the mixture components
  # alpha_beta: A vector of length 2
  size_f <- alpha_beta[1]  #alpha
  prob_f <- exp(log(alpha_beta[2]) - log(alpha_beta[2] + ei))  #beta / (beta + E)
  prob_f[prob_f < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f[prob_f > max_prob] <- max_prob
  log_f <- dnbinom(ni, size = size_f, prob = prob_f, log = TRUE)
  log_comp_f_cumul <- pnbinom(sum_limit, size = size_f, prob = prob_f,
    lower.tail = FALSE, log.p = TRUE
  )
  log_f - log_comp_f_cumul  #log(f / (1-f_cumul))
}
