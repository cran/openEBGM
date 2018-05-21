#' Calculate Qn
#'
#' \code{Qn} calculates \eqn{Q_n}, the posterior probability that \eqn{\lambda}
#' came from the first component of the mixture, given \emph{N = n} (Eq. 6,
#' DuMouchel 1999). \eqn{Q_n} is the mixture fraction for the posterior
#' distribution.
#'
#' @param theta_hat A numeric vector of hyperparameter estimates (likely from
#'   \code{\link{autoHyper}} or from directly minimizing
#'   \code{\link{negLLsquash}}) ordered as:
#'   \eqn{\alpha_1, \beta_1, \alpha_2, \beta_2, P}.
#' @param N A whole number vector of actual counts from
#'   \code{\link{processRaw}}.
#' @param E A numeric vector of expected counts from \code{\link{processRaw}}.
#' @return A numeric vector of probabilities.
#'
#' @details The hyperparameter estimates (\code{theta_hat}) are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameter estimates of the first
#'       component of the prior distribution}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameter estimates of the second
#'       component}
#'     \item{\eqn{P}: }{Mixture fraction estimate of the prior distribution}
#' }
#'
#' @examples
#' theta_init <- data.frame(
#'   alpha1 = c(0.2, 0.1),
#'   beta1  = c(0.1, 0.1),
#'   alpha2 = c(2,   10),
#'   beta2  = c(4,   10),
#'   p      = c(1/3, 0.2)
#' )
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' suppressWarnings(
#'   theta_hat <- autoHyper(data = squashed, theta_init = theta_init)$estimates
#' )
#' qn <- Qn(theta_hat, N = proc$N, E = proc$E)
#' head(qn)
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family posterior distribution functions
#' @keywords openEBGM
#' @seealso \code{\link{autoHyper}},  \code{\link{exploreHypers}},
#'   \code{\link{negLLsquash}}, \code{\link{negLL}},
#'   \code{\link{negLLzero}}, and \code{\link{negLLzeroSquash}} for
#'   hyperparameter estimation.
#' @seealso \code{\link{processRaw}} for finding counts.
#' @importFrom stats dnbinom
#' @export
Qn <- function(theta_hat, N, E) {

  .checkInputs_Qn(theta_hat, N, E)

  prob_f1 <- theta_hat[2] / (theta_hat[2] + E)  #beta1 / (beta1 + E)
  prob_f2 <- theta_hat[4] / (theta_hat[4] + E)  #beta2 / (beta2 + E)
  f1_NB   <- dnbinom(N, size = theta_hat[1], prob = prob_f1)  #Eq.5 (1999)
  f2_NB   <- dnbinom(N, size = theta_hat[3], prob = prob_f2)

  P_hat  <- theta_hat[5]
  Qn_num <- P_hat * f1_NB
  Qn_den <- (P_hat * f1_NB) + ((1 - P_hat) * f2_NB)

  Qn_num / Qn_den  #Eq. 6 (1999)
}


#' Calculate EBGM scores
#'
#' \code{ebgm} calculates the Empirical Bayes Geometric Mean (\emph{EBGM}),
#' which is \sQuote{the geometric mean of the empirical Bayes posterior
#' distribution of the \dQuote{true} \emph{RR}} (DuMouchel 1999, see Eq.11). The
#' \emph{EBGM} is essentially a version of the relative reporting ratio
#' (\emph{RR}) that uses Bayesian shrinkage.
#'
#' @inheritParams Qn
#' @param qn A numeric vector of posterior probabilities that \eqn{\lambda} came
#'   from the first component of the mixture, given \emph{N = n} (i.e., the
#'   mixture fraction). See function \code{\link{Qn}}.
#' @param digits A scalar whole number that determines the number of decimal
#'   places used when rounding the results.
#' @return A numeric vector of EBGM scores.
#'
#' @details The hyperparameter estimates (\code{theta_hat}) are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameter estimates of the first
#'       component of the prior distribution}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameter estimates of the second
#'       component}
#'     \item{\eqn{P}: }{Mixture fraction estimate of the prior distribution}
#' }
#'
#' @examples
#' theta_init <- data.frame(
#'   alpha1 = c(0.2, 0.1),
#'   beta1  = c(0.1, 0.1),
#'   alpha2 = c(2,   10),
#'   beta2  = c(4,   10),
#'   p      = c(1/3, 0.2)
#' )
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' suppressWarnings(
#'   theta_hat <- autoHyper(data = squashed, theta_init = theta_init)$estimates
#' )
#' qn <- Qn(theta_hat, N = proc$N, E = proc$E)
#' proc$EBGM <- ebgm(theta_hat, N = proc$N, E = proc$E, qn = qn)
#' head(proc)
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family posterior distribution functions
#' @keywords openEBGM
#' @seealso \code{\link{autoHyper}},  \code{\link{exploreHypers}},
#'   \code{\link{negLLsquash}}, \code{\link{negLL}},
#'   \code{\link{negLLzero}}, and \code{\link{negLLzeroSquash}} for
#'   hyperparameter estimation.
#' @seealso \code{\link{processRaw}} for finding counts.
#' @seealso \code{\link{Qn}} for finding mixture fractions.
#' @export
ebgm <- function(theta_hat, N, E, qn, digits = 2) {

  .checkInputs_ebgm(theta_hat, N, E, qn, digits)

  #If X ~ gamma(alpha, beta), then E[X] = alpha / beta
  #Also, E[ln(X)] = digamma(alpha) - ln(beta)
  expectation1 <- digamma(theta_hat[1] + N) - log(theta_hat[2] + E)  #Eq.9
  expectation2 <- digamma(theta_hat[3] + N) - log(theta_hat[4] + E)

  exp_log_lamda <- (qn * expectation1) + ((1 - qn) * expectation2)  #Eq.9 (1999)
  EBlog2        <- exp_log_lamda / log(2)  #Eq.10 (1999)

  EBGM <- 2 ^ EBlog2  #Eq.11 (1999)
  round(EBGM, digits)
}


#' Find quantiles of the posterior distribution
#'
#' \code{quantBisect} finds the desired quantile of the posterior distribution
#' using the bisection method. Used to create credibility limits.
#'
#' @inheritParams Qn
#' @inheritParams ebgm
#' @param percent A numeric scalar between 1 and 99 for the desired
#'   percentile (e.g., 5 for 5th percentile).
#' @param limits A whole number vector of length 2 for the upper and lower
#'   bounds of the search space.
#' @param max_iter A whole number scalar for the maximum number of iterations.
#'   Used to prevent infinite loops.
#' @return A numeric vector of quantile estimates.
#'
#' @details The hyperparameter estimates (\code{theta_hat}) are:
#'   \itemize{
#'     \item{\eqn{\alpha_1, \beta_1}: }{Parameter estimates of the first
#'       component of the prior distribution}
#'     \item{\eqn{\alpha_2, \beta_2}: }{Parameter estimates of the second
#'       component}
#'     \item{\eqn{P}: }{Mixture fraction estimate of the prior distribution}
#' }
#' @details Although this function can find any quantile of the posterior
#'   distribution, it will often be used to calculate the 5th and 95th
#'   percentiles to create a 90\% credibility interval.
#' @details The quantile is calculated by solving for \eqn{x} in the general
#'   equation \eqn{F(x) = cutoff}, or equivalently, \eqn{F(x) - cutoff = 0},
#'   where \eqn{F(x)} is the cumulative distribution function of the posterior
#'   distribution and \eqn{cutoff} is the appropriate cutoff level (e.g., 0.05
#'   for the 5th percentile).
#' @section Warning:
#'   The \code{digits} argument determines the tolerance for the bisection
#'   algorithm. The more decimal places you want returned, the longer the run
#'   time.
#'
#' @examples
#' theta_init <- data.frame(
#'   alpha1 = c(0.2, 0.1),
#'   beta1  = c(0.1, 0.1),
#'   alpha2 = c(2,   10),
#'   beta2  = c(4,   10),
#'   p      = c(1/3, 0.2)
#' )
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 100)
#' squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 20)
#' suppressWarnings(
#'   theta_hat <- autoHyper(data = squashed, theta_init = theta_init)$estimates
#' )
#' qn <- Qn(theta_hat, N = proc$N, E = proc$E)
#' proc$QUANT_05 <- quantBisect(percent = 5, theta = theta_hat, N = proc$N,
#'                              E = proc$E, qn = qn)
#' proc$QUANT_95 <- quantBisect(percent = 95, theta = theta_hat, N = proc$N,
#'                              E = proc$E, qn = qn)
#' head(proc)
#'
#' @family posterior distribution functions
#' @seealso \url{https://en.wikipedia.org/wiki/Bisection_method}
#' @keywords openEBGM
#' @seealso \code{\link{autoHyper}},  \code{\link{exploreHypers}},
#'   \code{\link{negLLsquash}}, \code{\link{negLL}},
#'   \code{\link{negLLzero}}, and \code{\link{negLLzeroSquash}} for
#'   hyperparameter estimation.
#' @seealso \code{\link{processRaw}} for finding counts.
#' @seealso \code{\link{Qn}} for finding mixture fractions.
#' @importFrom stats pgamma
#' @export
quantBisect <- function(percent, theta_hat, N, E, qn, digits = 2,
                        limits = c(-100000, 100000), max_iter = 2000) {

  .checkInputs_quantBisect(percent, theta_hat, N, E, qn, digits, limits,
                           max_iter)

  cdf_error <- function(guess, cut_point, qn, theta_hat, N, E) {
    #Difference between F(x) & the "cut point" (e.g. 0.05 for 5th percentile)
    cdf1 <- pgamma(guess, shape = theta_hat[1] + N, rate = theta_hat[2] + E)
    cdf2 <- pgamma(guess, shape = theta_hat[3] + N, rate = theta_hat[4] + E)
    post_cdf  <- (qn * cdf1) + ((1 - qn) * cdf2)
    post_cdf - cut_point
  }

  tol        <- 0.5 * (10 ^ (-digits))
  cut_point  <- percent / 100
  num_pts    <- length(N)
  guess_init <- rep(1, num_pts)  #null hypothesis: RR = 1
  error_init <- cdf_error(guess_init, cut_point, qn, theta_hat, N, E)
  is_too_big <- error_init > 0
  left       <- rep(limits[1], num_pts)
  right      <- rep(limits[2], num_pts)
  left       <- ifelse(is_too_big, left, guess_init)
  right      <- ifelse(is_too_big, guess_init, right)

  iter_count <- 1L
  while (iter_count <= max_iter) {
      mid_pt <- (left + right) / 2
      error  <- cdf_error(mid_pt, cut_point, qn, theta_hat, N, E)
      if (max((right - left) / 2) < tol) {
          quantiles <- round(mid_pt, digits)
          if (max(quantiles) == limits[2]) {
            stop("increase maximum for 'limits'")
          }
          return(quantiles)
      } else {
          iter_count <- iter_count + 1
          error_left <- cdf_error(left, cut_point, qn, theta_hat, N, E)
          left  <- ifelse(sign(error_left) == sign(error), mid_pt, left)
          right <- ifelse(sign(error_left) == sign(error), right, mid_pt)
      }
  }
  stop("failed to converge -- try adjusting 'limits' or 'max_iter'")
}
