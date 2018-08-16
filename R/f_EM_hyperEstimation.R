#' Estimate hyperparameters using an EM algorithm
#'
#' \code{hyperEM} finds hyperparameter estimates using a variation on the
#' Expectation-Maximization (EM) algorithm known as the Expectation/Conditional
#' Maximization (ECM) algorithm (Meng et al, 1993). The algorithm estimates each
#' element of the hyperparameter vector, \eqn{\theta}, while holding fixed
#' (conditioning on) the other parameters in the vector. Alternatively, it can
#' estimate both parameters for each distribution in the mixture while holding
#' the parameters from the other distribution and the mixing fraction fixed.
#'
#' @param data A data frame from \code{\link{processRaw}} or
#'   \code{\link{squashData}} containing columns named \emph{N}, \emph{E}, and
#'   (if squashed) \emph{weight}.
#' @param theta_init_vec A numeric vector of initial hyperparameter guesses
#'   ordered as: \eqn{\alpha_1, \beta_1, \alpha_2, \beta_2, P}.
#' @param squashed A scalar logical (\code{TRUE} or \code{FALSE}) indicating
#'   whether or not data squashing was used.
#' @param zeroes A scalar logical specifying if zero counts are included.
#' @param N_star A positive scalar whole number value for the minimum count
#'   size to be used for hyperparameter estimation. If zeroes are used, set
#'   \code{N_star} to \code{NULL}.
#' @param method A scalar string indicating which method (i.e. score functions
#'   or log-likelihood function) to use for the maximization steps. Possible
#'   values are \code{"score"} and \code{"nlminb"}.
#' @param profile A scalar string indicating which method to use to optimize the
#'   log-likelihood function if \code{method = "nlminb"} (ignored if
#'   \code{method = "score"}).  \code{profile = "parameter"} optimizes one
#'   parameter (\eqn{\alpha} or \eqn{\beta}) from the log-likelihood function at
#'   a time.  \code{profile = "distribution"} optimizes one distribution from
#'   the mixture at a time (\eqn{\alpha} and \eqn{\beta} simultaneously).
#' @param LL_tol A scalar numeric value for the tolerance used for determining
#'   when the change in log-likelihood at each iteration is sufficiently small.
#'   Used for convergence assessment.
#' @param consecutive A positive scalar whole number value for the number of
#'   consecutive iterations the change in log-likelihood must be below
#'   \code{LL_tol} in order to reach convergence. Larger values reduce the
#'   chance of getting stuck in a flat region of the curve.
#' @param param_lower A scalar numeric value for the smallest acceptable value
#'   for each \eqn{\alpha} and \eqn{\beta} estimate.
#' @param param_upper A scalar numeric value for the largest acceptable value
#'   for each \eqn{\alpha} and \eqn{\beta} estimate.
#' @param print_level A value that determines how much information is printed
#'   during execution. Possible value are \code{0} for no printing, \code{1} for
#'   minimal information, and \code{2} for maximal information.
#' @param max_iter A positive scalar whole number value for the maximum number
#'   of iterations to use.
#' @param conf_ints A scalar logical indicating if confidence intervals and
#'   standard errors should be returned.
#' @param conf_level A scalar string for the confidence level used if confidence
#'   intervals are requested.
#' @param track A scalar logical indicating whether or not to retain the
#'   hyperparameter estimates and log-likelihood value at each iteration. Can be
#'   used for plotting to better understand the algorithm's behavior.
#'
#' @return A list including the following:
#'   \itemize{
#'     \item{\emph{estimates}: }{The maximum likelihood estimate (MLE) of the
#'       hyperparameter, \eqn{\theta}.}
#'     \item{\emph{conf_int}: }{A data frame including the standard errors and
#'       confidence limits for \code{estimates}. Only included if
#'       \code{conf_ints = TRUE}.}
#'     \item{\emph{maximum}: }{The log-likelihood function evaluated at
#'       \code{estimates}.}
#'     \item{\emph{method}: }{The method used in the maximization step.}
#'     \item{\emph{elapsed}: }{The elapsed function execution time in seconds.}
#'     \item{\emph{iters}: }{The number of iterations used.}
#'     \item{\emph{score}: }{The score functions (i.e. score vector) evaluated
#'       at \code{estimates}. All elements should be close to zero.}
#'     \item{\emph{score_norm}: }{The Euclidean norm of the score vector
#'       evaluated at \code{estimates}. Should be close to zero.}
#'     \item{\emph{tracking}: }{The estimates of \eqn{\theta} at each iteration
#'       and the log-likelihood function evaluated at those estimates. Unless
#'       \code{track = TRUE}, only shows the starting point of the algorithm.}
#'   }
#' @details If \code{method = "score"}, the maximization step finds a root
#'   of the score function. If \code{method = "nlminb"},
#'   \code{\link[stats]{nlminb}} is used to find a minimum of the negative
#'   log-likelihood function.
#' @details If \code{method = "score"} and \code{zeroes = FALSE}, then
#'   \code{'N_star'} must equal \code{1}.
#' @details If \code{method = "score"}, the model is reparameterized. The
#'   parameters are transformed to force the parameter space to include all real
#'   numbers. This approach addresses numerical issues at the edge of the
#'   parameter space. The reparameterization follows:
#'   \eqn{\alpha_{prime} = log(\alpha)}, \eqn{\beta_{prime} = log(\beta)}, and
#'   \eqn{P_{prime} = tan(pi * P - pi / 2)}. However, the values returned in
#'   \code{estimates} are on the original scale (back-transformed).
#' @details On every 100th iteration, the procedure described in Millar (2011)
#'   is used to accelerate the estimate of \eqn{\theta}.
#' @details The score vector and its Euclidean norm should be close to zero at
#'   a local maximum and can therefore be used to help assess the reliability of
#'   the results. A local maximum might not be the global MLE, however.
#' @details Asymptotic normal confidence intervals, if requested, use standard
#'   errors calculated from the observed Fisher information matrix as discussed
#'   in DuMouchel (1999).
#'
#' @examples
#' data(caers)
#' proc <- processRaw(caers)
#' squashed <- squashData(proc, bin_size = 100, keep_pts = 0)
#' squashed <- squashData(squashed, count = 2, bin_size = 12, keep_pts = 24)
#' hyperEM(squashed, theta_init_vec = c(1, 1, 2, 2, .3), consecutive = 10)
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @references Meng X-L, Rubin D (1993). "Maximum likelihood estimation via the
#'   ECM algorithm: A general framework", \emph{Biometrika}, 80(2), 267-278.
#' @references Millar, Russell B (2011). "Maximum Likelihood Estimation and
#'   Inference", \emph{John Wiley & Sons, Ltd}, 111-112.
#' @family hyperparameter estimation functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{uniroot}} for finding a zero of the score
#'   function and \code{\link[stats]{nlminb}} for minimizing the negative
#'   log-likelihood function
#' @importFrom stats uniroot
#' @importFrom stats nlminb
#' @importFrom stats optim
#' @importFrom stats optimHess
#' @export

hyperEM <-
  function(data, theta_init_vec, squashed = TRUE, zeroes = FALSE, N_star = 1,
           method = c("score", "nlminb"),
           profile = c("parameter", "distribution"), LL_tol = 1e-04,
           consecutive = 100, param_lower = 1e-05, param_upper = 20,
           print_level = 2, max_iter = 5000, conf_ints = FALSE,
           conf_level = c("95", "80", "90", "99"), track = FALSE) {

  method  <- match.arg(method)
  profile <- match.arg(profile)
  conf_level <- match.arg(conf_level)
  .checkInputs_hyperEM(data, theta_init_vec, squashed, zeroes, N_star, method,
                       LL_tol, consecutive, param_lower, param_upper,
                       print_level, max_iter, conf_ints, track)

  #A lot of initialization
  start_time <- proc.time()
  N <- data$N; E <- data$E; W <- NULL
  if (squashed) W <- data$weight
  track_LL <- .deltaLL(theta_init_vec, theta_init_vec, N, E, W,
                       squashed, zeroes, N_star)$LL
  track_a1 <- theta_init_vec[1]; track_b1 <- theta_init_vec[2]
  track_a2 <- theta_init_vec[3]; track_b2 <- theta_init_vec[4]
  track_P  <- theta_init_vec[5]
  theta       <- theta_init_vec
  count_iter  <- count_converge <- count_iter_acc <- count_stuck_all <- 0L
  count_stuck <- rep.int(0L, 5); conf_int <- NULL

  #Responsibilities (E-step)
  if (method == "score") {
    marg_dens <- .marginalDensity(theta, N, E, zeroes)
    PQ        <- .updateProbs(theta, marg_dens$f1, marg_dens$f2)
  }

  while(count_converge <= consecutive &&
        count_stuck_all < 20L &&
        count_iter < max_iter) {

    count_iter     <- count_iter + 1
    count_iter_acc <- count_iter_acc + 1
    old_theta      <- theta

    #M-step
    if (method == "score") {
      theta <- .updateTheta(old_theta, PQ$P_vec, PQ$Q_vec, N, E, W, squashed,
                            zeroes, param_lower, param_upper)
      count_stuck <- ifelse(theta == old_theta, count_stuck + 1, 0L)
      if (any(count_stuck > 10L)) {
        theta <- .nudgeTheta(theta, count_stuck, 10L, param_lower)
      }
    } else {
      if (profile == "parameter") {
        theta <- .updateThetaLL(old_theta, N, E, W, squashed, zeroes, N_star,
                                param_lower, param_upper)
      } else {
        theta <- .updateThetaLLD(old_theta, N, E, W, squashed, zeroes, N_star,
                                 param_lower, param_upper)
      }
    }

    #Check if all hyper estimates are "stuck"
    if (all(theta == old_theta)) {
      count_stuck_all <- count_stuck_all + 1
    } else {
      count_stuck_all <- 0L
    }

    #Accelerate estimation every 100 iterations
    if (count_iter_acc == 98L) theta_98 <- theta
    if (count_iter_acc == 99L) theta_99 <- theta
    if (count_iter_acc == 100L) {
      theta_100 <- theta
      theta <- .accelerateTheta(theta_98, theta_99, theta_100,
                                param_lower, param_upper)
      count_iter_acc <- 0L
    }

    #E-step
    if (method == "score") {
      marg_dens <- .marginalDensity(theta, N, E, zeroes)
      PQ        <- .updateProbs(theta, marg_dens$f1, marg_dens$f2)
    }

    #Convergence check
    delta_LL <- .deltaLL(theta, old_theta, N, E, W, squashed, zeroes, N_star)
    .hyperEMmessage("A", print_level, count_iter, delta_LL$delta, theta)
    if (delta_LL$delta < LL_tol) {
      count_converge <- count_converge + 1
    } else {
      count_converge <- 0L
    }

    #Track hyperparameter estimates for plotting
    if (track) {
      tracking_position <- count_iter + 1
      track_LL[tracking_position] <- delta_LL$LL
      track_a1[tracking_position] <- theta[1]
      track_b1[tracking_position] <- theta[2]
      track_a2[tracking_position] <- theta[3]
      track_b2[tracking_position] <- theta[4]
      track_P[tracking_position]  <- theta[5]
    }
  }

  if (count_iter >= max_iter) stop("exceeded maximum number of iterations")
  if (count_stuck_all >= 20) stop("estimate seems to be stuck")

  #Finish up
  elapsed <- proc.time() - start_time
  if (conf_ints) {
    conf_int <- .hyperConfInts(data, theta, zeroes, squashed, N_star,
                               conf_level)
  }
  score <- .checkScore(theta, N, E, W, squashed, zeroes)
  .hyperEMmessage("B", print_level, count_iter, elapsed = elapsed)
  track_iters <- 0:(length(track_LL) - 1)
  tracking <- data.frame(iter = track_iters, logL = track_LL, alpha1 = track_a1,
                         beta1 = track_b1, alpha2 = track_a2, beta2 = track_b2,
                         P = track_P)

  list(estimates = theta, conf_int = conf_int, maximum = delta_LL$LL,
       method = method, elapsed = elapsed[[3]], iters = count_iter,
       score = score$score_vec, score_norm = score$score_norm,
       tracking = tracking)
}
