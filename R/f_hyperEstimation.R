#' Explore various hyperparameter estimates
#'
#' \code{exploreHypers} finds hyperparameter estimates using a variety of
#' starting points to examine the consistency of the optimization procedure.
#'
#' @param data A data frame from \code{\link{processRaw}} containing columns
#'   named \emph{N}, \emph{E}, and (if squashed) \emph{weight}.
#' @param theta_init A data frame of initial hyperparameter guesses with
#'   columns ordered as:
#'   \eqn{\alpha_1, \beta_1, \alpha_2, \beta_2, P}.
#' @param squashed A scalar logical (\code{TRUE} or \code{FALSE}) indicating
#'   whether or not data squashing was used.
#' @param zeroes A scalar logical specifying if zero counts are included.
#' @param N_star A positive scalar whole number value for the minimum count
#'   size to be used for hyperparameter estimation. If zeroes are used, set
#'   \code{N_star} to \code{NULL}.
#' @param method A scalar string indicating which optimization procedure is to
#'   be used. Choices are \code{"nlminb"}, \code{"nlm"}, or \code{"bfgs"}.
#' @param param_limit A scalar numeric value for the largest acceptable value
#'   for the \eqn{\alpha} and \eqn{\beta} estimates. Used to help protect
#'   against unreasonable/erroneous estimates.
#' @param max_pts A scalar whole number for the largest number of data points
#'   allowed. Used to help prevent extremely long run times.
#' @param std_errors A scalar logical indicating if standard errors should be
#'   returned for the hyperparameter estimates.
#'
#' @return A list including the data frame \code{estimates} of hyperparameter
#'   estimates corresponding to the initial guesses from \code{theta_init} (plus
#'   convergence results):
#'   \itemize{
#'     \item{\emph{code}: }{The convergence code returned by the chosen
#'       optimization function (see \code{\link[stats]{nlminb}},
#'       \code{\link[stats]{nlm}}, and \code{\link[stats]{optim}} for details).}
#'     \item{\emph{converge}: }{A logical indicating whether or not convergence
#'       was reached. See "Details" section for more information.}
#'     \item{\emph{in_bounds}: }{A logical indicating whether or not the
#'       estimates were within the bounds of the parameter space (upper bound
#'       for \eqn{\alpha_1, \beta_1, \alpha_2, and \beta_2} was determined by
#'       the \code{param_limit} argument).}
#'     \item{\emph{minimum}: }{The negative log-likelihood value corresponding
#'       to the estimated optimal value of the hyperparameter.}
#'   }
#'   Also returns the data frame \code{std_errs} if standard errors are
#'   requested.
#'
#' @section Warning: Make sure to properly specify the \code{squashed},
#'   \code{zeroes}, and \code{N_star} arguments for your data set, since these
#'   will determine the appropriate likelihood function. Also, this function
#'   will not filter out data points. For instance, if you use \code{N_star = 2}
#'   you must filter out the ones and zeroes (if present) from \code{data} prior
#'   to using this function.
#' @details The \code{method} argument determines which optimization procedure
#'   is used. All the options use functions from the \code{\link{stats}}
#'   package:
#'   \itemize{
#'     \item{\code{"nlminb":}} \code{\link[stats]{nlminb}}
#'     \item{\code{"nlm":}} \code{\link[stats]{nlm}}
#'     \item{\code{"bfgs":}} \code{\link[stats]{optim}} (\emph{method = "BFGS"})
#'   }
#' @details Since this function runs multiple optimization procedures, it is
#'   best to start with 5 or less initial starting points (rows in
#'   \code{theta_init}). If the function runs in a reasonable amount of time,
#'   this number can be increased.
#' @details This function should not be used with very large data sets unless
#'   data squashing is used first since each optimization call will take a long
#'   time.
#' @details It is recommended to use \code{N_star = 1} when practical. Data
#'   squashing (see \code{\link{squashData}}) can be used to reduce the number
#'   of data points.
#' @details The \emph{converge} column in the resulting data frame was
#'   determined by examining the convergence \emph{code} of the chosen
#'   optimization method. In some instances, the code is somewhat ambiguous. The
#'   determination of \emph{converge} was intended to be conservative (leaning
#'   towards FALSE when questionable). See the documentation for the chosen
#'   method for details about \emph{code}.
#' @details Standard errors, if requested, are calculated using the observed
#'   Fisher information matrix as discussed in DuMouchel (1999).
#'
#' @examples
#' #Start with 2 or more guesses
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
#'   exploreHypers(squashed, theta_init = theta_init)
#' )
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family hyperparameter estimation functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlminb}}, \code{\link[stats]{nlm}}, and
#'   \code{\link[stats]{optim}} for optimization details
#' @seealso \code{\link{squashData}} for data preparation
#' @importFrom stats nlm
#' @importFrom stats nlminb
#' @importFrom stats optim
#' @importFrom stats optimHess
#' @export
exploreHypers <-
  function(data, theta_init, squashed = TRUE, zeroes = FALSE, N_star = 1,
           method = c("nlminb", "nlm", "bfgs"), param_limit = 100,
           max_pts = 20000, std_errors = FALSE) {

  .checkInputs_exploreHypers(data, theta_init, squashed, zeroes, N_star,
                             method, param_limit, max_pts, std_errors)
  method <- match.arg(method)

  results <- data.frame(guess_num = 1:nrow(theta_init), a1_hat = NA,
                        b1_hat = NA, a2_hat = NA, b2_hat = NA, p_hat = NA,
                        code = NA, converge = NA, in_bounds = NA, minimum = NA)
  std_errs <- NULL
  if (std_errors) {
    std_errs <- data.frame(guess_num = 1:nrow(theta_init), a1_se = NA,
                           b1_se = NA, a2_se = NA, b2_se = NA, p_se = NA)
  }

  #Choose the proper arguments for the optimization function
  arg_theta     <- "as.numeric(theta_init[i, ]), "
  arg_squ_pts   <- "ni = data$N, ei = data$E, wi = data$weight"
  arg_unsqu_pts <- "N = data$N, E = data$E"
  if (squashed) {
    if (zeroes) {
      arg_pts <- arg_squ_pts
      lik_fun <- "negLLzeroSquash, "
    } else {
      arg_pts <- paste0(arg_squ_pts, ", N_star = N_star")
      lik_fun <- "negLLsquash, "
    }
  } else {
    if (zeroes) {
      arg_pts <- arg_unsqu_pts
      lik_fun <- "negLLzero, "
    } else {
      arg_pts <- paste0(arg_unsqu_pts,  ", N_star = N_star")
      lik_fun <- "negLL, "
    }
  }

  #Run optimization and store results
  for (i in 1:nrow(theta_init)) {
    try({
      if (method == "nlminb") {
        ctl_args <- paste0(", control = list(eval.max = 500, iter.max = 500, ",
                           "abs.tol = 1e-20, xf.tol = 5e-14)")
        argum <- paste0("start = ", arg_theta, "objective = ", lik_fun, arg_pts,
                        ctl_args)
        guess_i <- eval(parse(text = paste0("nlminb(", argum, ")")))
        results[i, "minimum"] <- guess_i$objective
      } else if (method == "nlm") {
        argum <- paste0("f =", lik_fun,"p =", arg_theta, arg_pts,
                        ", iterlim = 250")
        guess_i <- eval(parse(text = paste0("nlm(", argum, ")")))
        results[i, 2:6]        <- guess_i$estimate
        results[i, "minimum"]  <- guess_i$minimum
        results[i, "code"]     <- guess_i$code
        results[i, "converge"] <- guess_i$code == 1 || guess_i$code == 2
        #Check if within low & high bounds
        ab_low <- all(guess_i$estimate[1:4] > 0)
        ab_hi  <- all(guess_i$estimate[1:4] < param_limit)
        p_low  <- guess_i$estimate[5] > 0
        p_hi   <- guess_i$estimate[5] < 1
        results[i, "in_bounds"] <- ab_low && ab_hi && p_low && p_hi
      } else {
        argum <- paste0("par = ", arg_theta, "fn = ", lik_fun, arg_pts,
                        ', method = "BFGS"')
        guess_i <- eval(parse(text = paste0("optim(", argum, ")")))
        results[i, "minimum"] <- guess_i$value
      }

      if (method %in% c("nlminb", "bfgs")) {
        results[i, 2:6]        <- guess_i$par
        results[i, "code"]     <- guess_i$convergence
        results[i, "converge"] <- guess_i$convergence == 0
        ab_low <- all(guess_i$par[1:4] > 0)
        ab_hi  <- all(guess_i$par[1:4] < param_limit)
        p_low  <- guess_i$par[5] > 0
        p_hi   <- guess_i$par[5] < 1
        results[i, "in_bounds"] <- ab_low && ab_hi && p_low && p_hi
      }

      if (std_errors) {
        theta_hat_i <- results[i, 2:6]
        hess_i <- .hessianMatrix(data, theta_hat_i, zeroes, squashed, N_star)
        std_errs[i, 2:6] <- .hyperStdErrs(hess_i)
      }
    }, silent = TRUE)
  }

  list(estimates = results, std_errs = std_errs)
}

#' Semi-automated hyperparameter estimation
#'
#' \code{autoHyper} finds a single hyperparameter estimate using an algorithm
#' that evaluates results from multiple starting points (see
#' \code{\link{exploreHypers}}). The algorithm verifies that the optimization
#' converges within the bounds of the parameter space and that the chosen
#' estimate (smallest negative log-likelihood) is similar to at least
#' one (see \code{min_conv} argument) of the other convergent solutions.
#'
#' @inheritParams exploreHypers
#' @param tol A numeric vector of tolerances for determining how close the
#'   chosen estimate must be to at least \code{min_conv} convergent solutions.
#'   Order is \eqn{\alpha_1}, \eqn{\beta_1}, \eqn{\alpha_2}, \eqn{\beta_2},
#'   \eqn{P}.
#' @param min_conv A scalar positive whole number for defining the minimum
#'   number of convergent solutions that must be close to the convergent
#'   solution with the smallest negative log-likelihood. Must be at least one
#'   and at most one less than the number of rows in \code{theta_init}.
#' @param conf_ints A scalar logical indicating if confidence intervals and
#'   standard errors should be returned.
#' @param conf_level A scalar string for the confidence level used if confidence
#'   intervals are requested.
#' @return A list containing the following elements:
#'   \itemize{
#'     \item{\emph{method}: }{A scalar character string for the method used to
#'       find the hyperparameter estimate (possibilities are
#'       \dQuote{\code{nlminb}}, \dQuote{\code{nlm}}, and
#'       \dQuote{\code{bfgs}}).}
#'     \item{\emph{estimates}: }{A named numeric vector of length 5 for the
#'       hyperparameter estimate corresponding to the smallest log-likelihood.}
#'     \item{\emph{conf_int}: }{A data frame including the standard errors and
#'       confidence limits. Only included if \code{conf_ints = TRUE}.}
#'     \item{\emph{num_close}: }{A scalar integer for the number of other
#'       convergent solutions that were close (within tolerance) to the chosen
#'       estimate.}
#'     \item{\emph{theta_hats}: }{A data frame for the estimates corresponding
#'       to the initial starting points defined by \code{theta_init}. See
#'       \code{\link{exploreHypers}}}.
#'   }
#'
#' @details The algorithm first attempts to find a consistently convergent
#'   solution using \code{\link[stats]{nlminb}}. If it fails, it will next try
#'   \code{\link[stats]{nlm}}. If it still fails, it will try
#'   \code{\link[stats]{optim}} (\emph{method = "BFGS"}). If all three
#'   approaches fail, the function returns an error message.
#' @details Since this function runs multiple optimization procedures, it is
#'   best to start with 5 or less initial starting points (rows in
#'   \code{theta_init}). If the function runs in a reasonable amount of time,
#'   this number can be increased.
#' @details This function should not be used with very large data sets since
#'   each optimization call will take a long time. \code{\link{squashData}} can
#'   be used first to reduce the size of the data.
#' @details It is recommended to use \code{N_star = 1} when practical. Data
#'   squashing (see \code{\link{squashData}}) can be used to further reduce the
#'   number of data points.
#' @details Asymptotic normal confidence intervals, if requested, use standard
#'   errors calculated from the observed Fisher information matrix as discussed
#'   in DuMouchel (1999).
#'
#' @examples
#' #Start with 2 or more guesses
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
#'   autoHyper(squashed, theta_init = theta_init)
#' )
#'
#' @references DuMouchel W (1999). "Bayesian Data Mining in Large Frequency
#'   Tables, With an Application to the FDA Spontaneous Reporting System."
#'   \emph{The American Statistician}, 53(3), 177-190.
#' @family hyperparameter estimation functions
#' @keywords openEBGM
#' @seealso \code{\link[stats]{nlminb}}, \code{\link[stats]{nlm}}, and
#'   \code{\link[stats]{optim}} for optimization details
#' @seealso \code{\link{squashData}} for data preparation
#' @importFrom stats optimHess
#' @importFrom stats qnorm
#' @export
autoHyper <-
  function(data, theta_init, squashed = TRUE, zeroes = FALSE, N_star = 1,
           tol = c(0.05, 0.05, 0.2, 0.2, 0.025), min_conv = 1,
           param_limit = 100, max_pts = 20000, conf_ints = FALSE,
           conf_level = c("95", "80", "90", "99")) {

  .checkInputs_autoHyper(tol, min_conv, theta_init, conf_ints)
  conf_level <- match.arg(conf_level)

  conf_int <- NULL
  estimate_names <- c("a1_hat", "b1_hat", "a2_hat", "b2_hat", "p_hat")

  for (i in c("nlminb", "nlm", "bfgs")) {
    theta_hats <-
      exploreHypers(data = data, theta_init = theta_init, squashed = squashed,
                    zeroes = zeroes, N_star = N_star, method = i,
                    param_limit = param_limit, max_pts = max_pts)$estimates

    #Only care about convergent results within parameter space
    conv          <- theta_hats$converge == TRUE & !is.na(theta_hats$converge)
    within_bounds <- theta_hats$in_bounds == TRUE & !is.na(theta_hats$in_bounds)
    theta_conv    <- theta_hats[conv & within_bounds, ]

    #Candidate must be "close" to at least 'min_conv' other estimates
    if (nrow(theta_conv) >= min_conv + 1) {
      candidate <- theta_conv[theta_conv$minimum == min(theta_conv$minimum), ]
      other     <- theta_hats[theta_hats$guess_num != candidate$guess_num, ]
      a1_close  <- abs(other$a1_hat - candidate$a1_hat) < tol[1]
      b1_close  <- abs(other$b1_hat - candidate$b1_hat) < tol[2]
      a2_close  <- abs(other$a2_hat - candidate$a2_hat) < tol[3]
      b2_close  <- abs(other$b2_hat - candidate$b2_hat) < tol[4]
      p_close   <- abs(other$p_hat  - candidate$p_hat)  < tol[5]
      is_close  <- a1_close & a2_close & b1_close & b2_close & p_close
      num_close <- sum(is_close, na.rm = TRUE)

      if (num_close >= min_conv) {
        theta_hat <- as.numeric(candidate[, estimate_names])
        names(theta_hat) <- c("alpha1", "beta1", "alpha2", "beta2", "P")
        if (conf_ints) {
          conf_int <- .hyperConfInts(data, theta_hat, zeroes, squashed, N_star,
                                     conf_level)
        }
        return(list(method = i, estimates = theta_hat, conf_int = conf_int,
                    num_close = num_close, theta_hats = theta_hats))
      }
    }
  }

  failed <- paste0("consistent convergence failed --",
                   "\n  try squashing data with another 'bin_size' value --",
                   "\n  if that fails, try using zeroes with data squashing --",
                   "\n  or, try using neither zeroes nor data squashing")
  stop(failed)
}
