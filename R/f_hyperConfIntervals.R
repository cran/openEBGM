#Helper functions to calculate hyperparameter standard errors and confidence
#intervals

.critValNormal <- function(conf_level) {
  #Choose critical value from normal dist'n for given confidence level
  switch(conf_level,
         "95" = qnorm(0.975),
         "80" = qnorm(0.9),
         "90" = qnorm(0.95),
         "99" = qnorm(0.995)
  )
}

.hessianMatrix <- function(data, theta_hat, zeroes, squashed, N_star) {
  #Hessian matrix used to find hyperparameter standard errors (also covariances)
  N <- data$N
  E <- data$E
  if (zeroes) {
    if (squashed) {
      W <- data$weight
      return(optimHess(theta_hat, fn = negLLzeroSquash, ni = N, ei = E, wi = W))
    } else {
      return(optimHess(theta_hat, fn = negLLzero, N = N, E = E))
    }
  } else {
    if (squashed) {
      W <- data$weight
      return(optimHess(theta_hat, fn = negLLsquash, ni = N, ei = E, wi = W,
                       N_star = N_star))
    } else {
      return(optimHess(theta_hat, fn = negLL, N = N, E = E, N_star = N_star))
    }
  }
}

.hyperStdErrs <- function(hessian_matrix) {
  #Hyperparameter standard errors
  hessian_inverse <- solve(hessian_matrix)
  variances <- diag(hessian_inverse)
  sqrt(variances)
}

.hyperConfInts <- function(data, theta_hat, zeroes, squashed, N_star,
                           conf_level) {
  #Hyperparameter asymptotic normal confidence intervals
  crit_val <- .critValNormal(conf_level)
  estimate_names <- c("a1_hat", "b1_hat", "a2_hat", "b2_hat", "p_hat")
  conf_int <- matrix(NA, nrow = 5, ncol = 4)
  conf_int <- as.data.frame(conf_int, row.names = estimate_names)
  LL_name <- paste0("LL_", conf_level)
  UL_name <- paste0("UL_", conf_level)
  colnames(conf_int) <- c("pt_est", "SE", LL_name, UL_name)

  hess <- .hessianMatrix(data, theta_hat, zeroes, squashed, N_star)
  std_errs <- .hyperStdErrs(hess)
  conf_int$pt_est <- theta_hat
  conf_int$SE <- std_errs
  ME <- crit_val * std_errs
  conf_int[LL_name] <- theta_hat - ME
  conf_int[UL_name] <- theta_hat + ME
  conf_int <- round(conf_int, 4)
  conf_int
}
