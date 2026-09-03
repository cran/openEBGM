# Helper functions for updating theta (EM algorithm approach)

# -----------------------  Update theta ----------------------------------------
# Purpose: Find updated values of theta for a given iteration.
# ------------------------------------------------------------------------------

.updateTheta <- function(theta, P_vec, Q_vec, N, E, W, squashed, zeroes) {
  # Update theta using score functions and conditional optimization.
  # Optimizes one hyperparameter at a time.
  # method = "score"
  reparamP <- function(P) log(P) - log(1 - P)  #param space (-Inf, Inf)
  searchSpace <- function(x, distance) c(x - distance, x + distance)
  findRoot <- function(f, last_est, search_dist, fixed, probs) {
    root <- tryCatch({
      uniroot(f, interval = searchSpace(last_est, search_dist),
        fixed, probs, N = N, E = E, W = W, squashed = squashed, zeroes = zeroes,
        extendInt = "yes", maxiter = 50, tol = .0001
      )$root
    }, error = function(e) {
      return(last_est)
    })
    if (root > 11.5) root <- stats::runif(1, min = 8.5, max = 9.5)  #to "rein in" the estimate
    if (root < -11.5) root <- stats::runif(1, min = -11.5, max = -10)  #log(1e-5)
    root
  }
  search_dist_ab <- .1; search_dist_P <- .1
  theta_prime <- theta
  theta_prime[1:4] <- log(theta_prime[1:4])  #parameter space (-Inf, Inf)
  theta_prime[5] <- reparamP(theta_prime[5])
  P_prime_vec <- reparamP(P_vec); Q_prime_vec <- reparamP(Q_vec)
  old_a1_prime <- theta_prime[1]; old_b1_prime <- theta_prime[2]
  old_a2_prime <- theta_prime[3]; old_b2_prime <- theta_prime[4]
  old_P_prime <- theta_prime[5]
  a1_prime <- findRoot(.scoreAlpha,
    old_a1_prime, search_dist_ab, fixed = old_b1_prime, P_prime_vec
  )
  b1_prime <- findRoot(.scoreBeta,
    old_b1_prime, search_dist_ab, fixed = a1_prime, P_prime_vec
  )
  a2_prime <- findRoot(.scoreAlpha,
    old_a2_prime, search_dist_ab, fixed = old_b2_prime, Q_prime_vec
  )
  b2_prime <- findRoot(.scoreBeta,
    old_b2_prime, search_dist_ab, fixed = a2_prime, Q_prime_vec
  )
  search_sp_P <- searchSpace(theta_prime[5], search_dist_P)
  P_prime <- tryCatch({
    uniroot(.scoreFrac,
      interval = search_sp_P, alpha1 = a1_prime, beta1 = b1_prime,
      alpha2 = a2_prime, beta2 = b2_prime, N = N, E = E, W = W,
      squashed = squashed, zeroes = zeroes,
      extendInt = "yes", maxiter = 50, tol = .0001
    )$root
  }, error = function(e) {
    return(old_P_prime)
  })
  if (P_prime > 20) P_prime <- stats::runif(1, min = 15, max = 20)
  if (P_prime < -20) P_prime <- stats::runif(1, min = -20, max = -15)
  theta <- c(a1_prime, b1_prime, a2_prime, b2_prime, P_prime)
  # Back to original scale
  theta[1:4] <- exp(theta[1:4])
  theta[5] <- stats::plogis(theta[5])
  theta
}

.updateThetaLL <- function(theta, N, E, W, squashed, zeroes, N_star) {
  # Update theta using log-likelihood functions and conditional optimization.
  # Optimizes one parameter at a time.
  # method = "nlminb", profile = "parameter"
  theta[1:4] <- log(theta[1:4])
  theta[5] <- log(theta[5]) - log(1 - theta[5])
  olda1 <- theta[1]; oldb1 <- theta[2]; olda2 <- theta[3]; oldb2 <- theta[4]
  oldP <- theta[5]
  optimizeParam <- function(start, param_str) {
    par <- tryCatch({
      nlminb(start, objective = .LLtheta, param = param_str,
        N = N, E = E, W = W, N_star = N_star,
        theta = theta, squashed = squashed, zeroes = zeroes
      )$par
    }, error = function(e) {
      return(start)
    })
    if (par > 11.5) par <- stats::runif(1, min = 8.5, max = 9.5)  #to "rein in" the estimate
    if (par < -11.5) par <- stats::runif(1, min = -11.5, max = -10)  #log(1e-5)
    par
  }
  alpha1 <- optimizeParam(olda1, param_str = "alpha1")
  theta[1] <- alpha1
  beta1 <- optimizeParam(oldb1, param_str = "beta1")
  theta[2] <- beta1
  alpha2 <- optimizeParam(olda2, param_str = "alpha2")
  theta[3] <- alpha2
  beta2 <- optimizeParam(oldb2, param_str = "beta2")
  theta[4] <- beta2
  P <- tryCatch({
    nlminb(start = oldP, objective = .LLtheta, param = "P",
      N = N, E = E, W = W, N_star = N_star,
      theta = theta, squashed = squashed, zeroes = zeroes
    )$par
  }, error = function(e) {
    return(oldP)
  })
  if (P > 20) P <- stats::runif(1, min = 15, max = 20)
  if (P < -20) P <- stats::runif(1, min = -20, max = -15)
  theta <- c(alpha1, beta1, alpha2, beta2, P)
  # Back to original scale
  theta[1:4] <- exp(theta[1:4])
  theta[5] <- stats::plogis(theta[5])
  theta
}

.updateThetaLLD <- function(theta, N, E, W, squashed, zeroes, N_star) {
  # Update theta using log-likelihood functions and conditional optimization.
  # Optimizes one distribution at a time.
  # method = "nlminb", profile = "distribution"
  theta[1:4] <- log(theta[1:4])
  theta[5] <- log(theta[5]) - log(1 - theta[5])
  oldD1 <- theta[1:2]
  oldD2 <- theta[3:4]
  oldP <- theta[5]
  optimizeParam <- function(start, param_str) {
    par <- tryCatch({
      nlminb(start, objective = .LLthetaD, param = param_str,
        N = N, E = E, W = W, N_star = N_star,
        theta = theta, squashed = squashed, zeroes = zeroes
      )$par
    }, error = function(e) {
      return(start)
    })
    if (par[1] > 11.5) par[1] <- stats::runif(1, min = 8.5, max = 9.5)  #to "rein in" the estimate
    if (par[1] < -11.5) par[1] <- stats::runif(1, min = -11.5, max = -10)  #log(1e-5)
    if (par[2] > 11.5) par[2] <- stats::runif(1, min = 8.5, max = 9.5)
    if (par[2] < -11.5) par[2] <- stats::runif(1, min = -11.5, max = -10)
    par
  }
  D1 <- optimizeParam(oldD1, param_str = "D1")
  theta[1:2] <- D1
  D2 <- optimizeParam(oldD2, param_str = "D2")
  theta[3:4] <- D2
  P <- tryCatch({
    nlminb(start = oldP, objective = .LLtheta, param = "P",
      N = N, E = E, W = W, N_star = N_star,
      theta = theta, squashed = squashed, zeroes = zeroes
    )$par
  }, error = function(e) {
    return(oldP)
  })
  if (P > 20) P <- stats::runif(1, min = 15, max = 20)
  if (P < -20) P <- stats::runif(1, min = -20, max = -15)
  theta <- c(D1, D2, P)
  # Back to original scale
  theta[1:4] <- exp(theta[1:4])
  theta[5] <- stats::plogis(theta[5])
  theta
}

.LLtheta <- function(value, param_str, N, E, W, theta, N_star, squashed, zeroes) {
  # Negative log-likelihood for optimizing a single element in theta.
  if (param_str == "alpha1") {
    theta[1] <- value
  } else if (param_str == "beta1") {
    theta[2] <- value
  } else if (param_str == "alpha2") {
    theta[3] <- value
  } else if (param_str == "beta2") {
    theta[4] <- value
  } else if (param_str == "P") {
    theta[5] <- value
  } else {
    stop("actual argument for 'param' not recognized")
  }
  if (zeroes) {
    if (squashed) {
      return(negLLzeroSquash(theta, N, E, W, transformed = TRUE))
    } else {
      return(negLLzero(theta, N, E, transformed = TRUE))
    }
  } else {
    if (squashed) {
      return(negLLsquash(theta, N, E, W, N_star, transformed = TRUE))
    } else {
      return(negLL(theta, N, E, N_star, transformed = TRUE))
    }
  }
}

.LLthetaD <- function(value, param_str, N, E, W, theta, N_star, squashed, zeroes) {
  # Negative log-likelihood for optimizing a single distribution in theta.
  if (param_str == "D1") {
    theta[1] <- value[1]
    theta[2] <- value[2]
  } else if (param_str == "D2") {
    theta[3] <- value[1]
    theta[4] <- value[2]
  } else {
    stop("actual argument for 'param' not recognized")
  }
  if(zeroes) {
    if(squashed) {
      return(negLLzeroSquash(theta, N, E, W, transformed = TRUE))
    } else {
      return(negLLzero(theta, N, E, transformed = TRUE))
    }
  } else {
    if(squashed) {
      return(negLLsquash(theta, N, E, W, N_star, transformed = TRUE))
    } else {
      return(negLL(theta, N, E, N_star, transformed = TRUE))
    }
  }
}


# --------------------  Find marginal density ----------------------------------
# Purpose: Find the marginal density of each component.
# ------------------------------------------------------------------------------

.marginalDensity <- function(theta, N, E, zeroes) {
  # Marginal density (log scale) of each component
  alpha1 <- theta[1]; beta1 <- theta[2]; alpha2 <- theta[3]; beta2 <- theta[4]
  prob_f1 <- exp(log(beta1) - log(beta1 + E))  #beta1 / (beta1 + E)
  prob_f2 <- exp(log(beta2) - log(beta2 + E))
  prob_f1[prob_f1 < 1e-10] <- 1e-10
  prob_f2[prob_f2 < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f1[prob_f1 > max_prob] <- max_prob
  prob_f2[prob_f2 > max_prob] <- max_prob
  log_f1 <- dnbinom(N, alpha1, prob_f1, log = TRUE)
  log_f2 <- dnbinom(N, alpha2, prob_f2, log = TRUE)
  if (!zeroes) {
    log_f1_adjust <- pnbinom(0, alpha1, prob_f1, lower.tail = FALSE, log.p = TRUE)
    log_f2_adjust <- pnbinom(0, alpha2, prob_f2, lower.tail = FALSE, log.p = TRUE)
    log_f1_star <- log_f1 - log_f1_adjust
    log_f2_star <- log_f2 - log_f2_adjust
    log_f1 <- log_f1_star
    log_f2 <- log_f2_star
  }
  list(log_f1 = log_f1, log_f2 = log_f2)
}

# ------------------  Update point probabilities -------------------------------
# Purpose: Find updated values of point probabilities for a given iteration.
#          Point probabilities are updated by multiplying by the marginal
#            density & normalizing to sum to 1.
#          i.e., P(in first component) + P(in second component) = 1
# ------------------------------------------------------------------------------

.updateProbs <- function(theta, log_marg_dens_f1, log_marg_dens_f2) {
  # Update point probabilities (responsibilities)
  P <- theta[5]
  P_vec <- exp(log(P) + log_marg_dens_f1)
  Q_vec <- exp(log(1 - P) + log_marg_dens_f2)
  PQ_sum <- P_vec + Q_vec
  P_vec <- exp(log(P_vec) - log(PQ_sum))
  Q_vec  <- exp(log(Q_vec) - log(PQ_sum))
  P_vec[P_vec < 1e-10] <- 1e-10
  Q_vec[Q_vec < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  P_vec[P_vec > max_prob] <- max_prob
  Q_vec[Q_vec > max_prob] <- max_prob
  list(P_vec = P_vec, Q_vec = Q_vec)
}

# ------------------  Aitken acceleration of theta  ----------------------------
# Purpose: Accelerate the updating process for theta using the method described
#            in Section 5.3.3 of the book "Maximum Likelihood Estimation and
#            Inference: With Examples in R, SAS and ADMB" by Russell B. Millar.
# ------------------------------------------------------------------------------
.accelerateTheta <- function(theta1, theta2, theta3) {
  # Use Aitken acceleration to speed up EM algorithm
  acc_ratio_num <- theta3 - theta2
  acc_ratio_den <- theta2 - theta1
  acc_ratio <- acc_ratio_num / acc_ratio_den  #called 'a' in referenced book
  acc_ratio[is.nan(acc_ratio)] <- 2  #to be filtered in next line
  estimate_unstable <- abs(acc_ratio) > .99  #(1 - a) is close to zero
  theta_acc <- theta2 + acc_ratio_num / (1 - acc_ratio)  #'theta hat'
  theta_acc[estimate_unstable] <- theta3[estimate_unstable]
  alphas_betas <- theta_acc[1:4]
  P <- theta_acc[5]
  alphas_betas[alphas_betas < 1e-10] <- 1e-10
  P_lower <- 1e-10
  P_upper <- 1 - P_lower
  P <- max(P, P_lower)
  P <- min(P, P_upper)
  c(alphas_betas, P)
}
