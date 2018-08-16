# Helper functions for updating theta (EM algorithm approach)

# -----------------------  Update theta ----------------------------------------
# Purpose: Find updated values of theta for a given iteration.
# ------------------------------------------------------------------------------

.updateTheta <- function(theta, P_vec, Q_vec, N, E, W, squashed, zeroes,
                         param_lower, param_upper) {
  #Update theta using score functions and conditional optimization.
  #Optimizes one hyperparameter at a time.

  reparamP    <- function(P) tan(pi * P - pi / 2)  #param space (-Inf, Inf)
  searchSpace <- function(x, distance) c(x - distance, x + distance)

  lower_lim <- log(param_lower); upper_lim <- log(param_upper)

  findRoot <- function(f, last_est, search_dist, fixed, probs) {
    root <- tryCatch({
      uniroot(f, interval = searchSpace(last_est, search_dist),
              fixed, probs, N = N, E = E, W = W,
              squashed = squashed, zeroes = zeroes,
              #extendInt = "yes", maxiter = 25, trace = 2)$root  #for debugging
              extendInt = "yes", maxiter = 25)$root
    }, error = function(e) {
      return(last_est)
    })
    if (root > upper_lim) root <- .5 * upper_lim  #to "rein in" the estimate
    if (root < lower_lim) root <- lower_lim
    root
  }

  search_dist_ab <- .0001; search_dist_P <- .01

  theta_prime      <- theta
  theta_prime[1:4] <- log(theta_prime[1:4])  #parameter space (-Inf, Inf)
  theta_prime[5]   <- reparamP(theta_prime[5])
  P_prime_vec <- reparamP(P_vec); Q_prime_vec <- reparamP(Q_vec)

  old_a1_prime <- theta_prime[1]; old_b1_prime <- theta_prime[2]
  old_a2_prime <- theta_prime[3]; old_b2_prime <- theta_prime[4]
  old_P_prime  <- theta_prime[5]

  a1_prime <- findRoot(.scoreAlpha, old_a1_prime, search_dist_ab,
                       fixed = old_b1_prime, P_prime_vec)
  b1_prime <- findRoot(.scoreBeta, old_b1_prime, search_dist_ab,
                       fixed = a1_prime, P_prime_vec)
  a2_prime <- findRoot(.scoreAlpha, old_a2_prime, search_dist_ab,
                       fixed = old_b2_prime, Q_prime_vec)
  b2_prime <- findRoot(.scoreBeta, old_b2_prime, search_dist_ab,
                       fixed = a2_prime, Q_prime_vec)

  P_prime <- tryCatch({
    search_sp_P <- searchSpace(theta_prime[5], search_dist_P)
    P_prime <- uniroot(.scoreFrac, interval = search_sp_P, alpha1 = a1_prime,
                       beta1 = b1_prime, alpha2 = a2_prime, beta2 = b2_prime,
                       N = N, E = E, W = W, squashed = squashed,
                       zeroes = zeroes, extendInt = "yes", maxiter = 25)$root
  }, error = function(e) {
    return(old_P_prime)
  })

  # Back to original scale (i.e., backtransform)
  alpha1 <- exp(a1_prime); beta1 <- exp(b1_prime)
  alpha2 <- exp(a2_prime); beta2 <- exp(b2_prime)
  P      <- atan(P_prime) / pi + 0.5

  c(alpha1, beta1, alpha2, beta2, P)
}


.updateThetaLL <- function(theta, N, E, W, squashed, zeroes, N_star,
                           param_lower, param_upper) {
  #Update theta using log-likelihood functions and conditional optimization.
  #Optimizes one parameter at a time.

  olda1 <- theta[1]; oldb1 <- theta[2]; olda2 <- theta[3]; oldb2 <- theta[4]
  oldP <- theta[5]

  optimizeParam <- function(start, param_str) {
    nlminb(start, objective = .LLtheta, lower = param_lower, upper = param_upper,
           param = param_str, N = N, E = E, W = W, N_star = N_star,
           theta = theta, squashed = squashed, zeroes = zeroes)$par
  }

  alpha1 <- optimizeParam(olda1, param_str = "alpha1")
  theta[1] <- alpha1
  beta1 <- optimizeParam(oldb1, param_str = "beta1")
  theta[2] <- beta1
  alpha2 <- optimizeParam(olda2, param_str = "alpha2")
  theta[3] <- alpha2
  beta2 <- optimizeParam(oldb2, param_str = "beta2")
  theta[4] <- beta2
  P <- nlminb(start = oldP, objective = .LLtheta, lower = 0.001, upper = 0.999,
              param = "P", N = N, E = E, W = W, N_star = N_star, theta = theta,
              squashed = squashed, zeroes = zeroes)$par

  c(alpha1, beta1, alpha2, beta2, P)
}


.updateThetaLLD <- function(theta, N, E, W, squashed, zeroes, N_star,
                            param_lower, param_upper) {
  #Update theta using log-likelihood functions and conditional optimization.
  #Optimizes one distribution at a time.

  oldD1 <- c(theta[1], theta[2])
  oldD2 <- c(theta[3], theta[4])
  oldP <- theta[5]

  optimizeParam <- function(start, param_str) {
    nlminb(start, objective = .LLthetaD, lower = param_lower, upper = param_upper,
           param = param_str, N = N, E = E, W = W, N_star = N_star,
           theta = theta, squashed = squashed, zeroes = zeroes)$par
  }

  D1 <- optimizeParam(oldD1, param_str = "D1")
  theta[1] <- D1[1]
  theta[2] <- D1[2]
  D2 <- optimizeParam(oldD2, param_str = "D2")
  theta[3] <- D2[1]
  theta[4] <- D2[2]
  P <- nlminb(start = oldP, objective = .LLtheta, lower = 0.001, upper = 0.999,
              param = "P", N = N, E = E, W = W, N_star = N_star, theta = theta,
              squashed = squashed, zeroes = zeroes)$par
              #, control = list(xf.tol = 1e-4, x.tol = 1e-4, rel.tol = 1e-4))$par

  c(D1, D2, P)
}


.LLtheta <- function(value, param_str, N, E, W, theta, N_star, squashed,
                     zeroes) {
  #Negative log-likelihood for optimizing a single element in theta.

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

  if(zeroes) {
    if(squashed) {
      return(negLLzeroSquash(theta, N, E, W))
    } else {
      return(negLLzero(theta, N, E))
    }
  } else {
    if(squashed) {
      return(negLLsquash(theta, N, E, W, N_star))
    } else {
      return(negLL(theta, N, E, N_star))
    }
  }
}


.LLthetaD <- function(value, param_str, N, E, W, theta, N_star, squashed,
                      zeroes) {
  #Negative log-likelihood for optimizing a single distribution in theta.

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
      return(negLLzeroSquash(theta, N, E, W))
    } else {
      return(negLLzero(theta, N, E))
    }
  } else {
    if(squashed) {
      return(negLLsquash(theta, N, E, W, N_star))
    } else {
      return(negLL(theta, N, E, N_star))
    }
  }
}


# --------------------  Find marginal density ----------------------------------
# Purpose: Find the marginal density of each component.
# ------------------------------------------------------------------------------

.marginalDensity <- function(theta, N, E, zeroes) {
  #Marginal density of each component

  alpha1 <- theta[1]; beta1 <- theta[2]; alpha2 <- theta[3]; beta2 <- theta[4]

  prob_f1 <- beta1 / (beta1 + E)
  prob_f2 <- beta2 / (beta2 + E)
  f1 <- dnbinom(N, alpha1, prob_f1)
  # if (any(is.na(f1))) stop(".marginalDensity(): f1 has NAs")  #debugging
  # if (min(f1) == 0) stop(".marginalDensity(): f1 has zeroes")
  f2 <- dnbinom(N, alpha2, prob_f2)
  # if (any(is.na(f2))) stop(".marginalDensity(): f2 has NAs")
  # if (min(f2) == 0) stop(".marginalDensity(): f2 has zeroes")

  if(!zeroes) {
    f1_adjust <- 1 - dnbinom(0, alpha1, prob_f1)
    f2_adjust <- 1 - dnbinom(0, alpha2, prob_f2)
    # if (any(is.na(f1_adjust))) stop(".marginalDensity(): f1_adjust has NAs")
    # if (any(is.na(f2_adjust))) stop(".marginalDensity(): f2_adjust has NAs")
    # if (min(f1_adjust) == 0) stop(".marginalDensity(): f1_adjust has zeroes")
    # if (min(f2_adjust) == 0) stop(".marginalDensity(): f2_adjust has zeroes")
    f1_star <- f1 / f1_adjust
    f2_star <- f2 / f2_adjust
    # if (any(is.na(f1_star))) stop(".marginalDensity(): f1_star has NAs")
    # if (any(is.na(f2_star))) stop(".marginalDensity(): f2_star has NAs")
    f1 <- f1_star
    f2 <- f2_star
  }

  list(f1 = f1, f2 = f2)
}


# ------------------  Update point probabilities -------------------------------
# Purpose: Find updated values of point probabilities for a given iteration.
#          Point probabilities are updated by multiplying by the marginal
#            density & normalizing to sum to 1.
#          i.e., P(in first component) + P(in second component) = 1
# ------------------------------------------------------------------------------

.updateProbs <- function(theta, marg_dens_f1, marg_dens_f2) {
  #Update point probabilities (responsibilities)

  P      <- theta[5]
  P_vec  <- P * marg_dens_f1
  Q_vec  <- (1 - P) * marg_dens_f2
  PQ_sum <- P_vec + Q_vec
  # if (min(PQ_sum) == 0) stop(".updateProbs(): 'PQ_sum' has zeroes")
  P_vec  <- P_vec / PQ_sum
  # if (any(is.na(P_vec))) stop(".updateProbs(): P_vec has NAs")
  # if (any(is.infinite(P_vec))) stop(".updateProbs(): P_vec has Infs")
  Q_vec  <- Q_vec / PQ_sum
  # if (any(is.na(Q_vec))) stop(".updateProbs(): Q_vec has NAs")
  # if (any(is.infinite(Q_vec))) stop(".updateProbs(): Q_vec has Infs")

  list(P_vec = P_vec, Q_vec = Q_vec)
}


# -------------------------  Nudge theta  --------------------------------------
# Purpose: Give theta a nudge when it gets stuck for several iterations.
# ------------------------------------------------------------------------------
.nudgeTheta <- function(theta, count_stuck, stuck_limit, param_lower) {
  #Nudge theta when it gets stuck for several iterations.

  theta_nudge <- theta - 0.01
  mid_way <- 0.5 * (theta + param_lower)
  theta_nudge <- ifelse(theta_nudge < 0, mid_way, theta_nudge)
  theta[count_stuck > stuck_limit] <- theta_nudge[count_stuck > stuck_limit]
  theta
}


# ------------------  Aitken acceleration of theta  ----------------------------
# Purpose: Accelerate the updating process for theta using the method described
#            in Section 5.3.3 of the book "Maximum Likelihood Estimation and
#            Inference: With Examples in R, SAS and ADMB" by Russell B. Millar.
# ------------------------------------------------------------------------------
.accelerateTheta <- function(theta1, theta2, theta3, param_lower, param_upper) {
  #Use Aitken acceleration to speed up EM algorithm

  acc_ratio_num <- theta3 - theta2
  acc_ratio_den <- theta2 - theta1
  acc_ratio <- acc_ratio_num / acc_ratio_den  #called 'a' in referenced book
  acc_ratio[is.nan(acc_ratio)] <- 2  #to be filtered in next line
  estimate_unstable <- abs(acc_ratio) > .99  #(1 - a) is close to zero
  theta_acc <- theta2 + acc_ratio_num / (1 - acc_ratio)  #'theta hat'
  theta_acc[estimate_unstable] <- theta3[estimate_unstable]
  theta_acc[theta_acc < param_lower] <- param_lower
  theta_acc[theta_acc > param_upper] <- param_upper
  P_lower <- 1e-10
  P_upper <- 1 - P_lower
  if (theta_acc[5] < P_lower) theta_acc[5] <- P_lower
  if (theta_acc[5] > P_upper) theta_acc[5] <- P_upper

  theta_acc
}
