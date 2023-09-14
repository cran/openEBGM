# ------------------------  Score functions  -----------------------------------
# Purpose:  Partial derivatives of log-likelihood.
# ------------------------------------------------------------------------------

.scoreAlpha <- function(alpha, beta, P_vec, N, E, W, squashed, zeroes,
                        reparam = TRUE) {
  # Partial derivative of logL w.r.t. alpha
  # Caution: Default is to input reparameterized hyperparameters.
  if (reparam) {
    # Input reparameterization:
    #   a_prime = log(alpha); b_prime = log(beta)
    #   P_prime = tan(pi * P - pi / 2)
    alpha <- exp(alpha)
    beta  <- exp(beta)
    P_vec <- atan(P_vec) / pi + 0.5
  }
  alpha_base <- 1 + E / beta
  alpha_base_log <- log(alpha_base)
  core <- digamma(alpha + N) - digamma(alpha) - alpha_base_log
  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_base_log
    zero_adjust2 <- alpha_base ^ alpha - 1
    if (any(is.infinite(zero_adjust2))) stop("skip this iteration")  # k / Inf = 0
    zero_adjust <- zero_adjust1 / zero_adjust2
    score <- core - zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  if (is.infinite(score) || is.na(score)) stop("skip this iteration")
  score
}

.scoreBeta <- function(beta, alpha, P_vec, N, E, W, squashed, zeroes,
                       reparam = TRUE) {
  # Partial derivative of logL w.r.t. beta
  # Caution: Default is to input reparameterized hyperparameters.
  if (reparam) {
    # Input reparameterization:
    #   a_prime = log(alpha); b_prime = log(beta)
    #   P_prime = tan(pi * P - pi / 2)
    alpha <- exp(alpha)
    beta  <- exp(beta)
    P_vec <- atan(P_vec) / pi + 0.5
  }
  alpha_base <- 1 + E / beta
  alpha_E <- alpha * E
  core <- (alpha_E - N * beta) / (beta * (beta + E))
  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_E
    zero_adjust2 <- beta ^ 2 * (alpha_base ^ (alpha + 1) - alpha_base)
    if (any(is.infinite(zero_adjust2))) stop("skip this iteration")  # k / Inf = 0
    zero_adjust  <- zero_adjust1 / zero_adjust2
    score <- core + zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  if (is.infinite(score) || is.na(score)) stop("skip this iteration")
  score
}

.scoreFrac <- function(P, alpha1, beta1, alpha2, beta2, N, E, W, squashed,
                       zeroes, reparam = TRUE) {
  # Partial derivative of logL w.r.t. mixing fraction (P)
  # Caution: Default is to input reparameterized hyperparameters.
  if (reparam) {
    # Input reparameterization:
    #   a_prime = log(alpha); b_prime = log(beta)
    #   P_prime = tan(pi * P - pi / 2)
    alpha1_log <- alpha1; beta1_log <- beta1
    alpha2_log <- alpha2; beta2_log <- beta2
    alpha1 <- exp(alpha1); beta1 <- exp(beta1)
    alpha2 <- exp(alpha2); beta2 <- exp(beta2)
    P      <- atan(P) / pi + 0.5
  }
  prob_f1 <- beta1 / (beta1 + E)
  prob_f2 <- beta2 / (beta2 + E)
  f1 <- dnbinom(N, alpha1, prob_f1)
  f2 <- dnbinom(N, alpha2, prob_f2)
  if (zeroes) {
    num <- f1 - f2
    den <- P * f1 + (1 - P) * f2
  } else {
    f1_adjust <- 1 - dnbinom(0, alpha1, prob_f1)
    f2_adjust <- 1 - dnbinom(0, alpha2, prob_f2)
    f1_star <- f1 / f1_adjust
    f2_star <- f2 / f2_adjust
    num <- f1_star - f2_star
    den <- P * f1_star + (1 - P) * f2_star
  }
  score <- num / den
  if (squashed) score <- W * score
  score <- sum(score)
  if (is.infinite(score) || is.na(score)) stop("skip this iteration")
  score
}


.checkScore <- function(theta, N, E, W, squashed, zeroes) {
  # Calculate score vector & its Euclidean norm at theta (should be close to 0).
  marg_dens  <- .marginalDensity(theta, N, E, zeroes)
  PQ         <- .updateProbs(theta, marg_dens$f1, marg_dens$f2)
  score_a1   <- .scoreAlpha(theta[1], theta[2], PQ$P_vec, N, E, W, squashed,
                            zeroes, reparam = FALSE)
  score_b1   <- .scoreBeta(theta[2], theta[1], PQ$P_vec, N, E, W, squashed,
                           zeroes, reparam = FALSE)
  score_a2   <- .scoreAlpha(theta[3], theta[4], PQ$Q_vec, N, E, W, squashed,
                            zeroes, reparam = FALSE)
  score_b2   <- .scoreBeta(theta[4], theta[3], PQ$Q_vec, N, E, W, squashed,
                           zeroes, reparam = FALSE)
  score_P    <- .scoreFrac(theta[5], theta[1], theta[2], theta[3], theta[4],
                           N, E, W, squashed, zeroes, reparam = FALSE)
  score_vec  <- c(score_a1, score_b1, score_a2, score_b2, score_P)
  score_norm <- sqrt(sum(score_vec ^ 2))
  list(score_vec = score_vec, score_norm = score_norm)
}
