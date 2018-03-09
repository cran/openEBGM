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
  # if (any(is.na(alpha_base))) stop(".scoreAlpha() - NAs in 'alpha_base'")  #debugging
  # if (any(is.infinite(alpha_base))) stop(".scoreAlpha() - Infs in 'alpha_base'")
  # if (min(alpha_base) == 0) stop(".scoreAlpha() - zeroes in 'alpha_base'")

  alpha_base_log <- log(alpha_base)
  # if (any(is.na(alpha_base_log))) stop(".scoreAlpha() - NAs in 'alpha_base_log'")
  # if (any(is.infinite(alpha_base_log))) stop(".scoreAlpha() - Infs in 'alpha_base_log'")
  # digamma_alpha <- digamma(alpha)
  # if (is.nan(digamma_alpha)) stop(".scoreAlpha() - digamma(alpha) is NaN")
  # digamma_alpha_N <- digamma(alpha + N)
  # if (any(is.nan(digamma_alpha_N))) stop(".scoreAlpha() - digamma(alpha + N) is NaN")
  # core <- digamma_alpha_N - digamma_alpha - alpha_base_log
  core <- digamma(alpha + N) - digamma(alpha) - alpha_base_log
  # if (any(is.na(core))) stop(".scoreAlpha() - NAs in 'core'")
  # if (any(is.infinite(core))) stop(".scoreAlpha() - Infs in 'core'")

  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_base_log
    # if (any(is.na(zero_adjust1))) stop(".scoreAlpha() - NAs in 'zero_adjust1'")
    # if (any(is.infinite(zero_adjust1))) stop(".scoreAlpha() - Infs in 'zero_adjust1'")
    zero_adjust2 <- alpha_base ^ alpha - 1
    if (any(is.infinite(zero_adjust2))) stop("skip this iteration")  # k / Inf = 0
    # if (any(is.na(zero_adjust2))) stop(".scoreAlpha() - NAs in 'zero_adjust2'")
    # if (any(is.infinite(zero_adjust2))) stop(".scoreAlpha() - Infs in 'zero_adjust2'")
    # if (min(zero_adjust2) == 0) stop(".scoreAlpha() - zeroes in 'zero_adjust2'")
    zero_adjust <- zero_adjust1 / zero_adjust2
    score <- core - zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  # if (is.na(score)) stop(".scoreAlpha() - 'score' is NA")
  # if (is.infinite(score)) stop(".scoreAlpha() - 'score' is Inf")
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
  # if (any(is.na(alpha_base))) stop(".scoreBeta() - NAs in 'alpha_base'")
  # if (any(is.infinite(alpha_base))) stop(".scoreBeta() - Infs in 'alpha_base'")
  # if (min(alpha_base) == 0) stop(".scoreBeta() - zeroes in 'alpha_base'")

  alpha_E <- alpha * E
  # if (any(is.na(alpha_E))) stop(".scoreBeta() - NAs in 'alpha_E'")
  # if (any(is.infinite(alpha_E))) stop(".scoreBeta() - Infs in 'alpha_E'")
  core <- (alpha_E - N * beta) / (beta * (beta + E))
  # if (any(is.na(core))) stop(".scoreBeta() - NAs in 'core'")
  # if (any(is.infinite(core))) stop(".scoreBeta() - Infs in 'core'")

  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_E
    zero_adjust2 <- beta ^ 2 * (alpha_base ^ (alpha + 1) - alpha_base)
    if (any(is.infinite(zero_adjust2))) stop("skip this iteration")  # k / Inf = 0
    # if (any(is.na(zero_adjust2))) stop(".scoreBeta() - NAs in 'zero_adjust2'")
    # if (any(is.infinite(zero_adjust2))) stop(".scoreBeta() - Infs in 'zero_adjust2'")
    # if (min(zero_adjust2) == 0) stop(".scoreBeta() - zeroes in 'zero_adjust2'")
    zero_adjust  <- zero_adjust1 / zero_adjust2
    score <- core + zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  # if (is.na(score)) stop(".scoreBeta() - 'score' is NA")
  # if (is.infinite(score)) stop(".scoreBeta() - 'score' is Inf")
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
  # if (any(is.na(prob_f1))) stop(".scoreFrac() - NAs in 'prob_f1'")
  # if (any(is.infinite(prob_f1))) stop(".scoreFrac() - Infs in 'prob_f1'")
  # if (min(prob_f1) == 0) stop(".scoreFrac() - zeroes in 'prob_f1'")
  # if (any(is.na(prob_f2))) stop(".scoreFrac() - NAs in 'prob_f2'")
  # if (any(is.infinite(prob_f2))) stop(".scoreFrac() - Infs in 'prob_f2'")
  # if (min(prob_f2) == 0) stop(".scoreFrac() - zeroes in 'prob_f2'")
  f1 <- dnbinom(N, alpha1, prob_f1)
  f2 <- dnbinom(N, alpha2, prob_f2)
  # if (any(is.na(f1))) stop(".scoreFrac() - NAs in 'f1'")
  # if (any(is.infinite(f1))) stop(".scoreFrac() - Infs in 'f1'")
  # if (any(is.na(f2))) stop(".scoreFrac() - NAs in 'f2'")
  # if (any(is.infinite(f2))) stop(".scoreFrac() - Infs in 'f2'")

  if (zeroes) {
    num <- f1 - f2
    den <- P * f1 + (1 - P) * f2
  } else {
    f1_adjust <- 1 - dnbinom(0, alpha1, prob_f1)
    f2_adjust <- 1 - dnbinom(0, alpha2, prob_f2)
    # if (any(is.na(f1_adjust))) stop(".scoreFrac() - NAs in 'f1_adjust'")
    # if (any(is.na(f2_adjust))) stop(".scoreFrac() - NAs in 'f2_adjust'")
    # if (any(is.infinite(f1_adjust))) stop(".scoreFrac() - Infs in 'f1_adjust'")
    # if (any(is.infinite(f2_adjust))) stop(".scoreFrac() - Infs in 'f2_adjust'")
    # if (min(f1_adjust) == 0) stop(".scoreFrac() - zeroes in 'f1_adjust'")
    # if (min(f2_adjust) == 0) stop(".scoreFrac() - zeroes in 'f2_adjust'")
    f1_star <- f1 / f1_adjust
    f2_star <- f2 / f2_adjust
    # if (any(is.na(f1_star))) stop(".scoreFrac() - NAs in 'f1_star'")
    # if (any(is.infinite(f1_star))) stop(".scoreFrac() - Infs in 'f1_star'")
    # if (any(is.na(f2_star))) stop(".scoreFrac() - NAs in 'f2_star'")
    # if (any(is.infinite(f2_star))) stop(".scoreFrac() - Infs in 'f2_star'")
    num <- f1_star - f2_star
    den <- P * f1_star + (1 - P) * f2_star
  }
  # if (any(is.na(den))) stop(".scoreFrac() - NAs in 'den'")
  # if (any(is.infinite(den))) stop(".scoreFrac() - Infs in 'den'")
  # if (min(den) == 0) stop(".scoreFrac() - zeroes in 'den'")
  score <- num / den
  if (squashed) score <- W * score
  score <- sum(score)
  # if (is.na(score)) stop(".scoreFrac() - 'score' is NA")
  # if (is.infinite(score)) stop(".scoreFrac() - 'score' is Inf")
  if (is.infinite(score) || is.na(score)) stop("skip this iteration")
  score
}


.checkScore <- function(theta, N, E, W, squashed, zeroes) {
  #Calculate score vector & its Euclidean norm at theta (should be close to 0).

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
