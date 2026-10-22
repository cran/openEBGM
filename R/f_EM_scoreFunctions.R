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
    #   P_prime = logit(P)
    alpha <- exp(alpha)
    beta  <- exp(beta)
    P_vec <- stats::plogis(P_vec)  #inverse logit
  }
  alpha <- max(alpha, 1e-100)
  beta <- max(beta, 1e-100)
  alpha_base_log <- log1p(exp(log(E) - log(beta)))  #log(1 + E / beta)
  core <- digamma(alpha + N) - digamma(alpha) - alpha_base_log
  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_base_log
    zero_adjust2 <- expm1(alpha * alpha_base_log)  #alpha_base ^ alpha - 1
    zero_adjust <- exp(log(zero_adjust1) - log(zero_adjust2))  #zero_adjust1/zero_adjust2
    score <- core - zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  if (.isMissing_num(score)) stop("skip this iteration")
  score
}

# .scoreAlpha(alpha = 1, beta = 2, P_vec = c(.1, .2),
#   N = 1:2, E = c(.1, 1.3), W = 3:4, squashed = TRUE, zeroes = FALSE
# )

.scoreBeta <- function(beta, alpha, P_vec, N, E, W, squashed, zeroes,
                       reparam = TRUE) {
  # Partial derivative of logL w.r.t. beta
  # Caution: Default is to input reparameterized hyperparameters.
  if (reparam) {
    # Input reparameterization:
    #   a_prime = log(alpha); b_prime = log(beta)
    #   P_prime = logit(P)
    alpha <- exp(alpha)
    beta  <- exp(beta)
    P_vec <- stats::plogis(P_vec)
  }
  alpha <- max(alpha, 1e-100)
  beta <- max(beta, 1e-100)
  log_beta <- log(beta)
  alpha_base_log <- log1p(exp(log(E) - log_beta))  #log(1 + E / beta)
  alpha_E <- alpha * E  #alpha * E
  # core = (alpha_E - N * beta) / (beta * (beta + E))
  core <- (alpha_E - N * beta) / exp(log_beta + log(beta + E))
  if (zeroes) {
    score <- core
  } else {
    zero_adjust1 <- alpha_E
    # zero_adjust2 = beta ^ 2 * (alpha_base ^ (alpha + 1) - alpha_base)
    zero_adjust2_log <- 2 * log_beta + alpha_base_log + log(expm1(alpha * alpha_base_log))
    zero_adjust  <- exp(log(zero_adjust1) - zero_adjust2_log)  #zero_adjust1/zero_adjust2
    score <- core + zero_adjust
  }
  if (squashed) score <- W * score
  score <- sum(P_vec * score)
  if (.isMissing_num(score)) stop("skip this iteration")
  score
}

# .scoreBeta(beta = 2, alpha = 1, P_vec = c(.1, .2),
#   N = 1:2, E = c(.1, 1.3), W = 3:4, squashed = TRUE, zeroes = FALSE
# )

.scoreFrac <- function(P, alpha1, beta1, alpha2, beta2, N, E, W, squashed,
                       zeroes, reparam = TRUE) {
  # Partial derivative of logL w.r.t. mixing fraction (P)
  # Caution: Default is to input reparameterized hyperparameters.
  if (reparam) {
    # Input reparameterization:
    #   a_prime = log(alpha); b_prime = log(beta)
    #   P_prime = logit(P)
    alpha1_log <- alpha1; beta1_log <- beta1
    alpha2_log <- alpha2; beta2_log <- beta2
    alpha1 <- exp(alpha1); beta1 <- exp(beta1)
    alpha2 <- exp(alpha2); beta2 <- exp(beta2)
    P <- stats::plogis(P)
  }
  alpha1 <- max(alpha1, 1e-100)
  beta1 <- max(beta1, 1e-100)
  alpha2 <- max(alpha2, 1e-100)
  beta2 <- max(beta2, 1e-100)
  P <- max(P, 1e-10)
  P <- min(P, 1 - 1e-10)
  prob_f1 <- exp(log(beta1) - log(beta1 + E))  #beta1 / (beta1 + E)
  prob_f2 <- exp(log(beta2) - log(beta2 + E))
  prob_f1[prob_f1 < 1e-10] <- 1e-10
  prob_f2[prob_f2 < 1e-10] <- 1e-10
  max_prob <- 1 - 1e-10
  prob_f1[prob_f1 > max_prob] <- max_prob
  prob_f2[prob_f2 > max_prob] <- max_prob
  log_f1 <- dnbinom(N, alpha1, prob_f1, log = TRUE)
  log_f2 <- dnbinom(N, alpha2, prob_f2, log = TRUE)
  if (zeroes) {
    num <- exp(log_f1) - exp(log_f2)  #f1 - f2
    den <- exp(log(P) + log_f1) + exp(log(1 - P) + log_f2)  #P*f1 + (1-P)*f2
  } else {
    log_f1_adjust <- pnbinom(0, alpha1, prob_f1, lower.tail = FALSE, log.p = TRUE)
    log_f2_adjust <- pnbinom(0, alpha2, prob_f2, lower.tail = FALSE, log.p = TRUE)
    log_f1_star <- log_f1 - log_f1_adjust
    log_f2_star <- log_f2 - log_f2_adjust
    num <- exp(log_f1_star) - exp(log_f2_star)
    den <- exp(log(P) + log_f1_star) + exp(log(1 - P) + log_f2_star)
  }
  score <- num / den
  den_not_zero <- den != 0
  score <- score[den_not_zero]
  if (squashed) {
    W <- W[den_not_zero]
    score <- W * score
  }
  score <- sum(score)
  if (.isMissing_num(score)) stop("skip this iteration")
  score
}

# .scoreFrac(P = .1, alpha1 = 1, beta1 = 1, alpha2 = 2, beta2 = 2,
#   N = 1:2, E = c(.1, 1.3), W = 3:4, squashed = TRUE, zeroes = FALSE
# )

.checkScore <- function(theta, N, E, W, squashed, zeroes) {
  # Calculate score vector & its Euclidean norm at theta (should be close to 0).
  log_marg_dens <- .marginalDensity(theta, N, E, zeroes)
  PQ <- .updateProbs(theta, log_marg_dens$log_f1, log_marg_dens$log_f2)
  score_a1 <- .scoreAlpha(
    theta[1], theta[2], PQ$P_vec, N, E, W, squashed, zeroes, reparam = FALSE
  )
  score_b1 <- .scoreBeta(
    theta[2], theta[1], PQ$P_vec, N, E, W, squashed, zeroes, reparam = FALSE
  )
  score_a2 <- .scoreAlpha(
    theta[3], theta[4], PQ$Q_vec, N, E, W, squashed, zeroes, reparam = FALSE
  )
  score_b2 <- .scoreBeta(
    theta[4], theta[3], PQ$Q_vec, N, E, W, squashed, zeroes, reparam = FALSE
  )
  score_P <- .scoreFrac(
    theta[5], theta[1], theta[2], theta[3], theta[4], N, E, W, squashed, zeroes,
    reparam = FALSE
  )
  score_vec <- c(score_a1, score_b1, score_a2, score_b2, score_P)
  score_norm <- sqrt(sum(score_vec ^ 2))
  list(score_vec = score_vec, score_norm = score_norm)
}
