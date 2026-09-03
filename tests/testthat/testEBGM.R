#For testing if the posterior functions work properly

data.table::setDTthreads(2)  #only needed for CRAN checks

# Check if estimated hyperparameters are getting "reined in" correctly
testthat::test_that("check 'theta hat' is 'reined in' correctly", {
  x <- .reinInTheta(c(0, 0, .01, .01, 0),
    min_theta = c(.25, .25, .25, .25, .00001),
    max_theta = c(10, 10, 10, 10, .99999)
  )
  expect_all_true(x == c(0.25, 0.25, 0.25, 0.25, 0.00001))
  y <- .reinInTheta(c(100, 0, .05, 20, .000000001),
    min_theta = c(.25, .25, .25, .25, .00001),
    max_theta = c(10, 10, 10, 10, .99999)
  )
  expect_all_true(y == c(0.5, 0.25, 5, 10, 0.00001))
  z <- .reinInTheta(c(10000, 15, 5, 1e-4, 1),
    min_theta = c(.25, .25, .25, .25, .00001),
    max_theta = c(10, 10, 10, 10, .99999)
  )
  expect_all_true(z == c(10, 10, 0.5, 0.25, 0.99999))
})

data(caers)
proc_dat_ustrat <- processRaw(data = caers, stratify = FALSE)
N <- proc_dat_ustrat$N
E <- proc_dat_ustrat$E

theta <- c(1.23, 2.34, 1.4, 2.6, .1)
names(theta) <- c("alpha1", "beta1", "alpha2", "beta2", "p")

Qn_previous <- function (theta_hat, N, E) {
  # v0.9.1
  .checkInputs_Qn(theta_hat, N, E)
  prob_f1 <- theta_hat[2]/(theta_hat[2] + E)
  prob_f2 <- theta_hat[4]/(theta_hat[4] + E)
  f1_NB <- dnbinom(N, size = theta_hat[1], prob = prob_f1)
  f2_NB <- dnbinom(N, size = theta_hat[3], prob = prob_f2)
  P_hat <- theta_hat[5]
  Qn_num <- P_hat * f1_NB
  Qn_den <- (P_hat * f1_NB) + ((1 - P_hat) * f2_NB)
  Qn_num/Qn_den
}

ebgm_previous <- function (theta_hat, N, E, qn, digits = 2) {
  # v0.9.1
  .checkInputs_ebgm(theta_hat, N, E, qn, digits)
  expectation1 <- digamma(theta_hat[1] + N) - log(theta_hat[2] + 
      E)
  expectation2 <- digamma(theta_hat[3] + N) - log(theta_hat[4] + 
      E)
  exp_log_lamda <- (qn * expectation1) + ((1 - qn) * expectation2)
  EBlog2 <- exp_log_lamda/log(2)
  EBGM <- 2^EBlog2
  round(EBGM, digits)
}

quantBisect_previous <- function (
  percent, theta_hat, N, E, qn, digits = 2, limits = c(-1e+05, 1e+05), max_iter = 2000) {
  # v0.9.1
  .checkInputs_quantBisect(percent, theta_hat, N, E, qn, digits, 
      limits, max_iter)
  cdf_error <- function(guess, cut_point, qn, theta_hat, N, 
      E) {
      cdf1 <- pgamma(guess, shape = theta_hat[1] + N, rate = theta_hat[2] + 
          E)
      cdf2 <- pgamma(guess, shape = theta_hat[3] + N, rate = theta_hat[4] + 
          E)
      post_cdf <- (qn * cdf1) + ((1 - qn) * cdf2)
      post_cdf - cut_point
  }
  tol <- 0.5 * (10^(-digits))
  cut_point <- percent/100
  num_pts <- length(N)
  guess_init <- rep(1, num_pts)
  error_init <- cdf_error(guess_init, cut_point, qn, theta_hat, 
      N, E)
  is_too_big <- error_init > 0
  left <- rep(limits[1], num_pts)
  right <- rep(limits[2], num_pts)
  left <- ifelse(is_too_big, left, guess_init)
  right <- ifelse(is_too_big, guess_init, right)
  iter_count <- 1L
  while (iter_count <= max_iter) {
      mid_pt <- (left + right)/2
      error <- cdf_error(mid_pt, cut_point, qn, theta_hat, 
          N, E)
      if (max((right - left)/2) < tol) {
          quantiles <- round(mid_pt, digits)
          if (max(quantiles) == limits[2]) {
              stop("increase maximum for 'limits'")
          }
          return(quantiles)
      }
      else {
          iter_count <- iter_count + 1
          error_left <- cdf_error(left, cut_point, qn, theta_hat, 
              N, E)
          left <- ifelse(sign(error_left) == sign(error), mid_pt, 
              left)
          right <- ifelse(sign(error_left) == sign(error), 
              right, mid_pt)
      }
  }
  stop("failed to converge -- try adjusting 'limits' or 'max_iter'")
}

q_calc <- Qn(theta_hat = theta, N = N, E = E)
q_calc_previous <- Qn_previous(theta_hat = theta, N = N, E = E)
eb_scores <- ebgm(theta_hat = theta, N = N, E = E, qn = q_calc)
eb_scores_previous <- ebgm_previous(theta_hat = theta, N = N, E = E, qn = q_calc_previous)
eb_05 <- quantBisect(percent = 5, theta_hat = theta, N = N, E = E, qn = q_calc)
eb_05_previous <- quantBisect_previous(percent = 5, theta_hat = theta, N = N, E = E, qn = q_calc_previous)
eb_95 <- quantBisect(percent = 95, theta_hat = theta, N = N, E = E, qn = q_calc)
eb_95_previous <- quantBisect_previous(percent = 95, theta_hat = theta, N = N, E = E, qn = q_calc_previous)

testthat::test_that("check if Qn() is returning correct values", {
  x <- head(q_calc)
  y <- head(q_calc_previous)
  expect_all_true(abs(x - y) < 1e-5)
})

testthat::test_that("check if ebgm() is returning correct values", {
  x <- head(eb_scores)
  y <- head(eb_scores_previous)
  expect_all_true(abs(x - y) < 1e-5)
})

testthat::test_that("check if quantBisect() is returning correct values", {
  x <- head(eb_05)
  y <- head(eb_05_previous)
  expect_all_true(abs(x - y) < 1e-5)

  x <- head(eb_95)
  y <- head(eb_95_previous)
  expect_all_true(abs(x - y) < 1e-5)
})

testthat::test_that("Checking that the calculated values are logical", {
  expect_true(max(q_calc) <= 1)
  expect_true(min(q_calc) >= 0)

  expect_true(max(eb_scores) < Inf)
  expect_true(min(eb_scores) >= 0)

  expect_true(max(eb_95) < Inf)
  expect_true(min(eb_95) >= 0)
  expect_true(max(eb_05) < Inf)
  expect_true(min(eb_05) >= 0)

  expect_true(all(eb_95 >= eb_scores))
  expect_true(all(eb_05 <= eb_scores))
})

N_miss1 <- N
N_miss1[5] <- NA
N_miss2 <- N
N_miss2[3] <- NaN
N_miss3 <- N
N_miss3[4] <- Inf

testthat::test_that("do errors for Qn() get correctly printed?", {
  expect_error(
    Qn(theta_hat = theta[1:4], N = N, E = E),
    "'theta_hat' must contain 5 values",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = c(0.2, 0.2, -0.001, 3, 0.1), N = N, E = E),
    "'theta_hat' must contain only positive values",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = c(0.2, 0.2, 0.001, 3, 1.1), N = N, E = E),
    "'theta_hat[5]' (i.e., 'P') must be less than 1",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = theta, N = N[1:10], E = E),
    "'N' and 'E' must have the same length",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = theta, N = N_miss1, E = E),
    "missing or infinite values for 'N' and 'E' are not allowed",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = theta, N = N_miss2, E = E),
    "missing or infinite values for 'N' and 'E' are not allowed",
    fixed = TRUE
  )
  expect_error(
    Qn(theta_hat = theta, N = N_miss3, E = E),
    "missing or infinite values for 'N' and 'E' are not allowed",
    fixed = TRUE
  )
})

q_miss1 <- q_calc
q_miss1[5] <- NA
q_miss2 <- q_calc
q_miss2[3] <- NaN
q_miss3 <- q_calc
q_miss3[4] <- Inf

testthat::test_that("do errors for ebgm() get correctly printed?", {
  expect_error(
    ebgm(theta_hat = theta[1:4], N = N, E = E, qn = q_calc),
    "'theta_hat' must contain 5 values",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = c(0.2, 0.2, -0.001, 3, 0.1), N = N, E = E, qn = q_calc),
    "'theta_hat' must contain only positive values",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = c(0.2, 0.2, 0.001, 3, 1.1), N = N, E = E, qn = q_calc),
    "'theta_hat[5]' (i.e., 'P') must be less than 1",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = theta, N = N, E = E, qn = q_calc[1:10]),
    "'N', 'E', and 'qn' must have the same length",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = theta, N = N_miss1, E = E, qn = q_miss1),
    "missing or infinite values for 'N', 'E', and 'qn' are not allowed",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = theta, N = N_miss2, E = E, qn = q_miss2),
    "missing or infinite values for 'N', 'E', and 'qn' are not allowed",
    fixed = TRUE
  )
  expect_error(
    ebgm(theta_hat = theta, N = N_miss3, E = E, qn = q_miss3),
    "missing or infinite values for 'N', 'E', and 'qn' are not allowed",
    fixed = TRUE
  )
})

testthat::test_that("do errors for quantBisect() get correctly printed?", {
  expect_error(
    quantBisect(percent = 0.05, theta_hat = theta, N = N, E = E, qn = q_calc),
    "'percent' must be a value between 1 and 99 (inclusive)",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = c(0.2, 0.2, -0.001, 3, 0.1), N = N, E = E, qn = q_calc),
    "'theta_hat' must contain only positive values",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = c(0.2, 0.2, 0.001, 3, 1.1), N = N, E = E, qn = q_calc),
    "'theta_hat[5]' (i.e., 'P') must be less than 1",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = theta, N = N, E = E, qn = q_calc[1:10]),
    "'N', 'E', and 'qn' must have the same length",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = theta, N = N_miss3, E = E, qn = q_miss3),
    "missing or infinite values for 'N', 'E', and 'qn' are not allowed",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = theta, N = N + 10000000, E = E, qn = q_calc),
    "increase maximum for 'limits'",
    fixed = TRUE
  )
  expect_error(
    quantBisect(percent = 5, theta_hat = theta, N = N, E = E, qn = q_calc, max_iter = 2),
    "failed to converge -- try adjusting 'limits' or 'max_iter'",
    fixed = TRUE
  )
})
