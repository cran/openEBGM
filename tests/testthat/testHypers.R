#For testing if the functions that explore the hyperparameter space are
#functioning correctly

data.table::setDTthreads(2)  #only needed for CRAN checks

data(caers)
proc_dat_ustrat <- processRaw(caers, stratify = FALSE)
# table(proc_dat_ustrat$N)
proc_dat_ustrat_squash <- squashData(proc_dat_ustrat,
  bin_size = 100, keep_pts = 10
)
# table(proc_dat_ustrat_squash$N)
proc_dat_ustrat_squash <- squashData(proc_dat_ustrat_squash,
  count = 2, bin_size = 10, keep_pts = 10
)
# table(proc_dat_ustrat_squash$N)
proc_dat_ustrat_squash <- squashData(proc_dat_ustrat_squash,
  count = 3, bin_size = 2, keep_pts = 5, min_pts = 5
)
# table(proc_dat_ustrat_squash$N)

proc_dat_zeroes <- processRaw(caers, stratify = FALSE, zeroes = TRUE)
# table(proc_dat_zeroes$N)
proc_dat_zeroes <- proc_dat_zeroes[order(proc_dat_zeroes$N), ]
proc_dat_zeroes <- proc_dat_zeroes[2940000:nrow(proc_dat_zeroes), ]
# table(proc_dat_zeroes$N)
proc_dat_zeroes_squash <- squashData(proc_dat_zeroes,
  count = 0, bin_size = 100, keep_pts = 10, min_pts = 5
)
# table(proc_dat_zeroes_squash$N)
proc_dat_zeroes_squash <- squashData(proc_dat_zeroes_squash,
  count = 1, bin_size = 100, keep_pts = 10, min_pts = 5
)
# table(proc_dat_zeroes_squash$N)
proc_dat_zeroes_squash <- squashData(proc_dat_zeroes_squash,
  count = 2, bin_size = 10, keep_pts = 10, min_pts = 5
)
# table(proc_dat_zeroes_squash$N)
proc_dat_zeroes_squash <- squashData(proc_dat_zeroes_squash,
  count = 3, bin_size = 3, keep_pts = 10, min_pts = 5
)
# table(proc_dat_zeroes_squash$N)

theta_init <- data.frame(
  alpha1 = c(0.5, 1, 2),
  beta1 = c(0.5, 1, 2),
  alpha2 = c(2, 4, 5),
  beta2 = c(2, 4, 5),
  p = c(.1, 0.2, 0.3)
)

#Begin sanity checks -----------------------------------------------------------

#autoHyper()
hyper_no_squash <- autoHyper(proc_dat_ustrat, theta_init = theta_init,
  squashed = FALSE, max_pts = nrow(proc_dat_ustrat), conf_ints = TRUE
)
hyper_squash <- autoHyper(proc_dat_ustrat_squash, theta_init = theta_init,
  conf_ints = TRUE
)

testthat::test_that("check if autoHyper() produces similar point estimates for (un)squashed data", {
  x <- hyper_no_squash$estimates
  y <- hyper_squash$estimates
  expect_all_true(abs(x[1:2] - y[1:2]) < .05)
  expect_all_true(abs(x[3:4] - y[3:4]) < .2)
  expect_all_true(abs(x[5] - y[5]) < .025)
})

testthat::test_that("check if autoHyper() produces similar SEs for (un)squashed data", {
  x <- hyper_no_squash$conf_int$SE
  y <- hyper_squash$conf_int$SE
  expect_all_true(abs(x[1:2] - y[1:2]) < .05)
  expect_all_true(abs(x[3:4] - y[3:4]) < .075)
  expect_all_true(abs(x[5] - y[5]) < .025)
})

#exploreHypers()
hyper_squash_nlm <- exploreHypers(proc_dat_ustrat_squash, theta_init = theta_init,
  squashed = TRUE, 
  method = "nlm"
)
hyper_squash_bfgs <- exploreHypers(proc_dat_ustrat_squash, theta_init = theta_init,
  squashed = TRUE, 
  method = "bfgs"
)

testthat::test_that("check if 'nlm' & 'bfgs' method produce similar point estimates", {
  x <- hyper_squash_nlm$estimates
  y <- hyper_squash$theta_hats
  z <- hyper_squash_bfgs$estimates

  expect_equal(y[, 2:6], x[, 2:6], tolerance = .05)

  expect_all_true(abs(y[, 2] - z[, 2]) < 0.5)
  expect_all_true(abs(y[, 3] - z[, 3]) < 0.1)
  expect_equal(y[, 4:5], z[, 4:5], tolerance = .05)
  expect_all_true(abs(y[, 6] - z[, 6]) < 0.025)
})

hyper_zeroes <- exploreHypers(proc_dat_zeroes, theta_init = theta_init,
  squashed = FALSE, zeroes = TRUE, N_star = NULL, max_pts = nrow(proc_dat_zeroes),
  method = "nlm"
)
hyper_zeroes_squash <- exploreHypers(proc_dat_zeroes_squash, theta_init = theta_init,
  squashed = TRUE, zeroes = TRUE, N_star = NULL,
  method = "nlm"
)

testthat::test_that("check if exploreHypers() produces similar point estimates for zeroes/(un)squashed data", {
  expect_equal(
    hyper_zeroes$estimates[2:3],
    hyper_zeroes_squash$estimates[2:3],
    tolerance = .05
  )
  expect_equal(
    hyper_zeroes$estimates[4],
    hyper_zeroes_squash$estimates[4],
    tolerance = .3
  )
  x <- hyper_zeroes$estimates[, 5]
  y <- hyper_zeroes_squash$estimates[, 5]
  expect_all_true(abs(x - y) < 0.05)
  x <- hyper_zeroes$estimates[, 6]
  y <- hyper_zeroes_squash$estimates[, 6]
  expect_all_true(abs(x - y) < 0.025)
})

#hyperEM()
hyperEM_ests <- hyperEM(proc_dat_ustrat_squash,
  theta_init_vec = as.numeric(theta_init[1, ]),
  conf_ints = TRUE,
  # print_level = 2
  print_level = 0
)

hyperEM_ests2 <- hyperEM(proc_dat_ustrat_squash,
  theta_init_vec = as.numeric(theta_init[1, ]),
  method = "nlminb",
  profile = "parameter",
  conf_ints = TRUE,
  # print_level = 2
  print_level = 0
)

hyperEM_ests3 <- hyperEM(proc_dat_ustrat_squash,
  theta_init_vec = as.numeric(theta_init[1, ]),
  method = "nlminb",
  profile = "distribution",
  conf_ints = TRUE,
  consecutive = 25,
  # print_level = 2
  print_level = 0
)

testthat::test_that("check if hyperEM() produces similar pt ests & std errors for different method or profile", {
  x <- hyperEM_ests$estimates
  y <- hyperEM_ests2$estimates
  z <- hyperEM_ests3$estimates
  expect_all_true(abs(x - y) < .05)
  expect_all_true(abs(x - z) < .05)
  expect_all_true(abs(y - z) < .05)
  x <- hyperEM_ests$conf_int$SE
  y <- hyperEM_ests2$conf_int$SE
  z <- hyperEM_ests3$conf_int$SE
  expect_all_true(abs(x - y) < .05)
  expect_all_true(abs(x - z) < .05)
  expect_all_true(abs(y - z) < .05)
})

testthat::test_that("check if hyperEM() point estimates similar to autoHyper()", {
  x <- as.numeric(hyper_squash$estimates)
  y <- hyperEM_ests$estimates
  expect_all_true(abs(x - y) < .05)
})

# # Do not run on CRAN (takes a while) -------------------------------------------

# testthat::test_that("check if hyperEM() produces similar point estimates if (un)squashed", {
#   testthat::skip_on_cran()
#   hyperEM_ests_no_squash <- hyperEM(proc_dat_ustrat,
#     theta_init_vec = as.numeric(theta_init[1, ]),
#     squashed = FALSE,
#     conf_ints = TRUE,
#     # print_level = 2
#     print_level = 0
#   )
#   x <- hyperEM_ests_no_squash$estimates
#   y <- hyperEM_ests$estimates
#   expect_all_true(abs(x[c(1, 2, 5)] - y[c(1, 2, 5)]) < .05)
#   expect_all_true(abs(x[3:4] - y[3:4]) < .2)

#   # check if similar pt ests if unsquashed & 'nlminb' method
#   hyperEM_ests_no_squash2 <- hyperEM(proc_dat_ustrat,
#     theta_init_vec = as.numeric(theta_init[1, ]),
#     method = "nlminb",
#     profile = "parameter",
#     squashed = FALSE,
#     conf_ints = TRUE,
#     # print_level = 2
#     print_level = 0
#   )
#   suppressWarnings(
#     hyperEM_ests_no_squash3 <- hyperEM(proc_dat_ustrat,
#       theta_init_vec = as.numeric(theta_init[1, ]),
#       method = "nlminb",
#       profile = "distribution",
#       squashed = FALSE,
#       conf_ints = TRUE,
#       consecutive = 50,
#       # print_level = 2
#       print_level = 0
#     )
#   )
#   z1 <- hyperEM_ests_no_squash2$estimates
#   z2 <- hyperEM_ests_no_squash3$estimates
#   expect_all_true(abs(x - z1) < .05)
#   expect_all_true(abs(x - z2) < .05)
#   expect_all_true(abs(z1 - z2) < .05)
# })

# testthat::test_that("check if hyperEM() produces similar point estimates with zeroes if (un)squashed", {
#   testthat::skip_on_cran()
#   hyperEM_ests_zeroes <- hyperEM(proc_dat_zeroes_squash,
#     theta_init_vec = as.numeric(theta_init[1, ]),
#     zeroes = TRUE, N_star = NULL,
#     squashed = TRUE,
#     LL_tol = 1e-10,
#     # print_level = 2
#     print_level = 0
#   )
#   # takes a really long time to run
#   hyperEM_ests_zeroes_no_squash <- hyperEM(proc_dat_zeroes,
#     theta_init_vec = as.numeric(theta_init[1, ]),
#     zeroes = TRUE, N_star = NULL,
#     squashed = FALSE,
#     LL_tol = 1e-10,
#     # print_level = 2
#     print_level = 0
#   )
#   x <- hyperEM_ests_zeroes_no_squash$estimates
#   y <- hyperEM_ests_zeroes$estimates
#   expect_all_true(abs(x[1] - y[1]) < .25)  #takes too long to converge (for testing)
#   expect_all_true(abs(x[2:4] - y[2:4]) < .05)
#   expect_all_true(abs(x[5] - y[5]) < .025)

#   # check if similar pt ests if zeroes, squashed, & 'nlminb' method
#   hyperEM_ests_zeroes2 <- hyperEM(proc_dat_zeroes_squash,
#     theta_init_vec = as.numeric(theta_init[1, ]),
#     method = "nlminb",
#     profile = "parameter",
#     zeroes = TRUE, N_star = NULL,
#     squashed = TRUE,
#     LL_tol = 1e-10,
#     # print_level = 2
#     print_level = 0
#   )
#   suppressWarnings(
#     hyperEM_ests_zeroes3 <- hyperEM(proc_dat_zeroes_squash,
#       theta_init_vec = as.numeric(theta_init[1, ]),
#       method = "nlminb",
#       profile = "distribution",
#       zeroes = TRUE, N_star = NULL,
#       squashed = TRUE,
#       LL_tol = 1e-10,
#       # print_level = 2
#       print_level = 0
#     )
#   )
#   z1 <- hyperEM_ests_zeroes2$estimates
#   z2 <- hyperEM_ests_zeroes3$estimates
#   expect_all_true(abs(y[1] - z1[1]) < .65)  #takes too long to converge (for testing)
#   expect_all_true(abs(y[2:5] - z1[2:5]) < .05)
#   expect_all_true(abs(y[1] - z2[1]) < .2)
#   expect_all_true(abs(y[2:4] - z2[2:4]) < .05)
#   expect_all_true(abs(y[5] - z2[5]) < .025)
#   expect_all_true(abs(z1[1] - z2[1]) < .10)
#   expect_all_true(abs(z1[2:4] - z2[2:4]) < .05)
#   expect_all_true(abs(z1[5] - z2[5]) < .025)
# })

# Exploring lower level optimization functions ---------------------------------

# exploreHypers(proc_dat_zeroes, theta_init = theta_init,
#   squashed = FALSE, zeroes = TRUE, N_star = NULL, max_pts = nrow(proc_dat_zeroes),
#   std_errors = TRUE,
#   # method = "nlminb"
#   method = "nlm"
#   # method = "bfgs"
# )

# theta_init_transformed <- as.numeric(theta_init[1, ])
# theta_init_transformed[1:4] <- log(theta_init_transformed[1:4])
# theta_init_transformed[5] <- log(theta_init_transformed[5] / (1 - theta_init_transformed[5]))

# #negLLzero()
# x <- nlminb(start = theta_init_transformed, objective = negLLzero, 
#   N = proc_dat_zeroes$N, E = proc_dat_zeroes$E, transformed = TRUE,
#   control = list(
#     eval.max = 1000, iter.max = 1000,
#     trace = 2,
#     abs.tol = 1e-20, xf.tol = 5e-14
#   )
# )
# c(exp(x$par[1:4]), plogis(x$par[5]))

# z <- nlm(f = negLLzero, p = theta_init_transformed,
#   N = proc_dat_zeroes$N, E = proc_dat_zeroes$E,
#   iterlim = 1000,
#   transformed = TRUE,
#   typsize = c(log(rep(1, 4)), log(.2/.8)),
#   fscale = 1000,
#   stepmax = 100,
#   print.level = 2
# )
# c(exp(z$estimate[1:4]), plogis(z$estimate[5]))

# w <- optim(par = theta_init_transformed, fn = negLLzero,
#   N = proc_dat_zeroes$N, E = proc_dat_zeroes$E,
#   transformed = TRUE,
#   method = "BFGS",
#   hessian = TRUE,
#   control = list(
#     trace = 10,
#     fnscale = 1000,
#     # parscale = c(rep(1000, 4), 1000),
#     parscale = rep(1, 5),
#     ndeps = rep(1e-5, 5),
#     maxit = 1000,
#     REPORT = 1
#   )
# )
# c(exp(w$par[1:4]), plogis(w$par[5]))

# #negLLzeroSquash()
# x1 <- nlminb(start = theta_init_transformed, objective = negLLzeroSquash,
#   ni = proc_dat_zeroes_squash$N, ei = proc_dat_zeroes_squash$E, wi = proc_dat_zeroes_squash$weight,
#   transformed = TRUE,
#   control = list(
#     eval.max = 1000, iter.max = 1000,
#     trace = 2,
#     abs.tol = 1e-20,
#     xf.tol = 5e-14
#   )
# )
# c(exp(x1$par[1:4]), plogis(x1$par[5]))

# z1 <- nlm(f = negLLzeroSquash, p = theta_init_transformed,
#   ni = proc_dat_zeroes_squash$N, ei = proc_dat_zeroes_squash$E, wi = proc_dat_zeroes_squash$weight,
#   iterlim = 1000,
#   transformed = TRUE,
#   typsize = c(log(rep(1, 4)), log(.2/.8)),
#   fscale = 1000,
#   stepmax = 100,
#   print.level = 2
# )
# c(exp(z1$estimate[1:4]), plogis(z1$estimate[5]))

# w1 <- optim(par = theta_init_transformed, fn = negLLzeroSquash,
#   ni = proc_dat_zeroes_squash$N, ei = proc_dat_zeroes_squash$E, wi = proc_dat_zeroes_squash$weight,
#   transformed = TRUE,
#   method = "BFGS",
#   hessian = TRUE,
#   control = list(
#     trace = 10,
#     fnscale = 1000,
#     # parscale = c(rep(1000, 4), 1000),
#     parscale = rep(1, 5),
#     ndeps = rep(1e-5, 5),
#     maxit = 1000,
#     REPORT = 1
#   )
# )
# c(exp(w1$par[1:4]), plogis(w1$par[5]))

# #negLL(theta, N, E, N_star = 1)
# x2 <- nlminb(start = theta_init_transformed, objective = negLL,
#   N = proc_dat_ustrat$N, E = proc_dat_ustrat$E, transformed = TRUE,
#   control = list(
#     eval.max = 1000, iter.max = 1000,
#     trace = 2,
#     abs.tol = 1e-20, xf.tol = 5e-14
#   )
# )
# c(exp(x2$par[1:4]), plogis(x2$par[5]))

# z2 <- nlm(f = negLL, p = theta_init_transformed,
#   N = proc_dat_ustrat$N, E = proc_dat_ustrat$E,
#   iterlim = 1000,
#   transformed = TRUE,
#   typsize = c(log(rep(1, 4)), log(.2/.8)),
#   fscale = 1000,
#   stepmax = 100,
#   print.level = 2
# )
# c(exp(z2$estimate[1:4]), plogis(z2$estimate[5]))

# w2 <- optim(par = theta_init_transformed, fn = negLL,
#   N = proc_dat_ustrat$N, E = proc_dat_ustrat$E,
#   transformed = TRUE,
#   method = "BFGS",
#   hessian = TRUE,
#   control = list(
#     trace = 10,
#     fnscale = 1000,
#     # parscale = c(rep(1000, 4), 1000),
#     parscale = rep(1, 5),
#     ndeps = rep(1e-5, 5),
#     maxit = 1000,
#     REPORT = 1
#   )
# )
# c(exp(w2$par[1:4]), plogis(w2$par[5]))

# #negLLSquash()
# x3 <- nlminb(start = theta_init_transformed, objective = negLLsquash,
#   ni = proc_dat_ustrat_squash$N, ei = proc_dat_ustrat_squash$E, wi = proc_dat_ustrat_squash$weight,
#   transformed = TRUE,
#   control = list(
#     eval.max = 1000, iter.max = 1000,
#     trace = 2,
#     abs.tol = 1e-20,
#     xf.tol = 5e-14
#   )
# )
# c(exp(x3$par[1:4]), plogis(x3$par[5]))

# z3 <- nlm(f = negLLsquash, p = theta_init_transformed,
#   ni = proc_dat_ustrat_squash$N, ei = proc_dat_ustrat_squash$E, wi = proc_dat_ustrat_squash$weight,
#   iterlim = 1000,
#   hessian = TRUE,
#   transformed = TRUE,
#   typsize = c(log(rep(1, 4)), log(.2/.8)),
#   fscale = 1000,
#   stepmax = 100,
#   print.level = 2
# )
# z3
# c(exp(z3$estimate[1:4]), plogis(z3$estimate[5]))

# w3 <- optim(par = theta_init_transformed, fn = negLLsquash,
#   ni = proc_dat_ustrat_squash$N, ei = proc_dat_ustrat_squash$E, wi = proc_dat_ustrat_squash$weight,
#   transformed = TRUE,
#   method = "BFGS",
#   hessian = TRUE,
#   control = list(
#     trace = 10,
#     fnscale = 1000,
#     # parscale = c(rep(1000, 4), 1000),
#     parscale = rep(1, 5),
#     ndeps = rep(1e-5, 5),
#     maxit = 1000,
#     REPORT = 1
#   )
# )
# c(exp(w3$par[1:4]), plogis(w3$par[5]))
