## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----dataSquash example-------------------------------------------------------
library(openEBGM)
data(caers)

processed <- processRaw(caers)
processed[1:4, 1:4]

squashed <- squashData(processed) #Using defaults
head(squashed)

nrow(processed)
nrow(squashed)

## ----iterative squashing example----------------------------------------------
squash1 <- squashData(processed)
squash2 <- squashData(squash1, count = 2, bin_size = 8)

## ----autoSquash example-------------------------------------------------------
squash3 <- autoSquash(processed)
ftable(squash3[, c("N", "weight")])

## ----negLLsquash example------------------------------------------------------
theta_init <- c(alpha1 = 1, beta1 = 1, alpha2 = 2, beta2 = 2, p = .2)
theta_hats <- stats::nlm(negLLsquash, p = theta_init,
  ni = squashed$N, ei = squashed$E, wi = squashed$weight, N_star = 1
)$estimate
theta_hats

## ----negLLsquash transformed--------------------------------------------------
theta_init[1:4] <- log(theta_init[1:4])
theta_init[5] <- log(theta_init[5] / (1 - theta_init[5]))
theta_hats_transformed <- stats::nlm(negLLsquash, p = theta_init,
  ni = squashed$N, ei = squashed$E, wi = squashed$weight, N_star = 1,
  transformed = TRUE
)$estimate
theta_hats <- theta_hats_transformed
theta_hats[1:4] <- exp(theta_hats[1:4])
theta_hats[5] <- stats::plogis(theta_hats[5])
theta_hats

## -----------------------------------------------------------------------------
theta_init <- data.frame(
  alpha1 = c(1, 2, 3),
  beta1 = c(1, 2, 3),
  alpha2 = c(2, 4, 5),
  beta2 = c(2, 4, 5),
  p = c(.1, 0.2, 0.3)
)

## ----autoHyper example--------------------------------------------------------
system.time(
  hyper_estimates_full <- autoHyper(data = processed, theta_init = theta_init, 
    squashed = FALSE
  )
)

squashed <- squashData(processed, count = 1, bin_size = 25, keep_pts = 10)
squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 10)
system.time(
  hyper_estimates_squashed <- autoHyper(data = squashed, theta_init = theta_init)
)

hyper_estimates_full
hyper_estimates_squashed

## ----autoHyper example with CIs-----------------------------------------------
autoHyper(squashed, theta_init = theta_init, conf_ints = TRUE)$conf_int

## ----hyperparameter estimation example----------------------------------------
exploreHypers(data = squashed, theta_init = theta_init, std_errors = TRUE)

## ----hyperEM example----------------------------------------------------------
data(caers)
processed <- processRaw(caers)
squashed <- squashData(processed, count = 1, bin_size = 25, keep_pts = 10)
squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 10)
hyperEM_ests <- hyperEM(squashed, theta_init_vec = c(1, 1, 2, 2, .1),
  conf_ints = TRUE, LL_tol = 1e-6, track = TRUE
)
str(hyperEM_ests)

## ----hyperEM plotting example, fig.width = 7, fig.height = 10-----------------
library(ggplot2)
library(tidyr)
pdat <- pivot_longer(hyperEM_ests$tracking,
  cols = -iter, names_to = "metric", values_to = "value"
)
pdat$metric <- factor(pdat$metric, levels = unique(pdat$metric), ordered = TRUE)
ggplot(pdat, aes(x = iter, y = value)) +
  geom_line(linewidth = 1.1, col = "blue") +
  facet_grid(metric ~ ., scales = "free") +
  ggtitle("Convergence Assessment",
    subtitle = "Dashed red line indicates accelerated estimate"
  ) +
  labs(x = "Iteration Count", y = "Estimate") +
  geom_vline(
    xintercept = c(100, 200), linewidth = 1, linetype = 2, col = "red"
  )

## ----DEoptim example, message = FALSE-----------------------------------------
library(DEoptim)
set.seed(123456)
theta_hat <- DEoptim(negLLsquash,
  lower = rep(1e-03, 5),
  upper = c(rep(5, 4), .999),
  control = DEoptim.control(
    itermax = 2000,
    reltol = 1e-04,
    steptol = 200,
    NP = 100,
    CR = 0.85,
    F = 0.75,
    trace = 25
  ),
  ni = squashed$N, ei = squashed$E, wi = squashed$weight
)
(theta_hat <- as.numeric(theta_hat$optim$bestmem))

