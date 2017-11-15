## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----dataSquash example--------------------------------------------------
library(openEBGM)
data(caers)

processed <- processRaw(data = caers)
processed[1:4, 1:4]

squashed <- squashData(data = processed) #Using defaults
head(squashed)

nrow(processed)
nrow(squashed)

## ----negLLsquash example, warning = FALSE--------------------------------
theta_init <- c(alpha1 = 0.2, beta1 = 0.1, alpha2 = 2, beta2 = 4, p = 1/3)
stats::nlm(negLLsquash, p = theta_init,
           ni = squashed$N, ei = squashed$E, wi = squashed$weight, N_star = 1)

## ------------------------------------------------------------------------
theta_init <- data.frame(alpha1 = c(0.2, 0.1, 0.3),
                         beta1  = c(0.1, 0.2, 0.5),
                         alpha2 = c(2,   10,  6),
                         beta2  = c(4,   10,  6),
                         p      = c(1/3, 0.2, 0.5)
                         )

## ----autoHyper example, warning = FALSE----------------------------------
system.time(
hyper_estimates_full <- autoHyper(data = processed, theta_init = theta_init, 
                                  squashed = FALSE)
)

squashed2 <- squashData(squashed, count = 2, bin_size = 10, keep_bins = 5)
system.time(
hyper_estimates_squashed <- autoHyper(data = squashed2, theta_init = theta_init)
)

hyper_estimates_full

hyper_estimates_squashed

## ----autoHyper example with CIs, warning = FALSE-------------------------
autoHyper(squashed2, theta_init = theta_init, conf_ints = TRUE)$conf_int

## ----hyperparameter estimation example, warning = FALSE------------------
exploreHypers(data = squashed2, theta_init = theta_init, std_errors = TRUE)

