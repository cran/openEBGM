## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----hyperparameter estimation------------------------------------------------
library(openEBGM)
data(caers)  #subset of publicly available HFCS (formerly CAERS) data

processed <- processRaw(caers, stratify = FALSE, zeroes = FALSE)
squashed <- squashData(processed)
squashed <- squashData(squashed, count = 2, bin_size = 10, keep_pts = 50)
theta_init <- data.frame(
  alpha1 = c(1, 2, 3),
  beta1 = c(1, 2, 3),
  alpha2 = c(2, 4, 5),
  beta2 = c(2, 4, 5),
  p = c(.1, 0.2, 0.3)
)
theta_init
hyper_estimates <- autoHyper(squashed, theta_init = theta_init)
(theta_hat <- hyper_estimates$estimates)

## -----------------------------------------------------------------------------
qn <- Qn(theta_hat, N = processed$N, E = processed$E)
head(qn)
identical(length(qn), nrow(processed))
summary(qn)

## -----------------------------------------------------------------------------
processed$ebgm <- ebgm(theta_hat, N = processed$N, E = processed$E, qn  = qn)
head(processed)

## -----------------------------------------------------------------------------
processed$QUANT_05 <- quantBisect(5, theta_hat = theta_hat,
  N = processed$N, E = processed$E, qn = qn
)
processed$QUANT_95 <- quantBisect(95, theta_hat = theta_hat,
  N = processed$N, E = processed$E, qn = qn
)
head(processed)

## -----------------------------------------------------------------------------
suspicious <- processed[processed$QUANT_05 >= 2, ]
nrow(suspicious); nrow(processed); nrow(suspicious)/nrow(processed)

## -----------------------------------------------------------------------------
suspicious <- suspicious[order(suspicious$QUANT_05, decreasing = TRUE),
                         c("var1", "var2", "N", "E", "QUANT_05", "ebgm", 
                           "QUANT_95")]
head(suspicious, 5)

tabbed <- table(suspicious$var1)
head(tabbed[order(tabbed, decreasing = TRUE)])

