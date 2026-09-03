## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----Load data----------------------------------------------------------------
library(openEBGM)
data(caers_raw)
head(caers_raw, 3)

## ----One adverse event per row------------------------------------------------
library(tidyr)
dat <- separate_longer_delim(caers_raw,
  cols = CASE_MEDDRA_PREFERRED_TERMS, delim = ","
)
dat <- as.data.frame(dat)
head(dat, 5)

## ----Rename columns-----------------------------------------------------------
dat$id <- dat$REPORT_ID
dat$var1 <- dat$PRODUCT
dat$var2 <- dat$CASE_MEDDRA_PREFERRED_TERMS

## ----Rename sex column--------------------------------------------------------
dat$strat_sex <- dat$SEX
dat$strat_sex[dat$strat_sex == ""] <- "Not reported"
table(dat$strat_sex, useNA = "always")

## ----Display prepared data----------------------------------------------------
vars <- c("id", "var1", "var2", "strat_sex")
dat <- dat[, vars]
head(dat, 5)

