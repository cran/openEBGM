context("Data Squashing")
#For testing if the function that squashes data works properly

#Squash some processed data
set.seed(483726)
dat <- data.frame(var1 = letters[1:26], var2 = LETTERS[1:26],
                  N = c(rep(0, 11), rep(1, 10), rep(2, 4), rep(3, 1)),
                  E = abs(c(rnorm(11, 0), rnorm(10, 1), rnorm(4, 2),
                            rnorm(1, 3))))

dat_no_N <- dat[, !grepl("N", names(dat))]
dat_missing1 <- dat
dat_missing1[3, "N"] <- NA
dat_missing2 <- dat
dat_missing2[1, "E"] <- ""

zeroes <- squashData(dat, count = 0, bin_size = 3, keep_bins = 1,
                     min_bin = 2, min_pts = 2)
ones   <- squashData(zeroes, count = 1, bin_size = 2, keep_bins = 1,
                     min_bin = 2, min_pts = 2)
twos   <- squashData(ones, count = 2, bin_size = 2, keep_bins = 1,
                     min_bin = 2, min_pts = 2)

testthat::test_that("correct number of rows and columns from squashData()", {
  expect_equal(nrow(zeroes), 21)
  expect_equal(nrow(ones), 17)
  expect_equal(nrow(twos), 16)
  expect_equal(ncol(twos), 3)
})

testthat::test_that("do the weights make sense?", {
  expect_equal(squashData(dat, count = 0, bin_size = 3, keep_bins = 0,
                          min_bin = 2, min_pts = 2)$weight,
               c(3,3,3,2, rep(1, 15)))
  expect_equal(squashData(dat, count = 0, bin_size = 3, keep_bins = 1,
                          min_bin = 2, min_pts = 2)$weight,
               c(3,3,2,1,1,1, rep(1, 15)))
  expect_equal(squashData(dat, count = 0, bin_size = 3, keep_bins = 2,
                          min_bin = 2, min_pts = 2)$weight,
               c(3,2,1,1,1,1,1,1, rep(1, 15)))
})

testthat::test_that("correct column names from squashData()", {
  expect_match(paste0(names(ones), collapse = ""), "NEweight")
})

testthat::test_that("values returned from squashData() make sense", {
  expect_true(max(twos$E) < Inf)
  expect_true(min(twos$E) >= 0)
  expect_equal(max(twos$N), 3)
  expect_equal(min(twos$N), 0)
  expect_equal(max(twos$weight), 3)
  expect_equal(min(twos$weight), 1)
})

testthat::test_that("do warnings/errors get correctly printed?", {
  expect_error(squashData(as.matrix(dat)),
               "'data' must be a data frame",
               fixed = TRUE)
  expect_error(squashData(dat_no_N),
               "missing the appropriate column names (need 'N' and 'E')",
               fixed = TRUE)
  expect_error(squashData(dat_missing1),
               "missing values for 'N' or 'E' are not allowed",
               fixed = TRUE)
  expect_error(squashData(dat_missing2),
               "missing values for 'N' or 'E' are not allowed",
               fixed = TRUE)
  expect_error(squashData(dat, count = -1),
               "'count' must be non-negative",
               fixed = TRUE)
  expect_error(squashData(dat, count = -1.4),
               "'count' must be non-negative",
               fixed = TRUE)
  expect_error(squashData(dat, bin_size = 1),
               "'bin_size' must be >= 2",
               fixed = TRUE)
  expect_error(squashData(dat, keep_bins = -1),
               "'keep_bins' must be non-negative",
               fixed = TRUE)
  expect_error(squashData(ones, count = 1),
               "this data set has already been squashed for this count size",
               fixed = TRUE)
  expect_error(squashData(dat, count = 3, min_pts = 2),
               "not enough points for squashing",
               fixed = TRUE)
  expect_error(squashData(dat, count = 0, bin_size = 50, keep_bins = 0,
                          min_bin = 2, min_pts = 2),
               "not enough bins -- reduce 'bin_size'",
               fixed = TRUE)
  expect_error(squashData(dat, count = 0, bin_size = 2, keep_bins = 6,
                          min_bin = 2, min_pts = 2),
               "reduce 'bin_size' or 'keep_bins'",
               fixed = TRUE)
})

dat_E     <- dat[order(dat$N, dat$E), "E"]
dat_E1    <- squashData(dat, count = 0, bin_size = 5, keep_bins = 0,
                        min_bin = 2, min_pts = 2)
dat_E1    <- dat_E1[dat_E1$N == 0, "E"]

dat_E2    <- squashData(dat, count = 0, bin_size = 5, keep_bins = 1,
                        min_bin = 2, min_pts = 2)
dat_E2    <- dat_E2[dat_E2$N == 0, "E"]

dat_E3    <- squashData(dat, count = 1, bin_size = 2, keep_bins = 0,
                        min_bin = 2, min_pts = 2)
dat_E3    <- dat_E3[dat_E3$N == 1, "E"]

testthat::test_that("are the Es averaged correctly?", {
  expect_equal(dat_E1,
               c(mean(dat_E[1:5]), mean(dat_E[6:10]), dat_E[11]))
  expect_equal(dat_E2,
               c(mean(dat_E[1:5]), dat_E[6:11]))
  expect_equal(dat_E3,
               c(mean(dat_E[12:13]), mean(dat_E[14:15]), mean(dat_E[16:17]),
                 mean(dat_E[18:19]), mean(dat_E[20:21])))
})
