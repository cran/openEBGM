#For testing if the functions that calculate the count and other reporting
#measures work properly

data.table::setDTthreads(2)  #only needed for CRAN checks

#Process some raw data
dat <- data.frame(
  var1 = c(
    "product_A", "product_B", "product_B", "product_B", "product_A", "product_A",
    "product_A", "product_A", "product_A", "product_B", "product_B", "product_A",
    "product_A", "product_B", "product_A", "product_A", "product_C", "product_C"
  ),
  var2 = c(
    "event_2", "event_1", "event_1", "event_2", "event_1", "event_1", "event_1",
    "event_1", "event_1", "event_1", "event_1", "event_1", "event_2", "event_2",
    "event_2", "event_1", "event_1", "event_3"
  ),
  strat1 = c(
    "M", "M", "F", "F", "F", "F", "M", "F", "F", "M", "M", "M", "M", "M", "F",
    "M", "M", "M"
  ),
  strat2 = c(
    "age_cat1", "age_cat2", "age_cat2", "age_cat2", "age_cat1", "age_cat1",
    "age_cat1", "age_cat2", "age_cat2", "age_cat2", "age_cat1", "age_cat2",
    "age_cat1", "age_cat1", "age_cat1", "age_cat2", "age_cat2", "age_cat1"
  ),
  stringsAsFactors = FALSE
)
dat$id <- 1:nrow(dat)

dat_no_strat <- dat[, !grepl("strat", colnames(dat))]
dat_no_id <- dat[, !grepl("id", colnames(dat))]
dat_missing1 <- dat
dat_missing1[3, "var1"] <- ""
dat_missing2 <- dat
dat_missing2[1, "var2"] <- NA
dat_missing3 <- dat
dat_missing3[2, "id"] <- NA
dat_missing4 <- dat
dat_missing4[2, "strat2"] <- ""

# From version 0.9.1 -----------------------------------------------------------

testthat::test_that("correct number of rows and columns from processRaw", {
  expect_equal(
    nrow(processRaw(dat)),
    6
  )
  expect_equal(
    nrow(suppressWarnings(processRaw(dat, stratify = TRUE))),
    6
  )
  expect_equal(
    nrow(processRaw(dat, zeroes = TRUE)),
    9
  )
  expect_equal(
    nrow(suppressWarnings(processRaw(dat, zeroes = TRUE, stratify = TRUE))),
    9
  )
  expect_equal(
    ncol(processRaw(dat)),
    6
  )
})

testthat::test_that("correct column names from processRaw", {
  expect_match(
    paste0(colnames(processRaw(dat)), collapse = ""),
    "var1var2NERRPRR"
  )
})

testthat::test_that("correct values returned from processRaw", {
  dat1a <- processRaw(dat)
  expect_equal(
    dat1a$N,
    c(7, 3, 4, 2, 1, 1)
  )
  x <- dat1a$E
  y <- c(6.6666667, 2.7777778, 4.0000000, 1.6666667, 1.3333333, 0.1111111)
  expect_all_true(
    abs(x - y) < 1e-5
  )
  expect_equal(
    dat1a$RR,
    c(1.05, 1.08, 1.00, 1.20, 0.75, 9.00)
  )
  expect_equal(
    dat1a$PRR,
    c(1.12, 1.20, 1.00, 1.33, 0.73, Inf)
  )
  dat1b <- suppressWarnings(processRaw(dat, stratify = TRUE))
  expect_equal(
    dat1b$N,
    c(7, 3, 4, 2, 1, 1)
  )
  x <- dat1b$E
  y <- c(6.5000000, 3.0000000, 4.1666667, 1.5000000, 1.3333333, 0.1666667)
  expect_all_true(
    abs(x - y) < 1e-5
  )
  expect_equal(
    dat1b$RR,
    c(1.08, 1.00, 0.96, 1.33, 0.75, 6.00)
  )
  expect_equal(
    dat1b$PRR,
    c(1.12, 1.20, 1.00, 1.33, 0.73, Inf)
  )
  dat1c <- processRaw(dat, zeroes = TRUE)
  expect_equal(
    dat1c$N,
    c(7, 3, 0, 4, 2, 0, 1, 0, 1)
  )
  x <- dat1c$E
  y <- c(6.6667, 2.7778, 0.5556, 4.0000, 1.6667, 0.3333, 1.3333, 0.5556, 0.1111)
  expect_all_true(
    abs(x - y) < 1e-3
  )
  expect_equal(
    dat1c$RR,
    c(1.05, 1.08, 0.00, 1.00, 1.20, 0.00, 0.75, 0.00, 9.00)
  )
  expect_equal(
    dat1c$PRR,
    c(1.12, 1.20, 0.00, 1.00, 1.33, 0.00, 0.73, 0.00, Inf)
  )
  dat1d <- suppressWarnings(processRaw(dat, zeroes = TRUE, stratify = TRUE))
  expect_equal(
    dat1d$N,
    c(7, 3, 0, 4, 2, 0, 1, 0, 1)
  )
  x <- dat1d$E
  y <- c(6.5000, 3.0000, 0.5000, 4.1667, 1.5000, 0.3333, 1.3333, 0.5000, 0.1667)
  expect_all_true(
    abs(x - y) < 1e-3
  )
  expect_equal(
    dat1d$RR,
    c(1.08, 1.00, 0.00, 0.96, 1.33, 0.00, 0.75, 0.00, 6.00)
  )
  expect_equal(
    dat1d$PRR,
    c(1.12, 1.20, 0.00, 1.00, 1.33, 0.00, 0.73, 0.00, Inf)
  )
})

strat_error <- paste0(
  "at least one stratification variable contains more than 1 categories --",
  "\n  did you remember to categorize stratification variables?",
  "\n  if you really need more categories, increase 'max_cats'"
)

testthat::test_that("do warnings/errors get correctly printed?", {
  expect_error(
    processRaw(as.matrix(dat)),
    "'data' must be a data frame",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat, stratify = 4),
    "'stratify', 'zeroes', & 'list_ids' must be logical values",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat, zeroes = 4),
    "'stratify', 'zeroes', & 'list_ids' must be logical values",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat, list_ids = 4),
    "'stratify', 'zeroes', & 'list_ids' must be logical values",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_no_id),
    "missing the appropriate column names (need 'id', 'var1', & 'var2')",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_missing1),
    "missing values are not allowed for 'id', 'var1', or 'var2'",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_missing2),
    "missing values are not allowed for 'id', 'var1', or 'var2'",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_missing3),
    "missing values are not allowed for 'id', 'var1', or 'var2'",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_no_strat, stratify = TRUE),
    "no stratification variables found",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat_missing4, stratify = TRUE),
    "missing values are not allowed for stratification variables",
    fixed = TRUE
  )
  expect_message(
    suppressWarnings(processRaw(dat, stratify = TRUE)),
    "stratification variables used: strat1, strat2",
    fixed = TRUE
  )
  expect_error(
    processRaw(dat, stratify = TRUE, max_cats = 1),
    strat_error,
    fixed = TRUE
  )
  expect_message(
    suppressWarnings(processRaw(dat, stratify = TRUE)),
    "there were 4 strata",
    fixed = TRUE
  )
  expect_warning(
    processRaw(dat, stratify = TRUE),
    "at least one stratum contains less than 50 unique IDs",
    fixed = TRUE
  )
})

testthat::test_that("works with a single stratification variable?", {
  expect_equal(
    suppressWarnings(max(processRaw(dat[, -4], stratify = TRUE)$E)),
    6.7532468,
    tolerance = 1e-5
  )
})

# Manual calculations ----------------------------------------------------------

results_unstr <- processRaw(dat, zeroes = TRUE, digits = 3)
results_str <- suppressWarnings(
  suppressMessages(processRaw(dat, stratify = TRUE, zeroes = TRUE, digits = 3))
)

#Manually build unstratified results
var1s <- sort(unique(dat$var1))
var2s <- sort(unique(dat$var2))
num_var1 <- length(var1s)
num_var2 <- length(var2s)
dat_man <- data.frame(
  matrix(data = NA, nrow = num_var1 * num_var2, ncol = ncol(results_unstr)),
  stringsAsFactors = FALSE
)
colnames(dat_man) <- colnames(results_unstr)
dat_man$var1 <- sort(rep(var1s, num_var2))
dat_man$var2 <- rep(sort(var2s), num_var1)
tot_sum <- length(unique(dat$id))
for (i in var1s) {
  row_sum <- nrow(dat[dat$var1 == i, ])
  for (j in var2s) {
    col_sum <- nrow(dat[dat$var2 == j, ])
    N <- nrow(dat[dat$var1 == i & dat$var2 == j, ])
    dat_man[dat_man$var1 == i & dat_man$var2 == j, "N"] <- N
    dat_man[dat_man$var1 == i & dat_man$var2 == j, "E"] <- (row_sum * col_sum) / tot_sum
    PRR_num <- N / row_sum
    PRR_den <- (col_sum - N) / (tot_sum - row_sum)
    dat_man[dat_man$var1 == i & dat_man$var2 == j, "PRR"] <- round(PRR_num / PRR_den, 3)
  }
}
dat_man$RR <- round(dat_man$N / dat_man$E, 3)
dat_man$RR <- ifelse(is.nan(dat_man$RR), 0, dat_man$RR)
dat_man$var1 <- factor(dat_man$var1, levels = sort(unique(dat_man$var1)))
dat_man$var2 <- factor(dat_man$var2, levels = sort(unique(dat_man$var2)))

dat_man_no_zero <- dat_man[dat_man$N !=0, ]
rownames(dat_man_no_zero) <- NULL

#Manually build stratified results
strat1s <- sort(unique(dat$strat1))
strat2s <- sort(unique(dat$strat2))
num_strat1 <- length(strat1s)
num_strat2 <- length(strat2s)
num_strats <- num_strat1 * num_strat2
num_rows <- num_var1 * num_var2 * num_strat1 * num_strat2
cols <- c("strat1", "strat2", "var1", "var2", "N", "E")
dat_man2 <- data.frame(
  matrix(data = NA, nrow = num_rows, ncol = length(cols)),
  stringsAsFactors = FALSE
)
colnames(dat_man2) <- cols
dat_man2$strat1 <- sort(rep(strat1s, num_strat2 * num_var1 * num_var2))
dat_man2$strat2 <- rep(sort(rep(strat2s, num_var1 * num_var2)), num_strat1)
dat_man2$var1 <- rep(sort(rep(var1s, num_var2)), num_strats)
dat_man2$var2 <- rep(var2s, num_var1 * num_strats)

for (si in strat1s) {
  for (sj in strat2s) {
    dat_strat <- dat[dat$strat1 == si & dat$strat2 == sj, ]
    strat_sum <- nrow(dat_strat)
    for (i in var1s) {
      row_sum <- nrow(dat_strat[dat_strat$var1 == i, ])
      for (j in var2s) {
        col_sum <- nrow(dat_strat[dat_strat$var2 == j, ])
        N <- nrow(dat_strat[dat_strat$var1 == i & dat_strat$var2 == j, ])
        dat_man2[dat_man2$strat1 == si & dat_man2$strat2 == sj &
                 dat_man2$var1 == i & dat_man2$var2 == j, "N"] <- N
        E <- (row_sum * col_sum) / strat_sum
        dat_man2[dat_man2$strat1 == si & dat_man2$strat2 == sj &
                 dat_man2$var1 == i & dat_man2$var2 == j, "E"] <- E
      }
    }
  }
}

dat_man2_N <- aggregate(N ~ var1 + var2, data = dat_man2, FUN = sum)
dat_man2_E <- aggregate(E ~ var1 + var2, data = dat_man2, FUN = sum)
dat_man2_all <- merge(dat_man2_N, dat_man2_E, by = c("var1", "var2"))
dat_man2_all$RR <- round(dat_man2_all$N / dat_man2_all$E, 3)
dat_man2_all$RR <- ifelse(is.nan(dat_man2_all$RR), 0, dat_man2_all$RR)
dat_man2_all$var1 <- factor(dat_man2_all$var1, levels = sort(unique(dat_man2_all$var1)))
dat_man2_all$var2 <- factor(dat_man2_all$var2, levels = sort(unique(dat_man2_all$var2)))

dat_man2_all_no_zero <- dat_man2_all[dat_man2_all$N !=0, ]
rownames(dat_man2_all_no_zero) <- NULL

testthat::test_that("does processRaw() give the expected results?", {
  expect_equal(
    dat_man,
    processRaw(dat, stratify = FALSE, zeroes = TRUE, digits = 3)
  )
  expect_equal(
    dat_man_no_zero,
    processRaw(dat, stratify = FALSE, zeroes = FALSE, digits = 3)
  )
  expect_equal(
    dat_man2_all,
    suppressWarnings(processRaw(dat, stratify = TRUE, zeroes = TRUE, digits = 3))[, 1:5]
  )
  expect_equal(
    dat_man2_all_no_zero,
    suppressWarnings(processRaw(dat, stratify = TRUE, zeroes = FALSE, digits = 3))[, 1:5]
  )
})

dat_table <- data.table::as.data.table(dat)
dat_tib <- dplyr::as_tibble(dat)

testthat::test_that("does processRaw() work with data.tables and tibble_dfs?", {
  expect_equal(
    processRaw(dat),
    processRaw(dat_table)
  )
  expect_equal(
    suppressWarnings(processRaw(dat, stratify = TRUE)),
    suppressWarnings(processRaw(dat_table, stratify = TRUE))
  )
  expect_equal(
    suppressWarnings(processRaw(dat, stratify = TRUE, zeroes = TRUE)),
    suppressWarnings(processRaw(dat_table, stratify = TRUE, zeroes = TRUE))
  )

  expect_equal(
    processRaw(dat),
    processRaw(dat_tib)
  )
  expect_equal(
    suppressWarnings(processRaw(dat, stratify = TRUE)),
    suppressWarnings(processRaw(dat_tib, stratify = TRUE))
  )
  expect_equal(
    suppressWarnings(processRaw(dat, stratify = TRUE, zeroes = TRUE)),
    suppressWarnings(processRaw(dat_tib, stratify = TRUE, zeroes = TRUE))
  )
})

#Test to see if simple example gives expected results from hand calculation
dat2 <- data.frame(
  var1 = c(
    "product_A", rep("product_B", 3), "product_C", rep("product_A", 2),
    rep("product_B", 2), "product_C"
  ),
  var2 = c(
    "event_1", rep("event_2", 2), rep("event_3", 2), "event_2", rep("event_3", 3),
    "event_1"
  ),
  strat1 = c(
    rep("Male", 5), rep("Female", 3), rep("Male", 2)
  ),
  strat2 = c(
    rep("age_cat1", 5), rep("age_cat1", 3), rep("age_cat2", 2)
  ),
  stringsAsFactors = FALSE
)
dat2$id <- 1:nrow(dat2)

var1_hand <- c(rep("product_A", 3), rep("product_B", 3), rep("product_C", 3))
var2_hand <- c(rep(c("event_1", "event_2", "event_3"), 3))
hand <- data.frame(
  var1 = factor(var1_hand, levels = sort(unique(var1_hand))),
  var2 = factor(var2_hand, levels = sort(unique(var2_hand))),
  N = c(1, 1, 1, 0, 2, 3, 1, 0, 1),
  E = c(0.6, 0.9, 1.5, 1, 1.5, 2.5, 0.4, 0.6, 1),
  RR = c(1.6667, 1.1111, 0.6667, 0, 1.3333, 1.2, 2.5, 0, 1),
  PRR = c(2.3333, 1.1667, 0.5833, 0, 2, 1.5, 4, 0, 1),
  stringsAsFactors = FALSE
)

testthat::test_that("results match hand calculations, unstratified", {
  expect_equal(
    processRaw(dat2, stratify = FALSE, zeroes = TRUE, digits = 4),
    hand
  )
  expect_equal(
    processRaw(dat2, stratify = FALSE, zeroes = FALSE, digits = 4),
    hand[hand$N != 0, ],
    ignore_attr = TRUE
  )
})

var1_hand_s <- c(rep("product_A", 3), rep("product_B", 3), rep("product_C", 3))
var2_hand_s <- c(rep(c("event_1", "event_2", "event_3"), 3))
hand_s <- data.frame(
  var1 = factor(var1_hand_s, levels = sort(unique(var1_hand_s))),
  var2 = factor(var2_hand_s, levels = sort(unique(var2_hand_s))),
  N = c(1, 1, 1, 0, 2, 3, 1, 0, 1),
  E = c(0.2, 0.4+2/3, 0.4+4/3, 1.1, 1.2+1/3, 1.7+2/3, 0.7, 0.4, 0.9),
  RR = c(5, 0.9375, 0.5769, 0, 1.3043, 1.2676, 1.4286, 0, 1.1111),
  PRR = c(2.3333, 1.1667, 0.5833, 0, 2, 1.5, 4, 0, 1),
  stringsAsFactors = FALSE
)

testthat::test_that("results match hand calculations, stratified", {
  suppressWarnings(suppressMessages(
    expect_equal(
      processRaw(dat2, stratify = TRUE, zeroes = TRUE, digits = 4),
      hand_s
    )
  ))
  suppressWarnings(suppressMessages(
    expect_equal(
      processRaw(dat2, stratify = TRUE, zeroes = FALSE, digits = 4),
      hand_s[hand_s$N != 0, ],
      ignore_attr = TRUE
    )
  ))
})
