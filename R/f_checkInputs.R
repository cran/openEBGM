#Helper functions to check if actual argument values make sense

#processRaw() ------------------------------------------------------------------
.checkInputs_processRaw <- function(data = data, stratify = stratify,
                                    zeroes = zeroes) {
  #Sanity checks for actual arguments

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.logical(stratify) || !is.logical(zeroes)) {
    stop("'stratify' and 'zeroes' must be logical values")
  }
  if (!all(c("id", "var1", "var2") %in% colnames(data))) {
    stop("missing the appropriate column names (need 'id', 'var1', & 'var2')")
  }
  if (any(apply(data[, c("id", "var1", "var2")],
                MARGIN = 2, FUN = .isMissing_str))) {
    stop("missing values are not allowed for 'id', 'var1', or 'var2'")
  }
}

.checkStrata_processRaw <- function(data = data, max_cats = max_cats) {
  #Sanity checks for actual arguments. Also prints messages and
  #adds 'stratum' column to 'data'

  strat_vars <- colnames(data)[grepl("strat", colnames(data))]
  if (length(strat_vars) == 0) {
    stop("no stratification variables found")
  }
  if (any(apply(data[, strat_vars, with = FALSE],
                MARGIN = 2, FUN = .isMissing_str))) {
    stop("missing values are not allowed for stratification variables")
  }

  message(paste("stratification variables used:",
                paste(strat_vars, collapse = ", ")), sep = " ")
  if (max(apply(data[, strat_vars, with = FALSE], MARGIN = 2,
                FUN = .countUnique)) > max_cats) {
    error1 <- "at least one stratification variable contains more than "
    error2 <- "\n  did you remember to categorize stratification variables?"
    error3 <- "\n  if you really need more categories, increase 'max_cats'"
    stop(paste0(error1, max_cats, " categories --", error2, error3))
  }

  data$stratum <- do.call(paste, c(data[, strat_vars, with = FALSE], sep = "-"))
  message(paste("there were", .countUnique(data$stratum), "strata\n"))
  if (min(data[, .countUnique(id), by = .(stratum)]$V1) < 50) {
    warning("at least one stratum contains less than 50 unique IDs")
  }

  data
}

#squashData() ------------------------------------------------------------------
.checkInputs_squashData <- function(data = data, count = count,
                                    bin_size = bin_size, keep_bins = keep_bins,
                                    min_bin = min_bin, min_pts = min_pts) {
  #Sanity checks for actual arguments. Also coerces df to data table.
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!all(c("N", "E") %in% colnames(data))) {
    stop("missing the appropriate column names (need 'N' and 'E')")
  }
  count     <- as.integer(count)
  bin_size  <- as.integer(bin_size)
  keep_bins <- as.integer(keep_bins)
  min_bin   <- as.integer(min_bin)
  min_pts   <- as.integer(min_pts)
  data <- data.table::as.data.table(data)
  data <- data[, N := as.integer(N)]
  data <- data[, E := as.numeric(E)]

  if (.isMissing_num(data[, N]) | .isMissing_num(data[, E])) {
    stop("missing values for 'N' or 'E' are not allowed")
  }
  if (count < 0) {
    stop("'count' must be non-negative")
  }
  if (bin_size < 2) {
    stop("'bin_size' must be >= 2")
  }
  if (keep_bins < 0) {
    stop("'keep_bins' must be non-negative")
  }

  #In case we want to run the function iteratively (count = 1, then 2, etc.)
  weight_exists <- any(grepl("weight", colnames(data)))
  if (weight_exists) {
    if (max(data[N == count, weight], na.rm = TRUE) > 1) {
      stop("this data set has already been squashed for this count size")
    }
  } else {
    data <- data[, weight := 1L]
  }

  num_pts <- nrow(data[N == count, ])
  if (num_pts < min_pts) {
    stop("not enough points for squashing")
  }
  if (num_pts / bin_size < min_bin) {
    stop("not enough bins -- reduce 'bin_size'")
  }

  data
}

#Likelihood functions ----------------------------------------------------------
.checkInputs_negLLzero <- function(theta = theta, N = N, E = E) {
  #Sanity checks for actual arguments
  if (length(N) != length(E)) {
    stop("'N' & 'E' must be the same length")
  }
  if (min(N) != 0) {
    stop("no zero counts found -- please use another likelihood function")
  }
}

.checkInputs_negLLzeroSquash <- function(theta = theta, ni = ni, ei = ei,
                                         wi = wi) {
  ni_len <- length(ni)
  if (ni_len != length(ei) || ni_len != length(wi)) {
    stop("'ni', 'ei', & 'wi' must be the same length")
  }
  if (min(ni) != 0) {
    stop("no zero counts found -- please use another likelihood function")
  }
}

.checkInputs_negLL <- function(theta = theta, N = N, E = E, N_star = N_star) {
  if (length(N) != length(E)) {
    stop("'N' & 'E' must be the same length")
  }
  if (min(N) != N_star) {
    stop("'N_star' does not agree with 'N'")
  }
}

.checkInputs_negLLsquash <- function(theta = theta, ni = ni, ei = ei, wi = wi,
                                     N_star = N_star) {
  ni_len <- length(ni)
  if (ni_len != length(ei) || ni_len != length(wi)) {
    stop("'ni', 'ei', & 'wi' must be the same length")
  }
  if (min(ni) != N_star) {
    stop("'N_star' does not agree with 'ni'")
  }
}

#exploreHypers() ---------------------------------------------------------------
.checkInputs_exploreHypers <- function(data = data, theta_init = theta_init,
                                       squashed = squashed, zeroes = zeroes,
                                       N_star = N_star, method = method,
                                       param_limit = param_limit,
                                       max_pts = max_pts) {
  #Sanity checks for actual arguments
  if (!all(c("N", "E") %in% colnames(data))) {
    stop("missing the appropriate column names (need 'N' and 'E')")
  }
  if (.isMissing_num(data$N) || .isMissing_num(data$E)){
    stop("missing or infinite values for 'N' and 'E' are not allowed")
  }
  if (!is.logical(squashed) || !is.logical(zeroes)) {
    stop("'squashed' and 'zeroes' must be logical values")
  }
  if (zeroes && !is.null(N_star)) {
    stop("if zeroes are used, 'N_star' should be NULL")
  }
  if (!zeroes && is.null(N_star)) {
    stop("if zeroes are not used, 'N_star' must be specified")
  }
  if (zeroes && min(data$N) != 0) {
    stop("no zero counts found")
  }
  if (!zeroes && min(data$N) == 0) {
    stop("zero counts found")
  }
  if (!is.null(N_star)) {
    if (!is.numeric(N_star) || N_star < 1) {
      stop("'N_star' must be >= 1 or NULL")
    }
    N_star <- as.integer(N_star)
    if (N_star != min(as.integer(names(table(data$N))))) {
      stop("'N_star' does not agree with the data set")
    }
  }
  if (nrow(data) > max_pts) {
    error1 <- "there are more than "
    error2 <- "\n  either squash data, run individual optimizations,"
    error3 <- "\n  or increase 'max_pts' (caution: long run times)"
    stop(paste0(error1, max_pts, " points --", error2, error3))
  }
  if (ncol(theta_init) != 5) {
    stop("'theta_init' must contain 5 columns")
  }
  if (min(theta_init) <= 0) {
    stop("'theta_init' must contain nonnegative values")
  }
  if (max(theta_init[5]) >= 1) {
    stop("'theta[5]' (i.e., 'P') must be <1")
  }
  if (squashed && !"weight" %in% colnames(data)) {
    stop("'weight' column is missing -- are these data really squashed?")
  }
  if (!squashed && "weight" %in% colnames(data)) {
    stop("'weight' column was found -- were these data squashed?")
  }
}

#autoHyper() -------------------------------------------------------------------
.checkInputs_autoHyper <- function(data = data, theta_init = theta_init,
                                   squashed = squashed, zeroes = zeroes,
                                   N_star = N_star, tol = tol,
                                   min_conv = min_conv,
                                   param_limit = param_limit,
                                   max_pts = max_pts) {
  #Sanity checks for actual arguments
  if (length(tol) != 5) {
    stop("'tol' must have a length of 5")
  }
  if (min(tol) <= 0) {
    stop("'tol' must have only positive values")
  }
  if (min_conv < 1 || min_conv > nrow(theta_init) - 1) {
    stop(paste0("'min_conv' must be a positive number not more than one ",
                "\n  less than the number of rows in 'theta_init'"))
  }
}

#Qn() --------------------------------------------------------------------------
.checkInputs_Qn <- function(theta_hat = theta_hat, N = N, E = E) {
  #Sanity checks for actual arguments
  if (length(theta_hat) != 5) {
    stop("'theta_hat' must contain 5 values")
  }
  if (min(theta_hat) <= 0) {
    stop("'theta_hat' must contain only positive values")
  }
  if (theta_hat[5] >= 1) {
    stop("'theta_hat[5]' (i.e., 'P') must be less than 1")
  }
  if (length(N) != length(E)) {
    stop("'N' and 'E' must have the same length")
  }
  if (.isMissing_num(N) | .isMissing_num(E)) {
    stop("missing or infinite values for 'N' and 'E' are not allowed")
  }
}

#ebgm() ------------------------------------------------------------------------
.checkInputs_ebgm <- function(theta_hat = theta_hat, N = N, E = E, qn = qn,
                  digits = digits) {
  #Sanity checks for actual arguments
  if (length(theta_hat) != 5) {
    stop("'theta_hat' must contain 5 values")
  }
  if (min(theta_hat) <= 0) {
    stop("'theta_hat' must contain only positive values")
  }
  if (theta_hat[5] >= 1) {
    stop("'theta_hat[5]' (i.e., 'P') must be less than 1")
  }
  if (length(N) != length(E) || length(N) != length(qn)) {
    stop("'N', 'E', and 'qn' must have the same length")
  }
  if (.isMissing_num(N) || .isMissing_num(E) || .isMissing_num(qn)) {
    stop("missing or infinite values for 'N', 'E', and 'qn' are not allowed")
  }
}

#quantBisect() -----------------------------------------------------------------
.checkInputs_quantBisect <- function(percent = percent, theta_hat = theta_hat,
                                     N = N, E = E, qn = qn, digits = digits,
                                     limits = limits, max_iter = max_iter) {
  #Sanity checks for actual arguments
  if (percent < 1 || percent > 99) {
    stop("'percent' must be a value between 1 and 99 (inclusive)")
  }
  if (length(theta_hat) != 5) {
    stop("'theta_hat' must contain 5 values")
  }
  if (min(theta_hat) <= 0) {
    stop("'theta_hat' must contain only positive values")
  }
  if (theta_hat[5] >= 1) {
    stop("'theta_hat[5]' (i.e., 'P') must be less than 1")
  }
  if (length(N) != length(E) || length(N) != length(qn)) {
    stop("'N', 'E', and 'qn' must have the same length")
  }
  if (.isMissing_num(N) || .isMissing_num(E) || .isMissing_num(qn)) {
    stop("missing or infinite values for 'N', 'E', and 'qn' are not allowed")
  }
}

#ebScores() --------------------------------------------------------------------
.checkInputs_ebScores <- function(processed = processed,
                                  hyper_estimate = hyper_estimate,
                                  quantiles = quantiles, digits = digits) {
  if(!is.null(quantiles) & !is.numeric(quantiles)) {
    stop("'quantiles' must be NULL or a numeric vector of quantiles")
  }
  if(!is.list(hyper_estimate)) {
    stop("'hyper_estimate' must be the list output by autoHyper()")
  }
  if(!any(grepl("estimates", names(hyper_estimate)))) {
    stop("'hyper_estimate' must be a list containing an element of hyperparameter
         estimates. Was it actually calculated by autoHyper()?")
  }
  if(!is.data.frame(processed)) {
    stop("'processed' must be a data frame from processRaw()")
  }
  if(!any(grepl("var", names(processed)))) {
    stop("'processed' dataframe does not have 'var' variables. Was this dataframe
         actually created by processRaw()?")
  }
}
