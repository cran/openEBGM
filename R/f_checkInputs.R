#Helper functions to check if actual argument values make sense

#*******************************************************************************
#************************  General Functions  **********************************
#*******************************************************************************

.checkInputs_processed_data <- function(data) {

  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!all(c("N", "E") %in% colnames(data))) {
    stop("missing the appropriate column names (need 'N' and 'E')")
  }
  if (.isMissing_num(data$N) || .isMissing_num(data$E)){
    stop("missing or infinite values for 'N' and 'E' are not allowed")
  }
}

.checkInputs_squashed_zeroes <- function(data, squashed, zeroes) {

  if (!is.logical(squashed) || !is.logical(zeroes)) {
    stop("'squashed' and 'zeroes' must be logical values")
  }
  if (squashed && !"weight" %in% colnames(data)) {
    stop("'weight' column is missing -- are these data really squashed?")
  }
  if (!squashed && "weight" %in% colnames(data)) {
    stop("'weight' column was found -- were these data squashed?")
  }
  if (zeroes && min(data$N) != 0) {
    stop("no zero counts found")
  }
  if (!zeroes && min(data$N) == 0) {
    stop("zero counts found")
  }
}

.checkInputs_Nstar <- function(data, zeroes, N_star) {

  if (zeroes && !is.null(N_star)) {
    stop("if zeroes are used, 'N_star' should be NULL")
  }
  if (!zeroes && is.null(N_star)) {
    stop("if zeroes are not used, 'N_star' must be specified")
  }
  if (!is.null(N_star)) {
    if (!is.numeric(N_star) || N_star < 1 || N_star %% 1 != 0) {
      stop("'N_star' must be NULL or a positive whole number")
    }
    N_star <- as.integer(N_star)
    if (N_star != min(as.integer(names(table(data$N))))) {
      stop("'N_star' does not agree with the data set")
    }
  }
}

.checkInputs_theta_init <- function(theta_init) {

  if (ncol(theta_init) != 5) {
    stop("'theta_init' must contain 5 columns")
  }
  if (min(theta_init) <= 0) {
    stop("'theta_init' must contain nonnegative values")
  }
  if (max(theta_init[5]) >= 1) {
    stop("'theta[5]' (i.e., 'P') must be <1")
  }
}

.checkInputs_theta_hat <- function(theta_hat) {

  if (length(theta_hat) != 5) {
    stop("'theta_hat' must contain 5 values")
  }
  if (any(is.na((theta_hat)))) {
    stop("'theta_hat' cannot contain missing values")
  }
  if (min(theta_hat) <= 0) {
    stop("'theta_hat' must contain only positive values")
  }
  if (theta_hat[5] >= 1) {
    stop("'theta_hat[5]' (i.e., 'P') must be less than 1")
  }
}

.checkInputs_theta_init_vec <- function(theta_init_vec) {

  if (length(theta_init_vec) != 5) {
    stop("'theta_init_vec' must contain 5 elements")
  }
  if (any(is.na((theta_init_vec)))) {
    stop("'theta_init_vec' cannot contain missing values")
  }
  if (!is.numeric(theta_init_vec) || min(theta_init_vec) <= 0) {
    stop("'theta_init_vec' must contain positive numeric values")
  }
  if (theta_init_vec[5] >= 1) {
    stop("'theta_init_vec[5]' (i.e., 'P') must be <1")
  }
}

.checkInputs_score <- function(N_star) {

  if (N_star != 1 && !is.null(N_star)) {
    stop("if 'method' is 'score', N_star must be 1 or NULL")
  }
}

.checkInputs_LL_converge <- function(LL_tol, consecutive, max_iter) {

  if (!is.numeric(LL_tol) || LL_tol <= 0) {
    stop("'LL_tol' must be a positive numeric value")
  }
  if (!is.numeric(consecutive) || consecutive < 0 || consecutive %% 1 != 0) {
    stop("'consecutive' must be a nonnegative whole number")
  }
  if (!is.numeric(max_iter) || max_iter <= 0 || max_iter %% 1 != 0) {
    stop("'max_iter' must be a positive whole number")
  }
}

.checkInputs_param_lower_upper <- function(param_lower, param_upper) {

  if (!is.numeric(param_lower) || param_lower <= 0) {
    stop("'param_lower' must be a positive numeric value")
  }
  if (!is.numeric(param_upper) || param_upper <= 0) {
    stop("'param_upper' must be a positive numeric value")
  }
  if (param_lower >= param_upper) {
    stop("'param_lower' must be less than 'param_upper'")
  }
}

.checkInputs_print_level <- function(print_level) {

  if (!print_level %in% 0:2) {
    stop("'print_level' must be 0, 1, or 2")
  }
}

.checkInputs_N_E <- function(N, E) {

  if (length(N) != length(E)) {
    stop("'N' and 'E' must have the same length")
  }
  if (.isMissing_num(N) | .isMissing_num(E)) {
    stop("missing or infinite values for 'N' and 'E' are not allowed")
  }
}

.checkInputs_N_E_qn <- function(N, E, qn) {

  if (length(N) != length(E) || length(N) != length(qn)) {
    stop("'N', 'E', and 'qn' must have the same length")
  }
  if (.isMissing_num(N) || .isMissing_num(E) || .isMissing_num(qn)) {
    stop("missing or infinite values for 'N', 'E', and 'qn' are not allowed")
  }
}

#*******************************************************************************
#*******************  Data Processing Functions  *******************************
#*******************************************************************************

#processRaw() ------------------------------------------------------------------
.checkInputs_processRaw <- function(data, stratify, zeroes) {

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

.checkStrata_processRaw <- function(data, max_cats) {

  #Also prints messages and adds 'stratum' column to 'data'
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
.checkInputs_squashData <- function(data, count, bin_size, keep_pts,
                                    min_bin, min_pts) {

  #Also coerces data frame to data table.
  .checkInputs_processed_data(data)
  count    <- as.integer(count)
  bin_size <- as.integer(bin_size)
  keep_pts <- as.integer(keep_pts)
  min_bin  <- as.integer(min_bin)
  min_pts  <- as.integer(min_pts)
  data <- data.table::as.data.table(data)
  data[, N := as.integer(N)]
  data[, E := as.numeric(E)]

  if (.isMissing_num(data[, N]) | .isMissing_num(data[, E])) {
    stop("missing values for 'N' or 'E' are not allowed")
  }
  if (count < 0) {
    stop("'count' must be non-negative")
  }
  if (bin_size < 2) {
    stop("'bin_size' must be >= 2")
  }
  if (keep_pts < 0) {
    stop("'keep_pts' must be non-negative")
  }

  #In case we want to run the function iteratively (count = 1, then 2, etc.)
  weight_exists <- any(grepl("weight", colnames(data)))
  if (weight_exists) {
    if (max(data[N == count, weight], na.rm = TRUE) > 1) {
      stop("this data set has already been squashed for this count size")
    }
  } else {
    data[, weight := 1L]
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

#autoSquash() ------------------------------------------------------------------
.checkInputs_autoSquash <- function(data, keep_pts, cut_offs, num_super_pts) {

  .checkInputs_processed_data(data)
  if ("weight" %in% colnames(data)) {
    stop("'data' has already been squashed")
  }
  check_miss <- c(keep_pts, cut_offs, num_super_pts)
  if (any(c(is.na(check_miss), is.nan(check_miss)))) {
    stop("missing values not allowed in 'keep_pts', 'cut_offs', or 'num_super_pts'")
  }
  if (length(keep_pts) < 1 || length(cut_offs) < 1 || length(num_super_pts) < 1) {
    stop("'keep_pts', 'cut_offs', & 'num_super_pts' must have a length of 1 or more")
  }
  if (length(num_super_pts) != length(cut_offs) + 1) {
    stop("'num_super_pts' must have length 1 more than length of 'cut_offs'")
  }
  if (min(diff(cut_offs)) < 1) {
    stop("elements in 'cut_offs' must be in increasing order")
  }
}


#*******************************************************************************
#************************  Likelihood Functions  *******************************
#*******************************************************************************

.checkInputs_negLLzero <- function(theta, N, E) {

  if (length(N) != length(E)) {
    stop("'N' & 'E' must be the same length")
  }
  if (min(N) != 0) {
    stop("no zero counts found -- please use another likelihood function")
  }
}

.checkInputs_negLLzeroSquash <- function(theta, ni, ei, wi) {

  ni_len <- length(ni)
  if (ni_len != length(ei) || ni_len != length(wi)) {
    stop("'ni', 'ei', & 'wi' must be the same length")
  }
  if (min(ni) != 0) {
    stop("no zero counts found -- please use another likelihood function")
  }
}

.checkInputs_negLL <- function(theta, N, E, N_star) {

  if (length(N) != length(E)) {
    stop("'N' & 'E' must be the same length")
  }
  if (min(N) != N_star) {
    stop("'N_star' does not agree with 'N'")
  }
}

.checkInputs_negLLsquash <- function(theta, ni, ei, wi, N_star) {

  ni_len <- length(ni)
  if (ni_len != length(ei) || ni_len != length(wi)) {
    stop("'ni', 'ei', & 'wi' must be the same length")
  }
  if (min(ni) != N_star) {
    stop("'N_star' does not agree with 'ni'")
  }
}


#*******************************************************************************
#************************  Prior Functions  ************************************
#*******************************************************************************

#exploreHypers() ---------------------------------------------------------------
.checkInputs_exploreHypers <- function(data, theta_init, squashed, zeroes,
                                       N_star, method, param_limit, max_pts,
                                       std_errors) {

  .checkInputs_processed_data(data)
  .checkInputs_theta_init(theta_init)
  .checkInputs_squashed_zeroes(data, squashed, zeroes)
  .checkInputs_Nstar(data, zeroes, N_star)
  if (nrow(data) > max_pts) {
    error1 <- "there are more than "
    error2 <- "\n  either squash data, run individual optimizations,"
    error3 <- "\n  or increase 'max_pts' (caution: long run times)"
    stop(paste0(error1, max_pts, " points --", error2, error3))
  }
  if (!is.logical(std_errors)) {
    stop("'std_errors' must be a logical value")
  }
}

#autoHyper() -------------------------------------------------------------------
.checkInputs_autoHyper <- function(tol, min_conv, theta_init, conf_ints) {

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
  if (!is.logical(conf_ints)) {
    stop("'conf_ints' must be a logical value")
  }
}

#hyperEM() ---------------------------------------------------------------------
.checkInputs_hyperEM <-
  function(data, theta_init_vec, squashed, zeroes, N_star, method, LL_tol,
           consecutive, param_lower, param_upper, print_level, max_iter,
           conf_ints, track) {

    .checkInputs_processed_data(data)
    .checkInputs_theta_init_vec(theta_init_vec)
    .checkInputs_squashed_zeroes(data, squashed, zeroes)
    if (method == "score") .checkInputs_score(N_star)
    .checkInputs_Nstar(data, zeroes, N_star)
    .checkInputs_LL_converge(LL_tol, consecutive, max_iter)
    .checkInputs_param_lower_upper(param_lower, param_upper)
    .checkInputs_print_level(print_level)
    if (!is.logical(conf_ints)) {
      stop("'conf_ints' must be a logical value")
    }
    if (!is.logical(track)) {
      stop("'track' must be a logical value")
    }
}

#*******************************************************************************
#************************  Posterior Functions  ********************************
#*******************************************************************************

#Qn() --------------------------------------------------------------------------
.checkInputs_Qn <- function(theta_hat, N, E) {

  .checkInputs_theta_hat(theta_hat)
  .checkInputs_N_E(N, E)
}

#ebgm() ------------------------------------------------------------------------
.checkInputs_ebgm <- function(theta_hat, N, E, qn, digits) {

  .checkInputs_theta_hat(theta_hat)
  .checkInputs_N_E_qn(N, E, qn)
}

#quantBisect() -----------------------------------------------------------------
.checkInputs_quantBisect <- function(percent, theta_hat, N, E, qn, digits,
                                     limits, max_iter) {

  if (percent < 1 || percent > 99) {
    stop("'percent' must be a value between 1 and 99 (inclusive)")
  }
  .checkInputs_theta_hat(theta_hat = theta_hat)
  .checkInputs_N_E_qn(N, E, qn)
}

#ebScores() --------------------------------------------------------------------
.checkInputs_ebScores <- function(processed, hyper_estimate, quantiles,
                                  digits) {

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
  .checkInputs_processed_data(processed)
}
