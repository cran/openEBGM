# Various helper functions

.isMissing_str <- function(x) {
  # Determines if any elements in a string vector are missing.
  any(is.na(x) | x == "")
}

.isMissing_num <- function(x) {
  # Determines if any elements in a numeric vector are missing.
  any(is.na(x) | is.nan(x) | is.infinite(x))
}

.countUnique <- function(x) {
  # Counts the number of unique elements in a vector.
  length(unique(x))
}

.reinInTheta <- function(theta_hat, min_theta, max_theta) {
  # Helps prevent unusual shrinkage properties by adjusting extreme
  #   hyperparameter estimates
  if (length(min_theta) != 5 || length(max_theta) != 5) {
    stop("'min_theta' & 'max_theta' must have length 5")
  }
  if (.isMissing_num(min_theta) || .isMissing_num(max_theta)) {
    stop("missing or infinite values found")
  }
  if (any(min_theta <= 0) || any(max_theta <= 0)) {
    stop("negative values found")
  }
  if (any(max_theta - min_theta <= 0)) {
    stop("'min_theta' must be less than 'max_theta'")
  }
  theta_hat <- pmax(theta_hat, min_theta)
  theta_hat <- pmin(theta_hat, max_theta)
  # 'alpha' cannot be more than twice 'beta'
  theta_hat[1] <- min(theta_hat[1], 2 * theta_hat[2])
  theta_hat[3] <- min(theta_hat[3], 2 * theta_hat[4])
  # 'alpha' cannot be less than half of 'beta'
  theta_hat[1] <- max(theta_hat[1], 0.5 * theta_hat[2])
  theta_hat[3] <- max(theta_hat[3], 0.5 * theta_hat[4])
  theta_hat
}
