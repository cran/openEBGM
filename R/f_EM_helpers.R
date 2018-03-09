# Various helpers for the EM algorithm approach to hyperparameter estimation

# -----------------------  Print messages  -------------------------------------
# Purpose: Print messages depending on various levels of desired output.
# ------------------------------------------------------------------------------
.hyperEMmessage <- function(message = c("A", "B", "C"), print_level,
                            iter = NULL, delta_LL = NULL, theta = NULL,
                            elapsed = NULL) {
  #Print theta & change in LL for every 10th iteration.
  #Print iteration count every 50th iteration.
  if (message == "A") {
    if (print_level == 2) {
      if (iter %% 10 == 0) {
        delta_LL <- format(delta_LL, digits = 3, scientific = TRUE, width = 8)
        theta    <- as.character(round(theta, 3))
        theta    <- format(theta, digits = 2, scientific = FALSE, width = 5)
        cat("change in LL:", delta_LL, "| ")
        cat("theta:", theta, "\n")
        if (iter %% 50 == 0) {
          cat("    Current iterations:", iter, "\n")
        }
      }
    }
  #Print iteration count & time elapsed in seconds.
  } else if (message == "B") {
    if (print_level %in% 1:2) {
      cat("\n   Iterations used:", iter, "\n")
      cat("\n   Timing:  \n")
      print(elapsed)
      cat("\n")
    }
  #Print initial value used for theta (row number).
  } else if (message == "C") {
    if (print_level %in% 1:2) {
      start_pt <- c("\n    ******  ", "Starting point:", iter, "  ******\n\n")
      if (print_level %in% 1:2) cat(start_pt)
    }
  } else {
    stop("'message' argument not recognized")
  }
}

# ------------------  Change in log-likelihood  --------------------------------
# Purpose: Find change in log-likelihood function for a particular iteration.
#          Used for convergence assessment.
# ------------------------------------------------------------------------------

.deltaLL <- function(theta, old_theta, N, E, W, squashed, zeroes, N_star = 1) {
  #Change in log-likelihood
  if (zeroes) {
    if (squashed) {
      old_LL <- -negLLzeroSquash(old_theta, N, E, W)
      LL     <- -negLLzeroSquash(theta, N, E, W)
    } else {
      old_LL <- -negLLzero(old_theta, N, E)
      LL <- -negLLzero(theta, N, E)
    }
  } else {
    if (squashed) {
      old_LL <- -negLLsquash(old_theta, N, E, W, N_star = N_star)
      LL     <- -negLLsquash(theta, N, E, W, N_star = N_star)
    } else {
      old_LL <- -negLL(old_theta, N, E, N_star = N_star)
      LL     <- -negLL(theta, N, E, N_star = N_star)
    }
  }
  delta <- abs(LL - old_LL)

  list(delta = delta, LL = LL)
}

# ------------------------  Check results  -------------------------------------
# Purpose: Check the consistency of the results by looking at all the pairwise
#            ratios for all the results. (alpha1-to-alpha1 for all results, etc.)
# ------------------------------------------------------------------------------

.checkRatios <- function(theta_eb_df, ratio_limit) {
  #Check if ratios are "good"
  theta_eb_df <- theta_eb_df[, 2:6]
  max_ratio <- function(estimate) {
    max(estimate, na.rm = TRUE) / min(estimate, na.rm = TRUE)
  }
  theta_eb_ratio <- apply(theta_eb_df, 2, max_ratio)
  if (max(theta_eb_ratio) > ratio_limit) {
    warning("at least one starting point led to a different estimate")
  }

  theta_eb_ratio
}
