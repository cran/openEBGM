**openEBGM** version history
============================

## v0.11.0
* More improvements to numerical stability of hyperparameter calculations.
* hyperEM() now uses proportion change in log-likelihood in the convergence
  criterion instead of raw change.

## v0.10.0
* Improved numerical stability of hyperparameter calculations.
* Now requires 'numDeriv' package for calculating Hessian matrices.
* Removed 'param_lower' and 'param_upper' arguments from hyperEM().
* Added 'transformed' argument to the likelihood functions to allow for
  optimization on the log & logit scales.
* Added 'min_theta' & 'max_theta' arguments to Qn(), ebgm(), quantBisect(), &
  ebScores(). Alpha values now must also be between one half & twice the value
  of the corresponding beta values. Values outside these ranges are replaced by
  the limits.
* Updates to vignettes, data, examples, etc.
* Fixed NOTES regarding documentation.

## v0.9.1
* Ensured that no more than two threads are used for CRAN checks.
* Optimized examples for speed.
* Minor changes to vignettes.

## v0.9.0
* Use the new special sentinel "_PACKAGE" for package documentation.
* Various minor updates (URLs, etc.).

## v0.8.3
* Fixed unit testing problems caused by R's future switch to
  stringsAsFactors = FALSE.

## v0.8.2
* Adjusted calculation for expected counts using suggestion from
  Piotr Świnarski. Previously, calculation failed when marginal counts became
  too large for integer multiplication.

## v0.8.1
* Corrected unit test failures for processRaw() resulting from base R changes to
  the sample() function.
* Added DEoptim::DEoptim() example to hyperparameter estimation vignette.

## v0.8.0
* processRaw() now lists all strata when stratification is used.
* Added argument 'list_ids' to processRaw().

## v0.7.0
* Added the autoSquash() function to automate data squashing.
* Changed exit condition for while loop in hyperEM(). hyperEM() now throws an 
  error if the number of "stuck" or repeated estimates of theta exceeds 20
  when using 'method = "nlminb"'.
* Changed upper limit from 1 to 0.999 in hidden functions .updateThetaLL() and
  .updateThetaLLD(), which are called by hyperEM().
    
## v0.6.0
* Changed 'keep_bins' formal argument in squashData() to 'keep_pts' for added
  flexibility.

## v0.5.0
* Efficiency and code hygiene improvements to processRaw() and squashData().

## v0.4.0
* Added the hyperEM() function to estimate hyperparameters using an
  implementation of the EM algorithm.

## v0.3.0
* Added confidence intervals to autoHyper() and standard errors to autoHyper()
  and exploreHypers().
* processRaw() now returns Inf instead of 99999 when PRR results in division by
  zero.
* Fixed minor bug in exploreHypers().

## v0.2.0
* Minor aesthetic changes to plot(), summary(), and print() methods.
* Relaxed convergence requirements for exploreHypers() and autoHyper().
