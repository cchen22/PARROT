#' Statistical Inference Functions for PARROT
#'
#' Compute p-values, adjust for multiple testing, and assess significance.
#' Supports Wilks' asymptotic test, permutation, and bootstrap methods.

#' Compute P-values for Detected Change Points
#'
#' Computes, for each change point of a \code{parrot_result}, the p-value of the
#' split held fixed inside the window between the neighbouring change points.
#' The resampled statistics use the fit settings stored in the result.
#'
#' @param result A \code{parrot_result} (or a single-CP \code{cp_detection} for
#'   \code{method = "wilks"})
#' @param method "permutation" (conditional time-label permutation, exact when
#'   feasible), "wilks", "stationary_bootstrap" or "bootstrap" (parametric)
#' @param n_perm Number of random permutations / bootstrap samples when the
#'   p-value is not exact
#' @param permutation_scope "global" (default) or "local"
#' @param exact_max Enumerate all configurations when at most this many
#' @param verbose Print progress
#' @return List with one entry per change point (raw, Bonferroni and BH
#'   adjusted p-values)
#' @export
#' @examples
#' \dontrun{
#' res <- parrot(data$Y, Q = 2, multiple = FALSE, min_segment = 2, compute_pvalues = FALSE)
#' compute_pvalue(res, method = "wilks")
#' }
compute_pvalue <- function(result,
                           method = c("permutation", "wilks", "stationary_bootstrap", "bootstrap"),
                           n_perm = 99,
                           permutation_scope = c("global", "local"),
                           exact_max = 5040,
                           verbose = FALSE) {
  method <- match.arg(method)
  permutation_scope <- match.arg(permutation_scope)
  if (inherits(result, "cp_detection")) {
    if (method != "wilks") stop("Resampling p-values need the data: use a parrot_result.")
    return(list(list(cp = result$cp_raw, lr_stat = result$lr_stat, df = result$df,
                     pvalue_raw = .pvalue_wilks(result$lr_stat, result$df),
                     pvalue_method = "wilks")))
  }
  if (!inherits(result, "parrot_result")) stop("result must be a parrot_result object")
  .compute_pvalues_for_parrot(result, method, verbose, n_perm = n_perm,
                              permutation_scope = permutation_scope, exact_max = exact_max)
}

#' Adjust P-values for Multiple Testing
#'
#' Adjust p-values using various methods (Bonferroni, BH, etc.)
#'
#' @param pvalues Vector of p-values
#' @param method Adjustment method: "bonferroni", "holm", "BH", "BY", "none"
#' @return Adjusted p-values
#' @export
#' @examples
#' pvals <- c(0.01, 0.04, 0.03, 0.5)
#' adjust_pvalues(pvals, method = "BH")
adjust_pvalues <- function(pvalues,
                           method = c("bonferroni", "holm", "BH", "BY", "none")) {
  method <- match.arg(method)

  if (method == "none") {
    return(pvalues)
  }

  p.adjust(pvalues, method = method)
}

#' Test Significance of Detected Change Points
#'
#' Tests every change point of a \code{parrot_result} (see
#' \code{\link{compute_pvalue}}) and adjusts for multiple testing.
#'
#' @param result A parrot_result object
#' @param alpha Significance level
#' @param adjust Adjustment method for multiple testing
#' @param method P-value method passed to \code{\link{compute_pvalue}}
#' @param n_perm Number of random permutations / bootstrap samples
#' @param verbose Print progress
#' @return Object of class "cp_significance_test"
#' @export
test_significance <- function(result,
                              alpha = 0.05,
                              adjust = c("bonferroni", "BH", "holm", "BY", "none"),
                              method = c("permutation", "wilks", "stationary_bootstrap", "bootstrap"),
                              n_perm = 99,
                              verbose = FALSE) {

  adjust <- match.arg(adjust)
  method <- match.arg(method)
  if (!inherits(result, "parrot_result")) stop("result must be a parrot_result object")

  if (result$n_cps == 0) {
    if (verbose) message("No change points to test")
    out <- list(cps = integer(0), significant = logical(0), pvalues = numeric(0),
                adjusted_pvalues = numeric(0), n_significant = 0L, alpha = alpha,
                adjust_method = adjust)
    class(out) <- c("cp_significance_test", "list")
    return(out)
  }

  pv <- compute_pvalue(result, method = method, n_perm = n_perm, verbose = verbose)
  pvalues <- vapply(pv, function(x) x$pvalue_raw, numeric(1))
  adjusted_pvalues <- adjust_pvalues(pvalues, method = adjust)
  significant <- adjusted_pvalues < alpha
  significant[is.na(significant)] <- FALSE

  test_result <- list(
    cps = vapply(pv, function(x) as.numeric(x$cp), numeric(1)),
    significant = significant,
    pvalues = pvalues,
    adjusted_pvalues = adjusted_pvalues,
    n_significant = sum(significant),
    alpha = alpha,
    adjust_method = adjust,
    pvalue_method = method
  )
  class(test_result) <- c("cp_significance_test", "list")

  if (verbose) {
    message(sprintf("Tested %d change points", length(test_result$cps)))
    message(sprintf("Significant at alpha=%.2f after %s correction: %d",
                    alpha, adjust, sum(significant)))
  }

  test_result
}

#' Profile Likelihood Confidence Interval
#'
#' Compute confidence interval for a change point using profile likelihood.
#'
#' @param result A cp_detection object
#' @param level Confidence level (default: 0.95)
#' @return List with CI bounds and the confidence set
#' @export
profile_likelihood_ci <- function(result, level = 0.95) {
  
  if (!inherits(result, "cp_detection")) {
    stop("result must be a cp_detection object")
  }
  
  if (!result$significant || is.na(result$cp_est)) {
    return(list(
      lower = NA,
      upper = NA,
      level = level,
      ci_set = integer(0)
    ))
  }
  
  # Chi-squared cutoff
  alpha <- 1 - level
  cutoff <- result$alt_elbo - 0.5 * qchisq(1 - alpha, df = 1)
  
  # Find confidence set
  ci_set <- which(result$profile_ll >= cutoff)
  
  list(
    lower = if (length(ci_set) > 0) min(ci_set) else NA,
    upper = if (length(ci_set) > 0) max(ci_set) else NA,
    level = level,
    ci_set = ci_set,
    cutoff = cutoff
  )
}

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Wilks p-value
#' @keywords internal
.pvalue_wilks <- function(lr_stat, df) {
  if (is.na(lr_stat) || is.na(df)) return(NA)
  pchisq(lr_stat, df = df, lower.tail = FALSE)
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' Print method for cp_significance_test
#' @param x A `cp_significance_test` object returned by `test_significance()`.
#' @param ... Additional arguments (currently unused).
#' @export
print.cp_significance_test <- function(x, ...) {
  cat("\n=== Change Point Significance Test ===\n\n")
  cat(sprintf("Adjustment method: %s\n", x$adjust_method))
  cat(sprintf("Significance level: %.3f\n\n", x$alpha))
  
  if (length(x$cps) == 0) {
    cat("No change points to test\n")
  } else {
    cat("Results:\n")
    for (i in seq_along(x$cps)) {
      sig_str <- if (x$significant[i]) "***" else ""
      cat(sprintf("  CP %d (t=%d): p=%.4f, adj.p=%.4f %s\n",
                  i, x$cps[i], x$pvalues[i], x$adjusted_pvalues[i], sig_str))
    }
    cat(sprintf("\nSignificant: %d / %d\n", x$n_significant, length(x$cps)))
  }
  
  invisible(x)
}
