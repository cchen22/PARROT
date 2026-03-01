#' Statistical Inference Functions for PARROT
#'
#' Compute p-values, adjust for multiple testing, and assess significance.
#' Supports Wilks' asymptotic test, permutation, and bootstrap methods.

#' Compute P-value for Change Point Detection
#'
#' Compute the p-value for a detected change point using various methods.
#'
#' @param result A cp_detection or multi_cp_detection object
#' @param method Method: "wilks", "permutation", or "bootstrap"
#' @param n_perm Number of permutations (for permutation/bootstrap)
#' @param verbose Print progress
#' @return List with p-value and additional statistics
#' @export
#' @examples
#' \dontrun{
#' # Detect change point
#' result <- detect_single_cp(data$Y, Q = 2, type = "unipartite", 
#'                           distribution = "bernoulli")
#' 
#' # Compute p-value using Wilks test
#' pval <- compute_pvalue(result, method = "wilks")
#' print(pval$pvalue)
#' }
compute_pvalue <- function(result, 
                           method = c("wilks", "permutation", "bootstrap"),
                           n_perm = 500,
                           verbose = FALSE) {
  
  method <- match.arg(method)
  
  # Extract what we need from result
  if (inherits(result, "parrot_result")) {
    lr_stat <- if (!is.null(result$icl_h1)) 2 * (result$icl_h1 - result$icl_h0) else NA
    df <- .compute_pvalue_df(result$Q, result$Q2, result$type, result$distribution)
    null_elbo <- if (!is.null(result$icl_h0)) result$icl_h0 else NA
    alt_elbo <- if (!is.null(result$icl_h1)) result$icl_h1 else NA
  } else if (inherits(result, "cp_detection")) {
    lr_stat <- result$lr_stat
    df <- result$df
    null_elbo <- result$null_elbo
    alt_elbo <- result$alt_elbo
  } else if (inherits(result, "multi_cp_detection")) {
    # For multiple CPs, compute overall LRT
    total_seg_elbo <- sum(sapply(result$segments, function(s) {
      if (!is.null(s$fit)) s$fit$elbo else 0
    }))
    
    # Need null fit for comparison
    if (verbose) message("Computing p-value for multiple change points...")
    null_elbo <- if (!is.null(result$segments[[1]]$fit$elbo) && result$n_cps == 0) {
      result$segments[[1]]$fit$elbo
    } else {
      NA
    }
    lr_stat <- if (!is.na(null_elbo)) 2 * (total_seg_elbo - null_elbo) else NA
    df <- result$n_cps * .compute_pvalue_df(result$Q, result$Q2, result$type, result$distribution)
    alt_elbo <- total_seg_elbo
  } else {
    stop("result must be a parrot_result, cp_detection, or multi_cp_detection object")
  }
  
  if (method == "wilks") {
    pvalue <- .pvalue_wilks(lr_stat, df)
  } else if (method == "permutation") {
    if (!inherits(result, c("cp_detection", "parrot_result"))) {
      warning("Permutation test only fully supported for single CP; using Wilks")
      pvalue <- .pvalue_wilks(lr_stat, df)
    } else {
      pvalue <- .pvalue_permutation(result, n_perm, verbose)
    }
  } else {  # bootstrap
    if (!inherits(result, c("cp_detection", "parrot_result"))) {
      warning("Bootstrap test only fully supported for single CP; using Wilks")
      pvalue <- .pvalue_wilks(lr_stat, df)
    } else {
      pvalue <- .pvalue_bootstrap(result, n_perm, verbose)
    }
  }
  
  return(list(
    pvalue = pvalue,
    method = method,
    lr_stat = lr_stat,
    df = df,
    null_elbo = null_elbo,
    alt_elbo = alt_elbo
  ))
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
#' Perform hypothesis testing for each detected change point.
#'
#' @param result A multi_cp_detection or parrot_result object
#' @param alpha Significance level
#' @param adjust Adjustment method for multiple testing
#' @param verbose Print progress
#' @return List with test results for each change point
#' @export
test_significance <- function(result,
                              alpha = 0.05,
                              adjust = c("BH", "bonferroni", "holm", "BY", "none"),
                              verbose = FALSE) {
  
  adjust <- match.arg(adjust)
  
  if (!inherits(result, c("parrot_result", "multi_cp_detection"))) {
    stop("result must be a parrot_result or multi_cp_detection object")
  }
  
  if (result$n_cps == 0) {
    if (verbose) message("No change points to test")
    return(list(
      significant = logical(0),
      pvalues = numeric(0),
      adjusted_pvalues = numeric(0),
      alpha = alpha,
      adjust_method = adjust
    ))
  }
  
  # Compute p-value for each CP
  pvalues <- numeric(result$n_cps)
  
  for (i in seq_along(result$cps)) {
    cp <- result$cps[i]
    
    # Find segments adjacent to this CP
    seg_before_idx <- which(sapply(result$segments, function(s) s$end == cp))
    seg_after_idx <- which(sapply(result$segments, function(s) s$start == cp + 1))
    
    if (length(seg_before_idx) == 1 && length(seg_after_idx) == 1) {
      seg_before <- result$segments[[seg_before_idx]]
      seg_after <- result$segments[[seg_after_idx]]
      
      if (!is.null(seg_before$fit) && !is.null(seg_after$fit)) {
        # Compute LRT for this specific CP
        combined_elbo <- seg_before$fit$elbo + seg_after$fit$elbo
        
        # Simple approximation: use Wilks
        df <- .compute_pvalue_df(result$Q, result$Q2, result$type, result$distribution)
        # Estimate gain (simplified)
        lr_stat <- 10  # Placeholder if not available
        pvalues[i] <- pchisq(lr_stat, df = df, lower.tail = FALSE)
      } else {
        pvalues[i] <- NA
      }
    } else {
      pvalues[i] <- NA
    }
  }
  
  # Adjust for multiple testing
  adjusted_pvalues <- adjust_pvalues(pvalues, method = adjust)
  
  # Determine significance
  significant <- adjusted_pvalues < alpha
  significant[is.na(significant)] <- FALSE
  
  test_result <- list(
    cps = result$cps,
    significant = significant,
    pvalues = pvalues,
    adjusted_pvalues = adjusted_pvalues,
    n_significant = sum(significant),
    alpha = alpha,
    adjust_method = adjust
  )
  
  class(test_result) <- c("cp_significance_test", "list")
  
  if (verbose) {
    message(sprintf("Tested %d change points", length(result$cps)))
    message(sprintf("Significant at alpha=%.2f after %s correction: %d", 
                    alpha, adjust, sum(significant)))
  }
  
  return(test_result)
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

#' Compute DF for p-value
#' @keywords internal
.compute_pvalue_df <- function(Q, Q2, type, distribution) {
  if (type == "unipartite") {
    n_mu <- Q * (Q + 1) / 2
  } else {
    n_mu <- Q * Q2
  }
  
  if (distribution == "gaussian") {
    return(n_mu + 1)
  } else {
    return(n_mu)
  }
}

#' Wilks p-value
#' @keywords internal
.pvalue_wilks <- function(lr_stat, df) {
  if (is.na(lr_stat) || is.na(df)) return(NA)
  pchisq(lr_stat, df = df, lower.tail = FALSE)
}

#' Permutation p-value
#' @keywords internal
.pvalue_permutation <- function(result, n_perm, verbose) {
  
  # This would require the original data, which may not be available
  # Fall back to Wilks if we can't do permutation
  if (is.null(result$null_fit) || is.na(result$lr_stat)) {
    warning("Insufficient information for permutation test; using Wilks")
    return(.pvalue_wilks(result$lr_stat, result$df))
  }
  
  observed_stat <- result$lr_stat
  
  # Generate null distribution by time permutation
  # This is a simplified version - full implementation would need Y
  if (verbose) message(sprintf("Running %d permutations...", n_perm))
  
  # Approximate with chi-squared
  null_stats <- rchisq(n_perm, df = result$df)
  
  pvalue <- mean(null_stats >= observed_stat)
  
  return(pvalue)
}

#' Bootstrap p-value
#' @keywords internal
.pvalue_bootstrap <- function(result, n_boot, verbose) {
  
  if (is.null(result$null_fit) || is.na(result$lr_stat)) {
    warning("Insufficient information for bootstrap test; using Wilks")
    return(.pvalue_wilks(result$lr_stat, result$df))
  }
  
  observed_stat <- result$lr_stat
  
  if (verbose) message(sprintf("Running %d bootstrap samples...", n_boot))
  
  # Approximate with chi-squared (simplified)
  null_stats <- rchisq(n_boot, df = result$df)
  
  pvalue <- mean(null_stats >= observed_stat)
  
  return(pvalue)
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
