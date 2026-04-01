#' PARROT: Phase-Altering Regulatory Rewiring Over Time
#'
#' Main function for change point detection in dynamic networks.
#' Provides a unified interface for detecting single or multiple change points
#' using optimal default settings.
#'
#' @param Y Network tensor (T x N x N for unipartite, T x N1 x N2 for bipartite)
#' @param Q Number of clusters (or Q1 for bipartite). Can be:
#'   - NULL: automatically selected via ICL
#'   - Integer: use this exact value
#'   - Vector (e.g., 2:4): select best Q from this range via ICL
#' @param Q2 Number of column clusters for bipartite (NULL for unipartite).
#'   Can also be a range like Q.
#' @param type Network type: "unipartite" or "bipartite" (auto-detected if NULL)
#' @param distribution Edge distribution: "gaussian", "bernoulli", or "auto"
#' @param multiple Detect multiple change points (TRUE) or single (FALSE)
#' @param method Search method: "wbs" (Wild Binary Segmentation), "bs" (Binary Segmentation), "pelt" (Pruned Exact Linear Time)
#' @param scan_method Candidate-split scanning method: "score" (default; fast
#'   approximation using cumulative sufficient statistics from a single global
#'   fit) or "profile" (full refit at every candidate split; more expensive)
#' @param block_sigma Use block-specific variance for Gaussian (default: TRUE)
#' @param min_segment Minimum segment length (default: 5)
#' @param max_cps Maximum number of change points to detect (default: 10)
#' @param alpha Significance level (default: 0.05)
#' @param pvalue_method Method for p-values: "wilks", "permutation", "bootstrap"
#' @param use_cpp Use C++ acceleration (default: TRUE, highly recommended)
#' @param n_init Number of random initializations (default: 5)
#' @param max_iter Maximum VEM iterations (default: 100)
#' @param n_intervals Number of random intervals for WBS (default: 100)
#' @param compute_pvalues Compute p-values (default: FALSE; adds extra fits)
#' @param compute_cis Compute confidence intervals (default: TRUE; adds extra fits)
#' @param compute_inference Deprecated: sets both compute_pvalues and compute_cis
#' @param verbose Print progress (default: FALSE)
#' @return Object of class "parrot_result" with detected change points and diagnostics
#' @export
#' @examples
#' \dontrun{
#' # Simulate unipartite binary network with change point
#' data <- simulate_sbm_cp(N = 30, T_len = 60, cp_time = 30, Q = 2,
#'                         type = "unipartite", distribution = "bernoulli",
#'                         theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
#'                         theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2)))
#' 
#' # Detect change point
#' result <- parrot(data$Y, Q = 2)
#' print(result)
#' plot(result)
#' }
parrot <- function(Y, Q = NULL, Q2 = NULL,
                   type = NULL,
                   distribution = "auto",
                   multiple = TRUE,
                   method = c("wbs", "bs", "pelt"),
                   scan_method = c("score", "profile"),
                   block_sigma = TRUE,
                   min_segment = 5,
                   max_cps = 10,
                   alpha = 0.05,
                   pvalue_method = c("wilks", "permutation", "bootstrap"),
                   use_cpp = TRUE,
                   n_init = 5,
                   max_iter = 100,
                   n_intervals = 100,
                   compute_pvalues = FALSE,
                   compute_cis = TRUE,
                   compute_inference = NULL,
                   verbose = FALSE) {
  
  method <- match.arg(method)
  scan_method <- match.arg(scan_method)
  pvalue_method <- match.arg(pvalue_method)
  
  # Validate input
  if (!is.array(Y) || length(dim(Y)) != 3) {
    stop("Y must be a 3-dimensional array (T x N1 x N2)")
  }
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  # Auto-detect network type
  if (is.null(type)) {
    type <- if (N1 == N2) "unipartite" else "bipartite"
    if (verbose) message(sprintf("Auto-detected network type: %s", type))
  }
  
  # Auto-detect distribution
  if (distribution == "auto") {
    unique_vals <- unique(as.vector(Y))
    is_binary <- all(unique_vals %in% c(0, 1, NA))
    distribution <- if (is_binary) "bernoulli" else "gaussian"
    if (verbose) message(sprintf("Auto-detected distribution: %s", distribution))
  }
  
  # Handle Q - can be NULL, single value, or range
  if (is.null(Q)) {
    # Auto-select from default range 2:4
    Q <- .auto_select_Q(Y, type, distribution, 2:min(4, floor(N1/3)), use_cpp, verbose)
    if (verbose) message(sprintf("Auto-selected Q = %d", Q))
  } else if (length(Q) > 1) {
    # Q is a range - select best via ICL
    Q <- .auto_select_Q(Y, type, distribution, Q, use_cpp, verbose)
    if (verbose) message(sprintf("Selected Q = %d from range", Q))
  }
  # else Q is a single value, use as-is
  
  # Handle Q2 for bipartite
  if (type == "unipartite") {
    Q2 <- Q
  } else if (is.null(Q2)) {
    Q2 <- Q
  } else if (length(Q2) > 1) {
    # Q2 is a range - select best via ICL
    Q2 <- .auto_select_Q2(Y, type, distribution, Q, Q2, use_cpp, verbose)
    if (verbose) message(sprintf("Selected Q2 = %d from range", Q2))
  }
  
  if (verbose) {
    message("\n=== PARROT Change Point Detection ===")
    message(sprintf("Network: %s %s (T=%d, N1=%d, N2=%d)", 
                    type, distribution, T_len, N1, N2))
    message(sprintf("Clusters: Q1=%d, Q2=%d", Q, Q2))
    message(sprintf("Method: %s, Scan: %s, Multiple CPs: %s", method, scan_method, multiple))
  }
  
  # =========================================================================
  # CALL DETECTION FUNCTIONS FROM detect_cp.R
  # =========================================================================
  
  if (multiple) {
    # Use detect_multiple_cp from detect_cp.R
    det_result <- detect_multiple_cp(
      Y = Y, Q = Q, Q2 = Q2, 
      type = type, distribution = distribution,
      method = method, scan_method = scan_method, block_sigma = block_sigma,
      use_cpp = use_cpp, min_segment = min_segment,
      max_cps = max_cps, alpha = alpha,
      max_iter = max_iter, tol = 1e-6,
      n_init = n_init, n_intervals = n_intervals, verbose = verbose
    )
  } else {
    # Use detect_single_cp from detect_cp.R
    det_result <- detect_single_cp(
      Y = Y, Q = Q, Q2 = Q2,
      type = type, distribution = distribution,
      scan_method = scan_method,
      block_sigma = block_sigma, use_cpp = use_cpp,
      min_segment = min_segment, alpha = alpha,
      max_iter = max_iter, tol = 1e-6,
      n_init = n_init, verbose = verbose
    )
  }
  if (!is.null(det_result$scan_method)) {
    scan_method <- det_result$scan_method
  }
  
  # =========================================================================
  # BUILD PARROT RESULT
  # =========================================================================
  
  # Extract CPs from detection result
  if (multiple) {
    cps <- det_result$cps
    segments <- det_result$segments
  } else {
    # Single CP result
    if (!is.null(det_result$cp_est) && det_result$significant) {
      cps <- det_result$cp_est
      segments <- list(
        list(start = 1, end = det_result$cp_est, fit = det_result$left_fit),
        list(start = det_result$cp_est + 1, end = T_len, fit = det_result$right_fit)
      )
    } else {
      cps <- integer(0)
      segments <- list(list(start = 1, end = T_len, fit = det_result$null_fit))
    }
  }
  
  # Build result object
  result <- list(
    cps = cps,
    n_cps = length(cps),
    segments = segments,
    call = match.call(),
    type = type,
    distribution = distribution,
    Q = Q,
    Q2 = Q2,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    multiple = multiple,
    method = method,
    scan_method = scan_method,
    block_sigma = block_sigma,
    min_segment = min_segment,
    alpha = alpha,
    Y = Y,  # Store Y for p-value/CI computation
    detection_result = det_result  # Store original detection result
  )
  
  # Copy over profile likelihood if available
  if (!is.null(det_result$profile_ll)) {
    result$profile_ll <- det_result$profile_ll
  }
  if (!is.null(det_result$score_gain_by_cp)) {
    result$score_gain_by_cp <- det_result$score_gain_by_cp
  }
  if (!is.null(det_result$score_gain_curve)) {
    result$score_gain_curve <- det_result$score_gain_curve
  }
  if (!is.null(det_result$score_gain_table)) {
    result$score_gain_table <- det_result$score_gain_table
  }
  
  # Copy ICL values if available
  if (!is.null(det_result$icl_h1)) result$icl_h1 <- det_result$icl_h1
  if (!is.null(det_result$icl_h0)) result$icl_h0 <- det_result$icl_h0
  if (!is.null(det_result$null_elbo)) result$icl_h0 <- det_result$null_elbo
  if (!is.null(det_result$alt_elbo)) result$icl_h1 <- det_result$alt_elbo
  
  # Deprecated switch: set both flags if provided
  if (!is.null(compute_inference)) {
    compute_pvalues <- compute_inference
    compute_cis <- compute_inference
  }
  
  # =========================================================================
  # COMPUTE P-VALUES AND CIs
  # =========================================================================
  
  if (compute_pvalues && length(cps) > 0) {
    if (verbose) message("Computing p-values...")
    # Compute p-values using internal function (supports permutation/bootstrap)
    result$pvalues <- .compute_pvalues_for_parrot(result, pvalue_method, verbose)
  }
  
  if (compute_cis && length(cps) > 0) {
    if (verbose) message("Computing confidence intervals...")
    # Compute CIs using profile likelihood
    result$cis <- .compute_cis_for_parrot(result, alpha, use_cpp, n_init, max_iter, verbose)
  }
  
  class(result) <- c("parrot_result", "list")
  
  if (verbose) {
    message("\n=== Results ===")
    if (length(result$cps) == 0) {
      message("No change points detected")
    } else {
      message(sprintf("Detected %d change point(s): %s", 
                      length(result$cps), paste(result$cps, collapse = ", ")))
    }
  }
  
  return(result)
}

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Auto-select number of clusters Q from a range
#' @param Y Network tensor
#' @param type Network type
#' @param distribution Edge distribution
#' @param Q_range Vector of Q values to try
#' @param use_cpp Use C++ acceleration
#' @param verbose Print progress
#' @return Best Q value based on ICL
#' @keywords internal
.auto_select_Q <- function(Y, type, distribution, Q_range, use_cpp, verbose) {
  dims <- dim(Y)
  N1 <- dims[2]
  
  # Filter Q_range to valid values
  max_Q <- floor(N1 / 3)
  Q_candidates <- Q_range[Q_range >= 2 & Q_range <= max_Q]
  
  if (length(Q_candidates) == 0) {
    Q_candidates <- 2
  }
  
  best_icl <- -Inf
  best_Q <- Q_candidates[1]
  
  if (verbose) message("Selecting Q via ICL...")
  
  for (Q in Q_candidates) {
    Q2 <- if (type == "unipartite") Q else Q
    
    fit <- tryCatch(
      fit_sbm(Y, Q = Q, Q2 = Q2, type = type, distribution = distribution, 
              block_sigma = FALSE, use_cpp = use_cpp,
              max_iter = 30, tol = 1e-4, n_init = 2, verbose = FALSE),
      error = function(e) list(icl = -Inf)
    )
    
    if (fit$icl > best_icl) {
      best_icl <- fit$icl
      best_Q <- Q
    }
    
    if (verbose) message(sprintf("  Q=%d: ICL=%.2f", Q, fit$icl))
  }
  
  return(best_Q)
}

#' Auto-select Q2 for bipartite networks from a range
#' @keywords internal
.auto_select_Q2 <- function(Y, type, distribution, Q1, Q2_range, use_cpp, verbose) {
  dims <- dim(Y)
  N2 <- dims[3]
  
  # Filter Q2_range to valid values
  max_Q2 <- floor(N2 / 3)
  Q2_candidates <- Q2_range[Q2_range >= 2 & Q2_range <= max_Q2]
  
  if (length(Q2_candidates) == 0) {
    Q2_candidates <- 2
  }
  
  best_icl <- -Inf
  best_Q2 <- Q2_candidates[1]
  
  if (verbose) message("Selecting Q2 via ICL...")
  
  for (Q2 in Q2_candidates) {
    fit <- tryCatch(
      fit_sbm(Y, Q = Q1, Q2 = Q2, type = type, distribution = distribution, 
              block_sigma = FALSE, use_cpp = use_cpp,
              max_iter = 30, tol = 1e-4, n_init = 2, verbose = FALSE),
      error = function(e) list(icl = -Inf)
    )
    
    if (fit$icl > best_icl) {
      best_icl <- fit$icl
      best_Q2 <- Q2
    }
    
    if (verbose) message(sprintf("  Q2=%d: ICL=%.2f", Q2, fit$icl))
  }
  
  return(best_Q2)
}

#' Compute degrees of freedom for penalty
#' @keywords internal
.penalty_params <- function(Q, Q2, type, distribution) {
  # Number of parameters in the model
  if (type == "unipartite") {
    n_params <- Q * (Q + 1) / 2  # Symmetric matrix
  } else {
    n_params <- Q * Q2
  }
  
  # Add variance parameter for Gaussian
  if (distribution == "gaussian") {
    n_params <- n_params + 1
  }
  
  return(n_params)
}

#' Compute p-values for parrot result
#' @keywords internal
.compute_pvalues_for_parrot <- function(result, method, verbose,
                                       n_perm = 100, permutation_scope = "local") {
  
  if (length(result$cps) == 0) return(list())
  n_perm <- as.integer(n_perm[1])
  if (!is.finite(n_perm) || n_perm < 1L) n_perm <- 100L
  permutation_scope <- tolower(trimws(as.character(permutation_scope[1])))
  if (!(permutation_scope %in% c("local", "global"))) permutation_scope <- "local"
  
  n_cps <- length(result$cps)
  sorted_cps <- sort(result$cps)
  all_boundaries <- c(0, sorted_cps, result$T_len)
  
  pvalues <- list()
  raw_pvals <- numeric(n_cps)
  
  # Fast path: if single-CP detection_result already has LRT stats, reuse them
  if (!is.null(result$detection_result) &&
      inherits(result$detection_result, "cp_detection") &&
      method == "wilks" &&
      !is.null(result$detection_result$lr_stat) &&
      !is.null(result$detection_result$df)) {
    
    det <- result$detection_result
    cp_val <- result$cps[1]
    pval_raw <- pchisq(det$lr_stat, df = det$df, lower.tail = FALSE)
    
    pvalues[[1]] <- list(
      cp = cp_val,
      lr_stat = det$lr_stat,
      df = det$df,
      pvalue_raw = pval_raw,
      pvalue_method = "wilks"
    )
    
    raw_pvals[1] <- pval_raw
    
    # Bonferroni and BH (single CP => same)
    pvalues[[1]]$pvalue_bonferroni <- pval_raw
    pvalues[[1]]$pvalue_bh <- pval_raw
    pvalues[[1]]$significant_bonf <- pval_raw < result$alpha
    pvalues[[1]]$significant_bh <- pval_raw < result$alpha
    
    return(pvalues)
  }
  
  # Degrees of freedom for LRT (used when refitting nulls)
  df <- .penalty_params(result$Q, result$Q2, result$type, result$distribution)
  min_seg_for_p <- as.integer(if (!is.null(result$min_segment)) result$min_segment else 3L)
  if (!is.finite(min_seg_for_p) || min_seg_for_p < 1L) min_seg_for_p <- 3L
  
  for (i in seq_len(n_cps)) {
    cp <- sorted_cps[i]
    
    # Get segment boundaries for this CP
    left_start <- all_boundaries[i] + 1
    right_end <- all_boundaries[i + 2]
    
    # Use stored segment fits if available
    seg_left <- result$segments[[i]]
    seg_right <- result$segments[[i + 1]]
    
    # Extract data for this CP region
    Y_combined <- result$Y[left_start:right_end, , , drop = FALSE]
    local_cp <- cp - left_start + 1
    
    if (!is.null(seg_left$fit) && !is.null(seg_right$fit)) {
      alt_ll <- seg_left$fit$elbo + seg_right$fit$elbo
      
      # Null: fit single model to combined segment
      null_fit <- tryCatch(
        fit_sbm(Y_combined, result$Q, result$Q2, result$type, result$distribution,
                result$block_sigma, TRUE, 50, 1e-4, 3, FALSE),
        error = function(e) list(elbo = alt_ll - 1)
      )
      null_ll <- null_fit$elbo
      
      lr_stat <- 2 * (alt_ll - null_ll)
      lr_stat <- max(0, lr_stat)
    } else {
      lr_stat <- 0
      null_fit <- NULL
    }
    
    # Compute p-value based on method
    perm_scope_used <- NA_character_
    perm_space_size <- NA_real_
    perm_n_used <- NA_integer_

    if (method == "wilks") {
      pval_raw <- pchisq(lr_stat, df = df, lower.tail = FALSE)
      pval_method <- "wilks"
      
    } else if (method == "permutation") {
      if (verbose) message(sprintf("  Computing permutation p-value for CP %d...", i))
      perm_info <- .permutation_pvalue(
        Y_combined, local_cp, lr_stat,
        result$Q, result$Q2, result$type,
        result$distribution, result$block_sigma,
        n_perm = n_perm,
        scope = permutation_scope,
        Y_full = result$Y,
        cp_global = as.integer(cp),
        left_start = as.integer(left_start),
        right_end = as.integer(right_end),
        min_seg = min_seg_for_p,
        verbose = verbose
      )
      pval_raw <- as.numeric(perm_info$pvalue)
      pval_method <- "permutation"
      perm_scope_used <- as.character(perm_info$scope)
      perm_space_size <- as.numeric(perm_info$perm_space_size)
      perm_n_used <- as.integer(perm_info$n_perm_used)
      
    } else if (method == "bootstrap") {
      if (verbose) message(sprintf("  Computing bootstrap p-value for CP %d...", i))
      pval_raw <- .bootstrap_pvalue(Y_combined, local_cp, lr_stat, null_fit,
                                     result$Q, result$Q2, result$type,
                                     result$distribution, result$block_sigma,
                                     n_boot = 100, verbose = verbose)
      pval_method <- "bootstrap"
      
    } else {
      pval_raw <- pchisq(lr_stat, df = df, lower.tail = FALSE)
      pval_method <- "wilks"
    }
    
    raw_pvals[i] <- pval_raw
    
    pvalues[[i]] <- list(
      cp = cp,
      lr_stat = lr_stat,
      df = df,
      pvalue_raw = pval_raw,
      pvalue_method = pval_method,
      permutation_scope = perm_scope_used,
      permutation_space_size = perm_space_size,
      permutation_n_perm = perm_n_used
    )
  }
  
  # Bonferroni adjustment
  pvals_bonf <- pmin(1, raw_pvals * n_cps)
  
  # BH (FDR) adjustment
  pvals_bh <- p.adjust(raw_pvals, method = "BH")
  
  # Add adjusted p-values to each entry
  for (i in seq_len(n_cps)) {
    pvalues[[i]]$pvalue_bonferroni <- pvals_bonf[i]
    pvalues[[i]]$pvalue_bh <- pvals_bh[i]
    pvalues[[i]]$significant_bonf <- pvals_bonf[i] < result$alpha
    pvalues[[i]]$significant_bh <- pvals_bh[i] < result$alpha
  }
  
  return(pvalues)
}

#' Permutation test for CP p-value
#' @keywords internal
.permutation_pvalue <- function(Y, cp, observed_lr, Q, Q2, type, distribution, 
                                 block_sigma, n_perm = 100, scope = "local",
                                 Y_full = NULL, cp_global = NULL,
                                 left_start = NULL, right_end = NULL,
                                 min_seg = 3L,
                                 verbose = FALSE) {
  min_seg <- as.integer(min_seg[1])
  if (!is.finite(min_seg) || min_seg < 1L) min_seg <- 3L
  n_perm <- as.integer(n_perm[1])
  if (!is.finite(n_perm) || n_perm < 1L) n_perm <- 100L
  scope <- tolower(trimws(as.character(scope[1])))
  if (!(scope %in% c("local", "global"))) scope <- "local"

  factorial_size <- function(n) {
    n <- as.integer(n[1])
    if (!is.finite(n) || n < 0L) return(NA_real_)
    if (n > 170L) return(Inf)
    exp(lgamma(n + 1))
  }

  T_local <- dim(Y)[1]
  cp_local <- as.integer(cp[1])
  perm_space_size <- factorial_size(T_local)

  if (scope == "global") {
    if (is.null(Y_full)) {
      if (verbose) message("    Global permutation requested but Y_full is NULL; falling back to local scope.")
      scope <- "local"
    } else {
      T_full <- dim(Y_full)[1]
      cp_g <- as.integer(cp_global[1])
      l_start <- as.integer(left_start[1])
      r_end <- as.integer(right_end[1])

      if (!is.finite(cp_g)) cp_g <- cp_local
      if (!is.finite(l_start)) l_start <- cp_g - cp_local + 1L
      if (!is.finite(r_end)) r_end <- l_start + T_local - 1L

      l_start <- max(1L, l_start)
      r_end <- min(T_full, r_end)
      if (l_start > r_end || cp_g < l_start || cp_g >= r_end) {
        if (verbose) message("    Invalid global permutation bounds; falling back to local scope.")
        scope <- "local"
      } else {
        cp_local <- cp_g - l_start + 1L
        T_local <- r_end - l_start + 1L
        perm_space_size <- factorial_size(T_full)
      }
    }
  }

  if (!is.finite(cp_local) || cp_local < min_seg || T_local - cp_local < min_seg) {
    return(list(
      pvalue = NA_real_,
      scope = scope,
      perm_space_size = perm_space_size,
      n_perm_used = n_perm,
      count_extreme = NA_integer_
    ))
  }

  count_extreme <- 0L

  for (b in seq_len(n_perm)) {
    if (scope == "global") {
      perm_idx <- sample(dim(Y_full)[1])
      Y_perm_full <- Y_full[perm_idx, , , drop = FALSE]
      Y_perm <- Y_perm_full[l_start:r_end, , , drop = FALSE]
    } else {
      perm_idx <- sample(T_local)
      Y_perm <- Y[perm_idx, , , drop = FALSE]
    }
    
    null_fit <- tryCatch(
      fit_sbm(Y_perm, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
      error = function(e) list(elbo = -Inf)
    )
    
    Y_l <- Y_perm[1:cp_local, , , drop = FALSE]
    Y_r <- Y_perm[(cp_local + 1):T_local, , , drop = FALSE]
    
    fit_l <- tryCatch(
      fit_sbm(Y_l, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
      error = function(e) list(elbo = -Inf)
    )
    fit_r <- tryCatch(
      fit_sbm(Y_r, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
      error = function(e) list(elbo = -Inf)
    )
    
    lr_perm <- 2 * (fit_l$elbo + fit_r$elbo - null_fit$elbo)
    lr_perm <- max(0, lr_perm)
    if (lr_perm >= observed_lr) count_extreme <- count_extreme + 1L
    
    if (verbose && b %% 20 == 0) message(sprintf("    Permutation %d/%d (%s scope)", b, n_perm, scope))
  }
  
  pval <- (count_extreme + 1) / (n_perm + 1)
  list(
    pvalue = pval,
    scope = scope,
    perm_space_size = perm_space_size,
    n_perm_used = n_perm,
    count_extreme = count_extreme
  )
}

#' Parametric bootstrap for CP p-value
#' @keywords internal  
.bootstrap_pvalue <- function(Y, cp, observed_lr, null_fit, Q, Q2, type, 
                               distribution, block_sigma, n_boot = 100, verbose = FALSE) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  min_seg <- 3
  
  if (is.null(null_fit) || is.null(null_fit$prob) && is.null(null_fit$mu)) {
    df <- if (type == "unipartite") Q * (Q + 1) / 2 else Q * Q2
    if (distribution == "gaussian") df <- df + 1
    return(pchisq(observed_lr, df = df, lower.tail = FALSE))
  }
  
  count_extreme <- 0
  
  for (b in 1:n_boot) {
    Y_boot <- .simulate_from_fit(null_fit, T_len, N1, N2, type, distribution)
    
    if (is.null(Y_boot)) next
    
    null_boot <- tryCatch(
      fit_sbm(Y_boot, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
      error = function(e) list(elbo = -Inf)
    )
    
    if (cp >= min_seg && T_len - cp >= min_seg) {
      Y_l <- Y_boot[1:cp, , , drop = FALSE]
      Y_r <- Y_boot[(cp+1):T_len, , , drop = FALSE]
      
      fit_l <- tryCatch(
        fit_sbm(Y_l, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
        error = function(e) list(elbo = -Inf)
      )
      fit_r <- tryCatch(
        fit_sbm(Y_r, Q, Q2, type, distribution, block_sigma, TRUE, 30, 1e-3, 2, FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      lr_boot <- 2 * (fit_l$elbo + fit_r$elbo - null_boot$elbo)
      lr_boot <- max(0, lr_boot)
      
      if (lr_boot >= observed_lr) count_extreme <- count_extreme + 1
    }
    
    if (verbose && b %% 20 == 0) message(sprintf("    Bootstrap %d/%d", b, n_boot))
  }
  
  pval <- (count_extreme + 1) / (n_boot + 1)
  return(pval)
}

#' Simulate network from fitted model
#' @keywords internal
.simulate_from_fit <- function(fit, T_len, N1, N2, type, distribution) {
  if (is.null(fit)) return(NULL)
  
  if (!is.null(fit$tau_row)) {
    Z_row <- apply(fit$tau_row, 1, which.max)
  } else if (!is.null(fit$tau1)) {
    Z_row <- apply(fit$tau1, 1, which.max)
  } else if (!is.null(fit$tau)) {
    Z_row <- apply(fit$tau, 1, which.max)
  } else {
    return(NULL)
  }
  
  if (type == "bipartite") {
    if (!is.null(fit$tau_col)) {
      Z_col <- apply(fit$tau_col, 1, which.max)
    } else if (!is.null(fit$tau2)) {
      Z_col <- apply(fit$tau2, 1, which.max)
    } else {
      Z_col <- Z_row
    }
  } else {
    Z_col <- Z_row
  }
  
  Y <- array(0, dim = c(T_len, N1, N2))
  
  for (t in 1:T_len) {
    for (i in 1:N1) {
      j_start <- if (type == "unipartite") i + 1 else 1
      if (j_start > N2) next
      
      for (j in j_start:N2) {
        q <- Z_row[i]
        r <- Z_col[j]
        
        if (distribution == "bernoulli") {
          prob <- if (!is.null(fit$prob)) fit$prob[q, r] else 0.5
          Y[t, i, j] <- rbinom(1, 1, prob)
        } else {
          mu <- if (!is.null(fit$mu)) fit$mu[q, r] else 0
          sigma <- if (!is.null(fit$sigma)) {
            if (is.matrix(fit$sigma) && nrow(fit$sigma) > 1) {
              fit$sigma[q, r]
            } else if (is.matrix(fit$sigma)) {
              fit$sigma[1, 1]
            } else {
              fit$sigma
            }
          } else if (!is.null(fit$sigma_global)) {
            fit$sigma_global
          } else {
            1
          }
          Y[t, i, j] <- rnorm(1, mu, sigma)
        }
        
        if (type == "unipartite") Y[t, j, i] <- Y[t, i, j]
      }
    }
  }
  
  return(Y)
}

#' Compute CIs for parrot result
#' @keywords internal
.compute_cis_for_parrot <- function(result, alpha, use_cpp, n_init, max_iter, verbose) {
  
  if (length(result$cps) == 0) return(list())
  min_seg <- as.integer(if (!is.null(result$min_segment)) result$min_segment else 3L)
  
  # Ensure CI always contains the reported CP.
  .ensure_cp_in_ci <- function(cp, ci_lower, ci_upper, ci_set = integer(0)) {
    lo <- as.integer(min(ci_lower, cp))
    hi <- as.integer(max(ci_upper, cp))
    set <- as.integer(ci_set)
    if (length(set) == 0) set <- seq.int(lo, hi)
    if (!(cp %in% set)) set <- sort(unique(c(set, cp)))
    list(lo = lo, hi = hi, set = set)
  }
  
  # For single-CP runs, reuse the CI from detect_single_cp() to avoid
  # stochastic re-fit drift between cp estimate and CI.
  if (!isTRUE(result$multiple) &&
      !is.null(result$detection_result) &&
      inherits(result$detection_result, "cp_detection")) {
    det <- result$detection_result
    cp <- as.integer(result$cps[1])
    
    ci_lower <- cp
    ci_upper <- cp
    if (!is.null(det$ci) && length(det$ci) >= 2 &&
        !is.na(det$ci[1]) && !is.na(det$ci[2])) {
      ci_lower <- as.integer(det$ci[1])
      ci_upper <- as.integer(det$ci[2])
    }
    
    ci_set <- det$ci_set
    if (is.null(ci_set) || length(ci_set) == 0) {
      ci_set <- seq.int(ci_lower, ci_upper)
    }
    
    fixed_ci <- .ensure_cp_in_ci(cp, ci_lower, ci_upper, ci_set)
    ci_method <- if (!is.null(det$scan_method)) {
      paste0(det$scan_method, "_reuse")
    } else {
      "profile_reuse"
    }
    
    return(list(list(
      cp = cp,
      ci_lower = fixed_ci$lo,
      ci_upper = fixed_ci$hi,
      profile_ll = det$profile_ll,
      ci_set = fixed_ci$set,
      ci_method = ci_method
    )))
  }
  
  # Fast path for score-scan mode: derive CI from score profile without per-split refits
  use_score_ci <- !is.null(result$scan_method) && identical(result$scan_method, "score")
  if (use_score_ci) {
    base_fit <- NULL
    if (!is.null(result$detection_result$null_fit)) {
      base_fit <- result$detection_result$null_fit
    } else if (!is.null(result$detection_result$global_fit)) {
      base_fit <- result$detection_result$global_fit
    }
    
    if (!is.null(base_fit)) {
      prep <- tryCatch(
        .prepare_score_scan(result$Y, base_fit, result$type, result$distribution, result$block_sigma),
        error = function(e) NULL
      )
      if (!is.null(prep)) {
        n_cps <- length(result$cps)
        sorted_cps <- sort(result$cps)
        all_boundaries <- c(0, sorted_cps, result$T_len)
        cis <- vector("list", n_cps)
        
        for (i in seq_len(n_cps)) {
          cp <- sorted_cps[i]
          left_start <- all_boundaries[i] + 1L
          right_end <- all_boundaries[i + 2L]
          local_cp <- cp - left_start + 1L
          T_local <- right_end - left_start + 1L
          scan_range <- min_seg:(T_local - min_seg)
          
          if (length(scan_range) < 1) {
            cis[[i]] <- list(cp = cp, ci_lower = cp, ci_upper = cp, profile_ll = NULL, ci_method = "score")
            next
          }
          
          profile_global <- .score_profile_interval(prep, start = left_start, end = right_end, min_segment = min_seg)
          profile_ll <- profile_global[left_start:right_end]
          valid_idx <- which(is.finite(profile_ll))
          
          if (length(valid_idx) == 0) {
            cis[[i]] <- list(
              cp = cp,
              ci_lower = cp,
              ci_upper = cp,
              profile_ll = profile_ll,
              ci_set = cp,
              ci_method = "score"
            )
            next
          }
          
          max_ll <- max(profile_ll[valid_idx])
          cutoff <- max_ll - 0.5 * qchisq(1 - alpha, df = 1)
          ci_set_local <- which(profile_ll >= cutoff)
          
          ci_lower_local <- if (length(ci_set_local) > 0) min(ci_set_local) else local_cp
          ci_upper_local <- if (length(ci_set_local) > 0) max(ci_set_local) else local_cp
          ci_lower_global <- left_start + ci_lower_local - 1L
          ci_upper_global <- left_start + ci_upper_local - 1L
          ci_set_global <- ci_set_local + left_start - 1L
          fixed_ci <- .ensure_cp_in_ci(cp, ci_lower_global, ci_upper_global, ci_set_global)
          
          cis[[i]] <- list(
            cp = cp,
            ci_lower = fixed_ci$lo,
            ci_upper = fixed_ci$hi,
            profile_ll = profile_ll,
            ci_set = fixed_ci$set,
            ci_method = "score"
          )
        }
        
        return(cis)
      }
    }
    
    if (verbose) {
      message("Score-based CI preparation failed; reverting to profile-refit CI.")
    }
  }
  
  n_cps <- length(result$cps)
  sorted_cps <- sort(result$cps)
  all_boundaries <- c(0, sorted_cps, result$T_len)
  
  cis <- list()
  
  for (i in seq_len(n_cps)) {
    cp <- sorted_cps[i]
    
    left_start <- all_boundaries[i] + 1
    right_end <- all_boundaries[i + 2]
    
    Y_region <- result$Y[left_start:right_end, , , drop = FALSE]
    local_cp <- cp - left_start + 1
    T_local <- dim(Y_region)[1]
    
    scan_range <- min_seg:(T_local - min_seg)
    
    if (length(scan_range) < 1) {
      cis[[i]] <- list(cp = cp, ci_lower = cp, ci_upper = cp, profile_ll = NULL)
      next
    }
    
    profile_ll <- rep(-Inf, T_local)
    
    for (t in scan_range) {
      Y_l <- Y_region[1:t, , , drop = FALSE]
      Y_r <- Y_region[(t+1):T_local, , , drop = FALSE]
      
      fit_l <- tryCatch(
        fit_sbm(Y_l, result$Q, result$Q2, result$type, result$distribution,
                result$block_sigma, use_cpp, max_iter %/% 2, 1e-4, 
                max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      fit_r <- tryCatch(
        fit_sbm(Y_r, result$Q, result$Q2, result$type, result$distribution,
                result$block_sigma, use_cpp, max_iter %/% 2, 1e-4,
                max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      profile_ll[t] <- fit_l$elbo + fit_r$elbo
    }
    
    max_ll <- max(profile_ll[is.finite(profile_ll)])
    cutoff <- max_ll - 0.5 * qchisq(1 - alpha, df = 1)
    ci_set <- which(profile_ll >= cutoff)
    
    ci_lower_local <- if (length(ci_set) > 0) min(ci_set) else local_cp
    ci_upper_local <- if (length(ci_set) > 0) max(ci_set) else local_cp
    
    ci_lower <- left_start + ci_lower_local - 1
    ci_upper <- left_start + ci_upper_local - 1
    ci_set_global <- ci_set + left_start - 1
    fixed_ci <- .ensure_cp_in_ci(cp, ci_lower, ci_upper, ci_set_global)
    
    cis[[i]] <- list(
      cp = cp,
      ci_lower = fixed_ci$lo,
      ci_upper = fixed_ci$hi,
      profile_ll = profile_ll,
      ci_set = fixed_ci$set,
      ci_method = "profile"
    )
  }
  
  return(cis)
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' Print method for parrot_result
#' @param x A `parrot_result` object returned by `parrot()`.
#' @param ... Additional arguments (currently unused).
#' @export
print.parrot_result <- function(x, ...) {
  cat("\n=== PARROT Change Point Detection Results ===\n\n")
  
  cat(sprintf("Network: %s %s\n", x$type, x$distribution))
  cat(sprintf("Dimensions: T=%d, N1=%d, N2=%d\n", x$T_len, x$N1, x$N2))
  cat(sprintf("Clusters: Q=%d", x$Q))
  if (x$type == "bipartite") cat(sprintf(", Q2=%d", x$Q2))
  cat("\n\n")
  
  if (x$n_cps == 0) {
    cat("No change points detected.\n")
  } else {
    cat(sprintf("Detected %d change point(s):\n\n", x$n_cps))
    
    cat(sprintf("  %-4s  %-6s  %-10s  %-14s  %-14s  %-14s  %-12s  %-10s\n",
                "CP", "Time", "LR_stat", "p_raw", "p_BH", "p_Bonferroni", "95% CI", "Signif"))
    cat(paste0("  ", paste(rep("-", 98), collapse = ""), "\n"))
    
    for (i in seq_along(x$cps)) {
      cp <- x$cps[i]
      
      if (!is.null(x$pvalues) && i <= length(x$pvalues)) {
        pv <- x$pvalues[[i]]
        lr_stat <- sprintf("%.2f", pv$lr_stat)
        p_raw <- format(pv$pvalue_raw, scientific = TRUE, digits = 3)
        p_bh <- format(pv$pvalue_bh, scientific = TRUE, digits = 3)
        p_bonf <- format(pv$pvalue_bonferroni, scientific = TRUE, digits = 3)
        sig <- ifelse(isTRUE(pv$significant_bh), "YES", "NO")
      } else {
        lr_stat <- "NA"
        p_raw <- "NA"
        p_bh <- "NA"
        p_bonf <- "NA"
        sig <- "NA"
      }
      
      if (!is.null(x$cis) && i <= length(x$cis)) {
        ci <- x$cis[[i]]
        ci_str <- sprintf("[%d,%d]", ci$ci_lower, ci$ci_upper)
      } else {
        ci_str <- "NA"
      }
      
      cat(sprintf("  %-4d  %-6d  %-10s  %-14s  %-14s  %-14s  %-12s  %-10s\n",
                  i, cp, lr_stat, p_raw, p_bh, p_bonf, ci_str, sig))
    }
    cat("\n")
  }
  
  cat(sprintf("Segments: %d\n", length(x$segments)))
  for (s in seq_along(x$segments)) {
    seg <- x$segments[[s]]
    cat(sprintf("  Segment %d: [%d, %d] (length %d)\n", 
                s, seg$start, seg$end, seg$end - seg$start + 1))
  }
  
  invisible(x)
}

#' Summary method for parrot_result
#' @param object A `parrot_result` object returned by `parrot()`.
#' @param ... Additional arguments (currently unused).
#' @export
summary.parrot_result <- function(object, ...) {
  print.parrot_result(object)
  
  cat("\n--- Model Details ---\n")
  cat(sprintf("Method: %s\n", object$method))
  cat(sprintf("Block-specific variance: %s\n", object$block_sigma))
  cat(sprintf("Significance level: %.2f\n", object$alpha))
  
  if (!is.null(object$icl_h1)) {
    cat(sprintf("\nICL(H1): %.4f\n", object$icl_h1))
    cat(sprintf("ICL(H0): %.4f\n", object$icl_h0))
    cat(sprintf("Delta ICL: %.4f\n", object$icl_h1 - object$icl_h0))
  }
  
  invisible(object)
}

#' Plot method for parrot_result
#' @param x A `parrot_result` object returned by `parrot()`.
#' @param ... Additional plotting arguments passed to `plot_detection()`.
#' @export
plot.parrot_result <- function(x, ...) {
  plot_detection(x, ...)
}
