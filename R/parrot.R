#' PARROT: Phase-Altering Regulatory Rewiring Over Time
#'
#' Main function for change point detection in dynamic networks.
#' Provides a unified interface for detecting single or multiple change points.
#'
#' @param Y Network tensor (T x N x N for unipartite, T x N1 x N2 for bipartite)
#' @param Q Number of communities (Q1 for bipartite). NULL (default): selected
#'   with the ICL gain-ratio rule over \code{Q_range}; a vector: selected from
#'   these values; an integer: used as is.
#' @param Q2 Number of column communities for bipartite networks (default: Q).
#' @param Q_range Candidate values for the automatic selection of Q
#'   (default 1:10, truncated to at most N/2 nodes per community side).
#' @param q_stability Also report, for each candidate Q, the mean pairwise
#'   adjusted Rand index of the memberships fitted from \code{q_n_seeds}
#'   independent sets of initialisations (default FALSE; adds fits).
#' @param q_n_seeds Number of independent fits for \code{q_stability}.
#' @param type Network type: "unipartite" or "bipartite" (auto-detected if NULL)
#' @param distribution Edge distribution: "gaussian", "bernoulli", or "auto"
#' @param multiple Detect multiple change points (TRUE) or a single one (FALSE)
#' @param method Candidate search for multiple change points: "wbs", "bs" or "pelt"
#' @param stop_rule Multiple change points: "max_cps" (default; up to
#'   \code{max_cps} candidates from \code{method}, each then tested and
#'   Bonferroni/BH adjusted) or "test" (test-based binary segmentation: a
#'   segment is split while the p-value of its best split is below
#'   \code{alpha}, at most \code{max_cps} change points).
#' @param scan_method Candidate-split scan: "profile" (default; both segments
#'   refitted at every split, i.e. the profile likelihood) or "score" (fast
#'   approximation from cumulative sufficient statistics of one global fit)
#' @param block_sigma Block-specific variances for Gaussian edges (default
#'   FALSE: one variance shared by all blocks)
#' @param min_segment Minimum segment length (default: 5)
#' @param max_cps Maximum number of change points (default: 3)
#' @param alpha Significance level (default: 0.05)
#' @param pvalue_method "permutation" (default; time-label permutation test
#'   conditional on the window between the neighbouring change points, with the
#'   scan repeated on every permuted series; exact enumeration when feasible),
#'   "wilks" (chi-squared approximation; fastest), "stationary_bootstrap"
#'   (block resampling that preserves short-range temporal dependence) or
#'   "bootstrap" (parametric). In single-CP mode the change point is declared
#'   when this p-value is below \code{alpha}.
#' @param permutation_scope "global" (default; permute the time labels of the
#'   whole series) or "local" (permute only inside the window between the
#'   neighbouring change points)
#' @param rescan Repeat the scan within the window on every permuted /
#'   bootstrapped series (default TRUE; controls the type-I error). FALSE holds
#'   the split fixed at the detected position, which ignores that the split was
#'   selected by the scan and is anti-conservative (PARROT < 1.1.0).
#' @param exact_max Enumerate all distinct permutation configurations when
#'   there are at most this many (default 5040, i.e. all orderings of T = 7);
#'   otherwise Monte Carlo
#' @param fit_seed Seed of the segment fits (default 1). Every SBM fit of a
#'   segment uses its own seed derived from the segment's time points and
#'   \code{fit_seed}, so the result does not depend on the global random state
#'   and the scan and the tests use the same fits; change it to check the
#'   sensitivity to the random initialisations.
#' @param n_cores Number of cores for the segment fits of the permutation /
#'   bootstrap test (default 1)
#' @param n_perm Number of random permutations or bootstrap samples when the
#'   p-value is not exact (default: 99, enough for a decision at alpha = 0.05; the
#'   smallest attainable p-value is 1/(n_perm + 1), so use e.g. 999 to report a p-value)
#' @param sb_block Mean block length of the stationary bootstrap (default
#'   ceiling(T^(1/3)))
#' @param use_cpp Use C++ acceleration (default: TRUE)
#' @param n_init Number of random initialisations per fit (default: 20)
#' @param max_iter Maximum VEM iterations (default: 200)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param n_intervals Number of random intervals for WBS (default: 50)
#' @param compute_pvalues Compute p-values (default: TRUE)
#' @param compute_cis Compute confidence intervals (default: TRUE)
#' @param compute_inference Deprecated: sets both compute_pvalues and compute_cis
#' @param verbose Print progress (default: FALSE)
#' @section Speed versus accuracy:
#' The defaults favour accuracy. Fast options: \code{scan_method = "score"}
#' (one fit of the whole series instead of two segment fits per split; can miss
#' changes that leave the whole-series fit without block structure),
#' \code{pvalue_method = "wilks"} (no resampling; conservative), smaller
#' \code{n_init} / \code{max_iter}. Accurate options (defaults):
#' \code{scan_method = "profile"}, \code{pvalue_method = "permutation"} with
#' \code{rescan = TRUE} (the scan is repeated on every permuted series; the
#' cost grows with \code{n_perm} or with the number of enumerated
#' configurations), \code{n_init = 20}, \code{max_iter = 200}.
#' \code{rescan = FALSE} and \code{pvalue_method = "bootstrap"} are fast but
#' anti-conservative. \code{n_cores} parallelises the segment fits without
#' changing the result. A quick exploratory run is
#' \code{parrot(Y, scan_method = "score", pvalue_method = "wilks")}; a fixed
#' \code{Q} and \code{compute_cis = FALSE} also save time. The package README
#' ("Speed versus accuracy") orders the p-value methods by type-I error.
#' @return Object of class "parrot_result" with detected change points and diagnostics
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_sbm_cp(N = 30, T_len = 20, cp_time = 10, Q = 2,
#'                         type = "unipartite", distribution = "bernoulli",
#'                         theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
#'                         theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2)))
#' result <- parrot(data$Y, Q = 2, multiple = FALSE, min_segment = 2)
#' print(result)
#' # fast settings
#' fast <- parrot(data$Y, Q = 2, multiple = FALSE, min_segment = 2,
#'                scan_method = "score", pvalue_method = "wilks", n_init = 5)
#' }
parrot <- function(Y, Q = NULL, Q2 = NULL,
                   Q_range = 1:10,
                   q_stability = FALSE,
                   q_n_seeds = 5,
                   type = NULL,
                   distribution = "auto",
                   multiple = TRUE,
                   method = c("wbs", "bs", "pelt"),
                   stop_rule = c("max_cps", "test"),
                   scan_method = c("profile", "score"),
                   block_sigma = FALSE,
                   min_segment = 5,
                   max_cps = 3,
                   alpha = 0.05,
                   pvalue_method = c("permutation", "wilks", "stationary_bootstrap", "bootstrap"),
                   permutation_scope = c("global", "local"),
                   rescan = TRUE,
                   exact_max = 5040,
                   n_perm = 99,
                   n_cores = 1,
                   fit_seed = 1L,
                   sb_block = NULL,
                   use_cpp = TRUE,
                   n_init = 20,
                   max_iter = 200,
                   tol = 1e-6,
                   n_intervals = 50,
                   compute_pvalues = TRUE,
                   compute_cis = TRUE,
                   compute_inference = NULL,
                   verbose = FALSE) {

  method <- match.arg(method)
  stop_rule <- match.arg(stop_rule)
  scan_method <- match.arg(scan_method)
  pvalue_method <- match.arg(pvalue_method)
  permutation_scope <- match.arg(permutation_scope)

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

  # Q: NULL or a range -> ICL gain-ratio rule; a single value -> as is
  q_selection <- NULL
  if (is.null(Q) || length(Q) > 1) {
    cand <- if (is.null(Q)) Q_range else Q
    cand <- cand[cand >= 1 & cand <= max(1, floor(min(N1, N2) / 2))]
    q_selection <- .select_Q_icl(Y, type, distribution, cand, block_sigma, use_cpp,
                                 n_init, max_iter, tol, stability = q_stability,
                                 n_seeds = q_n_seeds, verbose = verbose)
    Q <- max(1L, q_selection$Q)
    if (verbose) message(sprintf("Selected Q = %d (ICL gain-ratio rule)", Q))
  }

  # Handle Q2 for bipartite
  if (type == "unipartite") {
    Q2 <- Q
  } else if (is.null(Q2)) {
    Q2 <- Q
  } else if (length(Q2) > 1) {
    Q2 <- .auto_select_Q2(Y, type, distribution, Q, Q2, use_cpp, verbose)
    if (verbose) message(sprintf("Selected Q2 = %d from range", Q2))
  }

  if (verbose) {
    message("\n=== PARROT Change Point Detection ===")
    message(sprintf("Network: %s %s (T=%d, N1=%d, N2=%d)",
                    type, distribution, T_len, N1, N2))
    message(sprintf("Clusters: Q1=%d, Q2=%d", Q, Q2))
    message(sprintf("Method: %s, Scan: %s, Multiple CPs: %s, stop rule: %s, p-value: %s",
                    method, scan_method, multiple, stop_rule, pvalue_method))
  }

  # =========================================================================
  # DETECTION
  # =========================================================================

  test_info <- NULL
  if (multiple && stop_rule == "test") {
    tb <- .test_based_bs(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                         min_segment, max_cps, alpha, max_iter, tol, n_init,
                         if (pvalue_method == "stationary_bootstrap") "stationary_bootstrap" else "permutation",
                         n_perm, permutation_scope, exact_max, sb_block, verbose,
                         rescan = rescan, n_cores = n_cores, fit_seed = fit_seed)
    det_result <- list(
      cps = tb$cps, n_cps = length(tb$cps),
      segments = .fit_final_segments(Y, tb$cps, Q, Q2, type, distribution,
                                     block_sigma, use_cpp, max_iter, tol, n_init),
      tests = tb$tests, method = "test_based_bs", scan_method = "profile",
      min_segment = min_segment, type = type, distribution = distribution,
      T_len = T_len, N1 = N1, N2 = N2, Q = Q, Q2 = Q2, block_sigma = block_sigma)
    class(det_result) <- c("multi_cp_detection", "list")
    test_info <- tb$tests
  } else if (multiple) {
    det_result <- detect_multiple_cp(
      Y = Y, Q = Q, Q2 = Q2,
      type = type, distribution = distribution,
      method = method, scan_method = scan_method, block_sigma = block_sigma,
      use_cpp = use_cpp, min_segment = min_segment,
      max_cps = max_cps, alpha = alpha,
      max_iter = max_iter, tol = tol,
      n_init = n_init, n_intervals = n_intervals, verbose = verbose,
      fit_seed = fit_seed
    )
  } else {
    det_result <- detect_single_cp(
      Y = Y, Q = Q, Q2 = Q2,
      type = type, distribution = distribution,
      scan_method = scan_method,
      block_sigma = block_sigma, use_cpp = use_cpp,
      min_segment = min_segment, alpha = alpha,
      max_iter = max_iter, tol = tol,
      n_init = n_init, verbose = verbose, fit_seed = fit_seed
    )
    # decision by the chosen test (Wilks is computed inside detect_single_cp)
    if (pvalue_method %in% c("permutation", "stationary_bootstrap") &&
        !is.null(det_result$cp_raw) && is.finite(det_result$cp_raw)) {
      cp_raw <- as.integer(det_result$cp_raw)
      te <- if (pvalue_method == "permutation") {
        .permutation_pvalue(Y, cp_raw, 1L, T_len, Q, Q2, type, distribution, block_sigma,
                            n_init, max_iter, tol, use_cpp, n_perm, permutation_scope,
                            exact_max, det_result$elbo_fun, verbose, min_segment = min_segment,
                            rescan = rescan, n_cores = n_cores, fit_seed = fit_seed)
      } else {
        .sb_pvalue(Y, cp_raw, 1L, T_len, Q, Q2, type, distribution, block_sigma,
                   n_init, max_iter, tol, use_cpp, n_perm, sb_block, det_result$elbo_fun, verbose,
                   min_segment = min_segment, rescan = rescan, n_cores = n_cores,
                   fit_seed = fit_seed)
      }
      det_result$pvalue_wilks <- det_result$pvalue
      det_result$pvalue <- te$pvalue
      det_result$pvalue_method <- pvalue_method
      det_result$test <- te
      det_result$significant <- te$pvalue < alpha
      det_result$cp_est <- if (det_result$significant) cp_raw else NA
      det_result$n_cps <- as.integer(det_result$significant)
    }
  }
  if (!is.null(det_result$scan_method)) {
    scan_method <- det_result$scan_method
  }

  # =========================================================================
  # BUILD PARROT RESULT
  # =========================================================================

  if (multiple) {
    cps <- det_result$cps
    segments <- det_result$segments
  } else {
    if (!is.null(det_result$cp_est) && isTRUE(det_result$significant)) {
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

  result <- list(
    cps = cps,
    n_cps = length(cps),
    segments = segments,
    call = match.call(),
    type = type,
    distribution = distribution,
    Q = Q,
    Q2 = Q2,
    q_selection = q_selection,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    multiple = multiple,
    method = method,
    stop_rule = if (multiple) stop_rule else NA_character_,
    scan_method = scan_method,
    block_sigma = block_sigma,
    min_segment = min_segment,
    alpha = alpha,
    pvalue_method = pvalue_method,
    rescan = rescan,
    n_cores = n_cores,
    fit_seed = fit_seed,
    n_init = n_init,
    max_iter = max_iter,
    tol = tol,
    use_cpp = use_cpp,
    Y = Y,
    detection_result = det_result
  )

  if (!is.null(det_result$profile_ll)) result$profile_ll <- det_result$profile_ll
  if (!is.null(det_result$score_gain_by_cp)) result$score_gain_by_cp <- det_result$score_gain_by_cp
  if (!is.null(det_result$score_gain_curve)) result$score_gain_curve <- det_result$score_gain_curve
  if (!is.null(det_result$score_gain_table)) result$score_gain_table <- det_result$score_gain_table
  if (!is.null(det_result$icl_h1)) result$icl_h1 <- det_result$icl_h1
  if (!is.null(det_result$icl_h0)) result$icl_h0 <- det_result$icl_h0
  if (!is.null(det_result$null_elbo)) result$icl_h0 <- det_result$null_elbo
  if (!is.null(det_result$alt_elbo)) result$icl_h1 <- det_result$alt_elbo

  if (!is.null(compute_inference)) {
    compute_pvalues <- compute_inference
    compute_cis <- compute_inference
  }

  # =========================================================================
  # P-VALUES AND CIs
  # =========================================================================

  if (compute_pvalues && length(cps) > 0) {
    if (!is.null(test_info)) {
      # test-based segmentation: report the tests of the accepted splits
      acc <- test_info[vapply(test_info, function(x) x$cp %in% cps && x$pvalue < alpha, logical(1))]
      result$pvalues <- lapply(acc[order(vapply(acc, `[[`, numeric(1), "cp"))], function(x) list(
        cp = x$cp, lr_stat = x$lr_stat, pvalue_raw = x$pvalue,
        pvalue_method = if (x$method == "stationary_bootstrap") "stationary_bootstrap" else "permutation",
        permutation_type = x$method, permutation_scope = x$scope,
        permutation_space_size = x$perm_space_size, permutation_n_perm = x$n_perm_used,
        window = x$segment))
    } else if (!multiple && !is.null(det_result$test)) {
      te <- det_result$test
      result$pvalues <- list(list(
        cp = cps[1], lr_stat = te$lr_stat, df = det_result$df, pvalue_raw = te$pvalue,
        pvalue_method = pvalue_method, permutation_type = te$method, rescan = rescan,
        permutation_scope = te$scope, permutation_space_size = te$perm_space_size,
        permutation_n_perm = te$n_perm_used, window = c(1L, T_len), null = te$null,
        pvalue_bonferroni = te$pvalue, pvalue_bh = te$pvalue,
        significant_bonf = te$pvalue < alpha, significant_bh = te$pvalue < alpha))
    } else {
      if (verbose) message("Computing p-values...")
      result$pvalues <- .compute_pvalues_for_parrot(result, pvalue_method, verbose,
                                                    n_perm = n_perm,
                                                    permutation_scope = permutation_scope,
                                                    exact_max = exact_max, sb_block = sb_block,
                                                    rescan = rescan, n_cores = n_cores)
    }
    if (multiple && stop_rule == "max_cps" && length(result$pvalues) > 0) {
      result$cps_significant <- sort(vapply(Filter(function(x) isTRUE(x$significant_bonf),
                                                   result$pvalues), `[[`, numeric(1), "cp"))
    }
  }

  if (compute_cis && length(cps) > 0) {
    if (verbose) message("Computing confidence intervals...")
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

#' Compute p-values for parrot result
#'
#' For every reported change point the split is held fixed and tested inside the
#' window bounded by the neighbouring change points (or the ends of the series).
#' All statistics are computed with the fit settings stored in \code{result}
#' (\code{n_init}, \code{max_iter}, \code{tol}), identical for observed and
#' resampled series. Bonferroni and BH adjustments are over the reported
#' change points.
#' @keywords internal
.compute_pvalues_for_parrot <- function(result, method = "permutation", verbose = FALSE,
                                        n_perm = 99, permutation_scope = "global",
                                        exact_max = 5040, sb_block = NULL,
                                        n_init = NULL, max_iter = NULL, tol = NULL,
                                        rescan = NULL, n_cores = NULL) {

  if (length(result$cps) == 0) return(list())
  or <- function(a, b) if (is.null(a)) b else a
  method <- match.arg(method, c("permutation", "wilks", "stationary_bootstrap", "bootstrap"))
  permutation_scope <- match.arg(tolower(as.character(permutation_scope[1])), c("global", "local"))
  n_perm <- as.integer(n_perm[1])
  if (!is.finite(n_perm) || n_perm < 1L) n_perm <- 199L
  n_init <- or(n_init, or(result$n_init, 20L))
  max_iter <- or(max_iter, or(result$max_iter, 200L))
  tol <- or(tol, or(result$tol, 1e-6))
  use_cpp <- or(result$use_cpp, TRUE)
  block_sigma <- isTRUE(result$block_sigma)
  rescan <- isTRUE(or(rescan, or(result$rescan, TRUE)))
  n_cores <- as.integer(or(n_cores, or(result$n_cores, 1L)))
  min_seg <- as.integer(or(result$min_segment, 2L))

  # Single-CP Wilks: reuse the statistic of the detection step
  if (method == "wilks" && !isTRUE(result$multiple) &&
      inherits(result$detection_result, "cp_detection") &&
      !is.null(result$detection_result$lr_stat)) {
    det <- result$detection_result
    p <- pchisq(det$lr_stat, df = det$df, lower.tail = FALSE)
    return(list(list(cp = result$cps[1], lr_stat = det$lr_stat, df = det$df,
                     pvalue_raw = p, pvalue_method = "wilks",
                     pvalue_bonferroni = p, pvalue_bh = p,
                     significant_bonf = p < result$alpha, significant_bh = p < result$alpha)))
  }

  sorted_cps <- sort(as.integer(result$cps))
  n_cps <- length(sorted_cps)
  bnd <- c(0L, sorted_cps, as.integer(result$T_len))
  df <- .compute_df(result$Q, result$Q2, result$type, result$distribution, block_sigma)
  same_settings <- identical(as.integer(n_init), as.integer(or(result$n_init, -1L))) &&
    identical(as.integer(max_iter), as.integer(or(result$max_iter, -1L))) &&
    isTRUE(all.equal(tol, or(result$tol, -1)))
  E <- if (same_settings && is.function(result$detection_result$elbo_fun)) {
    result$detection_result$elbo_fun
  } else {
    .segment_elbo_cache(result$Y, result$Q, result$Q2, result$type, result$distribution,
                        block_sigma, use_cpp, max_iter, tol, n_init,
                        seed = as.integer(or(result$fit_seed, 1L)))
  }

  pvalues <- vector("list", n_cps)
  raw <- numeric(n_cps)
  for (i in seq_len(n_cps)) {
    cp <- sorted_cps[i]; ls <- bnd[i] + 1L; re <- bnd[i + 2L]
    entry <- list(cp = cp, df = df, window = c(ls, re),
                  permutation_type = NA_character_, permutation_scope = NA_character_,
                  permutation_space_size = NA_real_, permutation_n_perm = NA_integer_)
    if (verbose) message(sprintf("  p-value for CP %d (window %d-%d, %s)...", cp, ls, re, method))
    if (method == "permutation") {
      te <- .permutation_pvalue(result$Y, cp, ls, re, result$Q, result$Q2, result$type,
                                result$distribution, block_sigma, n_init, max_iter, tol,
                                use_cpp, n_perm, permutation_scope, exact_max, E, verbose,
                                min_segment = min_seg, rescan = rescan, n_cores = n_cores)
      entry$lr_stat <- te$lr_stat; entry$pvalue_raw <- te$pvalue; entry$pvalue_method <- "permutation"
      entry$rescan <- rescan
      entry$permutation_type <- te$method; entry$permutation_scope <- te$scope
      entry$permutation_space_size <- te$perm_space_size; entry$permutation_n_perm <- te$n_perm_used
      entry$null <- te$null
    } else if (method == "stationary_bootstrap") {
      te <- .sb_pvalue(result$Y, cp, ls, re, result$Q, result$Q2, result$type,
                       result$distribution, block_sigma, n_init, max_iter, tol, use_cpp,
                       n_perm, sb_block, E, verbose, min_segment = min_seg, rescan = rescan,
                       n_cores = n_cores)
      entry$lr_stat <- te$lr_stat; entry$pvalue_raw <- te$pvalue; entry$rescan <- rescan
      entry$pvalue_method <- "stationary_bootstrap"; entry$permutation_n_perm <- te$n_perm_used
      entry$null <- te$null
    } else {
      lr <- .lr_from_cache(E, ls:cp, (cp + 1L):re)
      entry$lr_stat <- lr
      if (method == "wilks") {
        entry$pvalue_raw <- pchisq(lr, df = df, lower.tail = FALSE)
        entry$pvalue_method <- "wilks"
      } else {
        Yw <- result$Y[ls:re, , , drop = FALSE]
        null_fit <- tryCatch(fit_sbm(Yw, result$Q, result$Q2, result$type, result$distribution,
                                     block_sigma, use_cpp, max_iter, tol, n_init, FALSE),
                             error = function(e) NULL)
        entry$pvalue_raw <- .bootstrap_pvalue(Yw, cp - ls + 1L, lr, null_fit, result$Q, result$Q2,
                                              result$type, result$distribution, block_sigma,
                                              n_boot = n_perm, verbose = verbose)
        entry$pvalue_method <- "bootstrap"
      }
    }
    raw[i] <- entry$pvalue_raw
    pvalues[[i]] <- entry
  }

  bonf <- pmin(1, raw * n_cps)
  bh <- p.adjust(raw, method = "BH")
  for (i in seq_len(n_cps)) {
    pvalues[[i]]$pvalue_bonferroni <- bonf[i]
    pvalues[[i]]$pvalue_bh <- bh[i]
    pvalues[[i]]$significant_bonf <- bonf[i] < result$alpha
    pvalues[[i]]$significant_bh <- bh[i] < result$alpha
  }
  pvalues
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
