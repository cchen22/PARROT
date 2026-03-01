#' Change Point Detection Functions for PARROT
#'
#' Single and multiple change point detection for dynamic networks.
#' Uses Binary Segmentation (BS) or Wild Binary Segmentation (WBS).

#' Detect Single Change Point
#'
#' Scan for the optimal single change point location using profile likelihood.
#'
#' @param Y Network tensor (T x N1 x N2)
#' @param Q Number of clusters (or Q1 for bipartite)
#' @param Q2 Number of column clusters (NULL for unipartite)
#' @param type Network type: "unipartite" or "bipartite"
#' @param distribution Edge distribution: "gaussian" or "bernoulli"
#' @param scan_method Candidate-split scan method: "profile" (full segment refits)
#'   or "score" (fast approximation using a fixed global fit)
#' @param block_sigma Use block-specific variance (for Gaussian)
#' @param use_cpp Use C++ acceleration
#' @param min_segment Minimum segment length
#' @param alpha Significance level for CI
#' @param max_iter Maximum VEM iterations
#' @param tol Convergence tolerance
#' @param n_init Number of random initializations
#' @param verbose Print progress
#' @return Object of class "cp_detection" with results
#' @export
#' @examples
#' \dontrun{
#' # Simulate and detect
#' data <- simulate_sbm_cp(N = 30, T_len = 60, cp_time = 30, Q = 2,
#'                         type = "unipartite", distribution = "bernoulli",
#'                         theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
#'                         theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2)))
#' result <- detect_single_cp(data$Y, Q = 2, type = "unipartite", 
#'                           distribution = "bernoulli")
#' print(result$cp_est)
#' }
detect_single_cp <- function(Y, Q, Q2 = NULL,
                             type = c("unipartite", "bipartite"),
                             distribution = c("gaussian", "bernoulli"),
                             scan_method = c("profile", "score"),
                             block_sigma = FALSE,
                             use_cpp = TRUE,
                             min_segment = 5,
                             alpha = 0.05,
                             max_iter = 100,
                             tol = 1e-6,
                             n_init = 5,
                             verbose = FALSE) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  scan_method <- match.arg(scan_method)
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  if (type == "unipartite") {
    stopifnot(N1 == N2)
    Q2 <- Q
  } else {
    if (is.null(Q2)) Q2 <- Q
  }
  
  # Scan range
  scan_start <- min_segment
  scan_end <- T_len - min_segment
  
  if (scan_end <= scan_start) {
    warning("Time series too short for change point detection")
    return(.create_cp_result(
      cp = NA, significant = FALSE, n_cps = 0,
      type = type, distribution = distribution,
      T_len = T_len, N1 = N1, N2 = N2, Q = Q, Q2 = Q2
    ))
  }
  
  if (verbose) message("Fitting null model (no change point)...")
  
  # Fit null model
  null_fit <- fit_sbm(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                      max_iter, tol, n_init * 2, FALSE)
  
  # Scan candidate split points
  profile_ll <- rep(-Inf, T_len)
  best_cp <- NA_integer_
  best_scan_ll <- -Inf
  best_fits <- list(left = NULL, right = NULL)
  
  if (verbose) {
    message(sprintf("Scanning t = %d to %d using %s scan...", scan_start, scan_end, scan_method))
  }
  
  if (scan_method == "score") {
    prep <- tryCatch(
      .prepare_score_scan(Y, null_fit, type, distribution, block_sigma),
      error = function(e) NULL
    )
    
    if (!is.null(prep)) {
      profile_ll <- .score_profile_interval(prep, start = 1L, end = T_len, min_segment = min_segment)
      idx <- which(is.finite(profile_ll))
      if (length(idx) > 0) {
        best_cp <- as.integer(idx[which.max(profile_ll[idx])])
        best_scan_ll <- as.numeric(profile_ll[best_cp])
      }
    }
    
    if (!is.finite(best_scan_ll) || !is.finite(best_cp)) {
      warning("Score scan failed; falling back to profile scan refits.")
      scan_method <- "profile"
    }
  }
  
  if (scan_method == "profile") {
    for (t in scan_start:scan_end) {
      Y_left <- Y[1:t, , , drop = FALSE]
      Y_right <- Y[(t + 1):T_len, , , drop = FALSE]
      
      fit_l <- tryCatch(
        fit_sbm(Y_left, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter, tol, n_init, FALSE),
        error = function(e) list(elbo = -Inf)
      )
      fit_r <- tryCatch(
        fit_sbm(Y_right, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter, tol, n_init, FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      ll <- fit_l$elbo + fit_r$elbo
      profile_ll[t] <- ll
      
      if (ll > best_scan_ll) {
        best_scan_ll <- ll
        best_cp <- as.integer(t)
        best_fits$left <- fit_l
        best_fits$right <- fit_r
      }
      
      if (verbose && t %% 10 == 0) {
        message(sprintf("  t = %d: LL = %.2f (best = %d)", t, ll, best_cp))
      }
    }
  }
  
  if (!is.finite(best_cp)) {
    warning("No valid candidate split found.")
    return(.create_cp_result(
      cp = NA, significant = FALSE, n_cps = 0,
      type = type, distribution = distribution,
      T_len = T_len, N1 = N1, N2 = N2, Q = Q, Q2 = Q2
    ))
  }
  
  if (is.null(best_fits$left) || is.null(best_fits$right)) {
    Y_left <- Y[1:best_cp, , , drop = FALSE]
    Y_right <- Y[(best_cp + 1):T_len, , , drop = FALSE]
    
    best_fits$left <- tryCatch(
      fit_sbm(Y_left, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter, tol, n_init, FALSE),
      error = function(e) list(elbo = -Inf)
    )
    best_fits$right <- tryCatch(
      fit_sbm(Y_right, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter, tol, n_init, FALSE),
      error = function(e) list(elbo = -Inf)
    )
  }
  
  alt_ll <- as.numeric(best_fits$left$elbo + best_fits$right$elbo)
  
  # Compute test statistic and significance
  lr_stat <- 2 * (alt_ll - null_fit$elbo)
  lr_stat <- max(0, lr_stat)
  
  # Degrees of freedom for LRT
  df <- .compute_df(Q, Q2, type, distribution, block_sigma)
  pvalue <- pchisq(lr_stat, df = df, lower.tail = FALSE)
  significant <- pvalue < alpha
  
  # Confidence interval based on likelihood ratio
  cutoff <- best_scan_ll - 0.5 * qchisq(1 - alpha, df = 1)
  ci_set <- which(profile_ll >= cutoff)
  ci_lower <- if (length(ci_set) > 0) min(ci_set) else NA
  ci_upper <- if (length(ci_set) > 0) max(ci_set) else NA
  
  # Create result object
  result <- list(
    cp_est = if (significant) best_cp else NA,
    cp_raw = best_cp,
    significant = significant,
    n_cps = if (significant) 1 else 0,
    lr_stat = lr_stat,
    pvalue = pvalue,
    df = df,
    ci = c(ci_lower, ci_upper),
    ci_set = ci_set,
    profile_ll = profile_ll,
    scan_method = scan_method,
    null_elbo = null_fit$elbo,
    alt_elbo = alt_ll,
    null_fit = null_fit,
    left_fit = best_fits$left,
    right_fit = best_fits$right,
    type = type,
    distribution = distribution,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    Q = Q,
    Q2 = Q2,
    block_sigma = block_sigma,
    min_segment = min_segment,
    alpha = alpha
  )
  
  class(result) <- c("cp_detection", "list")
  
  if (verbose) {
    if (significant) {
      message(sprintf("\nChange point detected at t = %d (p = %.4f)", best_cp, pvalue))
      message(sprintf("95%% CI: [%d, %d]", ci_lower, ci_upper))
    } else {
      message(sprintf("\nNo significant change point (p = %.4f)", pvalue))
    }
  }
  
  return(result)
}

#' Detect Multiple Change Points
#'
#' Detect multiple change points using Wild Binary Segmentation (WBS),
#' standard Binary Segmentation (BS), or PELT (Pruned Exact Linear Time).
#'
#' @param Y Network tensor (T x N1 x N2)
#' @param Q Number of clusters (or Q1 for bipartite)
#' @param Q2 Number of column clusters (NULL for unipartite)
#' @param type Network type: "unipartite" or "bipartite"
#' @param distribution Edge distribution: "gaussian" or "bernoulli"
#' @param method Search method: "wbs" (default), "bs", or "pelt"
#' @param scan_method Candidate-split scan method: "profile" (full segment refits)
#'   or "score" (fast approximation using a fixed global fit)
#' @param block_sigma Use block-specific variance
#' @param use_cpp Use C++ acceleration
#' @param min_segment Minimum segment length
#' @param max_cps Maximum change points to detect
#' @param alpha Significance level
#' @param penalty Penalty for PELT: "BIC", "AIC", "MBIC", or numeric value
#' @param max_iter Maximum VEM iterations
#' @param tol Convergence tolerance
#' @param n_init Number of random initializations
#' @param n_intervals Number of random intervals for WBS
#' @param verbose Print progress
#' @return Object of class "multi_cp_detection" with results
#' @export
#' @examples
#' \dontrun{
#' # Simulate with 2 change points
#' theta1 <- list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2))
#' theta2 <- list(prob = matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2))
#' theta3 <- list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2))
#' 
#' data <- simulate_sbm_multi_cp(N = 30, T_len = 90, cp_times = c(30, 60),
#'                               Q = 2, type = "unipartite",
#'                               distribution = "bernoulli",
#'                               theta_list = list(theta1, theta2, theta3))
#' 
#' # Default: WBS
#' result <- detect_multiple_cp(data$Y, Q = 2, type = "unipartite",
#'                              distribution = "bernoulli", method = "wbs")
#' 
#' # PELT for optimal segmentation
#' result <- detect_multiple_cp(data$Y, Q = 2, type = "unipartite",
#'                              distribution = "bernoulli", method = "pelt")
#' print(result$cps)
#' }
detect_multiple_cp <- function(Y, Q, Q2 = NULL,
                               type = c("unipartite", "bipartite"),
                               distribution = c("gaussian", "bernoulli"),
                               method = c("wbs", "bs", "pelt"),
                               scan_method = c("profile", "score"),
                               block_sigma = FALSE,
                               use_cpp = TRUE,
                               min_segment = 5,
                               max_cps = 10,
                               alpha = 0.05,
                               penalty = "BIC",
                               max_iter = 100,
                               tol = 1e-6,
                               n_init = 5,
                               n_intervals = 100,
                               verbose = FALSE) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  method <- match.arg(method)
  scan_method <- match.arg(scan_method)
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  if (type == "unipartite") {
    stopifnot(N1 == N2)
    Q2 <- Q
  } else {
    if (is.null(Q2)) Q2 <- Q
  }
  
  if (verbose) {
    message(sprintf("Detecting multiple CPs using %s (scan=%s)...",
                    toupper(method), scan_method))
  }
  
  global_fit <- NULL
  score_prep <- NULL
  
  # PELT uses different algorithm path

  if (method == "pelt") {
    if (scan_method == "score") {
      warning("scan_method='score' is not implemented for method='pelt'; using profile scan path.")
      scan_method <- "profile"
    }
    pelt_result <- .detect_pelt(Y, Q, Q2, type, distribution, penalty, block_sigma,
                                use_cpp, min_segment, max_cps, max_iter, tol, n_init, verbose)
    pelt_result$scan_method <- scan_method
    pelt_result$min_segment <- min_segment
    return(pelt_result)
  }
  
  if (scan_method == "score") {
    global_fit <- tryCatch(
      fit_sbm(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter, tol, max(1, n_init), FALSE),
      error = function(e) NULL
    )
    score_prep <- if (!is.null(global_fit)) {
      tryCatch(.prepare_score_scan(Y, global_fit, type, distribution, block_sigma),
               error = function(e) NULL)
    } else {
      NULL
    }
    if (is.null(score_prep)) {
      warning("Score scan initialization failed; falling back to profile scan refits.")
      scan_method <- "profile"
    }
  }
  
  # Find candidate change points (WBS or BS)
  if (scan_method == "score") {
    if (method == "wbs") {
      candidates <- .wbs_candidates_score(score_prep, min_segment, n_intervals, verbose)
    } else {
      candidates <- .bs_candidates_score(score_prep, min_segment, verbose)
    }
  } else {
    if (method == "wbs") {
      candidates <- .wbs_candidates(Y, Q, Q2, type, distribution, block_sigma,
                                    use_cpp, min_segment, max_iter, tol, n_init,
                                    n_intervals, verbose)
    } else {
      candidates <- .bs_candidates(Y, Q, Q2, type, distribution, block_sigma,
                                   use_cpp, min_segment, max_iter, tol, n_init,
                                   verbose)
    }
  }
  
  if (length(candidates) == 0) {
    if (verbose) message("No change points detected")
    
    null_fit <- fit_sbm(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                        max_iter, tol, n_init, FALSE)
    
    result <- list(
      cps = integer(0),
      n_cps = 0,
      segments = list(list(start = 1, end = T_len, fit = null_fit)),
      candidates = integer(0),
      method = method,
      scan_method = scan_method,
      min_segment = min_segment,
      type = type,
      distribution = distribution,
      T_len = T_len,
      N1 = N1,
      N2 = N2,
      Q = Q,
      Q2 = Q2,
      global_fit = global_fit
    )
    class(result) <- c("multi_cp_detection", "list")
    return(result)
  }
  
  # Filter candidates to be well-separated
  final_cps <- .filter_cps(candidates, min_segment, max_cps)
  
  if (verbose) message(sprintf("Detected %d change point(s): %s", 
                               length(final_cps), paste(final_cps, collapse = ", ")))
  
  # Fit segments
  segments <- .fit_final_segments(Y, final_cps, Q, Q2, type, distribution,
                                   block_sigma, use_cpp, max_iter, tol, n_init)
  
  result <- list(
    cps = sort(final_cps),
    n_cps = length(final_cps),
    segments = segments,
    candidates = candidates,
    method = method,
    scan_method = scan_method,
    min_segment = min_segment,
    type = type,
    distribution = distribution,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    Q = Q,
    Q2 = Q2,
    block_sigma = block_sigma,
    global_fit = global_fit
  )
  
  class(result) <- c("multi_cp_detection", "list")
  
  return(result)
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Compute degrees of freedom for LRT
#' @keywords internal
.compute_df <- function(Q, Q2, type, distribution, block_sigma) {
  if (type == "unipartite") {
    n_mu <- Q * (Q + 1) / 2
  } else {
    n_mu <- Q * Q2
  }
  
  if (distribution == "gaussian") {
    if (block_sigma) {
      n_sigma <- n_mu
    } else {
      n_sigma <- 1
    }
    return(n_mu + n_sigma)
  } else {
    return(n_mu)
  }
}

#' Create empty CP result
#' @keywords internal
.create_cp_result <- function(cp, significant, n_cps, type, distribution,
                              T_len, N1, N2, Q, Q2) {
  result <- list(
    cp_est = cp,
    significant = significant,
    n_cps = n_cps,
    type = type,
    distribution = distribution,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    Q = Q,
    Q2 = Q2
  )
  class(result) <- c("cp_detection", "list")
  return(result)
}

#' Extract soft assignments for score scan
#' @keywords internal
.extract_score_tau <- function(fit, type) {
  if (!is.null(fit$tau1) && !is.null(fit$tau2)) {
    return(list(tau1 = fit$tau1, tau2 = fit$tau2))
  }
  if (!is.null(fit$tau_row) && !is.null(fit$tau_col)) {
    return(list(tau1 = fit$tau_row, tau2 = fit$tau_col))
  }
  if (!is.null(fit$tau) && type == "unipartite") {
    return(list(tau1 = fit$tau, tau2 = fit$tau))
  }
  stop("Unable to extract variational memberships from fitted model for score scan.")
}

#' Prepare cumulative sufficient statistics for score scan
#' @keywords internal
.prepare_score_scan <- function(Y, fit, type, distribution, block_sigma) {
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  tau_obj <- .extract_score_tau(fit, type)
  tau1 <- as.matrix(tau_obj$tau1)
  tau2 <- as.matrix(tau_obj$tau2)
  
  if (nrow(tau1) != N1 || nrow(tau2) != N2) {
    stop("Score-scan memberships have incompatible dimensions.")
  }
  
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  B <- Q1 * Q2
  
  if (type == "unipartite") {
    ones <- matrix(1, nrow = N1, ncol = N2)
    cross <- t(tau1) %*% ones %*% tau2
    self_term <- t(tau1) %*% tau2
    n_block_mat <- 0.5 * (cross - self_term)
  } else {
    ones <- matrix(1, nrow = N1, ncol = N2)
    n_block_mat <- t(tau1) %*% ones %*% tau2
  }
  n_block <- as.numeric(n_block_mat)
  n_block[!is.finite(n_block)] <- 0
  n_block <- pmax(n_block, 1e-12)
  
  sum1 <- matrix(0, nrow = T_len, ncol = B)
  sum2 <- if (distribution == "gaussian") matrix(0, nrow = T_len, ncol = B) else NULL
  
  for (tt in seq_len(T_len)) {
    Y_t <- Y[tt, , , drop = TRUE]
    if (!is.matrix(Y_t)) {
      Y_t <- matrix(Y_t, nrow = N1, ncol = N2)
    }
    Y_t[!is.finite(Y_t)] <- 0
    
    if (type == "unipartite") {
      diag(Y_t) <- 0
      S1 <- 0.5 * (t(tau1) %*% Y_t %*% tau2)
      if (distribution == "gaussian") {
        S2 <- 0.5 * (t(tau1) %*% (Y_t * Y_t) %*% tau2)
      }
    } else {
      S1 <- t(tau1) %*% Y_t %*% tau2
      if (distribution == "gaussian") {
        S2 <- t(tau1) %*% (Y_t * Y_t) %*% tau2
      }
    }
    
    sum1[tt, ] <- as.numeric(S1)
    if (distribution == "gaussian") sum2[tt, ] <- as.numeric(S2)
  }
  
  cum_sum1 <- rbind(rep(0, B), apply(sum1, 2, cumsum))
  if (B == 1L) cum_sum1 <- matrix(cum_sum1, ncol = 1)
  
  if (distribution == "gaussian") {
    cum_sum2 <- rbind(rep(0, B), apply(sum2, 2, cumsum))
    if (B == 1L) cum_sum2 <- matrix(cum_sum2, ncol = 1)
  } else {
    cum_sum2 <- matrix(0, nrow = T_len + 1L, ncol = B)
  }
  
  list(
    T_len = T_len,
    Q1 = Q1,
    Q2 = Q2,
    type = type,
    distribution = distribution,
    block_sigma = block_sigma,
    tau1 = tau1,
    tau2 = tau2,
    n_block = n_block,
    cum_sum1 = cum_sum1,
    cum_sum2 = cum_sum2
  )
}

#' Segment sufficient statistics from score prep
#' @keywords internal
.score_segment_stats <- function(prep, start, end) {
  if (!is.finite(start) || !is.finite(end) || start > end) {
    B <- length(prep$n_block)
    return(list(sum1 = rep(0, B), sum2 = rep(0, B), len = 0L))
  }
  sum1 <- prep$cum_sum1[end + 1L, ] - prep$cum_sum1[start, ]
  sum2 <- prep$cum_sum2[end + 1L, ] - prep$cum_sum2[start, ]
  list(sum1 = as.numeric(sum1), sum2 = as.numeric(sum2), len = as.integer(end - start + 1L))
}

#' Segment log-likelihood proxy under score scan
#' @keywords internal
.score_segment_ll <- function(sum1, sum2, len, n_block, distribution, block_sigma) {
  if (!is.finite(len) || len <= 0) return(-Inf)
  
  N <- pmax(len * as.numeric(n_block), 1e-12)
  S1 <- as.numeric(sum1)
  
  if (distribution == "bernoulli") {
    S1 <- pmin(pmax(S1, 0), N)
    p_hat <- pmin(pmax(S1 / N, 1e-8), 1 - 1e-8)
    ll <- sum(S1 * log(p_hat) + (N - S1) * log(1 - p_hat))
    return(as.numeric(ll))
  }
  
  S2 <- as.numeric(sum2)
  mu_hat <- S1 / N
  sse <- S2 - (S1 * S1) / N
  sse <- pmax(sse, 1e-8)
  
  if (block_sigma) {
    sigma2 <- pmax(sse / N, 1e-8)
    ll <- sum(-0.5 * N * (log(2 * pi * sigma2) + 1))
  } else {
    N_total <- sum(N)
    sse_total <- sum(sse)
    sigma2 <- pmax(sse_total / pmax(N_total, 1e-12), 1e-8)
    ll <- -0.5 * N_total * (log(2 * pi * sigma2) + 1)
  }
  as.numeric(ll)
}

#' Score profile over an interval
#' @keywords internal
.score_profile_interval <- function(prep, start, end, min_segment) {
  profile_ll <- rep(-Inf, prep$T_len)
  
  scan_start <- as.integer(start + min_segment - 1L)
  scan_end <- as.integer(end - min_segment)
  if (!is.finite(scan_start) || !is.finite(scan_end) || scan_end < scan_start) {
    return(profile_ll)
  }
  
  for (cp in scan_start:scan_end) {
    left <- .score_segment_stats(prep, start, cp)
    right <- .score_segment_stats(prep, cp + 1L, end)
    ll_left <- .score_segment_ll(left$sum1, left$sum2, left$len,
                                 prep$n_block, prep$distribution, prep$block_sigma)
    ll_right <- .score_segment_ll(right$sum1, right$sum2, right$len,
                                  prep$n_block, prep$distribution, prep$block_sigma)
    profile_ll[cp] <- ll_left + ll_right
  }
  
  profile_ll
}

#' Best score split in an interval
#' @keywords internal
.score_best_split_interval <- function(prep, start, end, min_segment) {
  profile_ll <- .score_profile_interval(prep, start, end, min_segment)
  idx <- which(is.finite(profile_ll))
  if (length(idx) == 0) {
    return(list(cp = NA_integer_, gain = 0, best_ll = -Inf, null_ll = -Inf))
  }
  
  seg <- .score_segment_stats(prep, start, end)
  null_ll <- .score_segment_ll(seg$sum1, seg$sum2, seg$len,
                               prep$n_block, prep$distribution, prep$block_sigma)
  best_cp <- as.integer(idx[which.max(profile_ll[idx])])
  best_ll <- as.numeric(profile_ll[best_cp])
  gain <- best_ll - null_ll
  if (!is.finite(gain)) gain <- 0
  
  list(cp = best_cp, gain = gain, best_ll = best_ll, null_ll = null_ll)
}

#' WBS candidates using score scan
#' @keywords internal
.wbs_candidates_score <- function(prep, min_segment, n_intervals, verbose) {
  T_len <- prep$T_len
  min_interval <- 2L * as.integer(min_segment) + 1L
  if (T_len < min_interval) return(integer(0))
  
  safe_sample <- function(x, n = 1L) {
    if (length(x) == 1L) x else sample(x, n)
  }
  
  candidates <- list()
  for (m in seq_len(n_intervals)) {
    s <- safe_sample(seq_len(T_len - min_interval + 1L), 1L)
    e <- safe_sample(seq.int(s + min_interval - 1L, T_len), 1L)
    if (e < s || e > T_len || s < 1L) next
    
    best <- .score_best_split_interval(prep, s, e, min_segment)
    if (is.finite(best$cp) && best$gain > 0) {
      candidates[[length(candidates) + 1L]] <- list(cp = as.integer(best$cp), gain = as.numeric(best$gain))
    }
    
    if (verbose && m %% 20 == 0) {
      message(sprintf("  WBS interval %d/%d (score scan)", m, n_intervals))
    }
  }
  
  if (length(candidates) == 0) return(integer(0))
  
  cps <- vapply(candidates, `[[`, integer(1), "cp")
  gains <- vapply(candidates, `[[`, numeric(1), "gain")
  unique_cps <- unique(cps)
  max_gains <- vapply(unique_cps, function(cp) max(gains[cps == cp]), numeric(1))
  ord <- order(max_gains, decreasing = TRUE)
  unique_cps[ord]
}

#' BS candidates using score scan
#' @keywords internal
.bs_candidates_score <- function(prep, min_segment, verbose) {
  T_len <- prep$T_len
  candidates <- integer(0)
  
  .bs_recursive_score <- function(start, end) {
    len <- end - start + 1L
    if (len < 2L * min_segment + 1L) return(invisible(NULL))
    
    best <- .score_best_split_interval(prep, start, end, min_segment)
    if (is.finite(best$cp) && best$gain > 0) {
      cp <- as.integer(best$cp)
      candidates <<- c(candidates, cp)
      .bs_recursive_score(start, cp)
      .bs_recursive_score(cp + 1L, end)
    }
    invisible(NULL)
  }
  
  .bs_recursive_score(1L, T_len)
  if (verbose) message(sprintf("  BS score candidates found: %d", length(candidates)))
  candidates
}

#' Wild Binary Segmentation to find candidates
#' @keywords internal
.wbs_candidates <- function(Y, Q, Q2, type, distribution, block_sigma,
                            use_cpp, min_segment, max_iter, tol, n_init,
                            n_intervals, verbose) {
  
  T_len <- dim(Y)[1]
  min_interval <- 2 * min_segment + 1
  
  if (T_len < min_interval) return(integer(0))
  
  # Safe sample function (avoids R's sample(n,1) returning sample(1:n,1))
  safe_sample <- function(x, n = 1) {
    if (length(x) == 1) x else sample(x, n)
  }
  
  candidates <- list()
  
  for (m in 1:n_intervals) {
    # Random interval
    s <- safe_sample(1:(T_len - min_interval + 1), 1)
    e <- safe_sample((s + min_interval - 1):T_len, 1)
    
    # Validate interval (defensive check)
    if (e < s || e > T_len || s < 1) next
    
    Y_sub <- Y[s:e, , , drop = FALSE]
    T_sub <- e - s + 1
    scan_range <- min_segment:(T_sub - min_segment)
    
    if (length(scan_range) < 1) next
    
    # Null fit for this interval
    null_fit <- tryCatch(
      fit_sbm(Y_sub, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
      error = function(e) list(elbo = -Inf)
    )
    
    best_cp <- NA
    best_gain <- 0
    
    for (t in scan_range) {
      Y_l <- Y_sub[1:t, , , drop = FALSE]
      Y_r <- Y_sub[(t+1):T_sub, , , drop = FALSE]
      
      fit_l <- tryCatch(
        fit_sbm(Y_l, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      fit_r <- tryCatch(
        fit_sbm(Y_r, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      gain <- fit_l$elbo + fit_r$elbo - null_fit$elbo
      
      # Guard against NA/NaN gains (can happen with -Inf elbos)
      if (!is.na(gain) && is.finite(gain) && gain > best_gain) {
        best_gain <- gain
        best_cp <- t
      }
    }
    
    if (!is.na(best_cp) && best_gain > 0) {
      cp_global <- s + best_cp - 1
      candidates[[length(candidates) + 1]] <- list(cp = cp_global, gain = best_gain)
    }
    
    if (verbose && m %% 20 == 0) {
      message(sprintf("  WBS interval %d/%d", m, n_intervals))
    }
  }
  
  if (length(candidates) == 0) return(integer(0))
  
  # Aggregate: keep unique positions with max gain
  cps <- sapply(candidates, `[[`, "cp")
  gains <- sapply(candidates, `[[`, "gain")
  
  unique_cps <- unique(cps)
  max_gains <- sapply(unique_cps, function(cp) max(gains[cps == cp]))
  
  # Sort by gain
  ord <- order(max_gains, decreasing = TRUE)
  return(unique_cps[ord])
}

#' Binary Segmentation to find candidates
#' @keywords internal
.bs_candidates <- function(Y, Q, Q2, type, distribution, block_sigma,
                           use_cpp, min_segment, max_iter, tol, n_init, verbose) {
  
  T_len <- dim(Y)[1]
  candidates <- integer(0)
  
  # Recursive helper
  .bs_recursive <- function(start, end) {
    len <- end - start + 1
    if (len < 2 * min_segment + 1) return()
    
    Y_sub <- Y[start:end, , , drop = FALSE]
    T_sub <- len
    scan_range <- min_segment:(T_sub - min_segment)
    
    if (length(scan_range) < 1) return()
    
    null_fit <- tryCatch(
      fit_sbm(Y_sub, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
      error = function(e) list(elbo = -Inf)
    )
    
    best_cp <- NA
    best_gain <- 0
    
    for (t in scan_range) {
      Y_l <- Y_sub[1:t, , , drop = FALSE]
      Y_r <- Y_sub[(t+1):T_sub, , , drop = FALSE]
      
      fit_l <- tryCatch(
        fit_sbm(Y_l, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      fit_r <- tryCatch(
        fit_sbm(Y_r, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      gain <- fit_l$elbo + fit_r$elbo - null_fit$elbo
      
      # Guard against NA/NaN gains (can happen with -Inf elbos)
      if (!is.na(gain) && is.finite(gain) && gain > best_gain) {
        best_gain <- gain
        best_cp <- t
      }
    }
    
    if (!is.na(best_cp) && best_gain > 0) {
      cp_global <- start + best_cp - 1
      candidates <<- c(candidates, cp_global)
      
      # Recurse on both sides
      .bs_recursive(start, cp_global)
      .bs_recursive(cp_global + 1, end)
    }
  }
  
  .bs_recursive(1, T_len)
  return(candidates)
}

#' Filter candidates to be well-separated
#' @keywords internal
.filter_cps <- function(candidates, min_segment, max_cps) {
  if (length(candidates) == 0) return(integer(0))

  # Preserve candidate priority order (e.g., WBS gain ranking) and remove duplicates.
  candidates <- as.integer(candidates[is.finite(candidates)])
  if (length(candidates) == 0) return(integer(0))
  candidates <- candidates[!duplicated(candidates)]

  selected <- integer(0)
  
  for (cp in candidates) {
    if (length(selected) == 0 || min(abs(selected - cp)) >= min_segment) {
      selected <- c(selected, cp)
    }
    if (length(selected) >= max_cps) break
  }
  
  return(sort(selected))
}

#' Fit segments after detecting CPs
#' @keywords internal
.fit_final_segments <- function(Y, cps, Q, Q2, type, distribution,
                                 block_sigma, use_cpp, max_iter, tol, n_init) {
  
  T_len <- dim(Y)[1]
  boundaries <- c(0, sort(cps), T_len)
  n_segments <- length(boundaries) - 1
  
  segments <- vector("list", n_segments)
  
  for (s in 1:n_segments) {
    start_t <- boundaries[s] + 1
    end_t <- boundaries[s + 1]
    
    Y_seg <- Y[start_t:end_t, , , drop = FALSE]
    
    fit <- tryCatch(
      fit_sbm(Y_seg, Q, Q2, type, distribution, block_sigma, use_cpp,
              max_iter, tol, n_init, FALSE),
      error = function(e) NULL
    )
    
    segments[[s]] <- list(start = start_t, end = end_t, fit = fit)
  }
  
  return(segments)
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' Print method for cp_detection
#' @param x A `cp_detection` object returned by `detect_single_cp()`.
#' @param ... Additional arguments (currently unused).
#' @export
print.cp_detection <- function(x, ...) {
  cat("\n=== Single Change Point Detection ===\n\n")
  cat(sprintf("Network: %s %s\n", x$type, x$distribution))
  cat(sprintf("Dimensions: T=%d, N1=%d, N2=%d\n", x$T_len, x$N1, x$N2))
  cat(sprintf("Clusters: Q=%d", x$Q))
  if (x$type == "bipartite") cat(sprintf(", Q2=%d", x$Q2))
  cat("\n\n")
  
  if (x$significant) {
    cat(sprintf("Change point detected: t = %d\n", x$cp_est))
    cat(sprintf("LR statistic: %.4f (df=%d)\n", x$lr_stat, x$df))
    cat(sprintf("P-value: %.4e\n", x$pvalue))
    if (!is.na(x$ci[1])) {
      cat(sprintf("%d%% CI: [%d, %d]\n", round((1 - x$alpha) * 100), x$ci[1], x$ci[2]))
    }
  } else {
    cat("No significant change point detected\n")
    cat(sprintf("P-value: %.4f (threshold: %.4f)\n", x$pvalue, x$alpha))
  }
  
  invisible(x)
}

# =============================================================================
# PELT (Pruned Exact Linear Time) Implementation
# =============================================================================

#' PELT algorithm for optimal change point detection
#' @keywords internal
.detect_pelt <- function(Y, Q, Q2, type, distribution, penalty, block_sigma,
                         use_cpp, min_segment, max_cps, max_iter, tol, 
                         n_init, verbose) {
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  # Compute penalty value (β per change point)
  beta <- .compute_pelt_penalty(penalty, Q, Q2, type, distribution, 
                                 block_sigma, N1, N2, T_len)
  
  if (verbose) {
    message(sprintf("PELT: T=%d, penalty=%.2f", T_len, beta))
    message("Precomputing segment costs...")
  }
  
  # Precompute segment costs
  cost_matrix <- .precompute_costs_pelt(Y, Q, Q2, type, distribution, 
                                         block_sigma, use_cpp, min_segment,
                                         max_iter, tol, n_init, verbose)
  
  if (verbose) message("Running PELT dynamic programming...")
  
  # Run PELT dynamic programming
  pelt_result <- .pelt_dp(cost_matrix, beta, min_segment, T_len, verbose)
  
  # Apply max_cps constraint if specified
  cps <- pelt_result$change_points
  if (!is.null(max_cps) && length(cps) > max_cps) {
    if (verbose) message(sprintf("Limiting to %d CPs (found %d)", max_cps, length(cps)))
    cps <- cps[1:max_cps]
  }
  
  if (verbose) {
    if (length(cps) > 0) {
      message(sprintf("Detected %d change point(s): %s", 
                      length(cps), paste(cps, collapse = ", ")))
    } else {
      message("No change points detected")
    }
  }
  
  # Fit segments with detected CPs
  segments <- .fit_final_segments(Y, cps, Q, Q2, type, distribution,
                                   block_sigma, use_cpp, max_iter, tol, n_init)
  
  result <- list(
    cps = sort(cps),
    n_cps = length(cps),
    segments = segments,
    candidates = cps,  # For PELT, candidates = final CPs
    method = "pelt",
    penalty = beta,
    total_cost = pelt_result$total_cost,
    type = type,
    distribution = distribution,
    T_len = T_len,
    N1 = N1,
    N2 = N2,
    Q = Q,
    Q2 = Q2,
    block_sigma = block_sigma
  )
  
  class(result) <- c("multi_cp_detection", "list")
  
  return(result)
}


#' Compute PELT penalty value
#' @keywords internal
.compute_pelt_penalty <- function(penalty, Q, Q2, type, distribution,
                                   block_sigma, N1, N2, T_len) {
  
  # Number of parameters per segment
  if (type == "unipartite") {
    n_mu <- Q * (Q + 1) / 2  # Symmetric matrix
    n_alpha <- Q - 1         # Mixture proportions
  } else {
    n_mu <- Q * Q2           # Full matrix
    n_alpha <- (Q - 1) + (Q2 - 1)
  }
  
  if (distribution == "gaussian") {
    if (block_sigma) {
      n_sigma <- if (type == "unipartite") Q * (Q + 1) / 2 else Q * Q2
    } else {
      n_sigma <- 1
    }
  } else {
    n_sigma <- 0
  }
  
  n_params <- n_mu + n_sigma + n_alpha
  
  if (is.numeric(penalty)) {
    return(penalty)
  }
  
  n_obs <- N1 * N2 * T_len
  
  if (penalty == "BIC") {
    beta <- 0.5 * n_params * log(n_obs)
  } else if (penalty == "AIC") {
    beta <- n_params
  } else if (penalty == "MBIC") {
    beta <- 0.5 * n_params * (log(n_obs) + log(T_len))
  } else {
    beta <- 0.5 * n_params * log(n_obs)  # Default to BIC
  }
  
  return(beta)
}


#' Precompute segment costs for PELT
#' @keywords internal
.precompute_costs_pelt <- function(Y, Q, Q2, type, distribution, block_sigma,
                                    use_cpp, min_segment, max_iter, tol, 
                                    n_init, verbose) {
  
  T_len <- dim(Y)[1]
  cost_matrix <- matrix(Inf, T_len, T_len)
  
  n_computed <- 0
  n_total <- sum(sapply(1:(T_len - min_segment + 1), function(s) T_len - s - min_segment + 2))
  
  for (s in 1:(T_len - min_segment + 1)) {
    for (e in (s + min_segment - 1):T_len) {
      
      Y_seg <- Y[s:e, , , drop = FALSE]
      
      fit <- tryCatch(
        fit_sbm(Y_seg, Q, Q2, type, distribution, block_sigma, use_cpp,
                max_iter %/% 2, tol * 10, max(1, n_init %/% 2), FALSE),
        error = function(e) list(elbo = -Inf)
      )
      
      cost_matrix[s, e] <- -fit$elbo
      n_computed <- n_computed + 1
      
      if (verbose && n_computed %% 50 == 0) {
        message(sprintf("  Segment costs: %d/%d (%.1f%%)", 
                        n_computed, n_total, 100 * n_computed / n_total))
      }
    }
  }
  
  return(cost_matrix)
}


#' PELT dynamic programming with pruning
#' @keywords internal
.pelt_dp <- function(cost_matrix, beta, min_segment, T_len, verbose) {
  
  F <- rep(Inf, T_len + 1)
  F[1] <- -beta
  cp <- rep(0, T_len + 1)
  R <- c(0)
  
  for (t in min_segment:T_len) {
    candidates <- R[R <= (t - min_segment)]
    
    if (length(candidates) == 0) {
      if (t >= min_segment && is.finite(cost_matrix[1, t])) {
        F[t + 1] <- cost_matrix[1, t] + beta
        cp[t + 1] <- 0
      }
      R <- unique(c(R, t))
      next
    }
    
    costs <- sapply(candidates, function(tau) {
      seg_cost <- cost_matrix[tau + 1, t]
      if (is.infinite(seg_cost)) return(Inf)
      F[tau + 1] + seg_cost + beta
    })
    
    if (all(is.infinite(costs))) {
      R <- unique(c(R, t))
      next
    }
    
    best_idx <- which.min(costs)
    F[t + 1] <- costs[best_idx]
    cp[t + 1] <- candidates[best_idx]
    
    # Pruning
    new_R <- c()
    for (tau in R) {
      if (tau + 1 <= t && t <= T_len) {
        seg_cost <- cost_matrix[tau + 1, t]
        if (is.finite(seg_cost) && F[tau + 1] + seg_cost <= F[t + 1]) {
          new_R <- c(new_R, tau)
        } else if (tau > t - min_segment) {
          new_R <- c(new_R, tau)
        }
      } else if (tau > t - min_segment) {
        new_R <- c(new_R, tau)
      }
    }
    R <- unique(c(new_R, t))
    
    if (verbose && t %% 20 == 0) {
      message(sprintf("  t=%d/%d, |R|=%d", t, T_len, length(R)))
    }
  }
  
  # Backtrack
  cps <- integer(0)
  current <- T_len
  while (cp[current + 1] > 0) {
    cps <- c(cp[current + 1], cps)
    current <- cp[current + 1]
  }
  
  return(list(change_points = cps, total_cost = F[T_len + 1]))
}

#' Print method for multi_cp_detection
#' @param x A `multi_cp_detection` object returned by `detect_multiple_cp()`.
#' @param ... Additional arguments (currently unused).
#' @export
print.multi_cp_detection <- function(x, ...) {
  cat("\n=== Multiple Change Point Detection ===\n\n")
  cat(sprintf("Method: %s\n", toupper(x$method)))
  cat(sprintf("Network: %s %s\n", x$type, x$distribution))
  cat(sprintf("Dimensions: T=%d, N1=%d, N2=%d\n", x$T_len, x$N1, x$N2))
  cat(sprintf("Clusters: Q=%d", x$Q))
  if (x$type == "bipartite") cat(sprintf(", Q2=%d", x$Q2))
  cat("\n\n")
  
  if (x$n_cps == 0) {
    cat("No change points detected\n")
  } else {
    cat(sprintf("Detected %d change point(s):\n", x$n_cps))
    for (i in seq_along(x$cps)) {
      cat(sprintf("  CP %d: t = %d\n", i, x$cps[i]))
    }
  }
  
  cat(sprintf("\nSegments: %d\n", length(x$segments)))
  for (s in seq_along(x$segments)) {
    seg <- x$segments[[s]]
    cat(sprintf("  Segment %d: [%d, %d] (length %d)\n",
                s, seg$start, seg$end, seg$end - seg$start + 1))
  }
  
  invisible(x)
}
