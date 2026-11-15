# =============================================================================
# Resampling tests and helpers: time-label permutation test (re-scanning,
# conditional on the window between neighbouring change points), stationary
# bootstrap, seeded segment fits, test-based binary segmentation and the ICL
# gain-ratio rule for Q
# =============================================================================

#' Seed of a segment fit
#'
#' Derived from the (sorted) time indices of the segment, so that every segment
#' fit is reproducible and independent of the order in which fits are computed
#' (and of the number of cores used).
#' @keywords internal
.subset_seed <- function(idx, seed = 1L) {
  idx <- sort(as.numeric(idx))
  h <- (as.numeric(seed) * 1000003 + sum(idx * 7919 * seq_along(idx)) + length(idx) * 104729) %% 2147483647
  as.integer(h)
}

#' Fit the SBM to a set of snapshots with a reproducible seed
#'
#' The global random-number state is restored afterwards. \code{.fit_segment_full}
#' returns the fitted model, \code{.fit_segment_seeded} its ELBO; both give the
#' same fit for the same time points and seed.
#' @keywords internal
.fit_segment_full <- function(Y, idx, Q, Q2, type, distribution, block_sigma,
                              use_cpp, max_iter, tol, n_init, seed = 1L) {
  idx <- sort(as.integer(idx))
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  } else assign(".Random.seed", old, envir = globalenv()))
  set.seed(.subset_seed(idx, seed))
  tryCatch(
    fit_sbm(Y[idx, , , drop = FALSE], Q, Q2, type, distribution, block_sigma,
            use_cpp, max_iter, tol, n_init, FALSE),
    error = function(e) list(elbo = -Inf))
}

#' @rdname dot-fit_segment_full
#' @keywords internal
.fit_segment_seeded <- function(Y, idx, Q, Q2, type, distribution, block_sigma,
                                use_cpp, max_iter, tol, n_init, seed = 1L) {
  as.numeric(.fit_segment_full(Y, idx, Q, Q2, type, distribution, block_sigma,
                               use_cpp, max_iter, tol, n_init, seed)$elbo)
}

#' Memoised segment ELBO
#'
#' The SBM likelihood of a segment does not depend on the order of its
#' snapshots, so each (multi)set of time indices is fitted once and cached.
#' Observed and resampled statistics therefore use identical fit settings.
#' @keywords internal
.segment_elbo_cache <- function(Y, Q, Q2, type, distribution, block_sigma,
                                use_cpp, max_iter, tol, n_init, seed = 1L) {
  cache <- new.env(hash = TRUE)
  spec <- list(Y = Y, Q = Q, Q2 = Q2, type = type, distribution = distribution,
               block_sigma = block_sigma, use_cpp = use_cpp, max_iter = max_iter,
               tol = tol, n_init = n_init, seed = seed)
  E <- function(idx) {
    idx <- sort(as.integer(idx))
    key <- paste(idx, collapse = ",")
    val <- cache[[key]]
    if (!is.null(val)) return(val)
    val <- .fit_segment_seeded(Y, idx, Q, Q2, type, distribution, block_sigma,
                               use_cpp, max_iter, tol, n_init, seed)
    assign(key, val, envir = cache)
    val
  }
  attr(E, "cache") <- cache
  attr(E, "spec") <- spec
  E
}

#' Fit the not yet cached segments in parallel
#' @keywords internal
.cache_fill <- function(E, idx_list, n_cores = 1L) {
  cache <- attr(E, "cache"); spec <- attr(E, "spec")
  if (is.null(cache) || is.null(spec)) {
    for (idx in idx_list) E(idx)
    return(invisible(NULL))
  }
  keys <- vapply(idx_list, function(i) paste(sort(as.integer(i)), collapse = ","), character(1))
  keep <- !duplicated(keys) & !vapply(keys, exists, logical(1), envir = cache, inherits = FALSE)
  todo <- idx_list[keep]; keys <- keys[keep]
  if (length(todo) == 0) return(invisible(NULL))
  if (n_cores > 1L && length(todo) > 1L) {
    cl <- parallel::makeCluster(min(as.integer(n_cores), length(todo)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(libs) { .libPaths(libs); suppressPackageStartupMessages(library(PARROT)); NULL }, .libPaths())
    parallel::clusterCall(cl, function(s) { assign(".parrot_spec", s, envir = globalenv()); NULL }, spec)
    worker <- function(idx) {
      s <- get(".parrot_spec", envir = globalenv())
      fit <- utils::getFromNamespace(".fit_segment_seeded", "PARROT")   # worker runs outside the namespace
      fit(s$Y, idx, s$Q, s$Q2, s$type, s$distribution, s$block_sigma, s$use_cpp, s$max_iter, s$tol, s$n_init, s$seed)
    }
    environment(worker) <- globalenv()
    vals <- parallel::parLapplyLB(cl, todo, worker)
    for (k in seq_along(todo)) assign(keys[k], vals[[k]], envir = cache)
  } else {
    for (idx in todo) E(idx)
  }
  invisible(NULL)
}

#' Likelihood-ratio statistic for a fixed split from a memoised ELBO function
#' @keywords internal
.lr_from_cache <- function(E, left, right) {
  max(0, 2 * (E(left) + E(right) - E(c(left, right))))
}

# segments needed for the statistic of one ordered sequence of snapshots
.window_segments <- function(sq, a, min_segment, rescan) {
  w <- length(sq)
  ts <- if (rescan) seq.int(min_segment, w - min_segment) else a
  ts <- ts[ts >= 1 & ts < w]
  c(list(sq), lapply(ts, function(t) sq[seq_len(t)]), lapply(ts, function(t) sq[(t + 1L):w]))
}

# statistic of one ordered sequence: LR at the split after position a (fixed
# split) or the maximum LR over all admissible splits (re-scan)
.window_stat <- function(E, sq, a, min_segment, rescan) {
  w <- length(sq)
  ts <- if (rescan) seq.int(min_segment, w - min_segment) else a
  ts <- ts[ts >= 1 & ts < w]
  if (length(ts) == 0) return(0)
  max(vapply(ts, function(t) .lr_from_cache(E, sq[seq_len(t)], sq[(t + 1L):w]), numeric(1)))
}

# all ordered sequences of length w drawn without replacement from pool
.ordered_sequences <- function(pool, w) {
  perms <- function(v) {
    if (length(v) <= 1) return(list(v))
    out <- vector("list", 0)
    for (i in seq_along(v)) for (p in perms(v[-i])) out[[length(out) + 1L]] <- c(v[i], p)
    out
  }
  combs <- if (length(pool) == w) list(pool) else utils::combn(pool, w, simplify = FALSE)
  unlist(lapply(combs, perms), recursive = FALSE)
}

#' Time-label permutation p-value for one change point
#'
#' The test is conditional on the window \code{left_start:right_end} bounded by
#' the neighbouring change points (or the ends of the series); \code{cp} is the
#' detected split. Under the null hypothesis the snapshots are exchangeable:
#' with \code{scope = "global"} the time labels of the whole series are permuted
#' and the snapshots falling on the window positions are used; with
#' \code{scope = "local"} only the snapshots inside the window are permuted.
#'
#' With \code{rescan = TRUE} (default) the statistic is the maximum
#' likelihood-ratio statistic over all admissible splits of the window
#' (segments of at least \code{min_segment} snapshots), for the observed and for
#' every permuted series, i.e. the scan that selected the change point is
#' repeated on the permuted data; the test then controls the type-I error at the
#' nominal level. With \code{rescan = FALSE} the split is held fixed at
#' \code{cp}; this ignores that the split was selected by the scan and is
#' anti-conservative (kept for comparison with PARROT < 1.1.0).
#'
#' Observed and permuted statistics are computed with identical fit settings.
#' If the number of distinct configurations is at most \code{exact_max}, all of
#' them are enumerated (exact p-value), otherwise \code{n_perm} random
#' permutations are used (Monte Carlo p-value, (1 + #exceedances)/(n_perm + 1)).
#' @keywords internal
.permutation_pvalue <- function(Y, cp, left_start, right_end, Q, Q2, type,
                                distribution, block_sigma, n_init = 20,
                                max_iter = 200, tol = 1e-6, use_cpp = TRUE,
                                n_perm = 99, scope = c("global", "local"),
                                exact_max = 5040, elbo_fun = NULL,
                                verbose = FALSE, min_segment = 2L, rescan = TRUE,
                                n_cores = 1L, fit_seed = 1L) {
  scope <- match.arg(scope)
  T_len <- dim(Y)[1]
  E <- if (is.null(elbo_fun)) {
    .segment_elbo_cache(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                        max_iter, tol, n_init, seed = fit_seed)
  } else elbo_fun
  obs_seq <- seq.int(left_start, right_end)
  a <- cp - left_start + 1L
  w <- length(obs_seq)
  pool <- if (scope == "global") seq_len(T_len) else obs_seq
  n <- length(pool)
  n_cfg <- round(if (rescan) exp(lfactorial(n) - lfactorial(n - w)) else choose(n, a) * choose(n - a, w - a))
  exact <- n_cfg <= exact_max
  seqs <- if (exact) {
    if (rescan) .ordered_sequences(pool, w) else {
      unlist(lapply(utils::combn(pool, a, simplify = FALSE), function(l) {
        rest <- setdiff(pool, l)
        rs <- if (length(rest) == w - a) list(rest) else utils::combn(rest, w - a, simplify = FALSE)
        lapply(rs, function(r) c(l, r))
      }), recursive = FALSE)
    }
  } else {
    lapply(seq_len(n_perm), function(k) pool[sample.int(n)][seq_len(w)])
  }
  if (n_cores > 1L) {
    need <- unlist(lapply(c(list(obs_seq), seqs), .window_segments, a = a,
                          min_segment = min_segment, rescan = rescan), recursive = FALSE)
    .cache_fill(E, need, n_cores)
  }
  obs <- .window_stat(E, obs_seq, a, min_segment, rescan)
  null <- vapply(seqs, function(sq) .window_stat(E, sq, a, min_segment, rescan), numeric(1))
  eps <- 1e-9 * max(1, abs(obs))
  pval <- if (exact) mean(null >= obs - eps) else (1 + sum(null >= obs - eps)) / (n_perm + 1)
  how <- if (exact) "exact" else "monte_carlo"
  if (verbose) message(sprintf("    CP %d: statistic = %.2f, %s %s permutation p = %.4g (%s, %g configurations)",
                               cp, obs, scope, if (rescan) "re-scanning" else "fixed-split",
                               pval, how, n_cfg))
  list(pvalue = pval, lr_stat = obs, method = how, scope = scope, rescan = rescan,
       perm_space_size = n_cfg, n_perm_used = length(null), null = null)
}

#' Stationary-bootstrap time indices (Politis and Romano 1994)
#'
#' Geometric block lengths with mean \code{mean_block}, circular wrap-around.
#' @keywords internal
.sb_index <- function(T_len, mean_block) {
  idx <- integer(T_len); t <- 1L
  while (t <= T_len) {
    start <- sample.int(T_len, 1L); len <- stats::rgeom(1L, 1 / mean_block) + 1L
    for (k in seq_len(len)) {
      if (t > T_len) break
      idx[t] <- ((start + k - 2L) %% T_len) + 1L; t <- t + 1L
    }
  }
  idx
}

#' Stationary-bootstrap p-value for one change point
#'
#' Resamples blocks of consecutive snapshots, which preserves short-range
#' temporal dependence under the null hypothesis. Same statistic as
#' \code{.permutation_pvalue} (re-scan within the window by default).
#' @keywords internal
.sb_pvalue <- function(Y, cp, left_start, right_end, Q, Q2, type, distribution,
                       block_sigma, n_init = 20, max_iter = 200, tol = 1e-6,
                       use_cpp = TRUE, n_boot = 99, mean_block = NULL,
                       elbo_fun = NULL, verbose = FALSE, min_segment = 2L,
                       rescan = TRUE, n_cores = 1L, fit_seed = 1L) {
  T_len <- dim(Y)[1]
  if (is.null(mean_block)) mean_block <- ceiling(T_len^(1 / 3))
  E <- if (is.null(elbo_fun)) {
    .segment_elbo_cache(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                        max_iter, tol, n_init, seed = fit_seed)
  } else elbo_fun
  obs_seq <- seq.int(left_start, right_end)
  a <- cp - left_start + 1L
  seqs <- lapply(seq_len(n_boot), function(k) .sb_index(T_len, mean_block)[obs_seq])
  if (n_cores > 1L) {
    need <- unlist(lapply(c(list(obs_seq), seqs), .window_segments, a = a,
                          min_segment = min_segment, rescan = rescan), recursive = FALSE)
    .cache_fill(E, need, n_cores)
  }
  obs <- .window_stat(E, obs_seq, a, min_segment, rescan)
  null <- vapply(seqs, function(sq) .window_stat(E, sq, a, min_segment, rescan), numeric(1))
  pval <- (1 + sum(null >= obs - 1e-9 * max(1, abs(obs)))) / (n_boot + 1)
  list(pvalue = pval, lr_stat = obs, method = "stationary_bootstrap", rescan = rescan,
       scope = NA_character_, perm_space_size = NA_real_, n_perm_used = n_boot,
       mean_block = mean_block, null = null)
}

#' Test-based binary segmentation
#'
#' Recursively splits a segment at its best profile split as long as the
#' permutation (or stationary-bootstrap) p-value of that segment is below
#' \code{alpha}; at most \code{max_cps} change points.
#' @keywords internal
.test_based_bs <- function(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                           min_segment, max_cps, alpha, max_iter, tol, n_init,
                           pvalue_method, n_perm, permutation_scope, exact_max,
                           sb_block, verbose, rescan = TRUE, n_cores = 1L,
                           fit_seed = 1L) {
  T_len <- dim(Y)[1]
  E <- .segment_elbo_cache(Y, Q, Q2, type, distribution, block_sigma, use_cpp,
                           max_iter, tol, n_init, seed = fit_seed)
  queue <- list(c(1L, T_len)); cps <- integer(0); tests <- list()
  while (length(queue) > 0 && length(cps) < max_cps) {
    seg <- queue[[1]]; queue <- queue[-1]
    s <- seg[1]; e <- seg[2]
    if (e - s + 1L < 2L * min_segment) next
    splits <- seq.int(s + min_segment - 1L, e - min_segment)
    lr <- vapply(splits, function(t) .lr_from_cache(E, s:t, (t + 1L):e), numeric(1))
    t_best <- splits[which.max(lr)]
    te <- if (pvalue_method == "stationary_bootstrap") {
      .sb_pvalue(Y, t_best, s, e, Q, Q2, type, distribution, block_sigma, n_init,
                 max_iter, tol, use_cpp, n_perm, sb_block, E, verbose,
                 min_segment = min_segment, rescan = rescan, n_cores = n_cores)
    } else {
      .permutation_pvalue(Y, t_best, s, e, Q, Q2, type, distribution, block_sigma,
                          n_init, max_iter, tol, use_cpp, n_perm, permutation_scope,
                          exact_max, E, verbose, min_segment = min_segment,
                          rescan = rescan, n_cores = n_cores)
    }
    te$null <- NULL; te$cp <- t_best; te$segment <- c(s, e)
    tests[[length(tests) + 1L]] <- te
    if (te$pvalue < alpha) {
      cps <- c(cps, t_best)
      queue <- c(queue, list(c(s, t_best)), list(c(t_best + 1L, e)))
    }
  }
  list(cps = sort(cps), tests = tests)
}

#' ICL gain-ratio rule for the number of communities
#'
#' With g(Q) = ICL(Q) - ICL(Q - 1), the rule selects the Q at which the marginal
#' ICL gain drops most, i.e. the Q maximising g(Q) / g(Q + 1), among the Q
#' values up to the ICL maximum.
#' @param Q Candidate numbers of communities (consecutive values)
#' @param icl ICL values for \code{Q}
#' @return Selected Q
#' @keywords internal
.q_gain_ratio_rule <- function(Q, icl) {
  o <- order(Q); Q <- Q[o]; icl <- icl[o]
  ok <- is.finite(icl); Q <- Q[ok]; icl <- icl[ok]
  if (length(Q) == 0) return(NA_integer_)
  qa <- which.max(icl)
  if (qa <= 2L) return(as.integer(Q[qa]))
  g <- c(NA, diff(icl))
  cand <- 2:(qa - 1L)
  ratio <- g[cand] / pmax(g[cand + 1L], .Machine$double.eps)
  as.integer(Q[cand][which.max(ratio)])
}

#' Adjusted Rand index (internal)
#' @keywords internal
.ari <- function(a, b) {
  tab <- table(a, b); n <- sum(tab)
  sc <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  ex <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == ex) return(1)
  (sc - ex) / (mx - ex)
}

#' Select Q with the ICL gain-ratio rule
#'
#' Fits the SBM to the whole series for each candidate Q (Q = 1 is always
#' included so that the first gain is defined) and applies
#' \code{.q_gain_ratio_rule}. Optionally checks the stability of the fitted
#' memberships across \code{n_seeds} independent sets of initialisations.
#' @keywords internal
.select_Q_icl <- function(Y, type, distribution, Q_range, block_sigma, use_cpp,
                          n_init, max_iter, tol, stability = FALSE, n_seeds = 5,
                          verbose = FALSE) {
  Q_range <- sort(unique(c(1L, as.integer(Q_range))))
  if (verbose) message("Selecting Q with the ICL gain-ratio rule...")
  fits <- lapply(Q_range, function(q) tryCatch(
    fit_sbm(Y, q, q, type, distribution, block_sigma, use_cpp, max_iter, tol, n_init, FALSE),
    error = function(e) NULL))
  icl <- vapply(fits, function(f) if (is.null(f)) -Inf else f$icl, numeric(1))
  elbo <- vapply(fits, function(f) if (is.null(f)) -Inf else f$elbo, numeric(1))
  Q_sel <- .q_gain_ratio_rule(Q_range, icl)
  gain <- c(NA, diff(icl))
  tab <- data.frame(Q = Q_range, elbo = elbo, icl = icl, gain = gain,
                    gain_ratio = gain / c(gain[-1], NA))
  if (verbose) for (i in seq_along(Q_range))
    message(sprintf("  Q=%d: ICL=%.2f%s", Q_range[i], icl[i],
                    if (identical(Q_range[i], Q_sel)) "  <- selected" else ""))
  stab <- NULL
  if (isTRUE(stability)) {
    stab <- do.call(rbind, lapply(Q_range[Q_range >= 2], function(q) {
      reps <- lapply(seq_len(n_seeds), function(s) tryCatch(
        fit_sbm(Y, q, q, type, distribution, block_sigma, use_cpp, max_iter, tol, n_init, FALSE),
        error = function(e) NULL))
      reps <- reps[!vapply(reps, is.null, logical(1))]
      z <- lapply(reps, function(f) if (!is.null(f$labels$z1)) f$labels$z1 else f$labels$z)
      ari <- if (length(z) >= 2) {
        pr <- utils::combn(length(z), 2)
        mean(apply(pr, 2, function(p) .ari(z[[p[1]]], z[[p[2]]])))
      } else NA_real_
      el <- vapply(reps, function(f) f$elbo, numeric(1))
      data.frame(Q = q, mean_pairwise_ari = ari, elbo_range = diff(range(el)))
    }))
    if (verbose && !is.na(Q_sel) && Q_sel >= 2)
      message(sprintf("  stability of the selected Q=%d: mean pairwise ARI %.3f",
                      Q_sel, stab$mean_pairwise_ari[stab$Q == Q_sel]))
  }
  list(Q = Q_sel, table = tab, stability = stab)
}
