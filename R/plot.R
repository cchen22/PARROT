#' Visualization Functions for PARROT
#'
#' Plot change point detection results, profile likelihoods, and network heatmaps.

#' Plot Profile Likelihood
#'
#' Plot the profile log-likelihood function for change point detection.
#'
#' @param result A cp_detection or parrot_result object
#' @param show_ci Show confidence interval region
#' @param show_cp Show detected change point
#' @param true_cp True change point (optional, for comparison)
#' @param main Plot title
#' @param ... Additional arguments passed to plot()
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' result <- detect_single_cp(data$Y, Q = 2, type = "unipartite",
#'                           distribution = "bernoulli")
#' plot_profile(result, true_cp = 30)
#' }
plot_profile <- function(result,
                         show_ci = TRUE,
                         show_cp = TRUE,
                         true_cp = NULL,
                         main = "Profile Log-Likelihood",
                         ...) {
  
  if (!inherits(result, c("cp_detection", "parrot_result"))) {
    stop("result must be a cp_detection or parrot_result object")
  }
  
  profile_ll <- result$profile_ll
  
  if (is.null(profile_ll)) {
    warning("No profile likelihood available to plot")
    return(invisible(NULL))
  }
  
  T_len <- length(profile_ll)
  t_vals <- 1:T_len
  
  # Remove -Inf for plotting
  valid <- is.finite(profile_ll)
  
  if (sum(valid) == 0) {
    warning("No valid profile likelihood values to plot")
    return(invisible(NULL))
  }
  
  y_range <- range(profile_ll[valid])
  y_range[1] <- y_range[1] - 0.1 * diff(y_range)
  y_range[2] <- y_range[2] + 0.1 * diff(y_range)
  
  plot(t_vals[valid], profile_ll[valid], type = "l", lwd = 2,
       xlab = "Time", ylab = "Profile Log-Likelihood",
       main = main, ylim = y_range, ...)
  
  # Confidence interval region
  if (show_ci && !is.null(result$ci_set) && length(result$ci_set) > 0) {
    ci_lower <- min(result$ci_set)
    ci_upper <- max(result$ci_set)
    
    rect(ci_lower, y_range[1], ci_upper, y_range[2],
         col = rgb(0, 0, 1, 0.1), border = NA)
  }
  
  # Detected change point
  if (show_cp) {
    cp_est <- result$cp_est
    if (!is.null(cp_est) && !is.na(cp_est)) {
      abline(v = cp_est, col = "blue", lwd = 2, lty = 1)
      points(cp_est, profile_ll[cp_est], pch = 19, cex = 1.5, col = "blue")
    }
  }
  
  # True change point
  if (!is.null(true_cp)) {
    abline(v = true_cp, col = "red", lwd = 2, lty = 2)
  }
  
  # Legend
  legend_items <- c("Profile LL")
  legend_cols <- c("black")
  legend_lty <- c(1)
  legend_pch <- c(NA)
  
  if (show_cp && !is.null(result$cp_est) && !is.na(result$cp_est)) {
    legend_items <- c(legend_items, sprintf("Detected CP (t=%d)", result$cp_est))
    legend_cols <- c(legend_cols, "blue")
    legend_lty <- c(legend_lty, 1)
    legend_pch <- c(legend_pch, 19)
  }
  
  if (!is.null(true_cp)) {
    legend_items <- c(legend_items, sprintf("True CP (t=%d)", true_cp))
    legend_cols <- c(legend_cols, "red")
    legend_lty <- c(legend_lty, 2)
    legend_pch <- c(legend_pch, NA)
  }
  
  legend("topright", legend = legend_items, col = legend_cols,
         lty = legend_lty, pch = legend_pch, lwd = 2, bty = "n")
  
  invisible(NULL)
}

#' Plot Detection Results
#'
#' Create a comprehensive visualization of change point detection results.
#'
#' @param result A parrot_result, cp_detection, or multi_cp_detection object
#' @param true_cps True change points (optional)
#' @param show_segments Show segment boundaries
#' @param main Plot title
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' result <- parrot(data$Y, Q = 2)
#' plot_detection(result, true_cps = c(30, 60))
#' }
plot_detection <- function(result,
                           true_cps = NULL,
                           show_segments = TRUE,
                           main = NULL,
                           ...) {
  
  if (is.null(main)) {
    main <- sprintf("Change Point Detection (%s %s)",
                    result$type, result$distribution)
  }
  
  # For single CP with profile likelihood
  if (inherits(result, "cp_detection") || 
      (inherits(result, "parrot_result") && !result$multiple)) {
    
    if (!is.null(result$profile_ll)) {
      plot_profile(result, show_ci = TRUE, show_cp = TRUE,
                   true_cp = if (length(true_cps) == 1) true_cps[1] else NULL,
                   main = main, ...)
    } else {
      .plot_timeline(result, true_cps, main, ...)
    }
    
  } else {
    # Multiple change points - timeline view
    .plot_timeline(result, true_cps, main, ...)
  }
  
  invisible(NULL)
}

#' Plot Network Heatmap
#'
#' Visualize the network adjacency matrix as a heatmap.
#'
#' @param Y Network matrix (N1 x N2) or tensor (T x N1 x N2)
#' @param time Time index to plot (if Y is tensor)
#' @param labels Optional cluster labels for ordering
#' @param main Plot title
#' @param col Color palette
#' @param ... Additional arguments passed to image()
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' data <- simulate_sbm(N = 30, T_len = 10, Q = 2,
#'                      type = "unipartite", distribution = "bernoulli",
#'                      theta = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)))
#' plot_network_heatmap(data$Y, time = 1, labels = data$labels$z)
#' }
plot_network_heatmap <- function(Y, time = 1, labels = NULL, 
                                  main = "Network Adjacency Matrix",
                                  col = NULL, ...) {
  
  # Handle tensor vs matrix
  if (length(dim(Y)) == 3) {
    Y_mat <- Y[time, , ]
  } else if (length(dim(Y)) == 2) {
    Y_mat <- Y
  } else {
    stop("Y must be a 2D or 3D array")
  }
  
  N1 <- nrow(Y_mat)
  N2 <- ncol(Y_mat)
  
  # Order by labels if provided
  if (!is.null(labels)) {
    if (is.list(labels)) {
      row_order <- order(labels$z1 %||% labels$z %||% 1:N1)
      col_order <- order(labels$z2 %||% labels$z %||% 1:N2)
    } else {
      row_order <- order(labels)
      col_order <- order(labels)
    }
    Y_mat <- Y_mat[row_order, col_order]
  }
  
  # Color palette
  if (is.null(col)) {
    # Detect if binary or continuous
    unique_vals <- unique(as.vector(Y_mat))
    if (all(unique_vals %in% c(0, 1, NA))) {
      col <- c("white", "darkblue")
    } else {
      col <- colorRampPalette(c("blue", "white", "red"))(100)
    }
  }
  
  # Plot
  image(1:N2, 1:N1, t(Y_mat[N1:1, ]),
        col = col, xlab = "", ylab = "",
        main = main, axes = FALSE, ...)
  axis(1, at = seq(1, N2, length.out = min(5, N2)))
  axis(2, at = seq(1, N1, length.out = min(5, N1)))
  box()
  
  invisible(NULL)
}

#' Plot Network Evolution
#'
#' Visualize how the network changes over time around change points.
#'
#' @param Y Network tensor (T x N1 x N2)
#' @param times Vector of time points to plot
#' @param cps Detected change points (optional)
#' @param labels Cluster labels
#' @param nrow Number of rows in plot grid
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_network_evolution <- function(Y, times = NULL, cps = NULL, 
                                    labels = NULL, nrow = 1, ...) {
  
  T_len <- dim(Y)[1]
  
  if (is.null(times)) {
    if (!is.null(cps) && length(cps) > 0) {
      # Plot around change points
      times <- c(1, cps, cps + 1, T_len)
      times <- sort(unique(times))
      times <- times[times >= 1 & times <= T_len]
    } else {
      times <- seq(1, T_len, length.out = min(4, T_len))
    }
  }
  
  n_plots <- length(times)
  ncol <- ceiling(n_plots / nrow)
  
  old_par <- par(mfrow = c(nrow, ncol), mar = c(2, 2, 3, 1))
  on.exit(par(old_par))
  
  for (t in times) {
    title <- sprintf("t = %d", t)
    if (!is.null(cps) && t %in% cps) {
      title <- paste(title, "(CP)")
    }
    plot_network_heatmap(Y, time = t, labels = labels, main = title, ...)
  }
  
  invisible(NULL)
}

#' Plot Segment Parameters
#'
#' Visualize the estimated parameters for each segment.
#'
#' @param result A parrot_result or multi_cp_detection object
#' @param param Which parameter to plot: "mu" or "prob"
#' @param ... Additional arguments
#' @return Invisible NULL
#' @export
plot_segment_params <- function(result, param = c("mu", "prob"), ...) {
  
  param <- match.arg(param)
  
  if (!inherits(result, c("parrot_result", "multi_cp_detection"))) {
    stop("result must be a parrot_result or multi_cp_detection object")
  }
  
  n_seg <- length(result$segments)
  
  if (n_seg == 0) {
    warning("No segments to plot")
    return(invisible(NULL))
  }
  
  # Set up plot grid
  ncol <- min(n_seg, 4)
  nrow <- ceiling(n_seg / ncol)
  old_par <- par(mfrow = c(nrow, ncol), mar = c(3, 3, 3, 1))
  on.exit(par(old_par))
  
  for (s in 1:n_seg) {
    seg <- result$segments[[s]]
    
    if (is.null(seg$fit)) {
      plot.new()
      title(main = sprintf("Segment %d: No fit", s))
      next
    }
    
    mat <- if (param == "mu") seg$fit$mu else seg$fit$prob
    
    if (is.null(mat)) {
      plot.new()
      title(main = sprintf("Segment %d: No %s", s, param))
      next
    }
    
    # Heatmap
    image(mat, main = sprintf("Segment %d: [%d, %d]", s, seg$start, seg$end),
          col = colorRampPalette(c("blue", "white", "red"))(100),
          axes = FALSE, ...)
    box()
  }
  
  invisible(NULL)
}

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Null-coalescing operator
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Plot timeline for multiple CPs
#' @keywords internal
.plot_timeline <- function(result, true_cps = NULL, main = NULL, ...) {
  
  T_len <- result$T_len
  cps <- if (!is.null(result$cps)) result$cps else integer(0)
  n_cps <- length(cps)
  
  if (is.null(main)) {
    main <- sprintf("Change Points Timeline (n=%d)", n_cps)
  }
  
  # Create a simple timeline plot
  plot(c(0, T_len), c(0, 1), type = "n",
       xlab = "Time", ylab = "", yaxt = "n", main = main, ...)
  
  # Draw timeline
  segments(0, 0.5, T_len, 0.5, lwd = 2)
  
  # Draw segment boundaries
  boundaries <- c(0, cps, T_len)
  segments <- length(boundaries) - 1
  
  colors <- rainbow(segments, alpha = 0.5)
  
  for (s in 1:segments) {
    rect(boundaries[s], 0.3, boundaries[s + 1], 0.7,
         col = colors[s], border = NA)
  }
  
  # Draw detected CPs
  if (n_cps > 0) {
    for (cp in cps) {
      abline(v = cp, col = "blue", lwd = 2)
      points(cp, 0.5, pch = 18, cex = 2, col = "blue")
    }
  }
  
  # Draw true CPs
  if (!is.null(true_cps)) {
    for (tcp in true_cps) {
      abline(v = tcp, col = "red", lwd = 2, lty = 2)
      points(tcp, 0.5, pch = 17, cex = 1.5, col = "red")
    }
  }
  
  # Legend
  legend_items <- character(0)
  legend_cols <- character(0)
  legend_pch <- integer(0)
  
  if (n_cps > 0) {
    legend_items <- c(legend_items, sprintf("Detected (%d)", n_cps))
    legend_cols <- c(legend_cols, "blue")
    legend_pch <- c(legend_pch, 18)
  }
  
  if (!is.null(true_cps)) {
    legend_items <- c(legend_items, sprintf("True (%d)", length(true_cps)))
    legend_cols <- c(legend_cols, "red")
    legend_pch <- c(legend_pch, 17)
  }
  
  if (length(legend_items) > 0) {
    legend("topright", legend = legend_items, col = legend_cols,
           pch = legend_pch, bty = "n")
  }
  
  invisible(NULL)
}
