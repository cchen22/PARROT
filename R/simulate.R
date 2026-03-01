#' Simulation Functions for PARROT
#' 
#' Supports 4 network types:
#' - Unipartite Weighted (Gaussian)
#' - Unipartite Binary (Bernoulli)
#' - Bipartite Weighted (Gaussian)
#' - Bipartite Binary (Bernoulli)

#' Simulate a Network from SBM
#'
#' Generate synthetic network data from a Stochastic Block Model
#' with static community structure.
#'
#' @param N Number of nodes (for unipartite) or N1 for bipartite
#' @param N2 Number of column nodes (NULL for unipartite)
#' @param T_len Number of time points
#' @param Q Number of clusters (for unipartite) or Q1 for bipartite
#' @param Q2 Number of column clusters (NULL for unipartite)
#' @param type Network type: "unipartite" or "bipartite"
#' @param distribution Edge distribution: "gaussian" or "bernoulli"
#' @param theta Model parameters (list with mu/prob and optionally sigma)
#' @param labels Optional: fixed node labels
#' @return List containing Y (network tensor) and labels
#' @export
#' @examples
#' \dontrun{
#' # Simulate unipartite binary network
#' data <- simulate_sbm(N = 30, T_len = 50, Q = 2,
#'                      type = "unipartite", distribution = "bernoulli",
#'                      theta = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)))
#'
#' # Simulate bipartite weighted network
#' data <- simulate_sbm(N = 20, N2 = 30, T_len = 50, Q = 2, Q2 = 2,
#'                      type = "bipartite", distribution = "gaussian",
#'                      theta = list(mu = matrix(c(1, 0, 0, 1), 2, 2), sigma = 0.5))
#' }
simulate_sbm <- function(N, N2 = NULL, T_len, Q, Q2 = NULL, 
                         type = c("unipartite", "bipartite"),
                         distribution = c("gaussian", "bernoulli"),
                         theta, labels = NULL) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  
  if (type == "unipartite") {
    return(.simulate_unipartite(N, T_len, Q, distribution, theta, labels))
  } else {
    if (is.null(N2)) N2 <- N
    if (is.null(Q2)) Q2 <- Q
    return(.simulate_bipartite(N, N2, T_len, Q, Q2, distribution, theta, labels))
  }
}

#' Simulate a Network with Single Change Point
#'
#' Generate synthetic network data with a single structural change point.
#'
#' @param N Number of nodes (or N1 for bipartite)
#' @param N2 Number of column nodes (NULL for unipartite)
#' @param T_len Total time points
#' @param cp_time Change point location (time before which theta1 applies)
#' @param Q Number of clusters (or Q1 for bipartite)
#' @param Q2 Number of column clusters (NULL for unipartite)
#' @param type Network type: "unipartite" or "bipartite"
#' @param distribution Edge distribution: "gaussian" or "bernoulli"
#' @param theta1 Parameters before CP (list with mu/prob)
#' @param theta2 Parameters after CP
#' @param Q_changes Logical: does Q change at CP? (for future use)
#' @param labels1 Labels before CP (NULL = generate balanced)
#' @param labels2 Labels after CP (NULL = same as labels1)
#' @return List with Y (network tensor), cp_true (true CP location), and labels
#' @export
#' @examples
#' \dontrun{
#' # Simulate network with change point at t=30
#' data <- simulate_sbm_cp(N = 30, T_len = 60, cp_time = 30, Q = 2,
#'                         type = "unipartite", distribution = "bernoulli",
#'                         theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
#'                         theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2)))
#' print(data$cp_true)  # Should be 30
#' }
simulate_sbm_cp <- function(N, N2 = NULL, T_len, cp_time,
                            Q, Q2 = NULL,
                            type = c("unipartite", "bipartite"),
                            distribution = c("gaussian", "bernoulli"),
                            theta1, theta2,
                            Q_changes = FALSE,
                            labels1 = NULL, labels2 = NULL) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  
  stopifnot(cp_time > 0, cp_time < T_len)
  
  if (type == "unipartite") {
    N2 <- N
    Q2 <- Q
  } else {
    if (is.null(N2)) N2 <- N
    if (is.null(Q2)) Q2 <- Q
  }
  
  # Generate labels if not provided
  if (is.null(labels1)) {
    if (type == "unipartite") {
      labels1 <- list(z = .balanced_labels(N, Q))
    } else {
      labels1 <- list(z1 = .balanced_labels(N, Q), 
                      z2 = .balanced_labels(N2, Q2))
    }
  }
  if (is.null(labels2)) labels2 <- labels1
  
  # Simulate before CP
  Y1 <- .simulate_segment(N, N2, cp_time, Q, Q2, type, distribution, theta1, labels1)
  
  # Simulate after CP
  T_after <- T_len - cp_time
  Y2 <- .simulate_segment(N, N2, T_after, 
                          if(Q_changes) nrow(theta2$mu %||% theta2$prob) else Q,
                          if(Q_changes) ncol(theta2$mu %||% theta2$prob) else Q2,
                          type, distribution, theta2, labels2)
  
  # Combine
  if (type == "unipartite") {
    Y <- array(0, dim = c(T_len, N, N))
    Y[1:cp_time, , ] <- Y1
    Y[(cp_time + 1):T_len, , ] <- Y2
  } else {
    Y <- array(0, dim = c(T_len, N, N2))
    Y[1:cp_time, , ] <- Y1
    Y[(cp_time + 1):T_len, , ] <- Y2
  }
  
  return(list(
    Y = Y,
    cp_true = cp_time,
    type = type,
    distribution = distribution,
    labels1 = labels1,
    labels2 = labels2,
    theta1 = theta1,
    theta2 = theta2,
    Q_changes = Q_changes
  ))
}

#' Simulate a Network with Multiple Change Points
#'
#' Generate synthetic network data with multiple structural change points.
#'
#' @param N Number of nodes
#' @param N2 Number of column nodes (NULL for unipartite)
#' @param T_len Total time points
#' @param cp_times Vector of change point locations (sorted)
#' @param Q Number of clusters (or vector for varying Q)
#' @param Q2 Number of column clusters
#' @param type Network type
#' @param distribution Edge distribution
#' @param theta_list List of parameter sets for each segment (length = n_cp + 1)
#' @return List with Y (network tensor), cp_true (true CP locations), and parameters
#' @export
#' @examples
#' \dontrun{
#' # Simulate network with 2 change points
#' theta1 <- list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2))
#' theta2 <- list(prob = matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2))
#' theta3 <- list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2))
#' 
#' data <- simulate_sbm_multi_cp(N = 30, T_len = 90, cp_times = c(30, 60),
#'                               Q = 2, type = "unipartite", 
#'                               distribution = "bernoulli",
#'                               theta_list = list(theta1, theta2, theta3))
#' }
simulate_sbm_multi_cp <- function(N, N2 = NULL, T_len, cp_times,
                                  Q, Q2 = NULL,
                                  type = c("unipartite", "bipartite"),
                                  distribution = c("gaussian", "bernoulli"),
                                  theta_list) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  
  # Validate cp_times
  cp_times <- sort(unique(cp_times))
  stopifnot(all(cp_times > 0), all(cp_times < T_len))
  
  if (type == "unipartite") {
    N2 <- N
    Q2 <- Q
  } else {
    if (is.null(N2)) N2 <- N
    if (is.null(Q2)) Q2 <- Q
  }
  
  n_segments <- length(cp_times) + 1
  if (length(theta_list) != n_segments) {
    stop(sprintf("theta_list must have %d elements (n_cp + 1), got %d", 
                 n_segments, length(theta_list)))
  }
  
  # Generate labels
  labels <- if (type == "unipartite") {
    list(z = .balanced_labels(N, Q))
  } else {
    list(z1 = .balanced_labels(N, Q), z2 = .balanced_labels(N2, Q2))
  }
  
  # Segment boundaries
  boundaries <- c(0, cp_times, T_len)
  
  if (type == "unipartite") {
    Y <- array(0, dim = c(T_len, N, N))
  } else {
    Y <- array(0, dim = c(T_len, N, N2))
  }
  
  for (s in 1:n_segments) {
    t_start <- boundaries[s] + 1
    t_end <- boundaries[s + 1]
    seg_len <- t_end - t_start + 1
    
    Y_seg <- .simulate_segment(N, N2, seg_len, Q, Q2, type, distribution, 
                               theta_list[[s]], labels)
    Y[t_start:t_end, , ] <- Y_seg
  }
  
  return(list(
    Y = Y,
    cp_true = cp_times,
    n_cps = length(cp_times),
    type = type,
    distribution = distribution,
    labels = labels,
    theta_list = theta_list,
    segment_boundaries = boundaries
  ))
}

# =============================================================================
# INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Generate balanced cluster labels
#' @keywords internal
.balanced_labels <- function(N, Q) {
  sort(rep(1:Q, length.out = N))
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Simulate unipartite network segment
#' @keywords internal
.simulate_unipartite <- function(N, T_len, Q, distribution, theta, labels) {
  if (is.null(labels)) {
    z <- .balanced_labels(N, Q)
  } else {
    z <- labels$z %||% labels
  }
  
  Y <- array(0, dim = c(T_len, N, N))
  
  if (distribution == "gaussian") {
    mu <- theta$mu
    # Support both global sigma and block-specific sigma_block
    sigma <- theta$sigma %||% 1
    sigma_block <- theta$sigma_block  # Can be NULL or a Q x Q matrix
    
    # Vectorized generation
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        # Use block-specific sigma if provided, otherwise global sigma
        if (!is.null(sigma_block)) {
          sigma_ij <- sigma_block[z[i], z[j]]
        } else {
          sigma_ij <- sigma
        }
        vals <- rnorm(T_len, mu[z[i], z[j]], sigma_ij)
        Y[, i, j] <- vals
        Y[, j, i] <- vals  # Symmetric
      }
    }
  } else {  # Bernoulli
    prob <- theta$prob
    
    for (i in 1:(N-1)) {
      for (j in (i+1):N) {
        vals <- rbinom(T_len, 1, prob[z[i], z[j]])
        Y[, i, j] <- vals
        Y[, j, i] <- vals
      }
    }
  }
  
  return(list(Y = Y, labels = list(z = z)))
}

#' Simulate bipartite network segment
#' @keywords internal
.simulate_bipartite <- function(N1, N2, T_len, Q1, Q2, distribution, theta, labels) {
  if (is.null(labels)) {
    z1 <- .balanced_labels(N1, Q1)
    z2 <- .balanced_labels(N2, Q2)
  } else {
    z1 <- labels$z1 %||% .balanced_labels(N1, Q1)
    z2 <- labels$z2 %||% .balanced_labels(N2, Q2)
  }
  
  Y <- array(0, dim = c(T_len, N1, N2))
  
  if (distribution == "gaussian") {
    mu <- theta$mu
    # Support both global sigma and block-specific sigma_block
    sigma <- theta$sigma %||% 1
    sigma_block <- theta$sigma_block  # Can be NULL or a Q1 x Q2 matrix
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        # Use block-specific sigma if provided, otherwise global sigma
        if (!is.null(sigma_block)) {
          sigma_ij <- sigma_block[z1[i], z2[j]]
        } else {
          sigma_ij <- sigma
        }
        Y[, i, j] <- rnorm(T_len, mu[z1[i], z2[j]], sigma_ij)
      }
    }
  } else {  # Bernoulli
    prob <- theta$prob
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        Y[, i, j] <- rbinom(T_len, 1, prob[z1[i], z2[j]])
      }
    }
  }
  
  return(list(Y = Y, labels = list(z1 = z1, z2 = z2)))
}

#' Simulate a segment (dispatch to unipartite or bipartite)
#' @keywords internal
.simulate_segment <- function(N1, N2, T_len, Q1, Q2, type, distribution, theta, labels) {
  if (type == "unipartite") {
    result <- .simulate_unipartite(N1, T_len, Q1, distribution, theta, labels)
    return(result$Y)
  } else {
    result <- .simulate_bipartite(N1, N2, T_len, Q1, Q2, distribution, theta, labels)
    return(result$Y)
  }
}
