#' Model Fitting Functions for PARROT
#'
#' Unified interface for fitting Stochastic Block Models (SBM) to all 4 network types:
#' - Unipartite Weighted (Gaussian)
#' - Unipartite Binary (Bernoulli)
#' - Bipartite Weighted (Gaussian)
#' - Bipartite Binary (Bernoulli)
#'
#' Supports both global variance and block-specific variance for weighted networks.
#' Uses C++ acceleration by default for 14-46x speedup.

#' Fit SBM to Network Data
#'
#' @param Y Network tensor (T x N x N for unipartite, T x N1 x N2 for bipartite)
#' @param Q Number of clusters (or Q1 for bipartite)
#' @param Q2 Number of column clusters (NULL for unipartite)
#' @param type Network type: "unipartite" or "bipartite"
#' @param distribution Edge distribution: "gaussian" or "bernoulli"
#' @param block_sigma Use block-specific variance instead of global (for Gaussian only)
#' @param use_cpp Use C++ implementation (TRUE) or pure R (FALSE)
#' @param max_iter Maximum EM iterations
#' @param tol Convergence tolerance
#' @param n_init Number of random initializations
#' @param verbose Print progress
#' @return List with model parameters and fit statistics
#' @export
#' @examples
#' \dontrun{
#' # Simulate unipartite binary network
#' data <- simulate_sbm(N = 30, T_len = 50, Q = 2,
#'                      type = "unipartite", distribution = "bernoulli",
#'                      theta = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)))
#'
#' # Fit SBM model
#' fit <- fit_sbm(data$Y, Q = 2, type = "unipartite", distribution = "bernoulli")
#' print(fit$elbo)
#' print(fit$labels)
#' }
fit_sbm <- function(Y, Q, Q2 = NULL,
                    type = c("unipartite", "bipartite"),
                    distribution = c("gaussian", "bernoulli"),
                    block_sigma = FALSE,
                    use_cpp = TRUE,
                    max_iter = 100,
                    tol = 1e-6,
                    n_init = 5,
                    verbose = FALSE) {
  
  type <- match.arg(type)
  distribution <- match.arg(distribution)
  
  dims <- dim(Y)
  if (length(dims) != 3) {
    stop("Y must be a 3-dimensional array (T x N1 x N2)")
  }
  
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  if (type == "unipartite") {
    stopifnot(N1 == N2)
    Q2 <- Q  # Same for unipartite
  } else {
    if (is.null(Q2)) Q2 <- Q
  }
  
  if (use_cpp) {
    result <- .fit_cpp(Y, Q, Q2, type, distribution, block_sigma, 
                       max_iter, tol, n_init)
  } else {
    result <- .fit_r(Y, Q, Q2, type, distribution, block_sigma,
                     max_iter, tol, n_init, verbose)
  }
  
  # Add metadata
  result$type <- type
  result$distribution <- distribution
  result$Q1 <- Q
  result$Q2 <- Q2
  result$N1 <- N1
  result$N2 <- N2
  result$T_len <- T_len
  result$block_sigma <- block_sigma
  
  # Compute ICL
  result$icl <- .compute_icl(result, type, distribution)
  
  class(result) <- c("sbm_fit", "list")
  
  return(result)
}

#' Fit SBM to a Segment of Network Data
#'
#' @param Y Network tensor
#' @param start_time Start time (inclusive)
#' @param end_time End time (inclusive)
#' @param Q Number of clusters
#' @param Q2 Number of column clusters
#' @param type Network type
#' @param distribution Edge distribution
#' @param block_sigma Use block-specific variance
#' @param use_cpp Use C++ implementation
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param n_init Number of initializations
#' @param verbose Print progress
#' @return List with model parameters and fit statistics
#' @export
fit_sbm_segment <- function(Y, start_time, end_time,
                            Q, Q2 = NULL,
                            type = c("unipartite", "bipartite"),
                            distribution = c("gaussian", "bernoulli"),
                            block_sigma = FALSE,
                            use_cpp = TRUE,
                            max_iter = 100,
                            tol = 1e-6,
                            n_init = 5,
                            verbose = FALSE) {
  
  Y_segment <- Y[start_time:end_time, , , drop = FALSE]
  
  result <- fit_sbm(Y_segment, Q, Q2, type, distribution, block_sigma,
                    use_cpp, max_iter, tol, n_init, verbose)
  result$start_time <- start_time
  result$end_time <- end_time
  
  return(result)
}

# =============================================================================
# C++ WRAPPERS
# =============================================================================

#' Call C++ VEM fitting functions
#' @keywords internal
.fit_cpp <- function(Y, Q1, Q2, type, distribution, block_sigma, 
                     max_iter, tol, n_init) {
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  # Flatten Y for C++ (column-major order)
  Y_flat <- as.vector(Y)
  
  if (type == "bipartite") {
    if (distribution == "gaussian") {
      result <- fit_vem_bipartite_weighted_cpp(
        Y_flat, T_len, N1, N2, Q1, Q2,
        max_iter, tol, n_init, block_sigma
      )
      result$labels <- list(z1 = result$z1, z2 = result$z2)
    } else {
      result <- fit_vem_bipartite_binary_cpp(
        Y_flat, T_len, N1, N2, Q1, Q2,
        max_iter, tol, n_init
      )
      result$labels <- list(z1 = result$z1, z2 = result$z2)
    }
  } else {
    if (distribution == "gaussian") {
      result <- fit_vem_unipartite_weighted_cpp(
        Y_flat, T_len, N1, Q1,
        max_iter, tol, n_init, block_sigma
      )
      result$labels <- list(z = result$z)
      result$tau1 <- result$tau
      result$tau2 <- result$tau
    } else {
      result <- fit_vem_unipartite_binary_cpp(
        Y_flat, T_len, N1, Q1,
        max_iter, tol, n_init
      )
      result$labels <- list(z = result$z)
      result$tau1 <- result$tau
      result$tau2 <- result$tau
    }
  }
  
  return(result)
}

# =============================================================================
# PURE R IMPLEMENTATION (FALLBACK)
# =============================================================================

#' Pure R VEM implementation
#' @keywords internal
.fit_r <- function(Y, Q1, Q2, type, distribution, block_sigma,
                   max_iter, tol, n_init, verbose) {
  
  dims <- dim(Y)
  T_len <- dims[1]
  N1 <- dims[2]
  N2 <- dims[3]
  
  best_result <- NULL
  best_elbo <- -Inf
  
  for (init in 1:n_init) {
    if (verbose) message("Initialization ", init, "/", n_init)
    
    # Hard random assignment to break symmetry (prevents saddle-point traps)
    make_hard_tau <- function(n, q) {
      tau <- matrix(0.05 / max(1, q - 1), n, q)
      assigned <- sample.int(q, n, replace = TRUE)
      for (i in seq_len(n)) tau[i, assigned[i]] <- 1.0 - (q - 1) * 0.05 / max(1, q - 1)
      tau
    }
    if (type == "unipartite") {
      tau1 <- make_hard_tau(N1, Q1)
      tau2 <- tau1
    } else {
      tau1 <- make_hard_tau(N1, Q1)
      tau2 <- make_hard_tau(N2, Q2)
    }
    
    # Initialize parameters
    if (distribution == "gaussian") {
      mu <- matrix(rnorm(Q1 * Q2, 0, 1), Q1, Q2)
      if (type == "unipartite") mu <- (mu + t(mu)) / 2
      sigma_global <- 1
      sigma_block <- matrix(1, Q1, Q2)
    } else {
      prob <- matrix(runif(Q1 * Q2, 0.2, 0.8), Q1, Q2)
      if (type == "unipartite") prob <- (prob + t(prob)) / 2
    }
    
    pi1 <- rep(1/Q1, Q1)
    pi2 <- rep(1/Q2, Q2)
    
    elbo_prev <- -Inf
    converged <- FALSE
    elbo_trace <- numeric(0)
    
    for (iter in 1:max_iter) {
      # E-step: Update tau using vectorized operations
      if (type == "unipartite") {
        if (distribution == "gaussian") {
          sigma_use <- if (block_sigma) sigma_block else matrix(sigma_global, Q1, Q2)
          tau1 <- .e_step_uni_gauss_vec(Y, tau1, mu, sigma_use, pi1)
          tau2 <- tau1
        } else {
          tau1 <- .e_step_uni_bin_vec(Y, tau1, prob, pi1)
          tau2 <- tau1
        }
      } else {
        if (distribution == "gaussian") {
          sigma_use <- if (block_sigma) sigma_block else matrix(sigma_global, Q1, Q2)
          e_result <- .e_step_bip_gauss_vec(Y, tau1, tau2, mu, sigma_use, pi1, pi2)
          tau1 <- e_result$tau1
          tau2 <- e_result$tau2
        } else {
          e_result <- .e_step_bip_bin_vec(Y, tau1, tau2, prob, pi1, pi2)
          tau1 <- e_result$tau1
          tau2 <- e_result$tau2
        }
      }
      
      # M-step
      if (distribution == "gaussian") {
        m_result <- .m_step_gauss_vec(Y, tau1, tau2, type, block_sigma)
        mu <- m_result$mu
        sigma_global <- m_result$sigma_global
        sigma_block <- m_result$sigma_block
      } else {
        prob <- .m_step_bin_vec(Y, tau1, tau2, type)
      }
      
      # Update pi
      pi1 <- colMeans(tau1)
      pi1 <- pmax(pi1, 1e-10)
      pi1 <- pi1 / sum(pi1)
      
      pi2 <- colMeans(tau2)
      pi2 <- pmax(pi2, 1e-10)
      pi2 <- pi2 / sum(pi2)
      
      # Compute ELBO
      if (distribution == "gaussian") {
        sigma_use <- if (block_sigma) sigma_block else matrix(sigma_global, Q1, Q2)
        elbo <- .compute_elbo_gauss_vec(Y, tau1, tau2, mu, sigma_use, pi1, pi2, type)
      } else {
        elbo <- .compute_elbo_bin_vec(Y, tau1, tau2, prob, pi1, pi2, type)
      }
      
      elbo_trace <- c(elbo_trace, elbo)
      
      if (verbose && iter %% 10 == 0) {
        message("  Iter ", iter, ": ELBO = ", round(elbo, 4))
      }
      
      rel_change <- abs(elbo - elbo_prev) / max(1, abs(elbo))
      if (rel_change < tol) {
        converged <- TRUE
        break
      }
      elbo_prev <- elbo
    }
    
    if (elbo > best_elbo) {
      best_elbo <- elbo
      
      # Create result
      if (type == "unipartite") {
        labels <- list(z = apply(tau1, 1, which.max))
      } else {
        labels <- list(z1 = apply(tau1, 1, which.max),
                       z2 = apply(tau2, 1, which.max))
      }
      
      best_result <- list(
        tau1 = tau1,
        tau2 = tau2,
        labels = labels,
        elbo = elbo,
        converged = converged,
        iterations = iter,
        elbo_trace = elbo_trace,
        pi1 = pi1,
        pi2 = pi2
      )
      
      if (distribution == "gaussian") {
        best_result$mu <- mu
        best_result$sigma_global <- sigma_global
        best_result$sigma_block <- sigma_block
        best_result$sigma <- if (block_sigma) sigma_block else sigma_global
      } else {
        best_result$prob <- prob
      }
    }
  }
  
  return(best_result)
}

# =============================================================================
# VECTORIZED E-STEP FUNCTIONS
# =============================================================================

#' Vectorized E-step for unipartite Gaussian
#' @keywords internal
.e_step_uni_gauss_vec <- function(Y, tau, mu, sigma, pi) {
  T_len <- dim(Y)[1]
  N <- dim(Y)[2]
  Q <- ncol(tau)
  
  # Aggregate Y over time
  Y_sum <- apply(Y, c(2, 3), sum)
  Y_sq_sum <- apply(Y^2, c(2, 3), sum)
  
  # Compute log-likelihoods
  log_tau <- matrix(0, N, Q)
  
  for (q in 1:Q) {
    log_prior <- log(pi[q] + 1e-10)
    
    log_lik_contrib <- matrix(0, N, Q)  # Fixed: should be N x Q
    for (r in 1:Q) {
      # Contribution from cluster r neighbors
      tau_r <- tau[, r]
      s <- sigma[q, r]
      m <- mu[q, r]
      
      # For each node i, sum over neighbors j in cluster r
      # log p(Y[t,i,j] | mu[q,r], sigma[q,r]) = -0.5*log(2*pi*sigma^2) - (y-mu)^2/(2*sigma^2)
      const <- -0.5 * log(2 * base::pi * s^2) * T_len
      
      for (i in 1:N) {
        ll <- 0
        for (j in 1:N) {
          if (i != j) {
            ll <- ll + tau_r[j] * (const - (Y_sq_sum[i,j] - 2*m*Y_sum[i,j] + m^2*T_len) / (2*s^2))
          }
        }
        log_lik_contrib[i, r] <- ll
      }
    }
    log_tau[, q] <- log_prior + rowSums(log_lik_contrib)
  }
  
  # Softmax normalization
  log_tau <- log_tau - apply(log_tau, 1, max)
  tau_new <- exp(log_tau)
  tau_new <- tau_new / rowSums(tau_new)
  tau_new <- pmax(tau_new, 1e-10)
  
  return(tau_new)
}

#' Vectorized E-step for unipartite Binary
#' @keywords internal
.e_step_uni_bin_vec <- function(Y, tau, prob, pi) {
  T_len <- dim(Y)[1]
  N <- dim(Y)[2]
  Q <- ncol(tau)
  
  # Count edges over time
  edge_sum <- apply(Y, c(2, 3), sum)
  
  prob <- pmax(pmin(prob, 1 - 1e-10), 1e-10)
  log_prob <- log(prob)
  log_1_prob <- log(1 - prob)
  
  log_tau <- matrix(0, N, Q)
  
  for (q in 1:Q) {
    log_prior <- log(pi[q] + 1e-10)
    
    log_lik <- rep(log_prior, N)
    for (r in 1:Q) {
      tau_r <- tau[, r]
      lp <- log_prob[q, r]
      l1p <- log_1_prob[q, r]
      
      for (i in 1:N) {
        for (j in 1:N) {
          if (i != j) {
            y_sum <- edge_sum[i, j]
            log_lik[i] <- log_lik[i] + tau_r[j] * (y_sum * lp + (T_len - y_sum) * l1p)
          }
        }
      }
    }
    log_tau[, q] <- log_lik
  }
  
  log_tau <- log_tau - apply(log_tau, 1, max)
  tau_new <- exp(log_tau)
  tau_new <- tau_new / rowSums(tau_new)
  tau_new <- pmax(tau_new, 1e-10)
  
  return(tau_new)
}

#' Vectorized E-step for bipartite Gaussian
#' @keywords internal
.e_step_bip_gauss_vec <- function(Y, tau1, tau2, mu, sigma, pi1, pi2) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  Y_sum <- apply(Y, c(2, 3), sum)
  Y_sq_sum <- apply(Y^2, c(2, 3), sum)
  
  # Update tau1
  log_tau1 <- matrix(0, N1, Q1)
  for (q in 1:Q1) {
    log_prior <- log(pi1[q] + 1e-10)
    log_lik <- rep(log_prior, N1)
    
    for (l in 1:Q2) {
      tau_l <- tau2[, l]
      s <- sigma[q, l]
      m <- mu[q, l]
      const <- -0.5 * log(2 * pi * s^2) * T_len
      
      for (i in 1:N1) {
        ll <- 0
        for (j in 1:N2) {
          ll <- ll + tau_l[j] * (const - (Y_sq_sum[i,j] - 2*m*Y_sum[i,j] + m^2*T_len) / (2*s^2))
        }
        log_lik[i] <- log_lik[i] + ll
      }
    }
    log_tau1[, q] <- log_lik
  }
  log_tau1 <- log_tau1 - apply(log_tau1, 1, max)
  tau1_new <- exp(log_tau1)
  tau1_new <- tau1_new / rowSums(tau1_new)
  tau1_new <- pmax(tau1_new, 1e-10)
  
  # Update tau2
  log_tau2 <- matrix(0, N2, Q2)
  for (l in 1:Q2) {
    log_prior <- log(pi2[l] + 1e-10)
    log_lik <- rep(log_prior, N2)
    
    for (q in 1:Q1) {
      tau_q <- tau1_new[, q]
      s <- sigma[q, l]
      m <- mu[q, l]
      const <- -0.5 * log(2 * pi * s^2) * T_len
      
      for (j in 1:N2) {
        ll <- 0
        for (i in 1:N1) {
          ll <- ll + tau_q[i] * (const - (Y_sq_sum[i,j] - 2*m*Y_sum[i,j] + m^2*T_len) / (2*s^2))
        }
        log_lik[j] <- log_lik[j] + ll
      }
    }
    log_tau2[, l] <- log_lik
  }
  log_tau2 <- log_tau2 - apply(log_tau2, 1, max)
  tau2_new <- exp(log_tau2)
  tau2_new <- tau2_new / rowSums(tau2_new)
  tau2_new <- pmax(tau2_new, 1e-10)
  
  return(list(tau1 = tau1_new, tau2 = tau2_new))
}

#' Vectorized E-step for bipartite Binary
#' @keywords internal
.e_step_bip_bin_vec <- function(Y, tau1, tau2, prob, pi1, pi2) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  edge_sum <- apply(Y, c(2, 3), sum)
  
  prob <- pmax(pmin(prob, 1 - 1e-10), 1e-10)
  log_prob <- log(prob)
  log_1_prob <- log(1 - prob)
  
  # Update tau1
  log_tau1 <- matrix(0, N1, Q1)
  for (q in 1:Q1) {
    log_prior <- log(pi1[q] + 1e-10)
    log_lik <- rep(log_prior, N1)
    
    for (l in 1:Q2) {
      tau_l <- tau2[, l]
      lp <- log_prob[q, l]
      l1p <- log_1_prob[q, l]
      
      for (i in 1:N1) {
        for (j in 1:N2) {
          y_sum <- edge_sum[i, j]
          log_lik[i] <- log_lik[i] + tau_l[j] * (y_sum * lp + (T_len - y_sum) * l1p)
        }
      }
    }
    log_tau1[, q] <- log_lik
  }
  log_tau1 <- log_tau1 - apply(log_tau1, 1, max)
  tau1_new <- exp(log_tau1)
  tau1_new <- tau1_new / rowSums(tau1_new)
  tau1_new <- pmax(tau1_new, 1e-10)
  
  # Update tau2
  log_tau2 <- matrix(0, N2, Q2)
  for (l in 1:Q2) {
    log_prior <- log(pi2[l] + 1e-10)
    log_lik <- rep(log_prior, N2)
    
    for (q in 1:Q1) {
      tau_q <- tau1_new[, q]
      lp <- log_prob[q, l]
      l1p <- log_1_prob[q, l]
      
      for (j in 1:N2) {
        for (i in 1:N1) {
          y_sum <- edge_sum[i, j]
          log_lik[j] <- log_lik[j] + tau_q[i] * (y_sum * lp + (T_len - y_sum) * l1p)
        }
      }
    }
    log_tau2[, l] <- log_lik
  }
  log_tau2 <- log_tau2 - apply(log_tau2, 1, max)
  tau2_new <- exp(log_tau2)
  tau2_new <- tau2_new / rowSums(tau2_new)
  tau2_new <- pmax(tau2_new, 1e-10)
  
  return(list(tau1 = tau1_new, tau2 = tau2_new))
}

# =============================================================================
# VECTORIZED M-STEP FUNCTIONS
# =============================================================================

#' Vectorized M-step for Gaussian
#' @keywords internal
.m_step_gauss_vec <- function(Y, tau1, tau2, type, block_sigma) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  # Compute outer product of tau1 and tau2
  # weight[q,l] = sum_ij tau1[i,q] * tau2[j,l]
  if (type == "unipartite") {
    # Only upper triangle for symmetric
    weight_sum <- matrix(0, Q1, Q2)
    mu_sum <- matrix(0, Q1, Q2)
    var_sum <- matrix(0, Q1, Q2)
    
    for (i in 1:(N1-1)) {
      for (j in (i+1):N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_ij <- Y[, i, j]
        y_sum <- sum(y_ij)
        y_sq_sum <- sum(y_ij^2)
        
        weight_sum <- weight_sum + outer_tau * T_len
        mu_sum <- mu_sum + outer_tau * y_sum
      }
    }
    weight_sum <- weight_sum + t(weight_sum)
    mu_sum <- mu_sum + t(mu_sum)
    
    mu <- ifelse(weight_sum > 1e-10, mu_sum / weight_sum, 0)
    mu <- (mu + t(mu)) / 2
    
    # Variance
    var_sum <- matrix(0, Q1, Q2)
    weight_sum_var <- matrix(0, Q1, Q2)  # Track weights separately for variance
    sigma_num <- 0
    sigma_denom <- 0
    
    for (i in 1:(N1-1)) {
      for (j in (i+1):N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_ij <- Y[, i, j]
        
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- outer_tau[q, l]
            if (w > 1e-10) {
              diff_sq <- sum((y_ij - mu[q, l])^2)
              var_sum[q, l] <- var_sum[q, l] + w * diff_sq
              weight_sum_var[q, l] <- weight_sum_var[q, l] + w * T_len
              sigma_num <- sigma_num + w * diff_sq
              sigma_denom <- sigma_denom + w * T_len
            }
          }
        }
      }
    }
    # Symmetrize var_sum and weight_sum_var for unipartite
    var_sum <- var_sum + t(var_sum)
    weight_sum_var <- weight_sum_var + t(weight_sum_var)
  } else {
    weight_sum <- matrix(0, Q1, Q2)
    mu_sum <- matrix(0, Q1, Q2)
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_ij <- Y[, i, j]
        y_sum <- sum(y_ij)
        
        weight_sum <- weight_sum + outer_tau * T_len
        mu_sum <- mu_sum + outer_tau * y_sum
      }
    }
    
    mu <- ifelse(weight_sum > 1e-10, mu_sum / weight_sum, 0)
    
    # Variance
    var_sum <- matrix(0, Q1, Q2)
    weight_sum_var <- matrix(0, Q1, Q2)  # Track weights for variance
    sigma_num <- 0
    sigma_denom <- 0
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_ij <- Y[, i, j]
        
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- outer_tau[q, l]
            if (w > 1e-10) {
              diff_sq <- sum((y_ij - mu[q, l])^2)
              var_sum[q, l] <- var_sum[q, l] + w * diff_sq
              weight_sum_var[q, l] <- weight_sum_var[q, l] + w * T_len
              sigma_num <- sigma_num + w * diff_sq
              sigma_denom <- sigma_denom + w * T_len
            }
          }
        }
      }
    }
  }
  
  sigma_global <- sqrt(sigma_num / max(1, sigma_denom))
  sigma_global <- max(0.01, sigma_global)
  
  # Use weight_sum_var (which is correctly tracked) for block sigma
  sigma_block <- sqrt(ifelse(weight_sum_var > 1e-10, var_sum / weight_sum_var, 1))
  sigma_block <- pmax(0.01, sigma_block)
  if (type == "unipartite") sigma_block <- (sigma_block + t(sigma_block)) / 2
  
  return(list(mu = mu, sigma_global = sigma_global, sigma_block = sigma_block))
}

#' Vectorized M-step for Binary
#' @keywords internal
.m_step_bin_vec <- function(Y, tau1, tau2, type) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  num <- matrix(0, Q1, Q2)
  denom <- matrix(0, Q1, Q2)
  
  if (type == "unipartite") {
    for (i in 1:(N1-1)) {
      for (j in (i+1):N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_sum <- sum(Y[, i, j])
        
        num <- num + outer_tau * y_sum
        denom <- denom + outer_tau * T_len
      }
    }
    num <- num + t(num)
    denom <- denom + t(denom)
  } else {
    for (i in 1:N1) {
      for (j in 1:N2) {
        outer_tau <- outer(tau1[i,], tau2[j,])
        y_sum <- sum(Y[, i, j])
        
        num <- num + outer_tau * y_sum
        denom <- denom + outer_tau * T_len
      }
    }
  }
  
  prob <- ifelse(denom > 1e-10, num / denom, 0.5)
  prob <- pmax(pmin(prob, 1 - 1e-6), 1e-6)
  if (type == "unipartite") prob <- (prob + t(prob)) / 2
  
  return(prob)
}

# =============================================================================
# ELBO COMPUTATION
# =============================================================================

#' Compute ELBO for Gaussian
#' @keywords internal
.compute_elbo_gauss_vec <- function(Y, tau1, tau2, mu, sigma, pi1, pi2, type) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  elbo <- 0
  
  # Prior terms
  elbo <- elbo + sum(tau1 * log(matrix(pi1, N1, Q1, byrow = TRUE) + 1e-10))
  if (type == "bipartite") {
    elbo <- elbo + sum(tau2 * log(matrix(pi2, N2, Q2, byrow = TRUE) + 1e-10))
  }
  
  # Likelihood term
  if (type == "unipartite") {
    for (i in 1:(N1-1)) {
      for (j in (i+1):N2) {
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- tau1[i, q] * tau2[j, l]
            if (w > 1e-10) {
              for (t in 1:T_len) {
                elbo <- elbo + w * dnorm(Y[t, i, j], mu[q, l], sigma[q, l], log = TRUE)
              }
            }
          }
        }
      }
    }
    elbo <- elbo * 2  # Symmetry
  } else {
    for (i in 1:N1) {
      for (j in 1:N2) {
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- tau1[i, q] * tau2[j, l]
            if (w > 1e-10) {
              for (t in 1:T_len) {
                elbo <- elbo + w * dnorm(Y[t, i, j], mu[q, l], sigma[q, l], log = TRUE)
              }
            }
          }
        }
      }
    }
  }
  
  # Entropy term
  elbo <- elbo - sum(tau1 * log(tau1 + 1e-10))
  if (type == "bipartite") {
    elbo <- elbo - sum(tau2 * log(tau2 + 1e-10))
  }
  
  return(elbo)
}

#' Compute ELBO for Binary
#' @keywords internal
.compute_elbo_bin_vec <- function(Y, tau1, tau2, prob, pi1, pi2, type) {
  T_len <- dim(Y)[1]
  N1 <- dim(Y)[2]
  N2 <- dim(Y)[3]
  Q1 <- ncol(tau1)
  Q2 <- ncol(tau2)
  
  prob <- pmax(pmin(prob, 1 - 1e-10), 1e-10)
  
  elbo <- 0
  
  # Prior terms
  elbo <- elbo + sum(tau1 * log(matrix(pi1, N1, Q1, byrow = TRUE) + 1e-10))
  if (type == "bipartite") {
    elbo <- elbo + sum(tau2 * log(matrix(pi2, N2, Q2, byrow = TRUE) + 1e-10))
  }
  
  # Likelihood term
  if (type == "unipartite") {
    for (i in 1:(N1-1)) {
      for (j in (i+1):N2) {
        y_sum <- sum(Y[, i, j])
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- tau1[i, q] * tau2[j, l]
            if (w > 1e-10) {
              ll <- y_sum * log(prob[q, l]) + (T_len - y_sum) * log(1 - prob[q, l])
              elbo <- elbo + w * ll
            }
          }
        }
      }
    }
    elbo <- elbo * 2
  } else {
    for (i in 1:N1) {
      for (j in 1:N2) {
        y_sum <- sum(Y[, i, j])
        for (q in 1:Q1) {
          for (l in 1:Q2) {
            w <- tau1[i, q] * tau2[j, l]
            if (w > 1e-10) {
              ll <- y_sum * log(prob[q, l]) + (T_len - y_sum) * log(1 - prob[q, l])
              elbo <- elbo + w * ll
            }
          }
        }
      }
    }
  }
  
  # Entropy term
  elbo <- elbo - sum(tau1 * log(tau1 + 1e-10))
  if (type == "bipartite") {
    elbo <- elbo - sum(tau2 * log(tau2 + 1e-10))
  }
  
  return(elbo)
}

# =============================================================================
# ICL COMPUTATION
# =============================================================================

#' Compute ICL criterion
#' @keywords internal
.compute_icl <- function(fit, type, distribution) {
  N1 <- fit$N1
  N2 <- fit$N2
  T_len <- fit$T_len
  Q1 <- fit$Q1
  Q2 <- fit$Q2

  # Guard against malformed inputs from failed fits.
  if (!is.finite(fit$elbo) || !is.finite(N1) || !is.finite(N2) ||
      !is.finite(T_len) || !is.finite(Q1) || !is.finite(Q2)) {
    return(NA_real_)
  }

  N1 <- max(2L, as.integer(N1))
  N2 <- max(2L, as.integer(N2))
  T_len <- max(1L, as.integer(T_len))
  Q1 <- max(1L, as.integer(Q1))
  Q2 <- max(1L, as.integer(Q2))

  if (type == "unipartite") {
    n_block <- Q1 * (Q1 + 1) / 2
    n_mix <- Q1 - 1
    mix_scale <- log(max(2, N1))
    block_scale <- log(max(2, N1 * (N1 - 1) / 2 * T_len))
  } else {
    n_block <- Q1 * Q2
    n_mix <- (Q1 - 1) + (Q2 - 1)
    mix_scale <- log(max(2, max(N1, N2)))
    block_scale <- log(max(2, N1 * N2 * T_len))
  }

  if (distribution == "gaussian") {
    n_sigma <- if (isTRUE(fit$block_sigma)) n_block else 1L
  } else {
    n_sigma <- 0L
  }

  penalty <- 0.5 * n_mix * mix_scale +
    0.5 * (n_block + n_sigma) * block_scale

  as.numeric(fit$elbo) - penalty
}

# =============================================================================
# S3 METHODS
# =============================================================================

#' Print method for sbm_fit
#' @param x An `sbm_fit` object returned by `fit_sbm()`.
#' @param ... Additional arguments (currently unused).
#' @export
print.sbm_fit <- function(x, ...) {
  cat("\nSBM Model Fit\n")
  cat("=============\n\n")
  cat(sprintf("Type: %s %s\n", x$type, x$distribution))
  cat(sprintf("Dimensions: T=%d, N1=%d, N2=%d\n", x$T_len, x$N1, x$N2))
  cat(sprintf("Clusters: Q1=%d, Q2=%d\n", x$Q1, x$Q2))
  cat(sprintf("ELBO: %.4f\n", x$elbo))
  cat(sprintf("ICL: %.4f\n", x$icl))
  cat(sprintf("Converged: %s (iterations=%d)\n", x$converged, x$iterations))
  
  if (!is.null(x$mu)) {
    cat("\nMean matrix (mu):\n")
    print(round(x$mu, 4))
  }
  if (!is.null(x$prob)) {
    cat("\nProbability matrix (prob):\n")
    print(round(x$prob, 4))
  }
  
  invisible(x)
}
