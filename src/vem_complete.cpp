// Complete VEM Fitting Functions for All Network Types
// Includes block-specific variance option
// Author: NetworkChangePoint Package

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

inline int idx3d(int t, int i, int j, int T, int N1) {
  return t + i * T + j * T * N1;
}

inline double log_sum_exp(const std::vector<double>& x) {
  if (x.empty()) return -1e100;
  double max_val = x[0];
  for (size_t i = 1; i < x.size(); i++) {
    if (x[i] > max_val) max_val = x[i];
  }
  if (!std::isfinite(max_val)) return -1e100;
  double sum_exp = 0.0;
  for (size_t i = 0; i < x.size(); i++) {
    sum_exp += std::exp(x[i] - max_val);
  }
  return max_val + std::log(sum_exp);
}

inline double log_dens_gaussian(double y, double mu, double sigma) {
  if (sigma <= 0) sigma = 0.01;
  double z = (y - mu) / sigma;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sigma) - 0.5 * z * z;
}

inline double log_dens_bernoulli(double y, double p) {
  p = std::max(1e-10, std::min(1.0 - 1e-10, p));
  return (y > 0.5) ? std::log(p) : std::log(1.0 - p);
}

// ============================================================================
// BIPARTITE WEIGHTED - FULL VEM WITH BLOCK-SPECIFIC VARIANCE
// ============================================================================

// [[Rcpp::export]]
List fit_vem_bipartite_weighted_cpp(
    NumericVector Y_flat,
    int T, int N1, int N2, int Q1, int Q2,
    int max_iter, double tol, int n_init,
    bool block_sigma) {
  
  double best_elbo = -1e100;
  List best_result;
  
  for (int init = 0; init < n_init; init++) {
    NumericMatrix tau1(N1, Q1);
    NumericMatrix tau2(N2, Q2);
    
    // Hard random assignment to break symmetry (critical for bipartite Gaussian)
    double eps = 0.05 / std::max(1, Q1 - 1);
    for (int i = 0; i < N1; i++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q1);
      if (assigned >= Q1) assigned = Q1 - 1;
      for (int q = 0; q < Q1; q++) {
        tau1(i, q) = (q == assigned) ? (1.0 - (Q1 - 1) * eps) : eps;
      }
    }
    
    eps = 0.05 / std::max(1, Q2 - 1);
    for (int j = 0; j < N2; j++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q2);
      if (assigned >= Q2) assigned = Q2 - 1;
      for (int l = 0; l < Q2; l++) {
        tau2(j, l) = (l == assigned) ? (1.0 - (Q2 - 1) * eps) : eps;
      }
    }
    
    // Initialize parameters
    NumericMatrix mu(Q1, Q2);
    NumericMatrix sigma_block(Q1, Q2);
    double sigma_global = 1.0;
    
    for (int q = 0; q < Q1; q++) {
      for (int l = 0; l < Q2; l++) {
        mu(q, l) = R::rnorm(0, 1);
        sigma_block(q, l) = 1.0;
      }
    }
    
    NumericVector pi1(Q1, 1.0 / Q1);
    NumericVector pi2(Q2, 1.0 / Q2);
    
    double elbo_prev = -1e100;
    bool converged = false;
    int iter = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
      if (iter > 0) {
      // E-step: Update tau1
      NumericMatrix new_tau1(N1, Q1);
      for (int i = 0; i < N1; i++) {
        std::vector<double> log_prob(Q1);
        for (int q = 0; q < Q1; q++) {
          double lp = std::log(pi1[q] + 1e-10);
          for (int j = 0; j < N2; j++) {
            for (int l = 0; l < Q2; l++) {
              double term = 0.0;
              double sig = block_sigma ? sigma_block(q, l) : sigma_global;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N1)];
                term += log_dens_gaussian(y, mu(q, l), sig);
              }
              lp += tau2(j, l) * term;
            }
          }
          log_prob[q] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int q = 0; q < Q1; q++) {
          new_tau1(i, q) = std::max(1e-10, std::exp(log_prob[q] - lse));
        }
      }
      tau1 = new_tau1;
      
      // E-step: Update tau2
      NumericMatrix new_tau2(N2, Q2);
      for (int j = 0; j < N2; j++) {
        std::vector<double> log_prob(Q2);
        for (int l = 0; l < Q2; l++) {
          double lp = std::log(pi2[l] + 1e-10);
          for (int i = 0; i < N1; i++) {
            for (int q = 0; q < Q1; q++) {
              double term = 0.0;
              double sig = block_sigma ? sigma_block(q, l) : sigma_global;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N1)];
                term += log_dens_gaussian(y, mu(q, l), sig);
              }
              lp += tau1(i, q) * term;
            }
          }
          log_prob[l] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int l = 0; l < Q2; l++) {
          new_tau2(j, l) = std::max(1e-10, std::exp(log_prob[l] - lse));
        }
      }
      tau2 = new_tau2;
      
      }
      // M-step: Update mu and sigma
      NumericMatrix weight_sum(Q1, Q2);
      NumericMatrix mu_new(Q1, Q2);
      NumericMatrix var_sum(Q1, Q2);
      double sigma_num = 0, sigma_denom = 0;
      
      // First pass: compute mu
      for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
          for (int q = 0; q < Q1; q++) {
            for (int l = 0; l < Q2; l++) {
              double w = tau1(i, q) * tau2(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  mu_new(q, l) += w * Y_flat[idx3d(t, i, j, T, N1)];
                }
                weight_sum(q, l) += w * T;
              }
            }
          }
        }
      }
      
      for (int q = 0; q < Q1; q++) {
        for (int l = 0; l < Q2; l++) {
          if (weight_sum(q, l) > 1e-10) {
            mu(q, l) = mu_new(q, l) / weight_sum(q, l);
          }
        }
      }
      
      // Second pass: compute variance
      for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
          for (int q = 0; q < Q1; q++) {
            for (int l = 0; l < Q2; l++) {
              double w = tau1(i, q) * tau2(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  double diff = Y_flat[idx3d(t, i, j, T, N1)] - mu(q, l);
                  var_sum(q, l) += w * diff * diff;
                  sigma_num += w * diff * diff;
                }
                sigma_denom += w * T;
              }
            }
          }
        }
      }
      
      // Update sigma
      if (block_sigma) {
        for (int q = 0; q < Q1; q++) {
          for (int l = 0; l < Q2; l++) {
            if (weight_sum(q, l) > 1e-10) {
              sigma_block(q, l) = std::sqrt(var_sum(q, l) / weight_sum(q, l));
              sigma_block(q, l) = std::max(0.01, sigma_block(q, l));
            }
          }
        }
      } else {
        sigma_global = std::sqrt(sigma_num / std::max(1.0, sigma_denom));
        sigma_global = std::max(0.01, sigma_global);
      }
      
      // Update pi
      for (int q = 0; q < Q1; q++) {
        pi1[q] = 0;
        for (int i = 0; i < N1; i++) pi1[q] += tau1(i, q);
        pi1[q] /= N1;
      }
      for (int l = 0; l < Q2; l++) {
        pi2[l] = 0;
        for (int j = 0; j < N2; j++) pi2[l] += tau2(j, l);
        pi2[l] /= N2;
      }
      
      // Compute ELBO
      double elbo = 0.0;
      
      // Entropy terms
      for (int i = 0; i < N1; i++) {
        for (int q = 0; q < Q1; q++) {
          if (tau1(i, q) > 1e-10 && pi1[q] > 1e-10) {
            elbo += tau1(i, q) * (std::log(pi1[q]) - std::log(tau1(i, q)));
          }
        }
      }
      for (int j = 0; j < N2; j++) {
        for (int l = 0; l < Q2; l++) {
          if (tau2(j, l) > 1e-10 && pi2[l] > 1e-10) {
            elbo += tau2(j, l) * (std::log(pi2[l]) - std::log(tau2(j, l)));
          }
        }
      }
      
      // Likelihood term
      for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
          for (int q = 0; q < Q1; q++) {
            for (int l = 0; l < Q2; l++) {
              double w = tau1(i, q) * tau2(j, l);
              if (w > 1e-10) {
                double sig = block_sigma ? sigma_block(q, l) : sigma_global;
                for (int t = 0; t < T; t++) {
                  elbo += w * log_dens_gaussian(Y_flat[idx3d(t, i, j, T, N1)], mu(q, l), sig);
                }
              }
            }
          }
        }
      }
      
      // Check convergence
      if (std::abs(elbo - elbo_prev) < tol) {
        converged = true;
        break;
      }
      elbo_prev = elbo;
    }
    
    // Compute final ELBO
    double final_elbo = elbo_prev;
    
    if (final_elbo > best_elbo) {
      best_elbo = final_elbo;
      
      // Get cluster assignments
      IntegerVector z1(N1), z2(N2);
      for (int i = 0; i < N1; i++) {
        int best_q = 0;
        for (int q = 1; q < Q1; q++) {
          if (tau1(i, q) > tau1(i, best_q)) best_q = q;
        }
        z1[i] = best_q + 1;
      }
      for (int j = 0; j < N2; j++) {
        int best_l = 0;
        for (int l = 1; l < Q2; l++) {
          if (tau2(j, l) > tau2(j, best_l)) best_l = l;
        }
        z2[j] = best_l + 1;
      }
      
      // Create a 1x1 matrix for global sigma if needed
      NumericMatrix sigma_single(1, 1);
      sigma_single(0, 0) = sigma_global;
      
      best_result = List::create(
        Named("tau1") = tau1,
        Named("tau2") = tau2,
        Named("mu") = mu,
        Named("sigma") = block_sigma ? sigma_block : sigma_single,
        Named("sigma_global") = sigma_global,
        Named("sigma_block") = sigma_block,
        Named("pi1") = pi1,
        Named("pi2") = pi2,
        Named("elbo") = best_elbo,
        Named("converged") = converged,
        Named("iterations") = iter + 1,
        Named("z1") = z1,
        Named("z2") = z2,
        Named("block_sigma") = block_sigma
      );
    }
  }
  
  return best_result;
}

// ============================================================================
// UNIPARTITE WEIGHTED - FULL VEM WITH BLOCK-SPECIFIC VARIANCE
// ============================================================================

// [[Rcpp::export]]
List fit_vem_unipartite_weighted_cpp(
    NumericVector Y_flat,
    int T, int N, int Q,
    int max_iter, double tol, int n_init,
    bool block_sigma) {
  
  double best_elbo = -1e100;
  List best_result;
  
  for (int init = 0; init < n_init; init++) {
    NumericMatrix tau(N, Q);
    double eps = 0.05 / std::max(1, Q - 1);
    for (int i = 0; i < N; i++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q);
      if (assigned >= Q) assigned = Q - 1;
      for (int q = 0; q < Q; q++) {
        tau(i, q) = (q == assigned) ? (1.0 - (Q - 1) * eps) : eps;
      }
    }
    
    // Initialize parameters (symmetric)
    NumericMatrix mu(Q, Q);
    NumericMatrix sigma_block(Q, Q);
    double sigma_global = 1.0;
    
    for (int q = 0; q < Q; q++) {
      for (int l = q; l < Q; l++) {
        mu(q, l) = mu(l, q) = R::rnorm(0, 1);
        sigma_block(q, l) = sigma_block(l, q) = 1.0;
      }
    }
    
    NumericVector pi(Q, 1.0 / Q);
    
    double elbo_prev = -1e100;
    bool converged = false;
    int iter = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
      if (iter > 0) {
      // E-step: Update tau
      NumericMatrix new_tau(N, Q);
      for (int i = 0; i < N; i++) {
        std::vector<double> log_prob(Q);
        for (int q = 0; q < Q; q++) {
          double lp = std::log(pi[q] + 1e-10);
          for (int j = 0; j < N; j++) {
            if (i == j) continue;
            for (int l = 0; l < Q; l++) {
              double term = 0.0;
              double sig = block_sigma ? sigma_block(q, l) : sigma_global;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N)];
                term += log_dens_gaussian(y, mu(q, l), sig);
              }
              lp += tau(j, l) * term;
            }
          }
          log_prob[q] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int q = 0; q < Q; q++) {
          new_tau(i, q) = std::max(1e-10, std::exp(log_prob[q] - lse));
        }
      }
      tau = new_tau;
      
      }
      // M-step: Update mu and sigma (use upper triangle, then symmetrize)
      NumericMatrix weight_sum(Q, Q);
      NumericMatrix mu_new(Q, Q);
      NumericMatrix var_sum(Q, Q);
      double sigma_num = 0, sigma_denom = 0;
      
      // First pass: compute mu
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int q = 0; q < Q; q++) {
            for (int l = 0; l < Q; l++) {
              double w = tau(i, q) * tau(j, l) + tau(i, l) * tau(j, q);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  mu_new(q, l) += w * Y_flat[idx3d(t, i, j, T, N)];
                }
                weight_sum(q, l) += w * T;
              }
            }
          }
        }
      }
      
      for (int q = 0; q < Q; q++) {
        for (int l = q; l < Q; l++) {
          double ws = weight_sum(q, l) + weight_sum(l, q);
          double ms = mu_new(q, l) + mu_new(l, q);
          if (ws > 1e-10) {
            mu(q, l) = mu(l, q) = ms / ws;
          }
        }
      }
      
      // Second pass: compute variance
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int q = 0; q < Q; q++) {
            for (int l = 0; l < Q; l++) {
              // Must use same weight scheme as mu computation
              double w = tau(i, q) * tau(j, l) + tau(i, l) * tau(j, q);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  double diff = Y_flat[idx3d(t, i, j, T, N)] - mu(q, l);
                  var_sum(q, l) += w * diff * diff;
                  sigma_num += w * diff * diff;
                }
                sigma_denom += w * T;
              }
            }
          }
        }
      }
      
      // Update sigma (symmetrize)
      if (block_sigma) {
        for (int q = 0; q < Q; q++) {
          for (int l = q; l < Q; l++) {
            double vs = var_sum(q, l) + var_sum(l, q);
            double ws = weight_sum(q, l) + weight_sum(l, q);
            if (ws > 1e-10) {
              sigma_block(q, l) = sigma_block(l, q) = std::sqrt(vs / ws);
              sigma_block(q, l) = sigma_block(l, q) = std::max(0.01, sigma_block(q, l));
            }
          }
        }
      } else {
        sigma_global = std::sqrt(sigma_num / std::max(1.0, sigma_denom));
        sigma_global = std::max(0.01, sigma_global);
      }
      
      // Update pi
      for (int q = 0; q < Q; q++) {
        pi[q] = 0;
        for (int i = 0; i < N; i++) pi[q] += tau(i, q);
        pi[q] /= N;
      }
      
      // Compute ELBO
      double elbo = 0.0;
      for (int i = 0; i < N; i++) {
        for (int q = 0; q < Q; q++) {
          if (tau(i, q) > 1e-10 && pi[q] > 1e-10) {
            elbo += tau(i, q) * (std::log(pi[q]) - std::log(tau(i, q)));
          }
        }
      }
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int q = 0; q < Q; q++) {
            for (int l = 0; l < Q; l++) {
              double w = tau(i, q) * tau(j, l);
              if (w > 1e-10) {
                double sig = block_sigma ? sigma_block(q, l) : sigma_global;
                for (int t = 0; t < T; t++) {
                  elbo += w * log_dens_gaussian(Y_flat[idx3d(t, i, j, T, N)], mu(q, l), sig);
                }
              }
            }
          }
        }
      }
      
      if (std::abs(elbo - elbo_prev) < tol) {
        converged = true;
        break;
      }
      elbo_prev = elbo;
    }
    
    double final_elbo = elbo_prev;
    
    if (final_elbo > best_elbo) {
      best_elbo = final_elbo;
      
      IntegerVector z(N);
      for (int i = 0; i < N; i++) {
        int best_q = 0;
        for (int q = 1; q < Q; q++) {
          if (tau(i, q) > tau(i, best_q)) best_q = q;
        }
        z[i] = best_q + 1;
      }
      
      // Create a 1x1 matrix for global sigma if needed
      NumericMatrix sigma_single(1, 1);
      sigma_single(0, 0) = sigma_global;
      
      best_result = List::create(
        Named("tau") = tau,
        Named("mu") = mu,
        Named("sigma") = block_sigma ? sigma_block : sigma_single,
        Named("sigma_global") = sigma_global,
        Named("sigma_block") = sigma_block,
        Named("pi") = pi,
        Named("elbo") = best_elbo,
        Named("converged") = converged,
        Named("iterations") = iter + 1,
        Named("z") = z,
        Named("block_sigma") = block_sigma
      );
    }
  }
  
  return best_result;
}

// ============================================================================
// BIPARTITE BINARY - FULL VEM
// ============================================================================

// [[Rcpp::export]]
List fit_vem_bipartite_binary_cpp(
    NumericVector Y_flat,
    int T, int N1, int N2, int Q1, int Q2,
    int max_iter, double tol, int n_init) {
  
  double best_elbo = -1e100;
  List best_result;
  
  for (int init = 0; init < n_init; init++) {
    NumericMatrix tau1(N1, Q1);
    NumericMatrix tau2(N2, Q2);
    
    double eps = 0.05 / std::max(1, Q1 - 1);
    for (int i = 0; i < N1; i++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q1);
      if (assigned >= Q1) assigned = Q1 - 1;
      for (int q = 0; q < Q1; q++) {
        tau1(i, q) = (q == assigned) ? (1.0 - (Q1 - 1) * eps) : eps;
      }
    }
    
    eps = 0.05 / std::max(1, Q2 - 1);
    for (int j = 0; j < N2; j++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q2);
      if (assigned >= Q2) assigned = Q2 - 1;
      for (int l = 0; l < Q2; l++) {
        tau2(j, l) = (l == assigned) ? (1.0 - (Q2 - 1) * eps) : eps;
      }
    }
    
    NumericMatrix prob(Q1, Q2);
    for (int q = 0; q < Q1; q++) {
      for (int l = 0; l < Q2; l++) {
        prob(q, l) = R::runif(0.2, 0.8);
      }
    }
    
    NumericVector pi1(Q1, 1.0 / Q1);
    NumericVector pi2(Q2, 1.0 / Q2);
    
    double elbo_prev = -1e100;
    bool converged = false;
    int iter = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
      // E-step: tau1
      NumericMatrix new_tau1(N1, Q1);
      for (int i = 0; i < N1; i++) {
        std::vector<double> log_prob(Q1);
        for (int q = 0; q < Q1; q++) {
          double lp = std::log(pi1[q] + 1e-10);
          for (int j = 0; j < N2; j++) {
            for (int l = 0; l < Q2; l++) {
              double term = 0.0;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N1)];
                term += log_dens_bernoulli(y, prob(q, l));
              }
              lp += tau2(j, l) * term;
            }
          }
          log_prob[q] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int q = 0; q < Q1; q++) {
          new_tau1(i, q) = std::max(1e-10, std::exp(log_prob[q] - lse));
        }
      }
      tau1 = new_tau1;
      
      // E-step: tau2
      NumericMatrix new_tau2(N2, Q2);
      for (int j = 0; j < N2; j++) {
        std::vector<double> log_prob(Q2);
        for (int l = 0; l < Q2; l++) {
          double lp = std::log(pi2[l] + 1e-10);
          for (int i = 0; i < N1; i++) {
            for (int q = 0; q < Q1; q++) {
              double term = 0.0;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N1)];
                term += log_dens_bernoulli(y, prob(q, l));
              }
              lp += tau1(i, q) * term;
            }
          }
          log_prob[l] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int l = 0; l < Q2; l++) {
          new_tau2(j, l) = std::max(1e-10, std::exp(log_prob[l] - lse));
        }
      }
      tau2 = new_tau2;
      
      // M-step: prob
      NumericMatrix num(Q1, Q2);
      NumericMatrix denom(Q1, Q2);
      
      for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
          for (int q = 0; q < Q1; q++) {
            for (int l = 0; l < Q2; l++) {
              double w = tau1(i, q) * tau2(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  num(q, l) += w * Y_flat[idx3d(t, i, j, T, N1)];
                }
                denom(q, l) += w * T;
              }
            }
          }
        }
      }
      
      for (int q = 0; q < Q1; q++) {
        for (int l = 0; l < Q2; l++) {
          if (denom(q, l) > 1e-10) {
            prob(q, l) = num(q, l) / denom(q, l);
            prob(q, l) = std::max(1e-6, std::min(1.0 - 1e-6, prob(q, l)));
          }
        }
      }
      
      // Update pi
      for (int q = 0; q < Q1; q++) {
        pi1[q] = 0;
        for (int i = 0; i < N1; i++) pi1[q] += tau1(i, q);
        pi1[q] /= N1;
      }
      for (int l = 0; l < Q2; l++) {
        pi2[l] = 0;
        for (int j = 0; j < N2; j++) pi2[l] += tau2(j, l);
        pi2[l] /= N2;
      }
      
      // ELBO
      double elbo = 0.0;
      for (int i = 0; i < N1; i++) {
        for (int q = 0; q < Q1; q++) {
          if (tau1(i, q) > 1e-10 && pi1[q] > 1e-10) {
            elbo += tau1(i, q) * (std::log(pi1[q]) - std::log(tau1(i, q)));
          }
        }
      }
      for (int j = 0; j < N2; j++) {
        for (int l = 0; l < Q2; l++) {
          if (tau2(j, l) > 1e-10 && pi2[l] > 1e-10) {
            elbo += tau2(j, l) * (std::log(pi2[l]) - std::log(tau2(j, l)));
          }
        }
      }
      for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) {
          for (int q = 0; q < Q1; q++) {
            for (int l = 0; l < Q2; l++) {
              double w = tau1(i, q) * tau2(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  elbo += w * log_dens_bernoulli(Y_flat[idx3d(t, i, j, T, N1)], prob(q, l));
                }
              }
            }
          }
        }
      }
      
      if (std::abs(elbo - elbo_prev) < tol) {
        converged = true;
        break;
      }
      elbo_prev = elbo;
    }
    
    if (elbo_prev > best_elbo) {
      best_elbo = elbo_prev;
      
      IntegerVector z1(N1), z2(N2);
      for (int i = 0; i < N1; i++) {
        int best_q = 0;
        for (int q = 1; q < Q1; q++) {
          if (tau1(i, q) > tau1(i, best_q)) best_q = q;
        }
        z1[i] = best_q + 1;
      }
      for (int j = 0; j < N2; j++) {
        int best_l = 0;
        for (int l = 1; l < Q2; l++) {
          if (tau2(j, l) > tau2(j, best_l)) best_l = l;
        }
        z2[j] = best_l + 1;
      }
      
      best_result = List::create(
        Named("tau1") = tau1,
        Named("tau2") = tau2,
        Named("prob") = prob,
        Named("pi1") = pi1,
        Named("pi2") = pi2,
        Named("elbo") = best_elbo,
        Named("converged") = converged,
        Named("iterations") = iter + 1,
        Named("z1") = z1,
        Named("z2") = z2
      );
    }
  }
  
  return best_result;
}

// ============================================================================
// UNIPARTITE BINARY - FULL VEM
// ============================================================================

// [[Rcpp::export]]
List fit_vem_unipartite_binary_cpp(
    NumericVector Y_flat,
    int T, int N, int Q,
    int max_iter, double tol, int n_init) {
  
  double best_elbo = -1e100;
  List best_result;
  
  for (int init = 0; init < n_init; init++) {
    NumericMatrix tau(N, Q);
    double eps = 0.05 / std::max(1, Q - 1);
    for (int i = 0; i < N; i++) {
      int assigned = (int)(R::runif(0.0, 1.0) * Q);
      if (assigned >= Q) assigned = Q - 1;
      for (int q = 0; q < Q; q++) {
        tau(i, q) = (q == assigned) ? (1.0 - (Q - 1) * eps) : eps;
      }
    }
    
    NumericMatrix prob(Q, Q);
    for (int q = 0; q < Q; q++) {
      for (int l = q; l < Q; l++) {
        prob(q, l) = prob(l, q) = R::runif(0.2, 0.8);
      }
    }
    
    NumericVector pi(Q, 1.0 / Q);
    
    double elbo_prev = -1e100;
    bool converged = false;
    int iter = 0;
    
    for (iter = 0; iter < max_iter; iter++) {
      // E-step
      NumericMatrix new_tau(N, Q);
      for (int i = 0; i < N; i++) {
        std::vector<double> log_prob(Q);
        for (int q = 0; q < Q; q++) {
          double lp = std::log(pi[q] + 1e-10);
          for (int j = 0; j < N; j++) {
            if (i == j) continue;
            for (int l = 0; l < Q; l++) {
              double term = 0.0;
              for (int t = 0; t < T; t++) {
                double y = Y_flat[idx3d(t, i, j, T, N)];
                term += log_dens_bernoulli(y, prob(q, l));
              }
              lp += tau(j, l) * term;
            }
          }
          log_prob[q] = lp;
        }
        double lse = log_sum_exp(log_prob);
        for (int q = 0; q < Q; q++) {
          new_tau(i, q) = std::max(1e-10, std::exp(log_prob[q] - lse));
        }
      }
      tau = new_tau;
      
      // M-step
      NumericMatrix num(Q, Q);
      NumericMatrix denom(Q, Q);
      
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int q = 0; q < Q; q++) {
            for (int l = 0; l < Q; l++) {
              double w = tau(i, q) * tau(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  double y = Y_flat[idx3d(t, i, j, T, N)];
                  num(q, l) += w * y;
                  num(l, q) += w * y;
                }
                denom(q, l) += w * T;
                denom(l, q) += w * T;
              }
            }
          }
        }
      }
      
      for (int q = 0; q < Q; q++) {
        for (int l = q; l < Q; l++) {
          if (denom(q, l) > 1e-10) {
            prob(q, l) = prob(l, q) = num(q, l) / denom(q, l);
            prob(q, l) = prob(l, q) = std::max(1e-6, std::min(1.0 - 1e-6, prob(q, l)));
          }
        }
      }
      
      for (int q = 0; q < Q; q++) {
        pi[q] = 0;
        for (int i = 0; i < N; i++) pi[q] += tau(i, q);
        pi[q] /= N;
      }
      
      // ELBO
      double elbo = 0.0;
      for (int i = 0; i < N; i++) {
        for (int q = 0; q < Q; q++) {
          if (tau(i, q) > 1e-10 && pi[q] > 1e-10) {
            elbo += tau(i, q) * (std::log(pi[q]) - std::log(tau(i, q)));
          }
        }
      }
      for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
          for (int q = 0; q < Q; q++) {
            for (int l = 0; l < Q; l++) {
              double w = tau(i, q) * tau(j, l);
              if (w > 1e-10) {
                for (int t = 0; t < T; t++) {
                  elbo += w * log_dens_bernoulli(Y_flat[idx3d(t, i, j, T, N)], prob(q, l));
                }
              }
            }
          }
        }
      }
      
      if (std::abs(elbo - elbo_prev) < tol) {
        converged = true;
        break;
      }
      elbo_prev = elbo;
    }
    
    if (elbo_prev > best_elbo) {
      best_elbo = elbo_prev;
      
      IntegerVector z(N);
      for (int i = 0; i < N; i++) {
        int best_q = 0;
        for (int q = 1; q < Q; q++) {
          if (tau(i, q) > tau(i, best_q)) best_q = q;
        }
        z[i] = best_q + 1;
      }
      
      best_result = List::create(
        Named("tau") = tau,
        Named("prob") = prob,
        Named("pi") = pi,
        Named("elbo") = best_elbo,
        Named("converged") = converged,
        Named("iterations") = iter + 1,
        Named("z") = z
      );
    }
  }
  
  return best_result;
}
