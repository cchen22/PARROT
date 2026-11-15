ari <- function(a, b) PARROT:::.ari(a, b)

test_that("Gaussian bipartite fits recover an easy block structure (no collapse)", {
  set.seed(1)
  z1 <- rep(1:2, each = 4); z2 <- rep(1:2, each = 5)
  M <- matrix(c(2, 0, 0, 2), 2)
  Y <- array(0, c(6, 8, 10))
  for (t in 1:6) Y[t, , ] <- M[z1, z2] + matrix(rnorm(80, 0, 0.5), 8)
  res <- vapply(1:10, function(s) {
    set.seed(s)
    f <- fit_sbm(Y, 2, 2, "bipartite", "gaussian", FALSE, TRUE, 200, 1e-6, 5, FALSE)
    ari(f$labels$z1, z1)
  }, numeric(1))
  expect_true(all(res == 1))
})

test_that("permutation p-values: exact enumeration, validity and reproducibility", {
  set.seed(2)
  sim <- simulate_sbm_cp(N = 12, T_len = 7, cp_time = 4, Q = 2, type = "unipartite",
                         distribution = "bernoulli",
                         theta1 = list(prob = matrix(c(.8, .2, .2, .8), 2)),
                         theta2 = list(prob = matrix(c(.2, .8, .8, .2), 2)))
  run <- function(rescan, n_cores = 1L) {
    set.seed(3)
    PARROT:::.permutation_pvalue(sim$Y, 4L, 1L, 7L, 2, 2, "unipartite", "bernoulli", FALSE,
                                 n_init = 3, max_iter = 50, tol = 1e-5, min_segment = 2L,
                                 rescan = rescan, n_cores = n_cores)
  }
  fx <- run(FALSE); rs <- run(TRUE); rs2 <- run(TRUE, 2L)
  expect_equal(fx$method, "exact"); expect_equal(fx$perm_space_size, choose(7, 4))
  expect_equal(rs$method, "exact"); expect_equal(rs$perm_space_size, factorial(7))
  for (x in list(fx, rs)) expect_true(x$pvalue >= 1 / x$perm_space_size && x$pvalue <= 1)
  expect_true(rs$pvalue >= fx$pvalue)          # the re-scan accounts for the selection of the split
  expect_identical(rs$pvalue, rs2$pvalue)       # segment seeds: independent of the number of cores
})

test_that("ICL gain-ratio rule picks the Q after which the gain drops most", {
  expect_equal(PARROT:::.q_gain_ratio_rule(1:5, c(0, 10, 18, 19, 19.5)), 3L)
  # candidates beyond the ICL maximum are ignored
  expect_equal(PARROT:::.q_gain_ratio_rule(1:6, c(0, 10, 18, 19, 18, 17)), 3L)
})

test_that("parrot() runs with the default and the fast settings", {
  set.seed(4)
  sim <- simulate_sbm_cp(N = 12, T_len = 10, cp_time = 5, Q = 2, type = "unipartite",
                         distribution = "bernoulli",
                         theta1 = list(prob = matrix(c(.8, .2, .2, .8), 2)),
                         theta2 = list(prob = matrix(c(.2, .8, .8, .2), 2)))
  r1 <- parrot(sim$Y, Q = 2, multiple = FALSE, min_segment = 2, n_init = 10, compute_cis = FALSE)
  r2 <- parrot(sim$Y, Q = 2, multiple = FALSE, min_segment = 2, n_init = 10, compute_cis = FALSE,
               scan_method = "score", pvalue_method = "wilks")
  expect_true(length(r1$cps) == 1 && abs(r1$cps - 5) <= 1)
  expect_true(length(r2$cps) <= 1)   # score scan: fast approximation (may miss or shift the change)
  expect_equal(r1$pvalues[[1]]$permutation_type, "monte_carlo")   # 10! orderings > exact_max
  expect_true(isTRUE(r1$pvalues[[1]]$rescan))
})

test_that("single-CP profile scan does not depend on the global seed and matches the test", {
  set.seed(5)
  sim <- simulate_sbm_cp(N = 12, T_len = 8, cp_time = 4, Q = 2, type = "unipartite",
                         distribution = "bernoulli",
                         theta1 = list(prob = matrix(c(.8, .2, .2, .8), 2)),
                         theta2 = list(prob = matrix(c(.2, .8, .8, .2), 2)))
  run <- function(s) { set.seed(s)
    parrot(sim$Y, Q = 2, distribution = "bernoulli", multiple = FALSE, min_segment = 2,
           n_init = 3, max_iter = 50, pvalue_method = "wilks", compute_cis = FALSE,
           compute_pvalues = FALSE)$detection_result }
  a <- run(1); b <- run(99)
  expect_identical(a$cp_raw, b$cp_raw)
  expect_equal(a$profile_ll, b$profile_ll)
  E <- a$elbo_fun
  expect_equal(a$lr_stat, PARROT:::.lr_from_cache(E, 1:a$cp_raw, (a$cp_raw + 1):8))
})
