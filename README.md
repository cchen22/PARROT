# PARROT

PARROT (Phase-Altering Regulatory Rewiring Over Time) is an R package for
change point detection in dynamic networks using stochastic block models (SBM).

It supports:

- unipartite weighted (Gaussian) networks
- unipartite binary (Bernoulli) networks
- bipartite weighted (Gaussian) networks
- bipartite binary (Bernoulli) networks

## Installation

Install directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("cchen22/PARROT")
```

Then load the package:

```r
library(PARROT)
```

## Quick Start: Single Change Point

```r
library(PARROT)

set.seed(1)

# Simulate a unipartite Bernoulli network with one change point
sim <- simulate_sbm_cp(
  N = 30,
  T_len = 60,
  cp_time = 30,
  Q = 2,
  type = "unipartite",
  distribution = "bernoulli",
  theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
  theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2))
)

# Detect change point (single-CP mode)
fit <- parrot(
  Y = sim$Y,
  Q = 2,
  multiple = FALSE,
  compute_cis = TRUE,
  compute_pvalues = TRUE,
  pvalue_method = "wilks"
)

print(fit)
summary(fit)
plot(fit)
```

## Quick Start: Multiple Change Points

```r
library(PARROT)

set.seed(2)

theta1 <- list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2))
theta2 <- list(prob = matrix(c(0.5, 0.5, 0.5, 0.5), 2, 2))
theta3 <- list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2))

sim_multi <- simulate_sbm_multi_cp(
  N = 30,
  T_len = 90,
  cp_times = c(30, 60),
  Q = 2,
  type = "unipartite",
  distribution = "bernoulli",
  theta_list = list(theta1, theta2, theta3)
)

fit_multi <- parrot(
  Y = sim_multi$Y,
  Q = 2,
  multiple = TRUE,
  method = "wbs",          # alternatives: "bs", "pelt"
  scan_method = "profile", # faster approximation: "score"
  max_cps = 5,
  compute_cis = TRUE,
  compute_pvalues = FALSE
)

print(fit_multi)
plot(fit_multi)
```

## Common Arguments

- `multiple`: `FALSE` for single-CP, `TRUE` for multiple-CP.
- `method`: multiple-CP search method (`"wbs"`, `"bs"`, `"pelt"`).
- `scan_method`: candidate split scan (`"profile"` or faster `"score"`).
- `compute_cis`: compute confidence intervals.
- `compute_pvalues`: compute p-values for detected change points.
- `pvalue_method`: `"wilks"`, `"permutation"`, or `"bootstrap"`.

## Citation

If you use PARROT in your research, please cite the corresponding manuscript
and repository:

- Repository: https://github.com/cchen22/PARROT
- Manuscript: 
