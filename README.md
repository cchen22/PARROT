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

## Quick Start: Single Change Point (Permutation P-value)

```r
library(PARROT)

set.seed(1)

# Simulate a unipartite Bernoulli network with one change point
sim <- simulate_sbm_cp(
  N = 24,
  T_len = 40,
  cp_time = 20,
  Q = 2,
  type = "unipartite",
  distribution = "bernoulli",
  theta1 = list(prob = matrix(c(0.8, 0.2, 0.2, 0.8), 2, 2)),
  theta2 = list(prob = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2))
)

# Detect a single change point and compute permutation p-value
fit <- parrot(
  Y = sim$Y,
  Q = 2,
  multiple = FALSE,
  compute_cis = FALSE,
  compute_pvalues = TRUE,
  pvalue_method = "permutation"
)

print(fit)

# Inspect permutation p-value for the detected change point
if (length(fit$pvalues) > 0) {
  fit$pvalues[[1]][c("cp", "pvalue_method", "pvalue_raw", "pvalue_bonferroni")]
}

plot(fit)
```

## Common Arguments

- `multiple`: use `FALSE` for fast single-CP detection.
- `scan_method`: candidate split scan (`"profile"` or faster `"score"`).
- `compute_cis`: compute confidence intervals.
- `compute_pvalues`: compute p-values for detected change points.
- `pvalue_method`: `"wilks"`, `"permutation"`, or `"bootstrap"`.

## Citation

If you use PARROT in your research, please cite the corresponding manuscript
and repository:

- Repository: https://github.com/cchen22/PARROT
