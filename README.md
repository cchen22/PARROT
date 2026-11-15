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

set.seed(42)

# Simulate a unipartite Bernoulli network with one change point at t=25
sim <- simulate_sbm_cp(
  N = 30,
  T_len = 50,
  cp_time = 25,
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
  scan_method = "profile",
  compute_cis = FALSE,
  compute_pvalues = TRUE,
  pvalue_method = "permutation"
)

print(fit)

# Inspect permutation p-value for the detected change point
if (length(fit$pvalues) > 0) {
  fit$pvalues[[1]][c("cp", "pvalue_method", "pvalue_raw", "pvalue_bh", "pvalue_bonferroni")]
}

plot(fit)
```

## Common Arguments

- `Q`: number of communities; `Q = NULL` selects Q with the ICL gain-ratio rule
  over `Q_range` (`q_stability = TRUE` adds a seed-stability check).
- `multiple`: `FALSE` for one change point, `TRUE` for several (`method = "wbs"`
  or `"bs"`, at most `max_cps`; `stop_rule = "test"` for test-based splitting).
- `scan_method`: candidate split scan, `"profile"` (default) or `"score"`.
- `pvalue_method`: `"permutation"` (default), `"wilks"`, `"stationary_bootstrap"`
  or `"bootstrap"`.
- `n_init`, `max_iter`: random initialisations and VEM iterations per fit.
- `n_cores`: parallel segment fits for the resampling tests.

## Speed versus accuracy

The defaults favour accuracy. With T time points and minimum segment length m,
the cost is dominated by the number of SBM fits:

| Setting | SBM fits | Speed | Accuracy |
|---|---|---|---|
| `scan_method = "score"` | 1 (whole series) | fast | approximate; can miss changes that leave the whole-series fit without block structure (e.g. two communities swapping their connection pattern) |
| `scan_method = "profile"` (default) | 2(T - 2m + 1) segment fits | slower | exact profile likelihood of every split |
| `pvalue_method = "wilks"` | none | fast | chi-squared approximation; below the nominal level in our simulations, but not guaranteed (no finite-sample justification) |
| `pvalue_method = "permutation"` (default) | the scan is repeated on every permuted series | slowest | valid type-I error for independent time points; exact when all configurations (at most `exact_max`) are enumerated. The statistic is the maximum profile likelihood ratio, also with `scan_method = "score"` |
| `pvalue_method = "stationary_bootstrap"` | as permutation | slowest | for temporally dependent series |
| `rescan = FALSE` | 2 per permuted series | fast | anti-conservative (the split is held fixed although it was chosen by the scan); not recommended |
| `n_init = 20`, `max_iter = 200` (default) | per fit | slower | fewer local optima; reduce (e.g. 5 and 50) for quick exploration |

Segment fits are cached by their set of time points, so the multiple change
point search and the tests never refit the same segment, and `n_cores > 1`
parallelises them without changing the result. For long series or large
networks, a fast exploratory run is `scan_method = "score", pvalue_method =
"wilks"`; the reported result should use the defaults.

Type-I error at the 5% level in the simulations of the paper (independent time
points, 6,000 series without change, 99 permutations): re-scan permutation 0.9%
(at most 2.6% for every network type and baseline) and Wilks 1.2%; fixed-split
permutation (`rescan = FALSE`) 5.5% overall but about 10% for binary networks. The parametric bootstrap
(`"bootstrap"`) also holds the split fixed and should not be used for decisions.
The permutation test is the default because its validity does not rest on an
asymptotic approximation.

## Citation

If you use PARROT in your research, please cite the corresponding manuscript
and repository:

- Chen C, Padi M, Quackenbush J. PARROT: Phase-Altering Regulatory Rewiring Over
  Time. bioRxiv 2026.06.24.734262 (2026).
  doi: [10.64898/2026.06.24.734262](https://doi.org/10.64898/2026.06.24.734262)
- Repository: https://github.com/cchen22/PARROT

```bibtex
@article{chen2026parrot,
  title   = {{PARROT}: Phase-Altering Regulatory Rewiring Over Time},
  author  = {Chen, Chen and Padi, Megha and Quackenbush, John},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.06.24.734262},
  note    = {Preprint 2026.06.24.734262}
}
```
