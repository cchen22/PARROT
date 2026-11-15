# PARROT 1.1.0

## Bug fixes

* Gaussian VEM fitters (unipartite and bipartite): the first iteration ran the
  E-step with random parameters (means N(0, 1), sigma = 1), which overwrote the
  random hard initialisation and made many fits collapse to a single community.
  The first iteration now starts with the M-step.
* `compute_pvalue(method = "permutation" / "bootstrap")` drew chi-squared values
  instead of resampling, and `test_significance()` used a placeholder statistic.
  Both now call the resampling tests below.
* The permutation test used weaker fit settings for the permuted series (2
  initialisations, 30 iterations) than for the observed series. Observed and
  permuted statistics are now computed with identical settings.

* Single-change-point profile scan: the segments were refitted with the global
  random state, so the best split could depend on the seed when two splits were
  close, and the permutation test (which uses seeded segment fits) could evaluate
  a different split than the one reported. The scan, the Wilks statistic and the
  tests now use the same seeded segment fits (argument `fit_seed`, default 1);
  results no longer depend on `set.seed()`.

## Changes to defaults

* `scan_method = "profile"` (refit of both segments at every split, i.e. the
  profile likelihood); `"score"` remains available as the fast option.
* `pvalue_method = "permutation"`: time-label permutation test conditional on the
  window bounded by the neighbouring change points, with the scan repeated within
  the window on every permuted series (`rescan = TRUE`). The previous test held the
  split fixed at the detected position, which ignores that the split was selected
  by the scan and is anti-conservative (about 10% false positives at the 5% level
  for small binary networks); it remains available as `rescan = FALSE`.
  `permutation_scope = "global"` (time labels of the whole series are permuted).
  All distinct configurations are enumerated when there are at most `exact_max`
  (default 5040), otherwise `n_perm` (default 99) random permutations are used.
  Segment fits can be run in parallel (`n_cores`); every segment fit has its own
  seed derived from its time points, so results do not depend on `n_cores`.
  In single-change-point mode the change point is declared by this test.
* `Q = NULL`: Q is selected with the ICL gain-ratio rule over `Q_range = 1:10`
  (the Q at which the marginal ICL gain drops most, among the values up to the
  ICL maximum). `q_stability = TRUE` additionally reports, for every Q, the mean
  pairwise adjusted Rand index of the memberships from `q_n_seeds` independent
  sets of initialisations.
* `n_init = 20`, `max_iter = 200`, `block_sigma = FALSE`, `max_cps = 3`,
  `n_intervals = 50`, `compute_pvalues = TRUE`.

## New features

* `pvalue_method = "stationary_bootstrap"` (Politis and Romano, 1994; mean block
  length `sb_block`, default ceiling(T^(1/3))) for temporally dependent series.
* `stop_rule = "test"` for multiple change points: test-based binary
  segmentation (a segment is split while the p-value of its best split is below
  `alpha`, at most `max_cps` change points). The default `stop_rule = "max_cps"`
  reports up to `max_cps` candidates, each with raw, Bonferroni and BH adjusted
  p-values (`cps_significant`: Bonferroni-significant candidates).
* Segment fits are cached by their set of time points, so WBS/BS and the
  permutation test never refit the same segment; the profile-scan candidate
  search uses the same fit settings as the test.
