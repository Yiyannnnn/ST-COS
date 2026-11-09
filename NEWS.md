# STCOS 0.1.1

- Fixed top-level seed propagation in gene-pattern generation by removing the
  inner seed reset from `generate_pattern()`.
- Aligned positive-definite matrix repair with the documented uniform
  eigenvalue-shift rule.
- Clarified that `p0` performs binomial thinning rather than zero-inflated
  mixture sampling.
- Added regression tests for seed handling, matrix repair, Cholesky
  orientation, and thinning.
- Added manuscript analysis, provenance, and validation materials under
  `reproducibility/`.
