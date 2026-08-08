# Changelog

## v0.1.0

- Initial extraction: `exubercore::radf()`, the recursive least-squares
  ADF/SADF/GSADF/BSADF statistic (Phillips, Shi & Yu 2015), ported verbatim
  from `exuber`'s `rls_gsadf()` (`exuber` src, pre-refactor). Pure Armadillo,
  no Rcpp/R dependency.
- Scope note: this release contains only the costly numerical routine.
  Monte Carlo / wild bootstrap / sieve bootstrap critical-value generation,
  date-stamping, and the bubble DGP simulators remain in `exuber`'s R code
  for now — they are RNG-driven orchestration around repeated calls to this
  routine, not the routine itself, and porting them is deferred.
- Golden fixtures under `tests/fixtures/golden/` frozen from the current
  `exuber` R package (`rls_gsadf()` on fixed-seed inputs, lag 0/1/2 at two
  sample sizes). Test tolerance: `1e-12` for the closed-form `lag == 0`
  path, `1e-9` for `lag > 0` (the O(n^2) sequential Sherman-Morrison update
  path accumulates cross-toolchain floating-point drift beyond 1e-12 even
  for bit-identical source — verified by compiling the unmodified
  pre-extraction source with the same toolchain and observing identical
  drift against the same golden output).
