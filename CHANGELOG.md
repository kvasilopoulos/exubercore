# Changelog

## v0.2.0

- Add `exubercore::radf_nested()`: the `radf()` statistics for every sample
  size `n` in `[n_min, N]` of one path in a single O(N^2) sweep. Every
  window a smaller `n` uses is a prefix window of the full path, so one pass
  over the (start, end) triangle -- with running sufficient statistics
  (`SSR = y'y - b'X'y`) instead of re-formed residuals -- serves all `n`.
  Returns the `badf` row (`W(0, e)` for every end `e`) plus `gsadf` per `n`;
  `bsadf` sequences are not returned (exuber's critical values only use
  `cummax(badf)`). Built for simulating critical-value tables across a whole
  grid of `n` at once: N = 4000 costs ~0.2 s (lag 0) to ~4 s (lag 4) per
  path versus hours summed over per-`n` `radf()` calls.
- Test: `radf_nested()` vs `radf()` on every prefix of a fixed-seed random
  walk, lag 0/1/2/4, agrees to ~1e-11 (tolerance 1e-7).
- `radf()` itself is unchanged.

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
