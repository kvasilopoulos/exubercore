# Changelog

## v0.3.1

- This release fixes a precision problem in the closed-form SSR introduced in
  v0.3.0. `radf()` and `radf_nested()` now regress `dy` instead of `y` on
  the same regressors and take the t-statistic on `gamma = beta - 1`. The
  statistic is identical, but `SSR = dy'dy - b'X'dy` is a well-conditioned
  quantity. With `y` as the regressand, the R^2 of the level regression is
  close to 1, so the two O(y^2) sums nearly cancelled and the `lag > 0` path
  was sensitive, at first order, to Sherman-Morrison drift in `b`. At
  n = 2000 with levels around 500 it differed from `lm()` by about 1e-5,
  whereas v0.2.0 with explicit residuals differed by about 2e-8. The error
  is now about 3e-11 for lag 0 and lag 1 at that size, and the golden
  fixtures agree to 9e-14 (lag 0) and 2e-12 (lag 2). The speed is the same
  as in v0.3.0.

## v0.3.0

- `radf()` is now O(n^2) instead of O(n^3). Both branches keep running
  cross-products across the window grid and use the closed-form
  `SSR = y'y - b'X'y` that `radf_nested()` introduced. Before, the full
  residual vector was formed again for each of the roughly n^2/2 windows.
  With `lag == 0` the per-window pass `u = y - a - b x` is gone. With
  `lag > 0` the code still updates `(X'X)^-1` by Sherman-Morrison, but it no
  longer recomputes `X b` for each window, and it works on plain arrays
  instead of Armadillo temporaries. The output is unchanged: the golden
  fixtures agree to 2e-13 (lag 0) and 1e-11 (lag 2), and the `radf_nested()`
  cross-checks are as before. At n = 400 a path takes a few milliseconds
  instead of about 140 ms. Every Monte Carlo or bootstrap critical-value
  loop downstream pays this cost once per replication.

## v0.2.0

- Added `exubercore::radf_nested()`. It returns the `radf()` statistics for
  every sample size `n` in `[n_min, N]` of one path in a single O(N^2) sweep.
  Every window that a smaller `n` uses is a prefix window of the full path,
  so one pass over the (start, end) triangle serves all `n`. The pass keeps
  running sufficient statistics (`SSR = y'y - b'X'y`) and does not re-form
  residuals. The function returns the `badf` row (`W(0, e)` for every end
  `e`) and `gsadf` for each `n`. It does not return the `bsadf` sequences,
  because the critical values in exuber use only `cummax(badf)`. It was built
  to simulate critical-value tables for a whole grid of `n` at once. For
  N = 4000 a path costs about 0.2 s (lag 0) to 4 s (lag 4), where summing
  per-`n` calls to `radf()` takes hours.
- A new test compares `radf_nested()` with `radf()` on every prefix of a
  fixed-seed random walk for lags 0, 1, 2 and 4. The two agree to about
  1e-11, against a tolerance of 1e-7.
- `radf()` itself is unchanged.

## v0.1.0

- Initial extraction of `exubercore::radf()`, the recursive least-squares
  ADF, SADF, GSADF and BSADF statistic (Phillips, Shi & Yu 2015). It was
  ported unchanged from `rls_gsadf()` in `exuber` (the source before the
  refactor) and uses Armadillo only, with no Rcpp or R dependency.
- This release contains only the costly numerical routine. Critical-value
  generation by Monte Carlo, wild bootstrap and sieve bootstrap,
  date-stamping and the bubble DGP simulators stay in the R code of
  `exuber` for now. They are driven by random number generators and call
  this routine repeatedly, and porting them is deferred.
- The golden fixtures in `tests/fixtures/golden/` were frozen from the
  `exuber` R package, by running `rls_gsadf()` on fixed-seed inputs at lags
  0, 1 and 2 and two sample sizes. The test tolerance is `1e-12` for the
  closed-form `lag == 0` path and `1e-9` for `lag > 0`. The `lag > 0` path
  is an O(n^2) sequential Sherman-Morrison update that accumulates
  floating-point drift across toolchains beyond 1e-12, even for identical
  source. We checked this by compiling the unmodified pre-extraction source
  with the same toolchain: it showed the same drift against the same golden
  output.
