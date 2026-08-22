# exubercore

Standalone C++ library: `exubercore::radf()`, the recursive least-squares
ADF/SADF/GSADF/BSADF statistic (Phillips, Shi & Yu 2015). Pure Armadillo,
no R/Python dependency. Consumed by [exuber](../exuber) (Rcpp) and
[pyexuber](../pyexuber) (pybind11, via CMake `FetchContent` at a pinned
tag) — this is the shared numerical core, not a standalone user-facing tool.

## Build

Requires system Armadillo (+ BLAS/LAPACK): `apt install libarmadillo-dev`,
`brew install armadillo`, or vcpkg's `armadillo` port (see
`.github/workflows/ci.yml` for the exact per-OS invocations, including the
Windows vcpkg toolchain-file flags — match those locally rather than
guessing new ones).

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Numerics: don't tighten tolerances without reading this

`radf()`'s `lag > 0` path is a sequential Sherman-Morrison rank-1 update
across O(n^2) windows — cross-toolchain bit-reproducibility is **not**
guaranteed at 1e-12 (verified: compiling the unmodified pre-extraction
source with the same toolchain reproduces the identical drift, so it's
compiler-flag floating-point drift, not a bug). Tests use 1e-9 for
`lag > 0`, 1e-12 for the closed-form `lag == 0` path. Golden fixtures live
in `tests/fixtures/golden/`, frozen from `exuber`'s `rls_gsadf()` — see
CHANGELOG.md for provenance.

## Scope: what's deliberately NOT here

Monte Carlo/wild/sieve bootstrap critical values, date-stamping, and bubble
DGP simulators are RNG-driven orchestration around repeated calls to
`radf()`, not the routine itself — they stay in each downstream binding's
host language (R for `exuber`, eventually Python for `pyexuber`) rather
than being ported here. Don't add them to this library; that's a decision
already made, not an oversight (see CHANGELOG.md "Scope note").

## Release mechanism

No package registry (this isn't distributed via CRAN/PyPI/vcpkg registry
itself) — versioned by `CHANGELOG.md` + a git tag, consumed by pinning that
tag in the downstream binding's `CMakeLists.txt` (`FetchContent_Declare`).
Changing `radf()`'s signature or behavior requires bumping the pin in both
`exuber/src/` (Rcpp glue) and `pyexuber/` (pybind11 glue) — check both
before assuming a core change is done, per root `CLAUDE.md`. No
deprecation cycle exists yet (pre-1.0, single exported function, no
external consumers besides the two bindings in this same workspace) — a
breaking change just needs the pin bumped everywhere, not a compat shim.
