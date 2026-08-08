# exubercore

Standalone C++ library for the recursive least-squares ADF/SADF/GSADF/BSADF
test statistic (Phillips, Shi & Yu 2015) underlying the `exuber` R package.
No R or Python dependency — plain Armadillo in, `arma::vec` out.

Currently ships one routine: `exubercore::radf()` (see
[include/exubercore/radf.hpp](include/exubercore/radf.hpp)), the costly
numerical core extracted verbatim from `exuber`'s `rls_gsadf()`. Monte
Carlo/bootstrap critical values, date-stamping, and DGP simulators are
RNG-driven orchestration that calls this routine repeatedly; they stay in
each downstream binding's host language for now (see CHANGELOG).

## Build

Requires a system Armadillo (with BLAS/LAPACK) — `apt install
libarmadillo-dev`, `brew install armadillo`, or vcpkg's `armadillo` port.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

## Numerics

`radf()`'s `lag > 0` path performs a sequential Sherman-Morrison rank-1
update across O(n^2) windows, so exact cross-toolchain bit-reproducibility
isn't guaranteed at 1e-12 (see CHANGELOG for the verification that this is
compiler-flag floating-point drift, not a correctness issue). The `lag ==
0` closed-form path matches to 1e-12 across toolchains.
