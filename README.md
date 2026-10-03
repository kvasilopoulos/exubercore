# exubercore

exubercore is a standalone C++ library that computes the recursive
least-squares ADF, SADF, GSADF and BSADF test statistics of Phillips, Shi and
Yu (2015). These are the statistics behind the `exuber` R package. The
library has no R or Python dependency: it takes an Armadillo matrix and
returns an `arma::vec`.

At present it contains one routine, `exubercore::radf()` (see
[include/exubercore/radf.hpp](include/exubercore/radf.hpp)), together with
`radf_nested()` for simulating critical values. `radf()` is the costly
numerical core, extracted unchanged from `rls_gsadf()` in `exuber`. Monte
Carlo and bootstrap critical values, date-stamping and the DGP simulators
are driven by random number generators and call this routine many times.
They remain in the host language of each downstream binding (see the
CHANGELOG).

## Build

You need a system installation of Armadillo with BLAS and LAPACK. Install
it with `apt install libarmadillo-dev`, `brew install armadillo`, or the
`armadillo` port in vcpkg.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

## Numerics

For `lag > 0`, `radf()` applies a sequential Sherman-Morrison rank-1 update
across O(n^2) windows. Results therefore cannot be guaranteed bit-identical
across toolchains at 1e-12. The CHANGELOG describes the check showing that
this is floating-point drift from compiler flags and not a correctness
problem. The closed-form `lag == 0` path agrees to 1e-12 across toolchains.
