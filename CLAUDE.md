# exubercore

A standalone C++ library with one public routine, `exubercore::radf()`, the
recursive least-squares ADF, SADF, GSADF and BSADF statistic (Phillips, Shi
& Yu 2015). It uses Armadillo only and has no R or Python dependency.
[exuber](../exuber) calls it through Rcpp, and [pyexuber](../pyexuber) calls
it through pybind11, fetching it with CMake `FetchContent` at a pinned tag.
It is the shared numerical core of the two packages and is not meant to be
used on its own.

## Build

You need a system Armadillo with BLAS and LAPACK. Install it with
`apt install libarmadillo-dev`, `brew install armadillo`, or vcpkg's
`armadillo` port. `.github/workflows/ci.yml` has the exact invocation for
each operating system, including the Windows vcpkg toolchain-file flags.
Use those flags locally and do not invent new ones.

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Numerics: read this before tightening a tolerance

For `lag > 0`, `radf()` applies a sequential Sherman-Morrison rank-1 update
across O(n^2) windows, so results cannot be guaranteed bit-identical across
toolchains at 1e-12. We checked that this is not a bug. Compiling the
unmodified pre-extraction source with the same toolchain reproduces the
identical drift, so it comes from compiler-flag floating-point behaviour.
The tests use 1e-9 for `lag > 0` and 1e-12 for the closed-form `lag == 0`
path. The golden fixtures are in `tests/fixtures/golden/`, frozen from
`rls_gsadf()` in `exuber`, and CHANGELOG.md records where they came from.

## Scope: what is deliberately not here

Critical values from Monte Carlo, wild bootstrap and sieve bootstrap,
date-stamping and the bubble DGP simulators are driven by random number
generators and call `radf()` many times. They are not part of the routine
itself, so they stay in the host language of each binding (R for `exuber`,
and Python for `pyexuber`). Do not add them to this library. The decision
was made on purpose (see the "Scope note" in CHANGELOG.md).

## Methodology record: `../docs/`

`../docs/` holds the provenance of the statistic (Phillips, Shi & Yu 2015
and the literature around it), every downstream method built on `radf()`,
and the binding that ships each one. `README.md` there is the map, and
`parity.md` is the per-method table. The validation record of this library
is the golden-fixture suite described above. If `radf()` ever gains a second
routine, verify its numbers in `../docs/replication/` first and freeze them
into `tests/fixtures/golden/` afterwards.

## Release mechanism

The library is not distributed through a package registry such as CRAN, PyPI
or the vcpkg registry. A release is a version in `CHANGELOG.md` plus a git
tag, and each downstream binding consumes it by pinning that tag in its
`CMakeLists.txt` (`FetchContent_Declare`). A change to the signature or the
behaviour of `radf()` means bumping the pin in both `exuber/src/` (Rcpp glue)
and `pyexuber/` (pybind11 glue), so check both before treating a core change
as finished, as the root `CLAUDE.md` says. There is no deprecation cycle yet.
The library is pre-1.0, exports a single function and has no consumers
outside the two bindings in this workspace, so a breaking change needs only
the pin bumped everywhere and no compatibility shim.

## Writing style (all user-facing text)

Applies to READMEs, vignettes, the website, `docs/`, NEWS/CHANGELOG,
roxygen and docstrings, and any prose a reader sees. Code comments and
CLAUDE.md files follow it too.

**Voice.** An applied economist writing for colleagues who also want
ordinary readers to be able to run the test. Precise, sober, a little
plain-spoken. Define a term at first use (what "explosive" means, what a
critical value is for) and give the idea in words before the formula.

**Rewrite, do not substitute.** Swapping an em dash for a comma, colon or
hyphen keeps the machine-written rhythm and is not acceptable. If a
sentence needed a dash, it was carrying two thoughts: split it into two
sentences, or fold the aside into the grammar (a relative clause, a
parenthesis only for a true aside, or a separate sentence). No U+2014 and
no spaced hyphen standing in for one. En dashes stay for numeric ranges
and joint names (Phillips–Shi–Yu).

**Patterns to remove at the sentence level.**
- Fragments stacked for effect, and "X, not Y" or "not just X, but Y"
  framings. State the claim directly.
- Triplets used for rhythm, and sentences that announce what they are about
  to say ("Importantly,", "It is worth noting that", "In essence").
- Telegraphic notes (dropped articles, arrows, semicolon chains, "confirmed,
  zero new code"). Write full sentences with a subject and a verb.
- Status-report voice: "genuinely", "confirmed", "now done", "picked
  clean", "the most topical candidate". Say what is true and give the date
  if it matters.
- Marketing and filler words: seamlessly, robust (unless a statistical
  sense is stated), leverage, delve, comprehensive, powerful, crucial,
  landscape, journey, "under the hood", "a rich set of".
- Hedge stacks, and bold used as emphasis inside running prose.
- Self-reference to the writing process ("this resolves the question this
  file flagged", "an earlier pass"). Keep history in dated notes, not in
  the body of explanations.

**Do keep.** Formulas, numbers, citations, function names and every fact.
This is a change of language, not of content. Vary sentence length. Prefer
"we" or the imperative to the passive, and say what a function does and
