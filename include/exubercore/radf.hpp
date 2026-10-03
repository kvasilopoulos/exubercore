#pragma once

#include <armadillo>

namespace exubercore {

// Recursive least-squares ADF, SADF, GSADF and BSADF test statistics
// (Phillips, Shi & Yu 2015). The matrix inversion lemma updates the inverse
// of the regressor cross-product from one window to the next, so it is never
// inverted from scratch.
//
// yxmat: column 0 is the dependent variable (levels). The remaining columns
// are the regressors (constant, lag, lagged differences) as produced by
// exuber's `unroot()`. min_win: minimum regression window size. lag: ADF
// lag length. A lag of 0 selects the closed-form single-regressor path.
//
// Returns a vector of length 2*total + 3, where total = nrow(yxmat) -
// min_win + 1:
//   [0, total)            badf
//   [total]                adf
//   [total + 1]            sadf
//   [total + 2]            gsadf
//   [total + 3, 2*total+3) bsadf
//
// Throws std::invalid_argument if min_win is not a valid window size for
// yxmat.
arma::vec radf(const arma::mat& yxmat, int min_win, int lag = 0);

// The statistics of radf() for every sample size n in [n_min, N] of one
// path, computed in a single pass. yxmat is unroot() of the full length-N
// path (R = N - 1 - lag rows). minw[k] is the minimum window for
// n = n_min + k, so K = minw.n_elem and N = n_min + K - 1. Every window that
// a smaller n uses is a prefix window of the full path, so one sweep over the
// (start, end) triangle serves all n. The sweep keeps running sufficient
// statistics (SSR = y'y - b'X'y) and does not re-form residuals, which makes
// it O(R^2) where a loop over radf() is O(R^3).
//
// Returns a vector of length R + K:
//   [0, R)     W(0, e): the ADF t-statistic on rows [0, e]. This is the
//              badf row of every n at once (NaN for e < min(minw) - 1). For
//              sample size n with window m, badf = W(0, m-1 .. R_n-1),
//              adf = W(0, R_n-1) and sadf = the maximum of that badf, where
//              R_n = n - 1 - lag.
//   [R, R+K)   gsadf for n = n_min + k.
// The bsadf sequences are not returned. The critical values in exuber use
// only cummax(badf), and the per-n sequences would add O(N^2) output.
//
// Throws std::invalid_argument if the sizes are inconsistent or if min(minw)
// is too small to fit the regression.
arma::vec radf_nested(const arma::mat& yxmat, const arma::ivec& minw, int n_min, int lag = 0);

} // namespace exubercore
