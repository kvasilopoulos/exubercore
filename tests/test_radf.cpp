// Compares exubercore::radf() against golden fixtures frozen from the
// current (pre-extraction) exuber R package's rls_gsadf(). Fixture format:
// <name>_input.csv (yxmat, header row + numeric rows), <name>_params.csv
// (min_win,lag,n_rows,n_cols header + one row), <name>_output.csv (value
// header + one column of expected results).
#include "exubercore/radf.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::vector<double>> read_csv_matrix(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::string line;
  std::getline(f, line); // header
  std::vector<std::vector<double>> rows;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::stringstream ss(line);
    std::vector<double> row;
    std::string cell;
    while (std::getline(ss, cell, ',')) row.push_back(std::stod(cell));
    rows.push_back(row);
  }
  return rows;
}

arma::mat to_arma(const std::vector<std::vector<double>>& rows) {
  if (rows.empty()) return arma::mat();
  arma::mat m(rows.size(), rows[0].size());
  for (size_t i = 0; i < rows.size(); ++i)
    for (size_t j = 0; j < rows[i].size(); ++j) m(i, j) = rows[i][j];
  return m;
}

bool run_case(const std::string& fixtures_dir, const std::string& name) {
  auto input_rows = read_csv_matrix(fixtures_dir + "/" + name + "_input.csv");
  auto params_rows = read_csv_matrix(fixtures_dir + "/" + name + "_params.csv");
  auto output_rows = read_csv_matrix(fixtures_dir + "/" + name + "_output.csv");

  arma::mat yxmat = to_arma(input_rows);
  int min_win = static_cast<int>(params_rows[0][0]);
  int lag = static_cast<int>(params_rows[0][1]);

  arma::vec expected(output_rows.size());
  for (size_t i = 0; i < output_rows.size(); ++i) expected(i) = output_rows[i][0];

  arma::vec actual = exubercore::radf(yxmat, min_win, lag);

  if (actual.n_elem != expected.n_elem) {
    std::printf("[FAIL] %s: length mismatch actual=%zu expected=%zu\n",
                name.c_str(), (size_t)actual.n_elem, (size_t)expected.n_elem);
    return false;
  }

  double max_abs_diff = 0.0;
  for (arma::uword i = 0; i < actual.n_elem; ++i) {
    double a = actual(i), e = expected(i);
    double diff = (std::isnan(a) && std::isnan(e)) ? 0.0 : std::fabs(a - e);
    max_abs_diff = std::max(max_abs_diff, diff);
  }

  // The lag > 0 path performs O(total^2) sequential rank-1 (Sherman-Morrison)
  // updates; floating-point non-associativity across compilers/flags can
  // accumulate past 1e-12 even for bit-identical source (confirmed by
  // compiling the unmodified pre-extraction exuber source with this same
  // toolchain and diffing against the same golden output -- identical
  // drift). 1e-12 holds for the closed-form lag == 0 path; lag > 0 gets a
  // looser bound that still catches a real translation bug.
  double tol = lag == 0 ? 1e-12 : 1e-9;
  bool ok = max_abs_diff < tol;
  std::printf("[%s] %s: max_abs_diff=%.3e (n=%zu, min_win=%d, lag=%d)\n",
              ok ? "PASS" : "FAIL", name.c_str(), max_abs_diff,
              (size_t)actual.n_elem, min_win, lag);
  return ok;
}

} // namespace

// radf_nested() must agree with radf() run separately on every prefix of
// the same path. No fixture: a fixed-seed random walk, exuber's psy_minw()
// per n, both lag paths. Tolerance is loose-ish because the nested sweep
// uses running sufficient statistics rather than explicit residuals.
arma::mat unroot(const arma::vec& y, int lag) {
  const int n = static_cast<int>(y.n_elem);
  const int rows = n - 1 - lag;
  arma::mat m(rows, lag == 0 ? 2 : lag + 3);
  for (int r = 0; r < rows; ++r) {
    const int t = r + lag + 1;
    if (lag == 0) {
      m(r, 0) = y(t); m(r, 1) = y(t - 1);
    } else {
      m(r, 0) = y(t); m(r, 1) = 1.0; m(r, 2) = y(t - 1);
      for (int k = 1; k <= lag; ++k) m(r, 2 + k) = y(t - k) - y(t - k - 1);
    }
  }
  return m;
}

int psy_minw(int n) { return static_cast<int>(std::floor((0.01 + 1.8 / std::sqrt(n)) * n)); }

bool run_nested_case(int lag, int n_min, int N) {
  arma::vec y(N);
  unsigned long long s = 12345 + lag;  // tiny LCG, deterministic across platforms
  double acc = 0;
  for (int t = 0; t < N; ++t) {
    s = s * 6364136223846793005ULL + 1442695040888963407ULL;
    acc += ((s >> 11) * (1.0 / 9007199254740992.0)) - 0.5;
    y(t) = acc;
  }
  const int K = N - n_min + 1;
  arma::ivec minw(K);
  for (int k = 0; k < K; ++k) minw(k) = psy_minw(n_min + k);
  arma::vec nested = exubercore::radf_nested(unroot(y, lag), minw, n_min, lag);
  const int R = N - 1 - lag;

  double worst = 0;
  for (int k = 0; k < K; ++k) {
    const int n = n_min + k, m = minw(k), Rn = n - 1 - lag, total = Rn - m + 1;
    arma::vec ref = exubercore::radf(unroot(y.head(n), lag), m, lag);
    arma::vec badf = nested.subvec(m - 1, Rn - 1);
    double d = 0;
    d = std::max(d, arma::abs(badf - ref.head(total)).max());
    d = std::max(d, std::fabs(badf(total - 1) - ref(total)));          // adf
    d = std::max(d, std::fabs(badf.max() - ref(total + 1)));            // sadf
    d = std::max(d, std::fabs(nested(R + k) - ref(total + 2)));         // gsadf
    worst = std::max(worst, d);
  }
  const double tol = 1e-7;
  std::printf("[%s] nested lag=%d n=%d..%d: max |diff| = %.3g (tol %.0e)\n",
              worst <= tol ? "OK" : "FAIL", lag, n_min, N, worst, tol);
  return worst <= tol;
}

int main(int argc, char** argv) {
  std::string fixtures_dir = argc > 1 ? argv[1] : "tests/fixtures/golden";
  const char* cases[] = {"small_lag0", "small_lag1", "medium_lag0", "medium_lag2"};

  bool all_ok = true;
  for (const char* c : cases) {
    all_ok &= run_case(fixtures_dir, c);
  }

  all_ok &= run_nested_case(0, 6, 200);
  all_ok &= run_nested_case(1, 8, 200);
  all_ok &= run_nested_case(2, 10, 150);
  all_ok &= run_nested_case(4, 16, 150);

  return all_ok ? 0 : 1;
}
