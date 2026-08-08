// Compares exubercore::radf() against golden fixtures frozen from the
// current (pre-extraction) exuber R package's rls_gsadf(). Fixture format:
// <name>_input.csv (yxmat, header row + numeric rows), <name>_params.csv
// (min_win,lag,n_rows,n_cols header + one row), <name>_output.csv (value
// header + one column of expected results).
#include "exubercore/radf.hpp"

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

int main(int argc, char** argv) {
  std::string fixtures_dir = argc > 1 ? argv[1] : "tests/fixtures/golden";
  const char* cases[] = {"small_lag0", "small_lag1", "medium_lag0", "medium_lag2"};

  bool all_ok = true;
  for (const char* c : cases) {
    all_ok &= run_case(fixtures_dir, c);
  }

  return all_ok ? 0 : 1;
}
