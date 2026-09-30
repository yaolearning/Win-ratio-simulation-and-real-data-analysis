// ============================================================================
// DIG thresholded pairwise kernel for 2- or 3-endpoint adverse-outcome analyses.
//
// General/GitHub version:
//   - Supports 2 or 3 endpoints
//   - Time-to-event thresholds are supplied in DAYS
//   - Binary/count endpoints use "lower is better"
//   - Evaluates all threshold vectors and endpoint orders in one subject-pair pass
//   - OpenMP parallelization is enabled when supported by the compiler
//
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]
// ============================================================================

#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>

using namespace Rcpp;

inline int compare_time(
    double a,
    double b,
    int ea,
    int eb,
    double t
) {
  if (!R_finite(a) ||
      !R_finite(b) ||
      ea == NA_INTEGER ||
      eb == NA_INTEGER) {
    return 0;
  }

  // For censored subjects, resolve only comparisons that are provable
  // from the observed event/censoring times.
  if (ea == 1 && eb == 1) {
    if (a - b > t) return 1;
    if (b - a > t) return -1;

  } else if (ea == 0 && eb == 1) {
    if (a - b > t) return 1;

  } else if (ea == 1 && eb == 0) {
    if (b - a > t) return -1;
  }

  return 0;
}

inline int compare_lower(double a, double b) {
  if (!R_finite(a) || !R_finite(b)) {
    return 0;
  }

  return (a < b) ? 1 : ((a > b) ? -1 : 0);
}

// [[Rcpp::export]]
List dig_evaluate_candidates_cpp(
    IntegerVector arm,
    IntegerVector type_code,
    NumericMatrix time_mat,
    IntegerMatrix event_mat,
    NumericMatrix value_mat,
    IntegerMatrix orders,
    NumericMatrix threshold_vectors,
    int n_threads = 8
) {
  const int n = arm.size();
  const int m = type_code.size();
  const int no = orders.nrow();
  const int nt = threshold_vectors.nrow();

  if (m != 2 && m != 3) {
    stop("Only 2 or 3 endpoints supported.");
  }

  if (
    orders.ncol() != m ||
    threshold_vectors.ncol() != m ||
    time_mat.ncol() != m ||
    event_mat.ncol() != m ||
    value_mat.ncol() != m ||
    time_mat.nrow() != n ||
    event_mat.nrow() != n ||
    value_mat.nrow() != n
  ) {
    stop("Inconsistent dimensions.");
  }

  if (n_threads < 1) {
    stop("n_threads must be positive.");
  }

  std::vector<int> ia;
  std::vector<int> ib;

  ia.reserve(n);
  ib.reserve(n);

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) {
      ia.push_back(i);
    } else if (arm[i] == 0) {
      ib.push_back(i);
    } else {
      stop("Treatment must be 0/1 and nonmissing.");
    }
  }

  if (ia.empty() || ib.empty()) {
    stop("Both treatment groups required.");
  }

  const long long na = ia.size();
  const long long nb = ib.size();

  const int nrows = nt * no;
  const int slots = nrows * m;

  std::vector<double> allw(slots, 0.0);
  std::vector<double> alll(slots, 0.0);
  std::vector<double> allu(slots, 0.0);

#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel
  {
    std::vector<double> lw(slots, 0.0);
    std::vector<double> ll(slots, 0.0);
    std::vector<double> lu(slots, 0.0);

#pragma omp for schedule(static)
    for (long long a = 0; a < na; ++a) {
      const int x = ia[a];

      for (long long b = 0; b < nb; ++b) {
        const int y = ib[b];

        for (int t = 0; t < nt; ++t) {
          int cmp[3] = {0, 0, 0};

          for (int ep = 0; ep < m; ++ep) {
            if (type_code[ep] == 1) {
              cmp[ep] = compare_time(
                time_mat(x, ep),
                time_mat(y, ep),
                event_mat(x, ep),
                event_mat(y, ep),
                threshold_vectors(t, ep)
              );
            } else {
              cmp[ep] = compare_lower(
                value_mat(x, ep),
                value_mat(y, ep)
              );
            }
          }

          for (int o = 0; o < no; ++o) {
            bool resolved = false;
            const int base = (t * no + o) * m;

            for (int r = 0; r < m && !resolved; ++r) {
              const int ep = orders(o, r) - 1;
              const int c = cmp[ep];

              if (c == 1) {
                lw[base + r] += 1.0;
                resolved = true;

              } else if (c == -1) {
                ll[base + r] += 1.0;
                resolved = true;

              } else {
                lu[base + r] += 1.0;
              }
            }
          }
        }
      }
    }

#pragma omp critical
    {
      for (int i = 0; i < slots; ++i) {
        allw[i] += lw[i];
        alll[i] += ll[i];
        allu[i] += lu[i];
      }
    }
  }

  NumericMatrix wins(nrows, m);
  NumericMatrix losses(nrows, m);
  NumericMatrix unresolved(nrows, m);

  for (int row = 0; row < nrows; ++row) {
    for (int r = 0; r < m; ++r) {
      const int slot = row * m + r;

      wins(row, r) = allw[slot];
      losses(row, r) = alll[slot];
      unresolved(row, r) = allu[slot];
    }
  }

  return List::create(
    _["wins"] = wins,
    _["losses"] = losses,
    _["unresolved_after"] = unresolved,
    _["n_pairs"] =
      static_cast<double>(na) * static_cast<double>(nb),
    _["n_treatment"] = static_cast<int>(na),
    _["n_control"] = static_cast<int>(nb)
  );
}

// [[Rcpp::export]]
int dig_openmp_threads_cpp() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
