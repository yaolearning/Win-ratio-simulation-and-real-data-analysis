// dig3_openmp_engine.cpp
// OpenMP pairwise engine for DIG 3-endpoint WR/WO analysis

// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <vector>
#include <cmath>
using namespace Rcpp;

inline int compare_time_endpoint(double ta, double tb, int ea, int eb) {
  if (!R_finite(ta) || !R_finite(tb)) return 0;
  if (ea == 1 && eb == 1) {
    if (ta > tb) return 1;
    if (tb > ta) return -1;
    return 0;
  }
  if (ea == 0 && eb == 1) {
    if (ta > tb) return 1;
    return 0;
  }
  if (ea == 1 && eb == 0) {
    if (tb > ta) return -1;
    return 0;
  }
  return 0;
}

inline int compare_lower_better(double va, double vb) {
  if (!R_finite(va) || !R_finite(vb)) return 0;
  if (va < vb) return 1;
  if (va > vb) return -1;
  return 0;
}

// [[Rcpp::export]]
List evaluate_six_orders_openmp(
    IntegerVector arm,
    IntegerVector type_code,
    NumericMatrix time_mat,
    IntegerMatrix event_mat,
    NumericMatrix value_mat,
    IntegerMatrix orders,
    int n_threads = 8) {

  const int n = arm.size();
  const int n_orders = orders.nrow();
  const int m = orders.ncol();

  std::vector<int> idxA, idxB;
  idxA.reserve(n);
  idxB.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) idxA.push_back(i);
    else if (arm[i] == 0) idxB.push_back(i);
  }

  const long long nA = static_cast<long long>(idxA.size());
  const long long nB = static_cast<long long>(idxB.size());
  const double n_pairs = static_cast<double>(nA) * static_cast<double>(nB);

  NumericMatrix wins(n_orders, m);
  NumericMatrix losses(n_orders, m);
  NumericVector ties(n_orders);

#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#endif

#pragma omp parallel
  {
    std::vector<double> local_wins(n_orders * m, 0.0);
    std::vector<double> local_losses(n_orders * m, 0.0);
    std::vector<double> local_ties(n_orders, 0.0);

#pragma omp for schedule(static)
    for (long long aa = 0; aa < nA; ++aa) {
      const int ia = idxA[aa];

      for (long long bb = 0; bb < nB; ++bb) {
        const int ib = idxB[bb];

        int cmp[3] = {0, 0, 0};
        for (int ep = 0; ep < m; ++ep) {
          const int typ = type_code[ep];
          if (typ == 1) {
            cmp[ep] = compare_time_endpoint(
              time_mat(ia, ep), time_mat(ib, ep),
              event_mat(ia, ep), event_mat(ib, ep)
            );
          } else {
            cmp[ep] = compare_lower_better(value_mat(ia, ep), value_mat(ib, ep));
          }
        }

        for (int oo = 0; oo < n_orders; ++oo) {
          bool resolved = false;
          for (int r = 0; r < m; ++r) {
            const int ep = orders(oo, r) - 1;
            const int c = cmp[ep];
            if (c == 1) {
              local_wins[oo * m + r] += 1.0;
              resolved = true;
              break;
            }
            if (c == -1) {
              local_losses[oo * m + r] += 1.0;
              resolved = true;
              break;
            }
          }
          if (!resolved) local_ties[oo] += 1.0;
        }
      }
    }

#pragma omp critical
    {
      for (int oo = 0; oo < n_orders; ++oo) {
        ties[oo] += local_ties[oo];
        for (int r = 0; r < m; ++r) {
          wins(oo, r) += local_wins[oo * m + r];
          losses(oo, r) += local_losses[oo * m + r];
        }
      }
    }
  }

  return List::create(
    _["wins"] = wins,
    _["losses"] = losses,
    _["ties"] = ties,
    _["nA"] = static_cast<double>(nA),
    _["nB"] = static_cast<double>(nB),
    _["n_pairs"] = n_pairs
  );
}

// [[Rcpp::export]]
int openmp_max_threads_cpp() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}
