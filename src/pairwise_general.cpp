#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static double wr_ratio_safe_cpp(double num, double den) {
  if (NumericVector::is_na(num) || NumericVector::is_na(den)) return NA_REAL;
  if (den == 0.0 && num == 0.0) return 1.0;
  if (den == 0.0 && num > 0.0) return R_PosInf;
  if (num == 0.0 && den > 0.0) return 0.0;
  return num / den;
}

static double wr_abslog_safe_cpp(double x) {
  if (NumericVector::is_na(x)) return NA_REAL;
  if (!R_finite(x)) return R_PosInf;
  if (x <= 0.0) return R_PosInf;
  return std::fabs(std::log(x));
}

static bool wr_ok_num_cpp(double x) {
  return !NumericVector::is_na(x) && R_finite(x);
}

static int wr_count_until_cpp(const NumericVector& times,
                              const IntegerVector& starts,
                              const IntegerVector& lens,
                              int subject,
                              double t) {
  if (!R_finite(t)) return 0;

  int len = lens[subject];
  if (len <= 0) return 0;

  int start = starts[subject];
  int lo = 0;
  int hi = len;

  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;

    if (times[start + mid] <= t) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return lo;
}

// [[Rcpp::export]]
List cpp_counts_for_candidate_core(IntegerVector idx_treatment,
                                   IntegerVector idx_control,
                                   IntegerVector order_vec,
                                   NumericVector weights,
                                   NumericVector thresholds_by_endpoint,
                                   IntegerVector type_code,
                                   NumericMatrix time_mat,
                                   IntegerMatrix event_mat,
                                   NumericMatrix value_mat,
                                   NumericMatrix followup_mat,
                                   List recurrent_times,
                                   List recurrent_start,
                                   List recurrent_len) {
  int n_pairs = idx_treatment.size();
  int m = type_code.size();

  if (idx_control.size() != n_pairs) {
    stop("Treatment and control pair-index vectors must have the same length.");
  }

  NumericVector wins_by_rank(m);
  NumericVector losses_by_rank(m);
  int tie_count = 0;

  for (int i = 0; i < n_pairs; ++i) {
    int ia = idx_treatment[i] - 1;
    int ib = idx_control[i] - 1;
    bool resolved = false;

    for (int r = 0; r < order_vec.size(); ++r) {
      int ep = order_vec[r] - 1;
      int typ = type_code[ep];
      int cmp = 0;

      if (typ == 1) {
        double threshold = thresholds_by_endpoint[ep];

        if (NumericVector::is_na(threshold) || !R_finite(threshold)) {
          threshold = 0.0;
        }

        double ta = time_mat(ia, ep);
        double tb = time_mat(ib, ep);
        int ea = event_mat(ia, ep);
        int eb = event_mat(ib, ep);

        bool ta_ok = wr_ok_num_cpp(ta);
        bool tb_ok = wr_ok_num_cpp(tb);
        bool ea1 = (ea == 1);
        bool eb1 = (eb == 1);

        if (ta_ok && tb_ok) {
          if (ea1 && eb1) {
            if (ta - tb > threshold) cmp = 1;
            else if (tb - ta > threshold) cmp = -1;
          } else if (!ea1 && eb1) {
            if (ta - tb > threshold) cmp = 1;
          } else if (ea1 && !eb1) {
            if (tb - ta > threshold) cmp = -1;
          }
        }
      } else if (typ == 2 || typ == 3 || typ == 5) {
        double a = value_mat(ia, ep);
        double b = value_mat(ib, ep);

        if (wr_ok_num_cpp(a) && wr_ok_num_cpp(b)) {
          if (b - a > 0.0) cmp = 1;
          else if (a - b > 0.0) cmp = -1;
        }
      } else if (typ == 4 || typ == 6) {
        double a = value_mat(ia, ep);
        double b = value_mat(ib, ep);

        if (wr_ok_num_cpp(a) && wr_ok_num_cpp(b)) {
          if (a - b > 0.0) cmp = 1;
          else if (b - a > 0.0) cmp = -1;
        }
      } else if (typ == 7) {
        double fua = followup_mat(ia, ep);
        double fub = followup_mat(ib, ep);

        if (wr_ok_num_cpp(fua) && wr_ok_num_cpp(fub)) {
          double common_followup = fua < fub ? fua : fub;

          NumericVector times = recurrent_times[ep];
          IntegerVector starts = recurrent_start[ep];
          IntegerVector lens = recurrent_len[ep];

          int a_count = wr_count_until_cpp(
            times,
            starts,
            lens,
            ia,
            common_followup
          );

          int b_count = wr_count_until_cpp(
            times,
            starts,
            lens,
            ib,
            common_followup
          );

          if (b_count > a_count) cmp = 1;
          else if (a_count > b_count) cmp = -1;
        }
      }

      if (cmp == 1) {
        wins_by_rank[r] += 1.0;
        resolved = true;
        break;
      }

      if (cmp == -1) {
        losses_by_rank[r] += 1.0;
        resolved = true;
        break;
      }
    }

    if (!resolved) tie_count += 1;
  }

  double weighted_win = 0.0;
  double weighted_loss = 0.0;

  for (int r = 0; r < m; ++r) {
    double w = weights[r];

    if (NumericVector::is_na(w) || !R_finite(w)) {
      w = 0.0;
    }

    weighted_win += w * wins_by_rank[r];
    weighted_loss += w * losses_by_rank[r];
  }

  double wr = wr_ratio_safe_cpp(
    weighted_win,
    weighted_loss
  );

  double wo = wr_ratio_safe_cpp(
    weighted_win + 0.5 * tie_count,
    weighted_loss + 0.5 * tie_count
  );

  return List::create(
    _["WR_statistic"] = wr,
    _["WO_statistic"] = wo,
    _["WR_abslog"] = wr_abslog_safe_cpp(wr),
    _["WO_abslog"] = wr_abslog_safe_cpp(wo),
    _["weighted_win"] = weighted_win,
    _["weighted_loss"] = weighted_loss,
    _["tie_count"] = tie_count,
    _["tie_proportion"] = static_cast<double>(tie_count) / n_pairs,
    _["n_pairs"] = n_pairs,
    _["wins_by_rank"] = wins_by_rank,
    _["losses_by_rank"] = losses_by_rank
  );
}
