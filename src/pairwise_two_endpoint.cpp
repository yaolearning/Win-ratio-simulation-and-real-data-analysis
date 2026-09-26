#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static int wr_count_hosp_until_cpp(const NumericVector& hosp_times,
                                   const IntegerVector& hosp_start,
                                   const IntegerVector& hosp_len,
                                   int idx,
                                   double t) {
  int len = hosp_len[idx];
  if (len <= 0) return 0;

  int start = hosp_start[idx];
  int lo = 0;
  int hi = len;

  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (hosp_times[start + mid] <= t) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }

  return lo;
}

// [[Rcpp::export]]
List cpp_two_endpoint_base_counts(NumericVector futime,
                                  IntegerVector cnsr,
                                  IntegerVector arm,
                                  NumericVector hosp_times,
                                  IntegerVector hosp_start,
                                  IntegerVector hosp_len) {
  int n = futime.size();
  std::vector<int> treatment;
  std::vector<int> control;

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) treatment.push_back(i);
    if (arm[i] == 0) control.push_back(i);
  }

  int n_treatment = treatment.size();
  int n_control = control.size();
  double total_pairs = static_cast<double>(n_treatment) * static_cast<double>(n_control);

  if (total_pairs <= 0) {
    stop("Both treatment and control groups must contain at least one subject.");
  }

  double D1_win = 0.0;
  double D1_loss = 0.0;
  double H2_win = 0.0;
  double H2_loss = 0.0;
  double H1_win = 0.0;
  double H1_loss = 0.0;
  double D2_win = 0.0;
  double D2_loss = 0.0;

  for (int ii = 0; ii < n_treatment; ++ii) {
    int pt1 = treatment[ii];
    double fu1 = futime[pt1];
    int c1_original = cnsr[pt1];

    for (int jj = 0; jj < n_control; ++jj) {
      int pt2 = control[jj];
      double fu2 = futime[pt2];
      int c2_original = cnsr[pt2];

      double common_followup = fu1 < fu2 ? fu1 : fu2;

      int hosp1 = wr_count_hosp_until_cpp(
        hosp_times, hosp_start, hosp_len, pt1, common_followup
      );
      int hosp2 = wr_count_hosp_until_cpp(
        hosp_times, hosp_start, hosp_len, pt2, common_followup
      );

      int hosp_sign = 0;
      if (hosp2 > hosp1) hosp_sign = 1;
      else if (hosp2 < hosp1) hosp_sign = -1;

      int c1 = c1_original;
      int c2 = c2_original;

      if (fu1 < fu2) c2 = 0;
      if (fu2 < fu1) c1 = 0;

      int death_sign = 0;
      if (c1 == 0 && c2 == 1) death_sign = 1;
      else if (c1 == 1 && c2 == 0) death_sign = -1;

      if (death_sign > 0) {
        D1_win += 1.0;
      } else if (death_sign < 0) {
        D1_loss += 1.0;
      } else if (hosp_sign > 0) {
        H2_win += 1.0;
      } else if (hosp_sign < 0) {
        H2_loss += 1.0;
      }

      if (hosp_sign > 0) {
        H1_win += 1.0;
      } else if (hosp_sign < 0) {
        H1_loss += 1.0;
      } else if (death_sign > 0) {
        D2_win += 1.0;
      } else if (death_sign < 0) {
        D2_loss += 1.0;
      }
    }
  }

  return List::create(
    _["total_pairs"] = total_pairs,
    _["D1_win"] = D1_win,
    _["D1_loss"] = D1_loss,
    _["H2_win"] = H2_win,
    _["H2_loss"] = H2_loss,
    _["H1_win"] = H1_win,
    _["H1_loss"] = H1_loss,
    _["D2_win"] = D2_win,
    _["D2_loss"] = D2_loss
  );
}

// [[Rcpp::export]]
List cpp_two_endpoint_threshold_counts(NumericVector futime,
                                       IntegerVector cnsr,
                                       IntegerVector arm,
                                       NumericVector hosp_times,
                                       IntegerVector hosp_start,
                                       IntegerVector hosp_len,
                                       NumericVector t_grid) {
  int n = futime.size();
  std::vector<int> treatment;
  std::vector<int> control;

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) treatment.push_back(i);
    if (arm[i] == 0) control.push_back(i);
  }

  int n_treatment = treatment.size();
  int n_control = control.size();
  double total_pairs = static_cast<double>(n_treatment) * static_cast<double>(n_control);

  if (total_pairs <= 0) {
    stop("Both treatment and control groups must contain at least one subject.");
  }

  int K = t_grid.size();

  NumericVector df_first_win(K);
  NumericVector df_first_loss(K);
  NumericVector df_second_win(K);
  NumericVector df_second_loss(K);
  NumericVector hf_first_win(K);
  NumericVector hf_first_loss(K);
  NumericVector hf_second_win(K);
  NumericVector hf_second_loss(K);

  for (int ii = 0; ii < n_treatment; ++ii) {
    int pt1 = treatment[ii];
    double fu1 = futime[pt1];
    int c1_original = cnsr[pt1];

    for (int jj = 0; jj < n_control; ++jj) {
      int pt2 = control[jj];
      double fu2 = futime[pt2];
      int c2_original = cnsr[pt2];

      double common_followup = fu1 < fu2 ? fu1 : fu2;

      int hosp1 = wr_count_hosp_until_cpp(
        hosp_times, hosp_start, hosp_len, pt1, common_followup
      );
      int hosp2 = wr_count_hosp_until_cpp(
        hosp_times, hosp_start, hosp_len, pt2, common_followup
      );

      int hosp_sign = 0;
      if (hosp2 > hosp1) hosp_sign = 1;
      else if (hosp2 < hosp1) hosp_sign = -1;

      int c1 = c1_original;
      int c2 = c2_original;

      if (fu1 < fu2) c2 = 0;
      if (fu2 < fu1) c1 = 0;

      int death_sign = 0;
      if (c1 == 0 && c2 == 1) death_sign = 1;
      else if (c1 == 1 && c2 == 0) death_sign = -1;

      double abs_diff = std::fabs(fu1 - fu2);

      for (int kk = 0; kk < K; ++kk) {
        int death_threshold_sign = 0;

        if (death_sign != 0 && abs_diff > t_grid[kk]) {
          death_threshold_sign = death_sign;
        }

        if (death_threshold_sign > 0) {
          df_first_win[kk] += 1.0;
        } else if (death_threshold_sign < 0) {
          df_first_loss[kk] += 1.0;
        } else if (hosp_sign > 0) {
          df_second_win[kk] += 1.0;
        } else if (hosp_sign < 0) {
          df_second_loss[kk] += 1.0;
        }

        if (hosp_sign > 0) {
          hf_first_win[kk] += 1.0;
        } else if (hosp_sign < 0) {
          hf_first_loss[kk] += 1.0;
        } else if (death_threshold_sign > 0) {
          hf_second_win[kk] += 1.0;
        } else if (death_threshold_sign < 0) {
          hf_second_loss[kk] += 1.0;
        }
      }
    }
  }

  return List::create(
    _["total_pairs"] = total_pairs,
    _["t_grid"] = t_grid,
    _["df_first_win"] = df_first_win,
    _["df_first_loss"] = df_first_loss,
    _["df_second_win"] = df_second_win,
    _["df_second_loss"] = df_second_loss,
    _["hf_first_win"] = hf_first_win,
    _["hf_first_loss"] = hf_first_loss,
    _["hf_second_win"] = hf_second_win,
    _["hf_second_loss"] = hf_second_loss
  );
}
