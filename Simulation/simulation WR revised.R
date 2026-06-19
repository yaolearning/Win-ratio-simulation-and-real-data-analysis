# Packages and global settings
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_if_missing("survival")
install_if_missing("Rcpp")

suppressPackageStartupMessages({
  library(survival)
  library(Rcpp)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

NSIM <- 1000L
B_PERM <- 500L
MASTER_SEED <- 2026L
ALPHA <- 0.05
OUTDIR <- "simulation_WR_15scenarios_WR_and_LR_one_two_sided_abs_log_1000x500_STANDALONE"
CHECKPOINT_EVERY <- 25L
RESUME_IF_EXISTS <- FALSE
VERBOSE <- TRUE
SAVE_EXAMPLE_PLOTS <- FALSE

# Weighted p grids.
# p is the weight assigned to the first endpoint in the selected hierarchy.
# For death-first order: p weights death/survival; 1-p weights hospitalization.
# For hospitalization-first order: p weights hospitalization; 1-p weights death.
# No p = 0 value is used. Values below 0.5 are exploratory only.
P_GRID_EXPLORATORY <- seq(0.01, 0.49, by = 0.01)
P_GRID_LOW <- P_GRID_EXPLORATORY
P_GRID_PRIMARY <- seq(0.50, 1.00, by = 0.01)
P_GRID_ALL <- sort(unique(c(P_GRID_EXPLORATORY, P_GRID_PRIMARY)))

CLINICAL_T_MONTHS <- c(1, 3, 6, 12, 18, 24)
MAX_T_GRID_SIZE <- 10
EPS_WR <- 1e-8

# Small safe helpers
safe_mean <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- as.numeric(x)
  if (sum(!is.na(x)) <= 1) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

safe_median <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}

mode_string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

clean_order_label <- function(x) {
  x <- as.character(x)
  x <- ifelse(x == "death_first", "death first",
              ifelse(x == "hospitalization_first", "hospitalization first", x))
  x
}

col_or_na <- function(d, candidates) {
  for (nm in candidates) {
    if (nm %in% names(d)) return(d[[nm]])
  }
  rep(NA_real_, nrow(d))
}

wr_log <- function(wr) {
  wr <- as.numeric(wr)
  out <- rep(NA_real_, length(wr))
  ok <- is.finite(wr) & wr > 0
  out[ok] <- log(wr[ok])
  out
}

# Primary two-sided WR statistic. Because WR is a ratio with null value 1,
# abs(log(WR)) treats reciprocal effects (for example, 2 and 0.5) as equally
# far from the null.
wr_abs_log <- function(wr) {
  abs(wr_log(wr))
}

right_tail_perm_p <- function(perm.stat, obs.stat) {
  perm.stat <- as.numeric(perm.stat)
  obs.stat <- as.numeric(obs.stat)[1]
  ok <- is.finite(perm.stat)
  if (!is.finite(obs.stat) || sum(ok) == 0) return(NA_real_)
  (1 + sum(perm.stat[ok] >= obs.stat)) / (sum(ok) + 1)
}

# Retained only for directional diagnostic columns. It is not used to form the
# primary two-sided WR p-value, which is based on the right tail of abs(log WR).
left_tail_perm_p <- function(perm.stat, obs.stat) {
  perm.stat <- as.numeric(perm.stat)
  obs.stat <- as.numeric(obs.stat)[1]
  ok <- is.finite(perm.stat)
  if (!is.finite(obs.stat) || sum(ok) == 0) return(NA_real_)
  (1 + sum(perm.stat[ok] <= obs.stat)) / (sum(ok) + 1)
}

# The main two-sided test does not use 2 * min(tail p-values).
combine_two_tail_p <- function(p.upper, p.lower) {
  p.upper <- as.numeric(p.upper)
  p.lower <- as.numeric(p.lower)
  out <- rep(NA_real_, max(length(p.upper), length(p.lower)))
  p.upper <- rep(p.upper, length.out = length(out))
  p.lower <- rep(p.lower, length.out = length(out))
  ok <- is.finite(p.upper) | is.finite(p.lower)
  out[ok] <- pmin(1, 2 * pmin(p.upper[ok], p.lower[ok], na.rm = TRUE))
  out
}

two_tail_perm_p_signed <- function(perm.stat, obs.stat) {
  p.upper <- right_tail_perm_p(perm.stat, obs.stat)
  p.lower <- left_tail_perm_p(perm.stat, obs.stat)
  combine_two_tail_p(p.upper, p.lower)[1]
}

wr_perm_p_one_sided <- function(perm.wr, obs.wr) {
  right_tail_perm_p(perm.wr, obs.wr)
}

# Primary two-sided WR p-value: transform observed and permuted WR values in
# exactly the same way and compare abs(log(WR)) in the right tail.
wr_perm_p_two_sided <- function(perm.wr, obs.wr) {
  right_tail_perm_p(wr_abs_log(perm.wr), wr_abs_log(obs.wr))
}

# The following compatibility helpers are retained. 
# They are not used to calculate the
# primary abs(log WR) two-sided p-value.
choose_two_sided_direction <- function(p.upper, p.lower) {
  p.upper <- as.numeric(p.upper)
  p.lower <- as.numeric(p.lower)
  out <- rep(NA_character_, max(length(p.upper), length(p.lower)))
  p.upper <- rep(p.upper, length.out = length(out))
  p.lower <- rep(p.lower, length.out = length(out))
  for (i in seq_along(out)) {
    if (!is.finite(p.upper[i]) && !is.finite(p.lower[i])) {
      out[i] <- NA_character_
    } else if (!is.finite(p.upper[i])) {
      out[i] <- "lower_harm"
    } else if (!is.finite(p.lower[i])) {
      out[i] <- "upper_benefit"
    } else if (p.lower[i] < p.upper[i]) {
      out[i] <- "lower_harm"
    } else {
      out[i] <- "upper_benefit"
    }
  }
  out
}

choose_vector_by_direction <- function(upper.vec, lower.vec, direction, methods = names(upper.vec)) {
  upper.vec <- upper.vec[methods]
  lower.vec <- lower.vec[methods]
  direction <- direction[methods]
  out <- upper.vec
  use.lower <- !is.na(direction) & direction == "lower_harm"
  out[use.lower] <- lower.vec[use.lower]
  out
}

choose_matrix_by_direction <- function(upper.mat, lower.mat, direction, methods = colnames(upper.mat)) {
  upper.mat <- upper.mat[, methods, drop = FALSE]
  lower.mat <- lower.mat[, methods, drop = FALSE]
  direction <- direction[methods]
  out <- upper.mat
  for (nm in methods) {
    if (!is.na(direction[nm]) && direction[nm] == "lower_harm") {
      out[, nm] <- lower.mat[, nm]
    }
  }
  out
}

choose_rows_by_direction <- function(upper.rows, lower.rows, direction,
                                     p.upper = NULL, p.lower = NULL,
                                     methods = NULL) {
  if (is.null(methods)) methods <- upper.rows$method
  upper.rows <- upper.rows[match(methods, upper.rows$method), , drop = FALSE]
  lower.rows <- lower.rows[match(methods, lower.rows$method), , drop = FALSE]
  direction <- direction[methods]
  out <- upper.rows
  for (i in seq_along(methods)) {
    if (!is.na(direction[i]) && direction[i] == "lower_harm") {
      out[i, names(lower.rows)] <- lower.rows[i, , drop = FALSE]
    }
  }
  out$two_sided_tail_direction <- as.character(direction)
  if (!is.null(p.upper)) out$two_sided_upper_tail_p_value <- as.numeric(p.upper[methods])
  if (!is.null(p.lower)) out$two_sided_lower_tail_p_value <- as.numeric(p.lower[methods])
  rownames(out) <- NULL
  out
}

extract_nearest <- function(x, value) {
  which.min(abs(x - value))
}


rbind_fill_base <- function(...) {
  xs.raw <- list(...)
  
  # Accept data frames, lists of data frames, or mixed combinations.
  # This keeps all checkpoint/final outputs while avoiding list/data-frame
  # aggregation bugs.
  flatten_data_frames <- function(z) {
    out <- list()
    for (item in z) {
      if (is.null(item)) {
        next
      } else if (is.data.frame(item)) {
        out <- c(out, list(item))
      } else if (is.list(item)) {
        out <- c(out, flatten_data_frames(item))
      }
    }
    out
  }
  
  xs <- flatten_data_frames(xs.raw)
  xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0, logical(1))]
  if (length(xs) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs2 <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, xs2)
  rownames(out) <- NULL
  out
}

method_label_final <- function(method) {
  map <- c(
    ordinaryWR = "Traditional WR, death first, p = 0.50",
    traditionalWR_hosp_first = "Traditional WR, hospitalization first, p = 0.50",
    traditionalOrderWR = "Traditional order-selected WR, p = 0.50",
    maxWRp_primary = "Death-first weighted max WR, p >= 0.50",
    maxWRp_low = "Death-first weighted max WR, p < 0.50",
    maxWRp_full = "Death-first weighted max WR, full p grid",
    maxOrderWR_primary = "Maximum-order weighted WR, p >= 0.50",
    maxOrderWR_full = "Maximum-order weighted WR, full p grid",
    maxWRt = "Threshold max WR(t)",
    logrank_death = "Log-rank death",
    logrank_composite = "Log-rank composite"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out
}

scenario_display_label_final <- function(x) {
  map <- c(
    S01_equivalence_equal_arms = "S01 Equal",
    S02_similarity_small_benefit_both = "S02 Similar +",
    S03_similarity_small_harm_both = "S03 Similar -",
    S04_superiority_moderate_better_both = "S04 Sup mod",
    S05_superiority_strong_better_both = "S05 Sup strong",
    S06_inferiority_moderate_worse_both = "S06 Inf mod",
    S07_inferiority_strong_worse_both = "S07 Inf strong",
    S08_death_benefit_only = "S08 Death +",
    S09_hosp_benefit_only = "S09 Hosp +",
    S10_death_benefit_hosp_harm = "S10 Death + / Hosp -",
    S11_death_harm_hosp_benefit = "S11 Death - / Hosp +",
    S12_high_random_censoring_superiority = "S12 High censor",
    S13_low_censoring_long_followup_superiority = "S13 Long FU",
    S14_strong_hosp_benefit_only = "S14 Strong hosp +",
    S15_weak_death_benefit_strong_hosp_benefit = "S15 Weak death + strong hosp +"
  )
  out <- as.character(x)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}


# Data simulation
simulate.one.dataset <- function(N = 50,
                                 mort.rate.ctrl = -log(0.6),
                                 mort.rate.trt = -log(0.6) * 0.6,
                                 evt.rate.shape.param = 5,
                                 evt.rate.scale.param.ctr = 1,
                                 evt.rate.scale.param.trt = 1/2,
                                 max.followup = 1) {
  if (length(N) == 1) N <- rep(N, 2)
  N <- as.integer(N)
  n <- sum(N)
  
  table.output <- data.frame(
    SUBJID = seq_len(n),
    ARM = rep(0:1, N),
    FUTIME = rep(max.followup, n),
    CNSR = rep(0L, n),
    SURVTIME = c(stats::rexp(N[1], rate = mort.rate.ctrl),
                 stats::rexp(N[2], rate = mort.rate.trt)),
    CNSRTIME = rep(Inf, n),
    FREQHOSP = c(stats::rgamma(N[1], shape = evt.rate.shape.param, scale = evt.rate.scale.param.ctr),
                 stats::rgamma(N[2], shape = evt.rate.shape.param, scale = evt.rate.scale.param.trt)),
    NUMHOSP = rep(0L, n),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n)) {
    table.output$CNSR[i] <- as.integer(
      table.output$SURVTIME[i] < table.output$FUTIME[i] &&
        table.output$SURVTIME[i] < table.output$CNSRTIME[i]
    )
    table.output$FUTIME[i] <- min(table.output$FUTIME[i],
                                  table.output$SURVTIME[i],
                                  table.output$CNSRTIME[i])
  }
  
  hosp.times.list <- vector("list", n)
  for (i in seq_len(n)) {
    gaps <- numeric(0)
    cum.time <- 0
    followup.time <- table.output$FUTIME[i]
    hosp.rate <- table.output$FREQHOSP[i]
    
    while (is.finite(hosp.rate) && hosp.rate > 0) {
      new.hosp.time <- stats::rexp(1, rate = hosp.rate)
      cum.time <- cum.time + new.hosp.time
      if (cum.time < followup.time) {
        gaps <- c(gaps, new.hosp.time)
      } else {
        break
      }
    }
    hosp.times.list[[i]] <- gaps
    table.output$NUMHOSP[i] <- length(gaps)
  }
  
  list(table.output = table.output, hosp.times.list = hosp.times.list)
}

# Independent random right censoring wrapper.
apply.random.censoring <- function(ds, censor.rate = 0, seed = NULL) {
  if (is.null(censor.rate) || is.na(censor.rate) || censor.rate <= 0) return(ds)
  if (!is.null(seed)) set.seed(seed)
  
  tab <- ds$table.output
  n <- nrow(tab)
  ctime <- stats::rexp(n, rate = censor.rate)
  
  for (i in seq_len(n)) {
    if (ctime[i] < tab$FUTIME[i]) {
      tab$FUTIME[i] <- ctime[i]
      tab$CNSRTIME[i] <- ctime[i]
      tab$CNSR[i] <- 0L
      
      x <- ds$hosp.times.list[[i]]
      if (length(x) == 0 || all(is.na(x))) {
        ds$hosp.times.list[[i]] <- numeric(0)
        tab$NUMHOSP[i] <- 0L
      } else {
        abs.x <- cumsum(as.numeric(x[!is.na(x)]))
        keep <- abs.x <= tab$FUTIME[i]
        abs.keep <- abs.x[keep]
        if (length(abs.keep) == 0) {
          ds$hosp.times.list[[i]] <- numeric(0)
          tab$NUMHOSP[i] <- 0L
        } else {
          ds$hosp.times.list[[i]] <- diff(c(0, abs.keep))
          tab$NUMHOSP[i] <- length(abs.keep)
        }
      }
    }
  }
  
  ds$table.output <- tab
  ds
}

prepare.ds.fast <- function(ds) {
  n <- nrow(ds$table.output)
  hosp.abs.times.list <- vector("list", n)
  
  for (i in seq_len(n)) {
    x <- ds$hosp.times.list[[i]]
    if (length(x) == 0 || all(is.na(x))) {
      hosp.abs.times.list[[i]] <- numeric(0)
    } else {
      x <- as.numeric(x[!is.na(x)])
      hosp.abs.times.list[[i]] <- cumsum(x)
    }
  }
  
  lens <- as.integer(lengths(hosp.abs.times.list))
  starts <- if (n == 0) integer(0) else as.integer(cumsum(c(0L, lens[-n])))
  flat <- as.numeric(unlist(hosp.abs.times.list, use.names = FALSE))
  if (length(flat) == 0) flat <- numeric(0)
  
  ds$hosp.abs.times.list <- hosp.abs.times.list
  ds$hosp.flat <- flat
  ds$hosp.start <- starts
  ds$hosp.len <- lens
  ds
}


# Fast WR core in C++
Rcpp::sourceCpp(code = '
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

int count_hosp_until_cpp(const NumericVector& hosp_times,
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

double safe_wr_cpp(double win_score, double loss_score, double total_pairs, double eps) {
  return ((win_score / total_pairs) + eps) / ((loss_score / total_pairs) + eps);
}

// [[Rcpp::export]]
List fast_wr_core_revised_cpp(NumericVector futime,
                              IntegerVector cnsr,
                              IntegerVector arm,
                              NumericVector hosp_times,
                              IntegerVector hosp_start,
                              IntegerVector hosp_len,
                              NumericVector p_grid,
                              NumericVector t_grid,
                              double eps = 1e-8) {
  int n = futime.size();
  std::vector<int> trt;
  std::vector<int> ctrl;
  trt.reserve(n);
  ctrl.reserve(n);

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) trt.push_back(i);
    if (arm[i] == 0) ctrl.push_back(i);
  }

  int n_trt = trt.size();
  int n_ctrl = ctrl.size();
  double total_pairs = static_cast<double>(n_trt) * static_cast<double>(n_ctrl);
  if (total_pairs <= 0) stop("Need at least one treatment patient and one control patient.");

  int K = t_grid.size();
  int P = p_grid.size();

  double D1_win = 0.0, D1_loss = 0.0, H2_win = 0.0, H2_loss = 0.0;
  double H1_win = 0.0, H1_loss = 0.0, D2_win = 0.0, D2_loss = 0.0;
  double true_tie_pairs = 0.0;

  NumericVector t_wins(K), t_losses(K), pr_win_t(K), pr_loss_t(K), pr_tie_t(K), WRt(K);

  for (int ii = 0; ii < n_trt; ++ii) {
    int ti = trt[ii];
    for (int jj = 0; jj < n_ctrl; ++jj) {
      int cj = ctrl[jj];

      int death_sign = 0;
      double ft = futime[ti];
      double fc = futime[cj];
      int dt = cnsr[ti];
      int dc = cnsr[cj];

      if (dt == 1 && dc == 1) {
        if (ft > fc) death_sign = 1;
        else if (ft < fc) death_sign = -1;
      } else if (dt == 0 && dc == 1 && ft >= fc) {
        death_sign = 1;
      } else if (dt == 1 && dc == 0 && fc >= ft) {
        death_sign = -1;
      }

      double common_t = ft < fc ? ft : fc;
      int ht = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, ti, common_t);
      int hc = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, cj, common_t);
      int hosp_sign = 0;
      if (ht < hc) hosp_sign = 1;
      else if (ht > hc) hosp_sign = -1;

      if (death_sign > 0) D1_win += 1.0;
      else if (death_sign < 0) D1_loss += 1.0;
      else if (hosp_sign > 0) H2_win += 1.0;
      else if (hosp_sign < 0) H2_loss += 1.0;
      else true_tie_pairs += 1.0;

      if (hosp_sign > 0) H1_win += 1.0;
      else if (hosp_sign < 0) H1_loss += 1.0;
      else if (death_sign > 0) D2_win += 1.0;
      else if (death_sign < 0) D2_loss += 1.0;

      for (int kk = 0; kk < K; ++kk) {
        int sign_t = 0;
        if (death_sign != 0 && std::fabs(ft - fc) >= t_grid[kk]) {
          sign_t = death_sign;
        } else {
          sign_t = hosp_sign;
        }
        if (sign_t > 0) t_wins[kk] += 1.0;
        else if (sign_t < 0) t_losses[kk] += 1.0;
      }
    }
  }

  NumericVector WR_death_first(P), win_death_first(P), loss_death_first(P);
  NumericVector WR_hosp_first(P), win_hosp_first(P), loss_hosp_first(P);
  for (int pp = 0; pp < P; ++pp) {
    double p = p_grid[pp];
    win_death_first[pp] = p * D1_win + (1.0 - p) * H2_win;
    loss_death_first[pp] = p * D1_loss + (1.0 - p) * H2_loss;
    WR_death_first[pp] = safe_wr_cpp(win_death_first[pp], loss_death_first[pp], total_pairs, eps);

    win_hosp_first[pp] = p * H1_win + (1.0 - p) * D2_win;
    loss_hosp_first[pp] = p * H1_loss + (1.0 - p) * D2_loss;
    WR_hosp_first[pp] = safe_wr_cpp(win_hosp_first[pp], loss_hosp_first[pp], total_pairs, eps);
  }

  double ordinary_win_score = D1_win + H2_win;
  double ordinary_loss_score = D1_loss + H2_loss;
  double ordinaryWR = safe_wr_cpp(ordinary_win_score, ordinary_loss_score, total_pairs, eps);

  for (int kk = 0; kk < K; ++kk) {
    pr_win_t[kk] = t_wins[kk] / total_pairs;
    pr_loss_t[kk] = t_losses[kk] / total_pairs;
    pr_tie_t[kk] = 1.0 - pr_win_t[kk] - pr_loss_t[kk];
    WRt[kk] = safe_wr_cpp(t_wins[kk], t_losses[kk], total_pairs, eps);
  }

  return List::create(
    Named("ordinaryWR") = ordinaryWR,
    Named("ordinary_win_score") = ordinary_win_score,
    Named("ordinary_loss_score") = ordinary_loss_score,
    Named("total_pairs") = total_pairs,
    Named("D1_win") = D1_win,
    Named("D1_loss") = D1_loss,
    Named("H2_win") = H2_win,
    Named("H2_loss") = H2_loss,
    Named("H1_win") = H1_win,
    Named("H1_loss") = H1_loss,
    Named("D2_win") = D2_win,
    Named("D2_loss") = D2_loss,
    Named("true_tie_pairs") = true_tie_pairs,
    Named("p_grid") = p_grid,
    Named("WR_death_first") = WR_death_first,
    Named("win_death_first") = win_death_first,
    Named("loss_death_first") = loss_death_first,
    Named("WR_hosp_first") = WR_hosp_first,
    Named("win_hosp_first") = win_hosp_first,
    Named("loss_hosp_first") = loss_hosp_first,
    Named("t_grid") = t_grid,
    Named("WRt") = WRt,
    Named("t_wins") = t_wins,
    Named("t_losses") = t_losses,
    Named("pr_win_t") = pr_win_t,
    Named("pr_loss_t") = pr_loss_t,
    Named("pr_tie_t") = pr_tie_t
  );
}
')


# WR engine in R
max_from_curve <- function(curve, method, lower, upper) {
  ix <- which(curve$p >= lower & curve$p <= upper)
  if (length(ix) == 0) {
    return(data.frame(method = method, order = NA_character_, p = NA_real_,
                      WR = NA_real_, win.score = NA_real_, loss.score = NA_real_,
                      stringsAsFactors = FALSE))
  }
  k <- ix[which.max(curve$WR[ix])]
  data.frame(method = method,
             order = as.character(curve$order[k]),
             p = curve$p[k],
             WR = curve$WR[k],
             win.score = curve$win.score[k],
             loss.score = curve$loss.score[k],
             stringsAsFactors = FALSE)
}

fast.wr.engine.revised <- function(ds,
                                   p.grid = P_GRID_ALL,
                                   t.grid = CLINICAL_T_MONTHS / 12,
                                   eps = EPS_WR) {
  if (is.null(ds$hosp.flat)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  
  core <- fast_wr_core_revised_cpp(
    futime = as.numeric(tab$FUTIME),
    cnsr = as.integer(tab$CNSR),
    arm = as.integer(tab$ARM),
    hosp_times = as.numeric(ds$hosp.flat),
    hosp_start = as.integer(ds$hosp.start),
    hosp_len = as.integer(ds$hosp.len),
    p_grid = as.numeric(p.grid),
    t_grid = as.numeric(t.grid),
    eps = eps
  )
  
  total.pairs <- as.numeric(core$total_pairs)
  D1.win <- as.numeric(core$D1_win)
  D1.loss <- as.numeric(core$D1_loss)
  H2.win <- as.numeric(core$H2_win)
  H2.loss <- as.numeric(core$H2_loss)
  H1.win <- as.numeric(core$H1_win)
  H1.loss <- as.numeric(core$H1_loss)
  D2.win <- as.numeric(core$D2_win)
  D2.loss <- as.numeric(core$D2_loss)
  true.tie.pairs <- as.numeric(core$true_tie_pairs)
  
  weighted.death.first <- data.frame(
    order = "death_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_death_first),
    win.score = as.numeric(core$win_death_first),
    loss.score = as.numeric(core$loss_death_first),
    stringsAsFactors = FALSE
  )
  weighted.death.first$win.pairs <- ifelse(weighted.death.first$p >= 1, D1.win, D1.win + H2.win)
  weighted.death.first$loss.pairs <- ifelse(weighted.death.first$p >= 1, D1.loss, D1.loss + H2.loss)
  weighted.death.first$tie.count <- total.pairs - weighted.death.first$win.pairs - weighted.death.first$loss.pairs
  weighted.death.first$tie.pr <- weighted.death.first$tie.count / total.pairs
  
  weighted.hosp.first <- data.frame(
    order = "hospitalization_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_hosp_first),
    win.score = as.numeric(core$win_hosp_first),
    loss.score = as.numeric(core$loss_hosp_first),
    stringsAsFactors = FALSE
  )
  weighted.hosp.first$win.pairs <- ifelse(weighted.hosp.first$p >= 1, H1.win, H1.win + D2.win)
  weighted.hosp.first$loss.pairs <- ifelse(weighted.hosp.first$p >= 1, H1.loss, H1.loss + D2.loss)
  weighted.hosp.first$tie.count <- total.pairs - weighted.hosp.first$win.pairs - weighted.hosp.first$loss.pairs
  weighted.hosp.first$tie.pr <- weighted.hosp.first$tie.count / total.pairs
  
  weighted.all.orders <- rbind(weighted.death.first, weighted.hosp.first)
  
  threshold <- data.frame(
    t = as.numeric(core$t_grid),
    t.months = as.numeric(core$t_grid) * 12,
    WR = as.numeric(core$WRt),
    win.score = as.numeric(core$t_wins),
    loss.score = as.numeric(core$t_losses),
    win.pairs = as.numeric(core$t_wins),
    loss.pairs = as.numeric(core$t_losses),
    pr.win = as.numeric(core$pr_win_t),
    pr.loss = as.numeric(core$pr_loss_t),
    pr.tie = as.numeric(core$pr_tie_t),
    stringsAsFactors = FALSE
  )
  threshold$tie.count <- total.pairs - threshold$win.pairs - threshold$loss.pairs
  threshold$tie.pr <- threshold$tie.count / total.pairs
  
  max.weighted.primary <- max_from_curve(weighted.death.first, "maxWRp_primary", 0.50, 1.00)
  max.weighted.low <- max_from_curve(weighted.death.first, "maxWRp_low", min(p.grid, na.rm = TRUE), 0.499999)
  max.weighted.full <- max_from_curve(weighted.death.first, "maxWRp_full", min(p.grid, na.rm = TRUE), max(p.grid, na.rm = TRUE))
  max.order.primary <- max_from_curve(weighted.all.orders, "maxOrderWR_primary", 0.50, 1.00)
  max.order.full <- max_from_curve(weighted.all.orders, "maxOrderWR_full", min(p.grid, na.rm = TRUE), max(p.grid, na.rm = TRUE))
  max.traditional.order <- max_from_curve(weighted.all.orders[abs(weighted.all.orders$p - 0.5) < 1e-9, ], "traditionalOrderWR", 0.50, 0.50)
  
  k.t <- if (nrow(threshold) == 0) NA_integer_ else which.max(threshold$WR)
  max.threshold <- if (is.na(k.t)) {
    data.frame(method = "maxWRt", t = NA_real_, t.months = NA_real_, WR = NA_real_,
               win.score = NA_real_, loss.score = NA_real_, stringsAsFactors = FALSE)
  } else {
    cbind(data.frame(method = "maxWRt", stringsAsFactors = FALSE), threshold[k.t, , drop = FALSE])
  }
  
  list(
    ordinaryWR = as.numeric(core$ordinaryWR),
    ordinary.win.score = as.numeric(core$ordinary_win_score),
    ordinary.loss.score = as.numeric(core$ordinary_loss_score),
    total.pairs = total.pairs,
    true.hierarchical.tie.count = true.tie.pairs,
    true.hierarchical.tie.pr = true.tie.pairs / total.pairs,
    weighted.death.first = weighted.death.first,
    weighted.hosp.first = weighted.hosp.first,
    weighted.all.orders = weighted.all.orders,
    threshold = threshold,
    max.weighted.primary = max.weighted.primary,
    max.weighted.low = max.weighted.low,
    max.weighted.full = max.weighted.full,
    max.traditional.order = max.traditional.order,
    max.order.primary = max.order.primary,
    max.order.full = max.order.full,
    max.threshold = max.threshold
  )
}

wr_at_order_p <- function(engine.out, order = "death_first", p = 0.50) {
  d <- if (order == "hospitalization_first") engine.out$weighted.hosp.first else engine.out$weighted.death.first
  if (is.null(d) || nrow(d) == 0) return(NA_real_)
  d$WR[which.min(abs(d$p - p))]
}

threshold_at_t <- function(engine.out, t) {
  d <- engine.out$threshold
  if (is.null(d) || nrow(d) == 0 || !is.finite(t)) return(NA_real_)
  d$WR[which.min(abs(d$t - t))]
}


# Composite endpoint and log-rank helpers
first_hosp_time_vec <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  sapply(ds$hosp.abs.times.list, function(x) if (length(x) > 0) x[1] else Inf)
}

compute.composite.endpoint.data <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  first.hosp <- first_hosp_time_vec(ds)
  death.time <- ifelse(tab$CNSR == 1, tab$FUTIME, Inf)
  comp.time <- pmin(death.time, first.hosp)
  comp.event <- is.finite(comp.time) & comp.time <= tab$FUTIME
  
  comp.type <- rep("censored", nrow(tab))
  comp.type[comp.event & death.time <= first.hosp] <- "death"
  comp.type[comp.event & first.hosp < death.time] <- "hospitalization"
  
  data.frame(
    SUBJID = tab$SUBJID,
    ARM = tab$ARM,
    time = ifelse(comp.event, comp.time, tab$FUTIME),
    event = as.integer(comp.event),
    event.type = comp.type,
    first.hosp.time = first.hosp,
    death.time = death.time,
    stringsAsFactors = FALSE
  )
}

composite.statistics <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  comp <- compute.composite.endpoint.data(ds)
  
  dat <- data.frame(
    ARM = tab$ARM,
    FUTIME = tab$FUTIME,
    CNSR = tab$CNSR,
    NUMHOSP = tab$NUMHOSP,
    FREQHOSP = tab$FREQHOSP,
    comp.time = comp$time,
    comp.event = comp$event,
    comp.type = comp$event.type,
    stringsAsFactors = FALSE
  )
  
  out <- do.call(rbind, lapply(split(dat, dat$ARM), function(d) {
    data.frame(
      ARM = unique(d$ARM),
      arm.label = ifelse(unique(d$ARM) == 0, "control", "treatment"),
      n = nrow(d),
      death.events = sum(d$CNSR == 1, na.rm = TRUE),
      death.event.rate = mean(d$CNSR == 1, na.rm = TRUE),
      subjects.with.hosp = sum(d$NUMHOSP > 0, na.rm = TRUE),
      hosp.subject.rate = mean(d$NUMHOSP > 0, na.rm = TRUE),
      total.hosp.events = sum(d$NUMHOSP, na.rm = TRUE),
      mean.num.hosp = mean(d$NUMHOSP, na.rm = TRUE),
      median.num.hosp = stats::median(d$NUMHOSP, na.rm = TRUE),
      mean.followup = mean(d$FUTIME, na.rm = TRUE),
      composite.events = sum(d$comp.event == 1, na.rm = TRUE),
      composite.event.rate = mean(d$comp.event == 1, na.rm = TRUE),
      composite.death.events = sum(d$comp.type == "death", na.rm = TRUE),
      composite.hosp.events = sum(d$comp.type == "hospitalization", na.rm = TRUE),
      mean.composite.time = mean(d$comp.time, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

# Threshold-grid helper
choose.threshold.grid.primary <- function(ds,
                                          clinical.months = CLINICAL_T_MONTHS,
                                          max.K = MAX_T_GRID_SIZE) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  max.fu <- max(tab$FUTIME, na.rm = TRUE)
  
  clinical <- as.numeric(clinical.months) / 12
  clinical <- clinical[is.finite(clinical) & clinical > 0 & clinical <= max.fu]
  
  trt <- tab$FUTIME[tab$ARM == 1]
  ctrl <- tab$FUTIME[tab$ARM == 0]
  empirical <- numeric(0)
  if (length(trt) > 0 && length(ctrl) > 0) {
    diffs <- as.vector(abs(outer(trt, ctrl, "-")))
    diffs <- diffs[is.finite(diffs) & diffs > 0 & diffs <= max.fu]
    if (length(diffs) > 0) {
      empirical <- as.numeric(stats::quantile(diffs, probs = c(0.25, 0.50, 0.75, 0.90),
                                              na.rm = TRUE, names = FALSE))
    }
  }
  
  candidates <- sort(unique(round(c(clinical, empirical), 8)))
  candidates <- candidates[is.finite(candidates) & candidates > 0 & candidates <= max.fu]
  if (length(candidates) == 0) candidates <- min(max.fu, 1 / 12)
  if (length(candidates) > max.K) {
    idx <- unique(round(seq(1, length(candidates), length.out = max.K)))
    candidates <- candidates[idx]
  }
  
  threshold.table <- data.frame(
    t = candidates,
    t_months = candidates * 12,
    stringsAsFactors = FALSE
  )
  diagnostics <- data.frame(
    max_followup = max.fu,
    n_candidates = length(candidates),
    candidate_months = paste(sprintf("%.2f", candidates * 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
  source.table <- data.frame(
    source = c(rep("clinical", length(clinical)), rep("empirical", length(empirical))),
    t = c(clinical, empirical),
    t_months = c(clinical, empirical) * 12,
    stringsAsFactors = FALSE
  )
  list(t.grid = candidates, threshold.table = threshold.table,
       diagnostics = diagnostics, source.table = source.table)
}


#15-scenario 
build.final.15.scenario.grid <- function() {
  base.mort <- -log(0.6)
  scenario_id <- c(
    "S01_equivalence_equal_arms",
    "S02_similarity_small_benefit_both",
    "S03_similarity_small_harm_both",
    "S04_superiority_moderate_better_both",
    "S05_superiority_strong_better_both",
    "S06_inferiority_moderate_worse_both",
    "S07_inferiority_strong_worse_both",
    "S08_death_benefit_only",
    "S09_hosp_benefit_only",
    "S10_death_benefit_hosp_harm",
    "S11_death_harm_hosp_benefit",
    "S12_high_random_censoring_superiority",
    "S13_low_censoring_long_followup_superiority",
    "S14_strong_hosp_benefit_only",
    "S15_weak_death_benefit_strong_hosp_benefit"
  )
  HR <- c(1.00, 0.90, 1.10, 0.75, 0.60, 1.25, 1.50,
          0.60, 1.00, 0.60, 1.50, 0.60, 0.60, 1.00, 0.90)
  hosp.scale.trt <- c(1.00, 0.90, 1.10, 0.75, 0.50, 1.25, 1.50,
                      1.00, 0.50, 1.50, 0.50, 0.50, 0.50, 0.25, 0.25)
  max.FU <- c(rep(1, 12), 3, 1, 1)
  censor.rate <- c(rep(0, 11), 1.50, 0.05, 0, 0)
  
  scenario_type <- c(
    "null / equivalence",
    "similarity / small benefit",
    "similarity / small harm",
    "superiority / moderate benefit",
    "superiority / strong benefit",
    "inferiority / moderate harm",
    "inferiority / strong harm",
    "death benefit only",
    "hospitalization benefit only",
    "discordant: death benefit, hospitalization harm",
    "discordant: death harm, hospitalization benefit",
    "sensitivity: high random censoring",
    "sensitivity: longer follow-up with light censoring",
    "targeted non-terminal benefit",
    "targeted weak death plus strong non-terminal benefit"
  )
  description <- c(
    "No treatment effect on death or hospitalization.",
    "Small benefit on both death and hospitalization.",
    "Small harm on both death and hospitalization.",
    "Moderate treatment benefit on both endpoints.",
    "Strong treatment benefit on both endpoints.",
    "Moderate treatment harm on both endpoints.",
    "Strong treatment harm on both endpoints.",
    "Treatment improves death only; hospitalization burden is unchanged.",
    "Treatment improves hospitalization only; death is unchanged.",
    "Treatment improves death but worsens hospitalization burden.",
    "Treatment worsens death but improves hospitalization burden.",
    "Same strong benefit pattern as S05 with high independent exponential censoring.",
    "Same strong benefit pattern as S05 with longer follow-up and light censoring.",
    "Treatment has no death effect but strongly reduces hospitalization; intended to favor order-change and maximum-order WR.",
    "Treatment has weak death benefit and strong hospitalization benefit; intended to favor weighted and order-adaptive WR."
  )
  
  out <- data.frame(
    scenario_index = seq_along(scenario_id),
    scenario_id = scenario_id,
    scenario_label = scenario_display_label_final(scenario_id),
    scenario_type = scenario_type,
    description = description,
    N0 = rep(50, length(scenario_id)),
    N1 = rep(50, length(scenario_id)),
    mort.rate.ctrl = rep(base.mort, length(scenario_id)),
    HR = HR,
    evt.rate.shape.param = rep(5, length(scenario_id)),
    evt.rate.scale.param.ctr = rep(1.0, length(scenario_id)),
    evt.rate.scale.param.trt = hosp.scale.trt,
    max.FU = max.FU,
    censor.rate = censor.rate,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

# WR curve selection helpers
filter_p_range <- function(d, lower = -Inf, upper = Inf) {
  if (is.null(d) || nrow(d) == 0 || !("p" %in% names(d))) return(d[FALSE, , drop = FALSE])
  d[!is.na(d$p) & d$p >= lower & d$p <= upper, , drop = FALSE]
}

max_score_row <- function(d, side = c("one", "two", "two_upper", "two_lower")) {
  side <- match.arg(side)
  if (is.null(d) || nrow(d) == 0 || !("WR" %in% names(d))) {
    return(data.frame())
  }
  
  z <- wr_log(d$WR)
  score <- switch(
    side,
    one = as.numeric(d$WR),
    two = abs(z),
    two_upper = z,
    two_lower = -z
  )
  
  score[!is.finite(score)] <- -Inf
  if (length(score) == 0 || all(score == -Inf)) return(d[FALSE, , drop = FALSE])
  d[which.max(score), , drop = FALSE]
}

row_at_order_p <- function(engine.out, order, p) {
  d <- if (order == "hospitalization_first") {
    engine.out$weighted.hosp.first
  } else {
    engine.out$weighted.death.first
  }
  if (is.null(d) || nrow(d) == 0) return(data.frame())
  idx <- which.min(abs(d$p - p))
  d[idx, , drop = FALSE]
}

row_at_threshold_t <- function(engine.out, t) {
  d <- engine.out$threshold
  if (is.null(d) || nrow(d) == 0) return(data.frame())
  idx <- which.min(abs(d$t - t))
  d[idx, , drop = FALSE]
}

selected_row_to_output <- function(method, d,
                                   side = c("one", "two", "two_upper", "two_lower")) {
  side <- match.arg(side)
  
  if (is.null(d) || nrow(d) == 0) {
    return(data.frame(
      method = method,
      selected_order = NA_character_,
      selected_p = NA_real_,
      selected_t = NA_real_,
      selected_t_months = NA_real_,
      selected_WR = NA_real_,
      selected_statistic = NA_real_,
      selected_log_WR = NA_real_,
      selected_abs_log_WR = NA_real_,
      selected_direction = NA_character_,
      selected_tie_count = NA_real_,
      selected_tie_pr = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  wr <- as.numeric(d$WR[1])
  log.wr <- wr_log(wr)
  abs.log.wr <- abs(log.wr)
  stat <- switch(
    side,
    one = wr,
    two = abs.log.wr,
    two_upper = log.wr,
    two_lower = -log.wr
  )
  direction <- if (!is.finite(log.wr)) {
    NA_character_
  } else if (log.wr >= 0) {
    "upper_benefit"
  } else {
    "lower_harm"
  }
  
  selected_order <- if ("order" %in% names(d)) as.character(d$order[1]) else NA_character_
  selected_p <- if ("p" %in% names(d)) as.numeric(d$p[1]) else NA_real_
  selected_t <- if ("t" %in% names(d)) as.numeric(d$t[1]) else NA_real_
  selected_t_months <- if ("t.months" %in% names(d)) {
    as.numeric(d$t.months[1])
  } else if (is.finite(selected_t)) {
    selected_t * 12
  } else {
    NA_real_
  }
  
  tie.count <- col_or_na(d, c("tie.count", "tie_count"))[1]
  tie.pr <- col_or_na(d, c("tie.pr", "tie_proportion", "pr.tie"))[1]
  
  data.frame(
    method = method,
    selected_order = selected_order,
    selected_p = selected_p,
    selected_t = selected_t,
    selected_t_months = selected_t_months,
    selected_WR = wr,
    selected_statistic = as.numeric(stat),
    selected_log_WR = as.numeric(log.wr),
    selected_abs_log_WR = as.numeric(abs.log.wr),
    selected_direction = direction,
    selected_tie_count = as.numeric(tie.count),
    selected_tie_pr = as.numeric(tie.pr),
    stringsAsFactors = FALSE
  )
}

selected_wr_rows <- function(engine.out,
                             side = c("one", "two", "two_upper", "two_lower")) {
  side <- match.arg(side)
  
  d05 <- row_at_order_p(engine.out, "death_first", 0.50)
  h05 <- row_at_order_p(engine.out, "hospitalization_first", 0.50)
  
  p.primary <- filter_p_range(engine.out$weighted.death.first, 0.50, 1.00)
  p.low <- filter_p_range(
    engine.out$weighted.death.first,
    min(engine.out$weighted.death.first$p, na.rm = TRUE),
    0.499999
  )
  p.full <- engine.out$weighted.death.first
  
  all.primary <- rbind(
    filter_p_range(engine.out$weighted.death.first, 0.50, 1.00),
    filter_p_range(engine.out$weighted.hosp.first, 0.50, 1.00)
  )
  all.full <- rbind(engine.out$weighted.death.first, engine.out$weighted.hosp.first)
  
  out <- rbind(
    selected_row_to_output("ordinaryWR", d05, side),
    selected_row_to_output("traditionalWR_hosp_first", h05, side),
    selected_row_to_output("traditionalOrderWR", max_score_row(rbind(d05, h05), side), side),
    selected_row_to_output("maxWRp_primary", max_score_row(p.primary, side), side),
    selected_row_to_output("maxWRp_low", max_score_row(p.low, side), side),
    selected_row_to_output("maxWRp_full", max_score_row(p.full, side), side),
    selected_row_to_output("maxOrderWR_primary", max_score_row(all.primary, side), side),
    selected_row_to_output("maxOrderWR_full", max_score_row(all.full, side), side),
    selected_row_to_output("maxWRt", max_score_row(engine.out$threshold, side), side)
  )
  rownames(out) <- NULL
  out
}

selected_stats_vector <- function(sel.rows) {
  out <- as.numeric(sel.rows$selected_statistic)
  names(out) <- sel.rows$method
  out
}

selected_wr_vector <- function(sel.rows) {
  out <- as.numeric(sel.rows$selected_WR)
  names(out) <- sel.rows$method
  out
}

selected_log_vector <- function(sel.rows) {
  out <- as.numeric(sel.rows$selected_log_WR)
  names(out) <- sel.rows$method
  out
}

fixed_eval_row <- function(engine.out, selected.row,
                           side = c("one", "two", "two_upper", "two_lower")) {
  side <- match.arg(side)
  method <- as.character(selected.row$method[1])
  
  if (method == "maxWRt" || is.finite(selected.row$selected_t[1])) {
    d <- row_at_threshold_t(engine.out, selected.row$selected_t[1])
  } else {
    ord <- as.character(selected.row$selected_order[1])
    p <- as.numeric(selected.row$selected_p[1])
    if (is.na(ord) || !nzchar(ord)) ord <- "death_first"
    if (!is.finite(p)) p <- 0.50
    d <- row_at_order_p(engine.out, ord, p)
  }
  
  selected_row_to_output(method, d, side)
}

fixed_eval_rows <- function(engine.out, selected.rows,
                            side = c("one", "two", "two_upper", "two_lower")) {
  side <- match.arg(side)
  out <- do.call(rbind, lapply(seq_len(nrow(selected.rows)), function(i) {
    fixed_eval_row(engine.out, selected.rows[i, , drop = FALSE], side = side)
  }))
  rownames(out) <- NULL
  out
}


# WR permutation test: output one-sided and two-sided p-values.
# One-sided WR statistic:
#   T = WR or max WR
# Two-sided WR statistic:
#   T = abs(log(WR)) for fixed candidates
#   T = max abs(log(WR)) for adaptive candidates

perm.test.revised <- function(ds,
                              B = B_PERM,
                              seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                              p.grid = if (exists("P_GRID_ALL")) P_GRID_ALL else seq(0.01, 1, by = 0.01),
                              t.grid = if (exists("CLINICAL_T_MONTHS")) CLINICAL_T_MONTHS / 12 else c(1, 3, 6, 12, 18, 24) / 12,
                              verbose = FALSE) {
  ds <- prepare.ds.fast(ds)
  set.seed(seed)
  
  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)
  
  # One-sided selection maximizes raw WR.
  # Two-sided selection maximizes abs(log(WR)).
  selected.max.one <- selected_wr_rows(obs, side = "one")
  selected.max.two <- selected_wr_rows(obs, side = "two")
  
  max.names <- selected.max.one$method
  selected.max.one <- selected.max.one[match(max.names, selected.max.one$method), , drop = FALSE]
  selected.max.two <- selected.max.two[match(max.names, selected.max.two$method), , drop = FALSE]
  
  T.obs.max.one <- selected_stats_vector(selected.max.one)
  T.obs.max.two <- selected_stats_vector(selected.max.two)
  T.obs.max.two.log <- selected_log_vector(selected.max.two)
  
  # Fixed-selected permutation p-values use the parameter selected from the
  # observed trial, then evaluate that same parameter in every permutation.
  selected.fixed.one <- selected.max.one
  selected.fixed.two <- selected.max.two
  
  T.obs.fixed.one.rows <- fixed_eval_rows(obs, selected.fixed.one, side = "one")
  T.obs.fixed.two.rows <- fixed_eval_rows(obs, selected.fixed.two, side = "two")
  
  T.obs.fixed.one <- selected_stats_vector(T.obs.fixed.one.rows)
  T.obs.fixed.two <- selected_stats_vector(T.obs.fixed.two.rows)
  T.obs.fixed.two.log <- selected_log_vector(T.obs.fixed.two.rows)
  
  perm.max.one <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                         dimnames = list(NULL, max.names))
  perm.max.two <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                         dimnames = list(NULL, max.names))
  perm.max.two.log <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                             dimnames = list(NULL, max.names))
  
  perm.fixed.one <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                           dimnames = list(NULL, max.names))
  perm.fixed.two <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                           dimnames = list(NULL, max.names))
  perm.fixed.two.log <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                               dimnames = list(NULL, max.names))
  
  perm.tie.max.one.count <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                   dimnames = list(NULL, max.names))
  perm.tie.max.one.pr <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                dimnames = list(NULL, max.names))
  perm.tie.max.two.count <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                   dimnames = list(NULL, max.names))
  perm.tie.max.two.pr <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                dimnames = list(NULL, max.names))
  
  perm.fixed.tie.one.count <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                     dimnames = list(NULL, max.names))
  perm.fixed.tie.one.pr <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                  dimnames = list(NULL, max.names))
  perm.fixed.tie.two.count <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                     dimnames = list(NULL, max.names))
  perm.fixed.tie.two.pr <- matrix(NA_real_, nrow = B, ncol = length(max.names),
                                  dimnames = list(NULL, max.names))
  
  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))
  
  
  # inside the 500-permutation loop is slower.
  selected.perm.one.list <- vector("list", B)
  selected.perm.two.list <- vector("list", B)
  
  perm.seeds <- seed + seq_len(B) * 1009L
  
  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b == B || b %% 50 == 0)) {
      cat("  Permutation", b, "of", B, "\n")
    }
    
    set.seed(perm.seeds[b])
    ds.b <- ds
    ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)
    
    out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)
    
    perm.pointwise.death.first[b, ] <- out.b$weighted.death.first$WR
    perm.pointwise.hosp.first[b, ] <- out.b$weighted.hosp.first$WR
    perm.pointwise.threshold[b, ] <- out.b$threshold$WR
    
    # Repeat the corresponding adaptive selection rule inside every
    # treatment-label permutation.
    sel.b.one <- selected_wr_rows(out.b, side = "one")
    sel.b.two <- selected_wr_rows(out.b, side = "two")
    sel.b.one <- sel.b.one[match(max.names, sel.b.one$method), , drop = FALSE]
    sel.b.two <- sel.b.two[match(max.names, sel.b.two$method), , drop = FALSE]
    
    perm.max.one[b, ] <- selected_stats_vector(sel.b.one)[max.names]
    perm.max.two[b, ] <- selected_stats_vector(sel.b.two)[max.names]
    perm.max.two.log[b, ] <- selected_log_vector(sel.b.two)[max.names]
    
    perm.tie.max.one.count[b, ] <- sel.b.one$selected_tie_count
    perm.tie.max.one.pr[b, ] <- sel.b.one$selected_tie_pr
    perm.tie.max.two.count[b, ] <- sel.b.two$selected_tie_count
    perm.tie.max.two.pr[b, ] <- sel.b.two$selected_tie_pr
    
    fix.b.one <- fixed_eval_rows(out.b, selected.fixed.one, side = "one")
    fix.b.two <- fixed_eval_rows(out.b, selected.fixed.two, side = "two")
    fix.b.one <- fix.b.one[match(max.names, fix.b.one$method), , drop = FALSE]
    fix.b.two <- fix.b.two[match(max.names, fix.b.two$method), , drop = FALSE]
    
    perm.fixed.one[b, ] <- selected_stats_vector(fix.b.one)[max.names]
    perm.fixed.two[b, ] <- selected_stats_vector(fix.b.two)[max.names]
    perm.fixed.two.log[b, ] <- selected_log_vector(fix.b.two)[max.names]
    
    perm.fixed.tie.one.count[b, ] <- fix.b.one$selected_tie_count
    perm.fixed.tie.one.pr[b, ] <- fix.b.one$selected_tie_pr
    perm.fixed.tie.two.count[b, ] <- fix.b.two$selected_tie_count
    perm.fixed.tie.two.pr[b, ] <- fix.b.two$selected_tie_pr
    
    tmp.one <- sel.b.one
    tmp.two <- sel.b.two
    tmp.one$b <- b
    tmp.two$b <- b
    tmp.two$two_sided_tail_direction <- tmp.two$selected_direction
    selected.perm.one.list[[b]] <- tmp.one
    selected.perm.two.list[[b]] <- tmp.two
  }
  
  selected.perm.one <- rbind_fill_base(selected.perm.one.list)
  selected.perm.two <- rbind_fill_base(selected.perm.two.list)
  
  # Primary adaptive permutation p-values.
  p.value.max.one <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.max.one[, nm], T.obs.max.one[nm])
  })
  p.value.max.two <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.max.two[, nm], T.obs.max.two[nm])
  })
  
  # Primary fixed-selected permutation p-values.
  p.value.fixed.one <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.fixed.one[, nm], T.obs.fixed.one[nm])
  })
  p.value.fixed.two <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.fixed.two[, nm], T.obs.fixed.two[nm])
  })
  
  # Directional diagnostic tails at the abs(log WR)-selected candidate. These
  # columns are retained for output compatibility but do not define the primary
  # two-sided p-value above.
  p.value.max.two.upper <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.max.two.log[, nm], T.obs.max.two.log[nm])
  })
  p.value.max.two.lower <- sapply(max.names, function(nm) {
    left_tail_perm_p(perm.max.two.log[, nm], T.obs.max.two.log[nm])
  })
  p.value.fixed.two.upper <- sapply(max.names, function(nm) {
    right_tail_perm_p(perm.fixed.two.log[, nm], T.obs.fixed.two.log[nm])
  })
  p.value.fixed.two.lower <- sapply(max.names, function(nm) {
    left_tail_perm_p(perm.fixed.two.log[, nm], T.obs.fixed.two.log[nm])
  })
  
  two.tail.direction.max <- ifelse(
    is.finite(T.obs.max.two.log),
    ifelse(T.obs.max.two.log >= 0, "upper_benefit", "lower_harm"),
    NA_character_
  )
  names(two.tail.direction.max) <- max.names
  
  two.tail.direction.fixed <- ifelse(
    is.finite(T.obs.fixed.two.log),
    ifelse(T.obs.fixed.two.log >= 0, "upper_benefit", "lower_harm"),
    NA_character_
  )
  names(two.tail.direction.fixed) <- max.names
  
  selected.max.two$two_sided_tail_direction <- as.character(two.tail.direction.max[max.names])
  selected.max.two$two_sided_upper_tail_p_value <- as.numeric(p.value.max.two.upper[max.names])
  selected.max.two$two_sided_lower_tail_p_value <- as.numeric(p.value.max.two.lower[max.names])
  
  selected.fixed.two$two_sided_tail_direction <- as.character(two.tail.direction.fixed[max.names])
  selected.fixed.two$two_sided_upper_tail_p_value <- as.numeric(p.value.fixed.two.upper[max.names])
  selected.fixed.two$two_sided_lower_tail_p_value <- as.numeric(p.value.fixed.two.lower[max.names])
  
  # Pointwise one-sided and two-sided permutation p-values.
  p.value.pointwise.death.first.one <- sapply(seq_along(p.grid), function(k) {
    wr_perm_p_one_sided(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k])
  })
  p.value.pointwise.death.first.two <- sapply(seq_along(p.grid), function(k) {
    wr_perm_p_two_sided(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k])
  })
  p.value.pointwise.death.first.two.upper <- sapply(seq_along(p.grid), function(k) {
    right_tail_perm_p(wr_log(perm.pointwise.death.first[, k]),
                      wr_log(obs$weighted.death.first$WR[k]))
  })
  p.value.pointwise.death.first.two.lower <- sapply(seq_along(p.grid), function(k) {
    left_tail_perm_p(wr_log(perm.pointwise.death.first[, k]),
                     wr_log(obs$weighted.death.first$WR[k]))
  })
  
  p.value.pointwise.hosp.first.one <- sapply(seq_along(p.grid), function(k) {
    wr_perm_p_one_sided(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k])
  })
  p.value.pointwise.hosp.first.two <- sapply(seq_along(p.grid), function(k) {
    wr_perm_p_two_sided(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k])
  })
  p.value.pointwise.hosp.first.two.upper <- sapply(seq_along(p.grid), function(k) {
    right_tail_perm_p(wr_log(perm.pointwise.hosp.first[, k]),
                      wr_log(obs$weighted.hosp.first$WR[k]))
  })
  p.value.pointwise.hosp.first.two.lower <- sapply(seq_along(p.grid), function(k) {
    left_tail_perm_p(wr_log(perm.pointwise.hosp.first[, k]),
                     wr_log(obs$weighted.hosp.first$WR[k]))
  })
  
  p.value.pointwise.threshold.one <- sapply(seq_along(t.grid), function(k) {
    wr_perm_p_one_sided(perm.pointwise.threshold[, k], obs$threshold$WR[k])
  })
  p.value.pointwise.threshold.two <- sapply(seq_along(t.grid), function(k) {
    wr_perm_p_two_sided(perm.pointwise.threshold[, k], obs$threshold$WR[k])
  })
  p.value.pointwise.threshold.two.upper <- sapply(seq_along(t.grid), function(k) {
    right_tail_perm_p(wr_log(perm.pointwise.threshold[, k]),
                      wr_log(obs$threshold$WR[k]))
  })
  p.value.pointwise.threshold.two.lower <- sapply(seq_along(t.grid), function(k) {
    left_tail_perm_p(wr_log(perm.pointwise.threshold[, k]),
                     wr_log(obs$threshold$WR[k]))
  })
  
  # Backwards-compatible aliases for fields that existed in the signed-log
  # They now refer to the signed log of the same candidate selected by
  # max abs(log WR); the primary two-sided test uses T.obs/T.perm.max.two.
  selected.max.two.upper <- selected.max.two
  selected.max.two.lower <- selected.max.two
  selected.fixed.two.upper <- selected.fixed.two
  selected.fixed.two.lower <- selected.fixed.two
  
  selected.perm.two.upper <- selected.perm.two
  selected.perm.two.lower <- selected.perm.two
  if (nrow(selected.perm.two.upper) > 0) {
    selected.perm.two.upper$two_sided_tail_direction <- "upper_benefit"
    selected.perm.two.lower$two_sided_tail_direction <- "lower_harm"
  }
  
  list(
    observed = obs,
    
    selected.max.one = selected.max.one,
    selected.max.two = selected.max.two,
    selected.max.two.upper = selected.max.two.upper,
    selected.max.two.lower = selected.max.two.lower,
    selected.fixed.one = selected.fixed.one,
    selected.fixed.two = selected.fixed.two,
    selected.fixed.two.upper = selected.fixed.two.upper,
    selected.fixed.two.lower = selected.fixed.two.lower,
    
    T.obs.max.one = T.obs.max.one,
    T.perm.max.one = perm.max.one,
    p.value.max.one = p.value.max.one,
    
    T.obs.max.two = T.obs.max.two,
    T.obs.max.two.log = T.obs.max.two.log,
    T.obs.max.two.upper = T.obs.max.two.log,
    T.obs.max.two.lower = T.obs.max.two.log,
    T.perm.max.two = perm.max.two,
    T.perm.max.two.log = perm.max.two.log,
    T.perm.max.two.upper = perm.max.two.log,
    T.perm.max.two.lower = perm.max.two.log,
    p.value.max.two = p.value.max.two,
    p.value.max.two.upper = p.value.max.two.upper,
    p.value.max.two.lower = p.value.max.two.lower,
    two.sided.tail.direction.max = two.tail.direction.max,
    
    T.obs.fixed.one = T.obs.fixed.one,
    T.perm.fixed.one = perm.fixed.one,
    p.value.fixed.one = p.value.fixed.one,
    
    T.obs.fixed.two = T.obs.fixed.two,
    T.obs.fixed.two.log = T.obs.fixed.two.log,
    T.obs.fixed.two.upper = T.obs.fixed.two.log,
    T.obs.fixed.two.lower = T.obs.fixed.two.log,
    T.perm.fixed.two = perm.fixed.two,
    T.perm.fixed.two.log = perm.fixed.two.log,
    T.perm.fixed.two.upper = perm.fixed.two.log,
    T.perm.fixed.two.lower = perm.fixed.two.log,
    p.value.fixed.two = p.value.fixed.two,
    p.value.fixed.two.upper = p.value.fixed.two.upper,
    p.value.fixed.two.lower = p.value.fixed.two.lower,
    two.sided.tail.direction.fixed = two.tail.direction.fixed,
    
    T.perm.tie.max.one.count = perm.tie.max.one.count,
    T.perm.tie.max.one.pr = perm.tie.max.one.pr,
    T.perm.tie.max.two.count = perm.tie.max.two.count,
    T.perm.tie.max.two.pr = perm.tie.max.two.pr,
    T.perm.tie.max.two.upper.count = perm.tie.max.two.count,
    T.perm.tie.max.two.upper.pr = perm.tie.max.two.pr,
    T.perm.tie.max.two.lower.count = perm.tie.max.two.count,
    T.perm.tie.max.two.lower.pr = perm.tie.max.two.pr,
    
    T.perm.tie.fixed.one.count = perm.fixed.tie.one.count,
    T.perm.tie.fixed.one.pr = perm.fixed.tie.one.pr,
    T.perm.tie.fixed.two.count = perm.fixed.tie.two.count,
    T.perm.tie.fixed.two.pr = perm.fixed.tie.two.pr,
    T.perm.tie.fixed.two.upper.count = perm.fixed.tie.two.count,
    T.perm.tie.fixed.two.upper.pr = perm.fixed.tie.two.pr,
    T.perm.tie.fixed.two.lower.count = perm.fixed.tie.two.count,
    T.perm.tie.fixed.two.lower.pr = perm.fixed.tie.two.pr,
    
    T.perm.pointwise = list(
      death_first = perm.pointwise.death.first,
      hospitalization_first = perm.pointwise.hosp.first,
      threshold = perm.pointwise.threshold
    ),
    p.value.pointwise.one = list(
      death_first = p.value.pointwise.death.first.one,
      hospitalization_first = p.value.pointwise.hosp.first.one,
      threshold = p.value.pointwise.threshold.one
    ),
    p.value.pointwise.two = list(
      death_first = p.value.pointwise.death.first.two,
      hospitalization_first = p.value.pointwise.hosp.first.two,
      threshold = p.value.pointwise.threshold.two
    ),
    p.value.pointwise.two.upper = list(
      death_first = p.value.pointwise.death.first.two.upper,
      hospitalization_first = p.value.pointwise.hosp.first.two.upper,
      threshold = p.value.pointwise.threshold.two.upper
    ),
    p.value.pointwise.two.lower = list(
      death_first = p.value.pointwise.death.first.two.lower,
      hospitalization_first = p.value.pointwise.hosp.first.two.lower,
      threshold = p.value.pointwise.threshold.two.lower
    ),
    
    selected.perm.one = selected.perm.one,
    selected.perm.two = selected.perm.two,
    selected.perm.two.upper = selected.perm.two.upper,
    selected.perm.two.lower = selected.perm.two.lower,
    
    # Backwards-compatible aliases: old code sees one-sided by default.
    T.obs.max = T.obs.max.one,
    T.perm.max = perm.max.one,
    p.value.max = p.value.max.one,
    T.obs.fixed = T.obs.fixed.one,
    T.perm.fixed = perm.fixed.one,
    p.value.fixed = p.value.fixed.one,
    p.value.pointwise = list(
      death_first = p.value.pointwise.death.first.one,
      hospitalization_first = p.value.pointwise.hosp.first.one,
      threshold = p.value.pointwise.threshold.one
    ),
    
    B = B,
    p.grid = p.grid,
    t.grid = t.grid
  )
}

# Log-rank: output both one-sided and two-sided.
# One-sided direction = treatment benefit, i.e. Cox HR < 1 for ARM = 1.
# If ARM=1 is treatment, Cox z < 0 supports benefit, so p_one = Phi(z).
logrank_one_endpoint <- function(dat, time.col, event.col, method, endpoint, note) {
  dat <- dat[is.finite(dat[[time.col]]) & !is.na(dat[[event.col]]) & !is.na(dat$ARM), , drop = FALSE]
  
  if (nrow(dat) == 0 || length(unique(dat$ARM)) < 2) {
    return(data.frame(
      method = method,
      endpoint = endpoint,
      statistic = NA_real_,
      chisq = NA_real_,
      z = NA_real_,
      HR = NA_real_,
      p.value.one.sided = NA_real_,
      p.value.two.sided = NA_real_,
      p.value = NA_real_,
      note = note,
      stringsAsFactors = FALSE
    ))
  }
  
  f <- stats::as.formula(paste0("survival::Surv(", time.col, ", ", event.col, ") ~ ARM"))
  
  lr <- tryCatch(survival::survdiff(f, data = dat), error = function(e) NULL)
  chisq <- if (is.null(lr)) NA_real_ else as.numeric(lr$chisq)
  p.two.lr <- if (is.finite(chisq)) {
    stats::pchisq(chisq, df = length(lr$n) - 1, lower.tail = FALSE)
  } else {
    NA_real_
  }
  
  fit <- tryCatch(survival::coxph(f, data = dat), error = function(e) NULL)
  
  z <- NA_real_
  hr <- NA_real_
  p.one <- NA_real_
  p.two <- p.two.lr
  
  if (!is.null(fit)) {
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    if (!is.null(sm) && nrow(sm$coef) >= 1) {
      beta <- as.numeric(sm$coef[1, "coef"])
      z <- as.numeric(sm$coef[1, "z"])
      hr <- exp(beta)
      if (is.finite(z)) {
        p.one <- stats::pnorm(z)
        p.two <- 2 * stats::pnorm(-abs(z))
      }
    }
  }
  
  data.frame(
    method = method,
    endpoint = endpoint,
    statistic = chisq,
    chisq = chisq,
    z = z,
    HR = hr,
    p.value.one.sided = p.one,
    p.value.two.sided = p.two,
    p.value = p.two,
    note = note,
    stringsAsFactors = FALSE
  )
}

run.logrank.tests <- function(ds) {
  tab <- ds$table.output
  
  death.dat <- data.frame(
    ARM = tab$ARM,
    FUTIME = tab$FUTIME,
    CNSR = tab$CNSR,
    stringsAsFactors = FALSE
  )
  
  comp <- compute.composite.endpoint.data(ds)
  
  rbind(
    logrank_one_endpoint(
      dat = death.dat,
      time.col = "FUTIME",
      event.col = "CNSR",
      method = "Log-rank death endpoint",
      endpoint = "death",
      note = "Death / survival endpoint only; one-sided direction is treatment HR < 1."
    ),
    logrank_one_endpoint(
      dat = comp,
      time.col = "time",
      event.col = "event",
      method = "Log-rank composite endpoint",
      endpoint = "composite",
      note = "Time to earliest of death or first hospitalization; one-sided direction is treatment HR < 1."
    )
  )
}


# Extract method-level rows from one simulated trial.
extract.method.rows.final <- function(scenario.row,
                                      sim.index,
                                      test.out,
                                      logrank.results,
                                      comp.stats = NULL,
                                      threshold.info = NULL,
                                      alpha = ALPHA) {
  wr.methods <- test.out$selected.max.one$method
  
  one <- test.out$selected.max.one
  two <- test.out$selected.max.two
  
  one <- one[match(wr.methods, one$method), , drop = FALSE]
  two <- two[match(wr.methods, two$method), , drop = FALSE]
  
  wr.rows <- data.frame(
    scenario_id = scenario.row$scenario_id,
    scenario_label = scenario_display_label_final(scenario.row$scenario_id),
    scenario_type = scenario.row$scenario_type,
    description = scenario.row$description,
    sim_index = sim.index,
    method = wr.methods,
    method_label = method_label_final(wr.methods),
    statistic_type = "Win ratio",
    
    observed_WR_one_sided_selection = one$selected_WR,
    observed_statistic_one_sided = one$selected_statistic,
    selected_order_one_sided = one$selected_order,
    selected_p_one_sided = one$selected_p,
    selected_t_months_one_sided = one$selected_t_months,
    permutation_p_value_one_sided = as.numeric(test.out$p.value.max.one[wr.methods]),
    fixed_parameter_p_value_one_sided = as.numeric(test.out$p.value.fixed.one[wr.methods]),
    
    observed_WR_two_sided_selection = two$selected_WR,
    observed_statistic_two_sided = two$selected_statistic,
    observed_log_WR_two_sided_selection = two$selected_log_WR,
    observed_abs_log_WR_two_sided_selection = two$selected_abs_log_WR,
    selected_order_two_sided = two$selected_order,
    selected_p_two_sided = two$selected_p,
    selected_t_months_two_sided = two$selected_t_months,
    two_sided_tail_direction = two$two_sided_tail_direction,
    two_sided_upper_tail_p_value = as.numeric(test.out$p.value.max.two.upper[wr.methods]),
    two_sided_lower_tail_p_value = as.numeric(test.out$p.value.max.two.lower[wr.methods]),
    permutation_p_value_two_sided = as.numeric(test.out$p.value.max.two[wr.methods]),
    fixed_two_sided_tail_direction = as.character(test.out$two.sided.tail.direction.fixed[wr.methods]),
    fixed_two_sided_upper_tail_p_value = as.numeric(test.out$p.value.fixed.two.upper[wr.methods]),
    fixed_two_sided_lower_tail_p_value = as.numeric(test.out$p.value.fixed.two.lower[wr.methods]),
    fixed_parameter_p_value_two_sided = as.numeric(test.out$p.value.fixed.two[wr.methods]),
    
    rejected_one_sided = as.numeric(test.out$p.value.max.one[wr.methods]) < alpha,
    rejected_two_sided = as.numeric(test.out$p.value.max.two[wr.methods]) < alpha,
    fixed_rejected_one_sided = as.numeric(test.out$p.value.fixed.one[wr.methods]) < alpha,
    fixed_rejected_two_sided = as.numeric(test.out$p.value.fixed.two[wr.methods]) < alpha,
    
    tie_count_one_sided_selection = one$selected_tie_count,
    tie_proportion_one_sided_selection = one$selected_tie_pr,
    tie_count_two_sided_selection = two$selected_tie_count,
    tie_proportion_two_sided_selection = two$selected_tie_pr,
    
    mean_perm_tie_count_one_sided = sapply(wr.methods, function(nm) safe_mean(test.out$T.perm.tie.max.one.count[, nm])),
    mean_perm_tie_proportion_one_sided = sapply(wr.methods, function(nm) safe_mean(test.out$T.perm.tie.max.one.pr[, nm])),
    mean_perm_tie_count_two_sided = sapply(wr.methods, function(nm) safe_mean(test.out$T.perm.tie.max.two.count[, nm])),
    mean_perm_tie_proportion_two_sided = sapply(wr.methods, function(nm) safe_mean(test.out$T.perm.tie.max.two.pr[, nm])),
    
    B = test.out$B,
    stringsAsFactors = FALSE
  )
  
  lr.map <- c(
    "Log-rank death endpoint" = "logrank_death",
    "Log-rank composite endpoint" = "logrank_composite"
  )
  
  lr.rows <- data.frame()
  if (!is.null(logrank.results) && nrow(logrank.results) > 0) {
    lr.rows <- data.frame(
      scenario_id = scenario.row$scenario_id,
      scenario_label = scenario_display_label_final(scenario.row$scenario_id),
      scenario_type = scenario.row$scenario_type,
      description = scenario.row$description,
      sim_index = sim.index,
      method = unname(lr.map[logrank.results$method]),
      method_label = method_label_final(unname(lr.map[logrank.results$method])),
      statistic_type = "Log-rank",
      
      observed_WR_one_sided_selection = NA_real_,
      observed_statistic_one_sided = logrank.results$z,
      selected_order_one_sided = NA_character_,
      selected_p_one_sided = NA_real_,
      selected_t_months_one_sided = NA_real_,
      permutation_p_value_one_sided = logrank.results$p.value.one.sided,
      fixed_parameter_p_value_one_sided = NA_real_,
      
      observed_WR_two_sided_selection = NA_real_,
      observed_statistic_two_sided = logrank.results$chisq,
      selected_order_two_sided = NA_character_,
      selected_p_two_sided = NA_real_,
      selected_t_months_two_sided = NA_real_,
      permutation_p_value_two_sided = logrank.results$p.value.two.sided,
      fixed_parameter_p_value_two_sided = NA_real_,
      
      rejected_one_sided = logrank.results$p.value.one.sided < alpha,
      rejected_two_sided = logrank.results$p.value.two.sided < alpha,
      fixed_rejected_one_sided = NA,
      fixed_rejected_two_sided = NA,
      
      tie_count_one_sided_selection = NA_real_,
      tie_proportion_one_sided_selection = NA_real_,
      tie_count_two_sided_selection = NA_real_,
      tie_proportion_two_sided_selection = NA_real_,
      
      mean_perm_tie_count_one_sided = NA_real_,
      mean_perm_tie_proportion_one_sided = NA_real_,
      mean_perm_tie_count_two_sided = NA_real_,
      mean_perm_tie_proportion_two_sided = NA_real_,
      
      B = NA_integer_,
      logrank_chisq = logrank.results$chisq,
      logrank_z = logrank.results$z,
      logrank_HR = logrank.results$HR,
      stringsAsFactors = FALSE
    )
  }
  
  # Add missing logrank columns to WR rows so rbind is stable.
  wr.rows$logrank_chisq <- NA_real_
  wr.rows$logrank_z <- NA_real_
  wr.rows$logrank_HR <- NA_real_
  
  out <- rbind_fill_base(wr.rows, lr.rows)
  rownames(out) <- NULL
  out
}

extract.trial.row.final <- function(scenario.row,
                                    sim.index,
                                    method.rows,
                                    comp.stats = NULL,
                                    threshold.info = NULL) {
  data.frame(
    scenario_id = scenario.row$scenario_id,
    scenario_label = scenario_display_label_final(scenario.row$scenario_id),
    scenario_type = scenario.row$scenario_type,
    description = scenario.row$description,
    sim_index = sim.index,
    N0 = scenario.row$N0,
    N1 = scenario.row$N1,
    HR = scenario.row$HR,
    hosp_scale_trt = scenario.row$evt.rate.scale.param.trt,
    max_FU = scenario.row$max.FU,
    censor_rate = scenario.row$censor.rate,
    n_methods = length(unique(method.rows$method)),
    min_p_one_sided = suppressWarnings(min(method.rows$permutation_p_value_one_sided, na.rm = TRUE)),
    min_p_two_sided = suppressWarnings(min(method.rows$permutation_p_value_two_sided, na.rm = TRUE)),
    any_reject_one_sided = any(method.rows$rejected_one_sided, na.rm = TRUE),
    any_reject_two_sided = any(method.rows$rejected_two_sided, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

extract.pointwise.rows.final <- function(scenario.row, sim.index, test.out) {
  obs <- test.out$observed
  p.grid <- test.out$p.grid
  t.grid <- test.out$t.grid
  
  make_p_rows <- function(pathway, d, p.one, p.two, p.two.upper = NULL, p.two.lower = NULL) {
    data.frame(
      scenario_id = scenario.row$scenario_id,
      scenario_label = scenario_display_label_final(scenario.row$scenario_id),
      sim_index = sim.index,
      pathway = pathway,
      candidate_type = "p",
      order = if ("order" %in% names(d)) d$order else NA_character_,
      p = d$p,
      t_months = NA_real_,
      WR = d$WR,
      log_WR = wr_log(d$WR),
      abs_log_WR = wr_abs_log(d$WR),
      permutation_p_value_one_sided = as.numeric(p.one),
      permutation_p_value_two_sided = as.numeric(p.two),
      two_sided_upper_tail_p_value = as.numeric(if (is.null(p.two.upper)) rep(NA_real_, length(p.two)) else p.two.upper),
      two_sided_lower_tail_p_value = as.numeric(if (is.null(p.two.lower)) rep(NA_real_, length(p.two)) else p.two.lower),
      win_pairs = col_or_na(d, c("win.pairs", "wins")),
      loss_pairs = col_or_na(d, c("loss.pairs", "losses")),
      tie_count = col_or_na(d, c("tie.count", "tie_count")),
      tie_proportion = col_or_na(d, c("tie.pr", "tie_proportion", "pr.tie")),
      stringsAsFactors = FALSE
    )
  }
  
  make_t_rows <- function(d, p.one, p.two, p.two.upper = NULL, p.two.lower = NULL) {
    data.frame(
      scenario_id = scenario.row$scenario_id,
      scenario_label = scenario_display_label_final(scenario.row$scenario_id),
      sim_index = sim.index,
      pathway = "threshold",
      candidate_type = "t",
      order = NA_character_,
      p = NA_real_,
      t_months = if ("t.months" %in% names(d)) d$t.months else d$t * 12,
      WR = d$WR,
      log_WR = wr_log(d$WR),
      abs_log_WR = wr_abs_log(d$WR),
      permutation_p_value_one_sided = as.numeric(p.one),
      permutation_p_value_two_sided = as.numeric(p.two),
      two_sided_upper_tail_p_value = as.numeric(if (is.null(p.two.upper)) rep(NA_real_, length(p.two)) else p.two.upper),
      two_sided_lower_tail_p_value = as.numeric(if (is.null(p.two.lower)) rep(NA_real_, length(p.two)) else p.two.lower),
      win_pairs = col_or_na(d, c("wins", "win.pairs")),
      loss_pairs = col_or_na(d, c("losses", "loss.pairs")),
      tie_count = col_or_na(d, c("tie.count", "tie_count")),
      tie_proportion = col_or_na(d, c("tie.pr", "tie_proportion", "pr.tie")),
      stringsAsFactors = FALSE
    )
  }
  
  out <- rbind(
    make_p_rows(
      "death_first",
      obs$weighted.death.first,
      test.out$p.value.pointwise.one$death_first,
      test.out$p.value.pointwise.two$death_first,
      test.out$p.value.pointwise.two.upper$death_first,
      test.out$p.value.pointwise.two.lower$death_first
    ),
    make_p_rows(
      "hospitalization_first",
      obs$weighted.hosp.first,
      test.out$p.value.pointwise.one$hospitalization_first,
      test.out$p.value.pointwise.two$hospitalization_first,
      test.out$p.value.pointwise.two.upper$hospitalization_first,
      test.out$p.value.pointwise.two.lower$hospitalization_first
    ),
    make_t_rows(
      obs$threshold,
      test.out$p.value.pointwise.one$threshold,
      test.out$p.value.pointwise.two$threshold,
      test.out$p.value.pointwise.two.upper$threshold,
      test.out$p.value.pointwise.two.lower$threshold
    )
  )
  
  rownames(out) <- NULL
  out
}


# Run one simulation
run.one.simulation.final <- function(scenario.row,
                                     sim.index,
                                     seed,
                                     B = B_PERM,
                                     save.example = FALSE,
                                     scenario.outdir = OUTDIR) {
  set.seed(seed)
  
  ds <- simulate.one.dataset(
    N = c(scenario.row$N0, scenario.row$N1),
    mort.rate.ctrl = scenario.row$mort.rate.ctrl,
    mort.rate.trt = scenario.row$mort.rate.ctrl * scenario.row$HR,
    evt.rate.shape.param = scenario.row$evt.rate.shape.param,
    evt.rate.scale.param.ctr = scenario.row$evt.rate.scale.param.ctr,
    evt.rate.scale.param.trt = scenario.row$evt.rate.scale.param.trt,
    max.followup = scenario.row$max.FU
  )
  
  ds <- apply.random.censoring(ds, censor.rate = scenario.row$censor.rate, seed = seed + 17L)
  ds <- prepare.ds.fast(ds)
  
  threshold.info <- choose.threshold.grid.primary(ds)
  
  test.out <- perm.test.revised(
    ds = ds,
    B = B,
    seed = seed + 100000L,
    p.grid = if (exists("P_GRID_ALL")) P_GRID_ALL else seq(0.01, 1, by = 0.01),
    t.grid = threshold.info$t.grid,
    verbose = FALSE
  )
  
  logrank.results <- run.logrank.tests(ds)
  comp.stats <- composite.statistics(ds)
  
  if (save.example) {
    if (!dir.exists(scenario.outdir)) dir.create(scenario.outdir, recursive = TRUE)
    saveRDS(
      list(
        ds = ds,
        threshold.info = threshold.info,
        test.out = test.out,
        logrank.results = logrank.results,
        comp.stats = comp.stats
      ),
      file.path(scenario.outdir, paste0("example_sim", sim.index, "_full_object.rds"))
    )
  }
  
  method.rows <- extract.method.rows.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out,
    logrank.results = logrank.results,
    comp.stats = comp.stats,
    threshold.info = threshold.info
  )
  
  trial.row <- extract.trial.row.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    method.rows = method.rows,
    comp.stats = comp.stats,
    threshold.info = threshold.info
  )
  
  pointwise.rows <- extract.pointwise.rows.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out
  )
  
  list(
    trial.row = trial.row,
    method.rows = method.rows,
    pointwise.rows = pointwise.rows
  )
}


# Summaries
summarise.method.power.final <- function(method.rows, alpha = ALPHA) {
  if (is.null(method.rows) || nrow(method.rows) == 0) return(data.frame())
  
  keys <- unique(method.rows[, c("scenario_id", "scenario_label", "scenario_type",
                                 "description", "method", "method_label", "statistic_type"), drop = FALSE])
  
  out <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    d <- method.rows[
      method.rows$scenario_id == k$scenario_id &
        method.rows$method == k$method,
      ,
      drop = FALSE
    ]
    
    data.frame(
      scenario_id = k$scenario_id,
      scenario_label = k$scenario_label,
      scenario_type = k$scenario_type,
      description = k$description,
      method = k$method,
      method_label = k$method_label,
      statistic_type = k$statistic_type,
      nsim.available = nrow(d),
      
      mean_p_one_sided = safe_mean(d$permutation_p_value_one_sided),
      median_p_one_sided = safe_median(d$permutation_p_value_one_sided),
      rejection_proportion_one_sided = mean(d$permutation_p_value_one_sided < alpha, na.rm = TRUE),
      
      mean_p_two_sided = safe_mean(d$permutation_p_value_two_sided),
      median_p_two_sided = safe_median(d$permutation_p_value_two_sided),
      rejection_proportion_two_sided = mean(d$permutation_p_value_two_sided < alpha, na.rm = TRUE),
      
      mean_fixed_p_one_sided = safe_mean(d$fixed_parameter_p_value_one_sided),
      fixed_rejection_proportion_one_sided = mean(d$fixed_parameter_p_value_one_sided < alpha, na.rm = TRUE),
      
      mean_fixed_p_two_sided = safe_mean(d$fixed_parameter_p_value_two_sided),
      fixed_rejection_proportion_two_sided = mean(d$fixed_parameter_p_value_two_sided < alpha, na.rm = TRUE),
      
      mean_WR_one_sided_selection = safe_mean(d$observed_WR_one_sided_selection),
      mean_WR_two_sided_selection = safe_mean(d$observed_WR_two_sided_selection),
      mean_log_WR_two_sided_selection = safe_mean(d$observed_log_WR_two_sided_selection),
      mean_abs_log_WR_two_sided_selection = safe_mean(d$observed_abs_log_WR_two_sided_selection),
      
      mean_selected_p_one_sided = safe_mean(d$selected_p_one_sided),
      mean_selected_p_two_sided = safe_mean(d$selected_p_two_sided),
      mode_selected_order_one_sided = mode_string(d$selected_order_one_sided),
      mode_selected_order_two_sided = mode_string(d$selected_order_two_sided),
      mean_selected_t_months_one_sided = safe_mean(d$selected_t_months_one_sided),
      mean_selected_t_months_two_sided = safe_mean(d$selected_t_months_two_sided),
      
      mean_tie_proportion_one_sided_selection = safe_mean(d$tie_proportion_one_sided_selection),
      mean_tie_proportion_two_sided_selection = safe_mean(d$tie_proportion_two_sided_selection),
      mean_perm_tie_proportion_one_sided = safe_mean(d$mean_perm_tie_proportion_one_sided),
      mean_perm_tie_proportion_two_sided = safe_mean(d$mean_perm_tie_proportion_two_sided),
      
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summarise.pointwise.final <- function(pointwise.rows, alpha = ALPHA) {
  if (is.null(pointwise.rows) || nrow(pointwise.rows) == 0) return(data.frame())
  
  keys <- unique(pointwise.rows[, c("scenario_id", "scenario_label", "pathway",
                                    "candidate_type", "order", "p", "t_months"), drop = FALSE])
  
  out <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    
    same <- pointwise.rows$scenario_id == k$scenario_id &
      pointwise.rows$pathway == k$pathway &
      pointwise.rows$candidate_type == k$candidate_type
    
    if (!is.na(k$order)) same <- same & pointwise.rows$order == k$order
    if (!is.na(k$p)) same <- same & abs(pointwise.rows$p - k$p) < 1e-9
    if (!is.na(k$t_months)) same <- same & abs(pointwise.rows$t_months - k$t_months) < 1e-9
    
    d <- pointwise.rows[same, , drop = FALSE]
    
    data.frame(
      scenario_id = k$scenario_id,
      scenario_label = k$scenario_label,
      pathway = k$pathway,
      candidate_type = k$candidate_type,
      order = k$order,
      p = k$p,
      t_months = k$t_months,
      nsim.available = nrow(d),
      mean_WR = safe_mean(d$WR),
      mean_log_WR = safe_mean(d$log_WR),
      mean_abs_log_WR = safe_mean(d$abs_log_WR),
      mean_p_one_sided = safe_mean(d$permutation_p_value_one_sided),
      rejection_proportion_one_sided = mean(d$permutation_p_value_one_sided < alpha, na.rm = TRUE),
      mean_p_two_sided = safe_mean(d$permutation_p_value_two_sided),
      rejection_proportion_two_sided = mean(d$permutation_p_value_two_sided < alpha, na.rm = TRUE),
      mean_tie_proportion = safe_mean(d$tie_proportion),
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}


# Scenario and global runners
run.scenario.final <- function(scenario.row,
                               nsim = NSIM,
                               B = B_PERM,
                               master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                               outdir = OUTDIR) {
  scenario.label <- scenario_display_label_final(scenario.row$scenario_id)
  scenario.outdir <- file.path(outdir, paste0(sprintf("%02d", scenario.row$scenario_index), "_", scenario.row$scenario_id))
  if (!dir.exists(scenario.outdir)) dir.create(scenario.outdir, recursive = TRUE)
  
  cat("\n===== Scenario", scenario.row$scenario_index, ":", scenario.label, "=====\n")
  cat("HR =", scenario.row$HR,
      "; hosp scale =", scenario.row$evt.rate.scale.param.trt,
      "; FU =", scenario.row$max.FU,
      "; censor =", scenario.row$censor.rate, "\n")
  cat("NSIM =", nsim, "; B_PERM =", B, "\n")
  
  trial.list <- vector("list", nsim)
  method.list <- vector("list", nsim)
  pointwise.list <- vector("list", nsim)
  
  for (s in seq_len(nsim)) {
    if (s == 1 || s == nsim || s %% CHECKPOINT_EVERY == 0) {
      cat("  simulated trial", s, "of", nsim, "\n")
    }
    
    seed <- as.integer(master.seed + scenario.row$scenario_index * 1000000L + s * 104729L)
    
    res <- tryCatch(
      run.one.simulation.final(
        scenario.row = scenario.row,
        sim.index = s,
        seed = seed,
        B = B,
        save.example = (SAVE_EXAMPLE_PLOTS && s == 1),
        scenario.outdir = scenario.outdir
      ),
      error = function(e) {
        message("  ERROR in scenario ", scenario.row$scenario_id, ", sim ", s, ": ", conditionMessage(e))
        err.method <- data.frame(
          scenario_id = scenario.row$scenario_id,
          scenario_label = scenario.label,
          scenario_type = scenario.row$scenario_type,
          description = scenario.row$description,
          sim_index = s,
          method = NA_character_,
          method_label = NA_character_,
          statistic_type = NA_character_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
        list(
          trial.row = data.frame(
            scenario_id = scenario.row$scenario_id,
            scenario_label = scenario.label,
            scenario_type = scenario.row$scenario_type,
            description = scenario.row$description,
            sim_index = s,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          ),
          method.rows = err.method,
          pointwise.rows = data.frame()
        )
      }
    )
    
    trial.list[[s]] <- res$trial.row
    method.list[[s]] <- res$method.rows
    pointwise.list[[s]] <- res$pointwise.rows
    
    if (s %% CHECKPOINT_EVERY == 0 || s == nsim) {
      trial.tmp <- rbind_fill_base(trial.list[seq_len(s)])
      method.tmp <- rbind_fill_base(method.list[seq_len(s)])
      point.tmp <- rbind_fill_base(pointwise.list[seq_len(s)])
      
      write.csv(trial.tmp, file.path(scenario.outdir, "checkpoint_trial_summary.csv"), row.names = FALSE)
      write.csv(method.tmp, file.path(scenario.outdir, "checkpoint_method_trial_results.csv"), row.names = FALSE)
      write.csv(point.tmp, file.path(scenario.outdir, "checkpoint_pointwise_results.csv"), row.names = FALSE)
    }
  }
  
  trial.rows <- rbind_fill_base(trial.list)
  method.rows <- rbind_fill_base(method.list)
  pointwise.rows <- rbind_fill_base(pointwise.list)
  
  power.summary <- summarise.method.power.final(method.rows, alpha = ALPHA)
  pointwise.summary <- summarise.pointwise.final(pointwise.rows, alpha = ALPHA)
  
  write.csv(trial.rows, file.path(scenario.outdir, "trial_summary.csv"), row.names = FALSE)
  write.csv(method.rows, file.path(scenario.outdir, "method_trial_results.csv"), row.names = FALSE)
  write.csv(pointwise.rows, file.path(scenario.outdir, "pointwise_results.csv"), row.names = FALSE)
  write.csv(power.summary, file.path(scenario.outdir, "power_summary.csv"), row.names = FALSE)
  write.csv(pointwise.summary, file.path(scenario.outdir, "pointwise_summary.csv"), row.names = FALSE)
  
  list(
    trial = trial.rows,
    method = method.rows,
    pointwise = pointwise.rows,
    power = power.summary,
    pointwise.summary = pointwise.summary
  )
}

run.all.scenarios.final <- function(scenario.grid = build.final.15.scenario.grid(),
                                    nsim = NSIM,
                                    B = B_PERM,
                                    outdir = OUTDIR,
                                    master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  scenario.grid$scenario_index <- seq_len(nrow(scenario.grid))
  scenario.grid$scenario_label <- scenario_display_label_final(scenario.grid$scenario_id)
  
  write.csv(scenario.grid, file.path(outdir, "GLOBAL_scenario_config.csv"), row.names = FALSE)
  
  settings.table <- data.frame(
    setting = c(
      "NSIM",
      "B_PERM",
      "MASTER_SEED",
      "ALPHA",
      "number_of_scenarios",
      "WR_outputs",
      "WR_one_sided_statistic",
      "WR_two_sided_statistic",
      "WR_adaptive_one_sided_rule",
      "WR_adaptive_two_sided_rule",
      "WR_permutation_scheme",
      "logrank_outputs",
      "one_sided_logrank_direction",
      "p_grid_low_exploratory",
      "p_grid_primary",
      "p_grid_full",
      "clinical_t_months",
      "checkpoint_every",
      "resume_if_exists",
      "save_example_plots",
      "output_folder"
    ),
    value = c(
      as.character(nsim),
      as.character(B),
      as.character(master.seed),
      as.character(ALPHA),
      as.character(nrow(scenario.grid)),
      "one-sided and two-sided for every WR pathway",
      "WR; larger values favor treatment",
      "abs(log(WR)); right-tail permutation comparison on distance from the null value 1",
      "max WR, selection repeated inside each permutation",
      "max abs(log(WR)), with the same selection repeated inside each permutation",
      "treatment-label permutation",
      "one-sided and two-sided for death and composite endpoints",
      "treatment benefit, Cox HR < 1 for ARM = 1",
      if (exists("P_GRID_EXPLORATORY")) paste0(sprintf("%.2f", min(P_GRID_EXPLORATORY)), " to ", sprintf("%.2f", max(P_GRID_EXPLORATORY))) else NA_character_,
      if (exists("P_GRID_PRIMARY")) paste0(sprintf("%.2f", min(P_GRID_PRIMARY)), " to ", sprintf("%.2f", max(P_GRID_PRIMARY))) else NA_character_,
      if (exists("P_GRID_ALL")) paste0(sprintf("%.2f", min(P_GRID_ALL)), " to ", sprintf("%.2f", max(P_GRID_ALL)), "; no p = 0") else NA_character_,
      if (exists("CLINICAL_T_MONTHS")) paste(CLINICAL_T_MONTHS, collapse = ", ") else "1, 3, 6, 12, 18, 24",
      as.character(CHECKPOINT_EVERY),
      as.character(RESUME_IF_EXISTS),
      as.character(SAVE_EXAMPLE_PLOTS),
      outdir
    ),
    stringsAsFactors = FALSE
  )
  write.csv(settings.table, file.path(outdir, "GLOBAL_settings.csv"), row.names = FALSE)
  
  cat("\n===== FINAL Simulation settings =====\n")
  cat("Output folder:", outdir, "\n")
  cat("NSIM =", nsim, "; B_PERM =", B, "; scenarios =", nrow(scenario.grid), "\n")
  cat("WR outputs: one-sided and two-sided\n")
  cat("WR two-sided statistic: abs(log(WR)); adaptive rule: max abs(log(WR)) repeated inside each permutation\n")
  cat("Log-rank outputs: one-sided benefit and two-sided\n")
  
  all.trial <- data.frame()
  all.method <- data.frame()
  all.pointwise <- data.frame()
  all.power <- data.frame()
  all.pointwise.summary <- data.frame()
  
  for (i in seq_len(nrow(scenario.grid))) {
    scenario.row <- scenario.grid[i, , drop = FALSE]
    
    out.i <- run.scenario.final(
      scenario.row = scenario.row,
      nsim = nsim,
      B = B,
      master.seed = master.seed,
      outdir = outdir
    )
    
    all.trial <- rbind_fill_base(all.trial, out.i$trial)
    all.method <- rbind_fill_base(all.method, out.i$method)
    all.pointwise <- rbind_fill_base(all.pointwise, out.i$pointwise)
    all.power <- rbind_fill_base(all.power, out.i$power)
    all.pointwise.summary <- rbind_fill_base(all.pointwise.summary, out.i$pointwise.summary)
    
    write.csv(all.trial, file.path(outdir, "GLOBAL_trial_summary.csv"), row.names = FALSE)
    write.csv(all.method, file.path(outdir, "GLOBAL_method_trial_results.csv"), row.names = FALSE)
    write.csv(all.pointwise, file.path(outdir, "GLOBAL_pointwise_results.csv"), row.names = FALSE)
    write.csv(all.power, file.path(outdir, "GLOBAL_power_summary.csv"), row.names = FALSE)
    write.csv(all.pointwise.summary, file.path(outdir, "GLOBAL_pointwise_summary.csv"), row.names = FALSE)
  }
  
  manifest <- data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE)
  write.csv(manifest, file.path(outdir, "GLOBAL_output_manifest.csv"), row.names = FALSE)
  
  list(
    scenario.grid = scenario.grid,
    trial = all.trial,
    method = all.method,
    pointwise = all.pointwise,
    power = all.power,
    pointwise.summary = all.pointwise.summary,
    settings = settings.table,
    manifest = manifest
  )
}

run.all.scenarios <- function(scenario.grid = build.final.15.scenario.grid(),
                              nsim = NSIM,
                              B = B_PERM,
                              outdir = OUTDIR,
                              master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L) {
  run.all.scenarios.final(
    scenario.grid = scenario.grid,
    nsim = nsim,
    B = B,
    outdir = outdir,
    master.seed = master.seed
  )
}





# output:
# scenario-level raw/method/pointwise/power/average/pathway files
# global old-style aliases such as GLOBAL_all_scenarios_* files
# selected example_sim1 tables and simple figures
# GLOBAL_full_results_bundle.rds
# update:
# 15 scenarios
# NSIM = 1000, B_PERM = 500
# WR one-sided and two-sided outputs
# log-rank one-sided and two-sided outputs


OUTDIR <- "simulation_WR_15scenarios_WR_and_LR_one_two_sided_abs_log_1000x500_FULL_OUTPUT_PARALLEL_SAFE"
NSIM <- 1000L
B_PERM <- 500L
ALPHA <- 0.05
CHECKPOINT_EVERY <- 25L
RESUME_IF_EXISTS <- TRUE
SAVE_EXAMPLE_PLOTS <- TRUE
VERBOSE <- TRUE

# ---------------------------------------------------------------------
# Small output helpers
# ---------------------------------------------------------------------
safe_write_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(x, file, row.names = FALSE)
  invisible(file)
}

safe_save_rds <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file)
  invisible(file)
}

safe_png_open <- function(file, width = 2400, height = 1400, res = 160) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    grDevices::png(file, width = width, height = height, res = res)
    TRUE
  }, error = function(e) FALSE)
  ok
}

safe_dev_off <- function(ok) {
  if (isTRUE(ok)) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  invisible(NULL)
}

finite_or_na <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}


# Example-trial tables, restored.
make.example.weighted.curve <- function(test.out, order = c("death_first", "hospitalization_first")) {
  order <- match.arg(order)
  d <- if (order == "death_first") test.out$observed$weighted.death.first else test.out$observed$weighted.hosp.first
  if (is.null(d) || nrow(d) == 0) return(data.frame())
  p.one <- if (order == "death_first") test.out$p.value.pointwise.one$death_first else test.out$p.value.pointwise.one$hospitalization_first
  p.two <- if (order == "death_first") test.out$p.value.pointwise.two$death_first else test.out$p.value.pointwise.two$hospitalization_first
  p.two.upper <- if (order == "death_first") test.out$p.value.pointwise.two.upper$death_first else test.out$p.value.pointwise.two.upper$hospitalization_first
  p.two.lower <- if (order == "death_first") test.out$p.value.pointwise.two.lower$death_first else test.out$p.value.pointwise.two.lower$hospitalization_first
  out <- d
  out$log_WR <- wr_log(out$WR)
  out$abs_log_WR <- wr_abs_log(out$WR)
  out$pointwise_p_value_one_sided <- as.numeric(p.one)
  out$pointwise_p_value_two_sided <- as.numeric(p.two)
  out$pointwise_p_value_two_sided_upper_tail <- as.numeric(p.two.upper)
  out$pointwise_p_value_two_sided_lower_tail <- as.numeric(p.two.lower)
  out
}

make.example.threshold.curve <- function(test.out) {
  d <- test.out$observed$threshold
  if (is.null(d) || nrow(d) == 0) return(data.frame())
  out <- d
  out$log_WR <- wr_log(out$WR)
  out$abs_log_WR <- wr_abs_log(out$WR)
  out$pointwise_p_value_one_sided <- as.numeric(test.out$p.value.pointwise.one$threshold)
  out$pointwise_p_value_two_sided <- as.numeric(test.out$p.value.pointwise.two$threshold)
  out$pointwise_p_value_two_sided_upper_tail <- as.numeric(test.out$p.value.pointwise.two.upper$threshold)
  out$pointwise_p_value_two_sided_lower_tail <- as.numeric(test.out$p.value.pointwise.two.lower$threshold)
  out
}

make.example.pvalue.table <- function(test.out, logrank.results = NULL) {
  wr.methods <- names(test.out$p.value.max.one)
  wr <- data.frame(
    method = wr.methods,
    method_label = method_label_final(wr.methods),
    statistic_family = "WR",
    observed_statistic_one_sided = as.numeric(test.out$T.obs.max.one[wr.methods]),
    max_permutation_p_value_one_sided = as.numeric(test.out$p.value.max.one[wr.methods]),
    fixed_parameter_p_value_one_sided = as.numeric(test.out$p.value.fixed.one[wr.methods]),
    observed_statistic_two_sided = as.numeric(test.out$T.obs.max.two[wr.methods]),
    two_sided_tail_direction = as.character(test.out$two.sided.tail.direction.max[wr.methods]),
    two_sided_upper_tail_p_value = as.numeric(test.out$p.value.max.two.upper[wr.methods]),
    two_sided_lower_tail_p_value = as.numeric(test.out$p.value.max.two.lower[wr.methods]),
    max_permutation_p_value_two_sided = as.numeric(test.out$p.value.max.two[wr.methods]),
    fixed_two_sided_tail_direction = as.character(test.out$two.sided.tail.direction.fixed[wr.methods]),
    fixed_two_sided_upper_tail_p_value = as.numeric(test.out$p.value.fixed.two.upper[wr.methods]),
    fixed_two_sided_lower_tail_p_value = as.numeric(test.out$p.value.fixed.two.lower[wr.methods]),
    fixed_parameter_p_value_two_sided = as.numeric(test.out$p.value.fixed.two[wr.methods]),
    selected_order_one_sided = test.out$selected.max.one$selected_order[match(wr.methods, test.out$selected.max.one$method)],
    selected_p_one_sided = test.out$selected.max.one$selected_p[match(wr.methods, test.out$selected.max.one$method)],
    selected_t_months_one_sided = test.out$selected.max.one$selected_t_months[match(wr.methods, test.out$selected.max.one$method)],
    selected_order_two_sided = test.out$selected.max.two$selected_order[match(wr.methods, test.out$selected.max.two$method)],
    selected_p_two_sided = test.out$selected.max.two$selected_p[match(wr.methods, test.out$selected.max.two$method)],
    selected_t_months_two_sided = test.out$selected.max.two$selected_t_months[match(wr.methods, test.out$selected.max.two$method)],
    stringsAsFactors = FALSE
  )
  
  lr <- data.frame()
  if (!is.null(logrank.results) && nrow(logrank.results) > 0) {
    map <- c("Log-rank death endpoint" = "logrank_death",
             "Log-rank composite endpoint" = "logrank_composite")
    mids <- unname(map[logrank.results$method])
    lr <- data.frame(
      method = mids,
      method_label = method_label_final(mids),
      statistic_family = "Log-rank",
      observed_statistic_one_sided = logrank.results$z,
      max_permutation_p_value_one_sided = logrank.results$p.value.one.sided,
      fixed_parameter_p_value_one_sided = NA_real_,
      observed_statistic_two_sided = logrank.results$chisq,
      max_permutation_p_value_two_sided = logrank.results$p.value.two.sided,
      fixed_parameter_p_value_two_sided = NA_real_,
      selected_order_one_sided = NA_character_,
      selected_p_one_sided = NA_real_,
      selected_t_months_one_sided = NA_real_,
      selected_order_two_sided = NA_character_,
      selected_p_two_sided = NA_real_,
      selected_t_months_two_sided = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  rbind_fill_base(wr, lr)
}

make.example.selected.max.table <- function(test.out) {
  one <- test.out$selected.max.one
  two <- test.out$selected.max.two
  methods <- one$method
  out <- data.frame(
    method = methods,
    method_label = method_label_final(methods),
    one_sided_selected_order = one$selected_order,
    one_sided_selected_p = one$selected_p,
    one_sided_selected_t_months = one$selected_t_months,
    one_sided_selected_WR = one$selected_WR,
    one_sided_statistic = one$selected_statistic,
    one_sided_permutation_p_value = as.numeric(test.out$p.value.max.one[methods]),
    one_sided_fixed_parameter_p_value = as.numeric(test.out$p.value.fixed.one[methods]),
    one_sided_tie_count = one$selected_tie_count,
    one_sided_tie_proportion = one$selected_tie_pr,
    two_sided_selected_order = two$selected_order[match(methods, two$method)],
    two_sided_selected_p = two$selected_p[match(methods, two$method)],
    two_sided_selected_t_months = two$selected_t_months[match(methods, two$method)],
    two_sided_selected_WR = two$selected_WR[match(methods, two$method)],
    two_sided_statistic_log_WR = two$selected_log_WR[match(methods, two$method)],
    two_sided_statistic_abs_log_WR = two$selected_statistic[match(methods, two$method)],
    two_sided_tail_direction = two$two_sided_tail_direction[match(methods, two$method)],
    two_sided_upper_tail_p_value = as.numeric(test.out$p.value.max.two.upper[methods]),
    two_sided_lower_tail_p_value = as.numeric(test.out$p.value.max.two.lower[methods]),
    two_sided_permutation_p_value = as.numeric(test.out$p.value.max.two[methods]),
    two_sided_fixed_tail_direction = as.character(test.out$two.sided.tail.direction.fixed[methods]),
    two_sided_fixed_upper_tail_p_value = as.numeric(test.out$p.value.fixed.two.upper[methods]),
    two_sided_fixed_lower_tail_p_value = as.numeric(test.out$p.value.fixed.two.lower[methods]),
    two_sided_fixed_parameter_p_value = as.numeric(test.out$p.value.fixed.two[methods]),
    two_sided_tie_count = two$selected_tie_count[match(methods, two$method)],
    two_sided_tie_proportion = two$selected_tie_pr[match(methods, two$method)],
    mean_perm_tie_count_one_sided = sapply(methods, function(nm) safe_mean(test.out$T.perm.tie.max.one.count[, nm])),
    mean_perm_tie_proportion_one_sided = sapply(methods, function(nm) safe_mean(test.out$T.perm.tie.max.one.pr[, nm])),
    mean_perm_tie_count_two_sided = sapply(methods, function(nm) safe_mean(test.out$T.perm.tie.max.two.count[, nm])),
    mean_perm_tie_proportion_two_sided = sapply(methods, function(nm) safe_mean(test.out$T.perm.tie.max.two.pr[, nm])),
    stringsAsFactors = FALSE
  )
  out
}

make.example.pointwise.table <- function(test.out) {
  d1 <- make.example.weighted.curve(test.out, "death_first")
  d2 <- make.example.weighted.curve(test.out, "hospitalization_first")
  dt <- make.example.threshold.curve(test.out)
  
  if (nrow(d1) > 0) {
    d1$family <- "weighted_p"
    d1$endpoint_order <- "death_first"
    d1$parameter_name <- "p"
    d1$parameter <- d1$p
    d1$parameter_months <- NA_real_
  }
  if (nrow(d2) > 0) {
    d2$family <- "weighted_p"
    d2$endpoint_order <- "hospitalization_first"
    d2$parameter_name <- "p"
    d2$parameter <- d2$p
    d2$parameter_months <- NA_real_
  }
  if (nrow(dt) > 0) {
    dt$family <- "threshold_t"
    dt$endpoint_order <- NA_character_
    dt$parameter_name <- "t_months"
    dt$parameter <- dt$t
    dt$parameter_months <- if ("t.months" %in% names(dt)) dt$t.months else dt$t * 12
  }
  rbind_fill_base(d1, d2, dt)
}

make.example.max.vs.fixed.table <- function(test.out) {
  methods <- names(test.out$p.value.max.one)
  data.frame(
    method = methods,
    method_label = method_label_final(methods),
    max_p_value_one_sided = as.numeric(test.out$p.value.max.one[methods]),
    fixed_parameter_p_value_one_sided = as.numeric(test.out$p.value.fixed.one[methods]),
    max_minus_fixed_p_value_one_sided = as.numeric(test.out$p.value.max.one[methods]) - as.numeric(test.out$p.value.fixed.one[methods]),
    max_p_value_two_sided = as.numeric(test.out$p.value.max.two[methods]),
    fixed_parameter_p_value_two_sided = as.numeric(test.out$p.value.fixed.two[methods]),
    max_minus_fixed_p_value_two_sided = as.numeric(test.out$p.value.max.two[methods]) - as.numeric(test.out$p.value.fixed.two[methods]),
    stringsAsFactors = FALSE
  )
}

plot.example.weighted.curves <- function(test.out, outdir, prefix) {
  d1 <- make.example.weighted.curve(test.out, "death_first")
  d2 <- make.example.weighted.curve(test.out, "hospitalization_first")
  if (nrow(d1) == 0 && nrow(d2) == 0) return(invisible(NULL))
  ok <- safe_png_open(file.path(outdir, paste0(prefix, "_weighted_p_curve.png")))
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  y <- c(d1$WR, d2$WR)
  y <- y[is.finite(y)]
  ylim <- if (length(y) == 0) c(0, 2) else range(y, na.rm = TRUE)
  if (diff(ylim) == 0) ylim <- ylim + c(-0.2, 0.2)
  graphics::plot(d1$p, d1$WR, type = "l", lwd = 2, ylim = ylim,
                 xlab = "p", ylab = "Observed WR", main = paste0(prefix, ": weighted p curves"))
  if (nrow(d2) > 0) graphics::lines(d2$p, d2$WR, lwd = 2, lty = 2)
  graphics::abline(h = 1, lty = 3)
  graphics::legend("topright", legend = c("death first", "hospitalization first"), lty = c(1, 2), lwd = 2, bty = "n")
  invisible(NULL)
}

plot.example.threshold.curve <- function(test.out, outdir, prefix) {
  dt <- make.example.threshold.curve(test.out)
  if (nrow(dt) == 0) return(invisible(NULL))
  x <- if ("t.months" %in% names(dt)) dt$t.months else dt$t * 12
  ok <- safe_png_open(file.path(outdir, paste0(prefix, "_threshold_curve.png")))
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  graphics::plot(x, dt$WR, type = "b", lwd = 2,
                 xlab = "Threshold t, months", ylab = "Observed WR",
                 main = paste0(prefix, ": threshold WR curve"))
  graphics::abline(h = 1, lty = 3)
  invisible(NULL)
}

plot.example.pvalues <- function(pvalue.table, outdir, prefix) {
  if (is.null(pvalue.table) || nrow(pvalue.table) == 0) return(invisible(NULL))
  key <- pvalue.table[pvalue.table$method %in% c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite"), , drop = FALSE]
  if (nrow(key) == 0) return(invisible(NULL))
  ok <- safe_png_open(file.path(outdir, paste0(prefix, "_pvalue_key_methods.png")), width = 2600, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(
    one_sided = finite_or_na(key$max_permutation_p_value_one_sided),
    two_sided = finite_or_na(key$max_permutation_p_value_two_sided)
  )
  colnames(mat) <- method_label_final(key$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = "p-value", main = paste0(prefix, ": one- vs two-sided p-values"),
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

save.example.outputs.final <- function(ds, test.out, threshold.info, logrank.results, comp.stats,
                                       scenario.row, sim.index, method.rows, pointwise.rows,
                                       scenario.outdir) {
  prefix <- paste0(scenario.row$scenario_id, "_example_sim", sim.index)
  
  weighted.death <- make.example.weighted.curve(test.out, "death_first")
  weighted.hosp <- make.example.weighted.curve(test.out, "hospitalization_first")
  threshold.curve <- make.example.threshold.curve(test.out)
  pvalue.table <- make.example.pvalue.table(test.out, logrank.results)
  selected.table <- make.example.selected.max.table(test.out)
  pointwise.table <- make.example.pointwise.table(test.out)
  max.fixed.table <- make.example.max.vs.fixed.table(test.out)
  
  safe_write_csv(weighted.death, file.path(scenario.outdir, paste0(prefix, "_weighted_death_first_curve.csv")))
  safe_write_csv(weighted.hosp, file.path(scenario.outdir, paste0(prefix, "_weighted_hospitalization_first_curve.csv")))
  safe_write_csv(threshold.curve, file.path(scenario.outdir, paste0(prefix, "_threshold_curve.csv")))
  if (!is.null(threshold.info$diagnostics)) safe_write_csv(threshold.info$diagnostics, file.path(scenario.outdir, paste0(prefix, "_threshold_diagnostics_first_endpoint_only.csv")))
  if (!is.null(threshold.info$threshold.table)) safe_write_csv(threshold.info$threshold.table, file.path(scenario.outdir, paste0(prefix, "_threshold_grid.csv")))
  safe_write_csv(logrank.results, file.path(scenario.outdir, paste0(prefix, "_logrank_results.csv")))
  safe_write_csv(comp.stats, file.path(scenario.outdir, paste0(prefix, "_composite_statistics.csv")))
  safe_write_csv(pvalue.table, file.path(scenario.outdir, paste0(prefix, "_pvalue_table.csv")))
  safe_write_csv(selected.table, file.path(scenario.outdir, paste0(prefix, "_selected_max_parameter_table.csv")))
  safe_write_csv(pointwise.table, file.path(scenario.outdir, paste0(prefix, "_pointwise_pvalue_and_tie_table.csv")))
  safe_write_csv(method.rows, file.path(scenario.outdir, paste0(prefix, "_pathway_results_table.csv")))
  safe_write_csv(max.fixed.table, file.path(scenario.outdir, paste0(prefix, "_max_vs_fixed_parameter_pvalues.csv")))
  
  safe_save_rds(
    list(
      ds = ds,
      threshold.info = threshold.info,
      test.out = test.out,
      logrank.results = logrank.results,
      comp.stats = comp.stats,
      method.rows = method.rows,
      pointwise.rows = pointwise.rows
    ),
    file.path(scenario.outdir, paste0(prefix, "_full_object.rds"))
  )
  
  if (isTRUE(SAVE_EXAMPLE_PLOTS)) {
    plot.example.weighted.curves(test.out, scenario.outdir, prefix)
    plot.example.threshold.curve(test.out, scenario.outdir, prefix)
    plot.example.pvalues(pvalue.table, scenario.outdir, prefix)
  }
  invisible(NULL)
}


# Scenario-level summaries and compatibility names-
make.average.statistics.final <- function(method.rows, alpha = ALPHA) {
  # In the new long-format output, this is the method-level average table.
  # It intentionally contains both p-value power and descriptive averages.
  summarise.method.power.final(method.rows, alpha = alpha)
}

make.pathway.method.comparison.final <- function(power.summary) {
  if (is.null(power.summary) || nrow(power.summary) == 0) return(data.frame())
  out <- power.summary
  out$pathway_group <- ifelse(out$statistic_type == "Log-rank", "logrank",
                              ifelse(grepl("Order", out$method, ignore.case = TRUE), "order_adaptive_WR",
                                     ifelse(grepl("maxWRp", out$method), "weighted_p_WR",
                                            ifelse(grepl("maxWRt", out$method), "threshold_WR", "fixed_or_traditional_WR"))))
  keep.first <- c("scenario_id", "scenario_label", "scenario_type", "description", "pathway_group", "method", "method_label", "statistic_type")
  out <- out[, c(keep.first, setdiff(names(out), keep.first)), drop = FALSE]
  out
}

plot.scenario.power.final <- function(power.summary, scenario.outdir, scenario.id) {
  if (is.null(power.summary) || nrow(power.summary) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  d <- power.summary[power.summary$method %in% key.methods, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  ok <- safe_png_open(file.path(scenario.outdir, paste0(scenario.id, "_power_comparison.png")), width = 2600, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(
    one_sided = finite_or_na(d$rejection_proportion_one_sided),
    two_sided = finite_or_na(d$rejection_proportion_two_sided)
  )
  colnames(mat) <- method_label_final(d$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = paste0("Rejection proportion, alpha = ", ALPHA),
                    main = paste0(scenario.id, ": one- vs two-sided power"),
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

plot.scenario.pvalues.final <- function(method.rows, scenario.outdir, scenario.id) {
  if (is.null(method.rows) || nrow(method.rows) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  d <- method.rows[method.rows$method %in% key.methods, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  # Plot median p-values by method; this is light and works even for large nsim.
  agg <- summarise.method.power.final(d, alpha = ALPHA)
  ok <- safe_png_open(file.path(scenario.outdir, paste0(scenario.id, "_permutation_pvalue_comparison.png")), width = 2600, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(
    one_sided = finite_or_na(agg$median_p_one_sided),
    two_sided = finite_or_na(agg$median_p_two_sided)
  )
  colnames(mat) <- method_label_final(agg$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = "Median permutation/log-rank p-value",
                    main = paste0(scenario.id, ": median p-values"),
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

plot.scenario.max.fixed.final <- function(method.rows, scenario.outdir, scenario.id) {
  if (is.null(method.rows) || nrow(method.rows) == 0) return(invisible(NULL))
  d <- method.rows[method.rows$statistic_type == "Win ratio", , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  agg <- summarise.method.power.final(d, alpha = ALPHA)
  ok <- safe_png_open(file.path(scenario.outdir, paste0(scenario.id, "_max_vs_fixed_selected_pvalues.png")), width = 2600, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(
    max_one_sided = finite_or_na(agg$median_p_one_sided),
    fixed_one_sided = finite_or_na(agg$mean_fixed_p_one_sided),
    max_two_sided = finite_or_na(agg$median_p_two_sided),
    fixed_two_sided = finite_or_na(agg$mean_fixed_p_two_sided)
  )
  colnames(mat) <- method_label_final(agg$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = "p-value", main = paste0(scenario.id, ": max vs fixed selected p-values"),
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.75), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

write.scenario.outputs.final <- function(scenario.id, scenario.outdir, trial.rows, method.rows, pointwise.rows,
                                         power.summary, average.statistics, method.comparison,
                                         pointwise.summary) {
  # New explicit names
  safe_write_csv(trial.rows, file.path(scenario.outdir, "trial_summary.csv"))
  safe_write_csv(method.rows, file.path(scenario.outdir, "method_trial_results.csv"))
  safe_write_csv(pointwise.rows, file.path(scenario.outdir, "pointwise_results.csv"))
  safe_write_csv(power.summary, file.path(scenario.outdir, "power_summary.csv"))
  safe_write_csv(average.statistics, file.path(scenario.outdir, "average_statistics.csv"))
  safe_write_csv(method.comparison, file.path(scenario.outdir, "pathway_method_comparison.csv"))
  safe_write_csv(pointwise.summary, file.path(scenario.outdir, "pointwise_summary.csv"))
  
  # Old-style compatibility names
  safe_write_csv(method.rows, file.path(scenario.outdir, paste0(scenario.id, "_raw_results.csv")))
  safe_write_csv(pointwise.rows, file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.csv")))
  safe_write_csv(power.summary, file.path(scenario.outdir, paste0(scenario.id, "_power_summary.csv")))
  safe_write_csv(average.statistics, file.path(scenario.outdir, paste0(scenario.id, "_average_statistics.csv")))
  safe_write_csv(method.comparison, file.path(scenario.outdir, paste0(scenario.id, "_pathway_method_comparison.csv")))
  safe_write_csv(pointwise.summary, file.path(scenario.outdir, paste0(scenario.id, "_pointwise_average_pvalue_tie_summary.csv")))
  
  safe_save_rds(method.rows, file.path(scenario.outdir, paste0(scenario.id, "_raw_results.rds")))
  safe_save_rds(pointwise.rows, file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.rds")))
  invisible(NULL)
}


# Override run.one.simulation.final.
run.one.simulation.final <- function(scenario.row,
                                     sim.index,
                                     seed,
                                     B = B_PERM,
                                     save.example = FALSE,
                                     scenario.outdir = OUTDIR) {
  set.seed(seed)
  
  ds <- simulate.one.dataset(
    N = c(scenario.row$N0, scenario.row$N1),
    mort.rate.ctrl = scenario.row$mort.rate.ctrl,
    mort.rate.trt = scenario.row$mort.rate.ctrl * scenario.row$HR,
    evt.rate.shape.param = scenario.row$evt.rate.shape.param,
    evt.rate.scale.param.ctr = scenario.row$evt.rate.scale.param.ctr,
    evt.rate.scale.param.trt = scenario.row$evt.rate.scale.param.trt,
    max.followup = scenario.row$max.FU
  )
  
  ds <- apply.random.censoring(ds, censor.rate = scenario.row$censor.rate, seed = seed + 17L)
  ds <- prepare.ds.fast(ds)
  threshold.info <- choose.threshold.grid.primary(ds)
  
  test.out <- perm.test.revised(
    ds = ds,
    B = B,
    seed = seed + 100000L,
    p.grid = if (exists("P_GRID_ALL")) P_GRID_ALL else seq(0.01, 1, by = 0.01),
    t.grid = threshold.info$t.grid,
    verbose = FALSE
  )
  
  logrank.results <- run.logrank.tests(ds)
  comp.stats <- composite.statistics(ds)
  
  method.rows <- extract.method.rows.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out,
    logrank.results = logrank.results,
    comp.stats = comp.stats,
    threshold.info = threshold.info
  )
  
  trial.row <- extract.trial.row.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    method.rows = method.rows,
    comp.stats = comp.stats,
    threshold.info = threshold.info
  )
  
  pointwise.rows <- extract.pointwise.rows.final(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out
  )
  
  if (isTRUE(save.example)) {
    save.example.outputs.final(
      ds = ds,
      test.out = test.out,
      threshold.info = threshold.info,
      logrank.results = logrank.results,
      comp.stats = comp.stats,
      scenario.row = scenario.row,
      sim.index = sim.index,
      method.rows = method.rows,
      pointwise.rows = pointwise.rows,
      scenario.outdir = scenario.outdir
    )
  }
  
  list(
    trial.row = trial.row,
    method.rows = method.rows,
    pointwise.rows = pointwise.rows
  )
}


# Override scenario runner
run.scenario.final <- function(scenario.row,
                               nsim = NSIM,
                               B = B_PERM,
                               master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                               outdir = OUTDIR,
                               resume = RESUME_IF_EXISTS,
                               checkpoint.every = CHECKPOINT_EVERY,
                               verbose = VERBOSE) {
  scenario.label <- scenario_display_label_final(scenario.row$scenario_id)
  scenario.outdir <- file.path(outdir, paste0(sprintf("%02d", scenario.row$scenario_index), "_", scenario.row$scenario_id))
  dir.create(scenario.outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(verbose)) {
    cat("\n===== Scenario", scenario.row$scenario_index, ":", scenario.label, "=====\n")
    cat("HR =", scenario.row$HR,
        "; hosp scale =", scenario.row$evt.rate.scale.param.trt,
        "; FU =", scenario.row$max.FU,
        "; censor =", scenario.row$censor.rate, "\n")
    cat("NSIM =", nsim, "; B_PERM =", B, "\n")
  }
  
  scenario.id <- as.character(scenario.row$scenario_id)
  raw.rds <- file.path(scenario.outdir, paste0(scenario.id, "_raw_results.rds"))
  pointwise.rds <- file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.rds"))
  trial.rds <- file.path(scenario.outdir, paste0(scenario.id, "_trial_summary.rds"))
  
  existing.trial <- data.frame()
  existing.method <- data.frame()
  existing.pointwise <- data.frame()
  completed <- 0L
  
  if (isTRUE(resume) && file.exists(raw.rds)) {
    existing.method <- tryCatch(readRDS(raw.rds), error = function(e) data.frame())
    existing.pointwise <- if (file.exists(pointwise.rds)) tryCatch(readRDS(pointwise.rds), error = function(e) data.frame()) else data.frame()
    existing.trial <- if (file.exists(trial.rds)) tryCatch(readRDS(trial.rds), error = function(e) data.frame()) else data.frame()
    if (nrow(existing.method) > 0 && "sim_index" %in% names(existing.method)) {
      completed <- suppressWarnings(max(as.integer(existing.method$sim_index), na.rm = TRUE))
      if (!is.finite(completed)) completed <- 0L
    }
    if (isTRUE(verbose) && completed > 0L) cat("Resuming", scenario.id, "from simulated trial", completed + 1L, "\n")
  }
  
  if (completed >= nsim) {
    trial.rows <- existing.trial
    method.rows <- existing.method
    pointwise.rows <- existing.pointwise
  } else {
    nnew <- nsim - completed
    trial.list <- vector("list", nnew)
    method.list <- vector("list", nnew)
    pointwise.list <- vector("list", nnew)
    
    for (ii in seq_len(nnew)) {
      s <- completed + ii
      if (isTRUE(verbose) && (s == 1 || s == nsim || s %% checkpoint.every == 0)) {
        cat("  simulated trial", s, "of", nsim, "\n")
      }
      seed <- as.integer(master.seed + scenario.row$scenario_index * 1000000L + s * 104729L)
      
      res <- tryCatch(
        run.one.simulation.final(
          scenario.row = scenario.row,
          sim.index = s,
          seed = seed,
          B = B,
          save.example = (isTRUE(SAVE_EXAMPLE_PLOTS) && s == 1),
          scenario.outdir = scenario.outdir
        ),
        error = function(e) {
          message("  ERROR in scenario ", scenario.id, ", sim ", s, ": ", conditionMessage(e))
          err.method <- data.frame(
            scenario_id = scenario.id,
            scenario_label = scenario.label,
            scenario_type = scenario.row$scenario_type,
            description = scenario.row$description,
            sim_index = s,
            method = NA_character_,
            method_label = NA_character_,
            statistic_type = NA_character_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
          list(
            trial.row = data.frame(
              scenario_id = scenario.id,
              scenario_label = scenario.label,
              scenario_type = scenario.row$scenario_type,
              description = scenario.row$description,
              sim_index = s,
              error = conditionMessage(e),
              stringsAsFactors = FALSE
            ),
            method.rows = err.method,
            pointwise.rows = data.frame()
          )
        }
      )
      
      trial.list[[ii]] <- res$trial.row
      method.list[[ii]] <- res$method.rows
      pointwise.list[[ii]] <- res$pointwise.rows
      
      if (s %% checkpoint.every == 0 || s == nsim) {
        trial.tmp <- rbind_fill_base(existing.trial, trial.list[seq_len(ii)])
        method.tmp <- rbind_fill_base(existing.method, method.list[seq_len(ii)])
        point.tmp <- rbind_fill_base(existing.pointwise, pointwise.list[seq_len(ii)])
        safe_write_csv(trial.tmp, file.path(scenario.outdir, "checkpoint_trial_summary.csv"))
        safe_write_csv(method.tmp, file.path(scenario.outdir, "checkpoint_method_trial_results.csv"))
        safe_write_csv(point.tmp, file.path(scenario.outdir, "checkpoint_pointwise_results.csv"))
        safe_save_rds(trial.tmp, trial.rds)
        safe_save_rds(method.tmp, raw.rds)
        safe_save_rds(point.tmp, pointwise.rds)
      }
    }
    
    trial.rows <- rbind_fill_base(existing.trial, trial.list)
    method.rows <- rbind_fill_base(existing.method, method.list)
    pointwise.rows <- rbind_fill_base(existing.pointwise, pointwise.list)
  }
  
  power.summary <- summarise.method.power.final(method.rows, alpha = ALPHA)
  average.statistics <- make.average.statistics.final(method.rows, alpha = ALPHA)
  method.comparison <- make.pathway.method.comparison.final(power.summary)
  pointwise.summary <- summarise.pointwise.final(pointwise.rows, alpha = ALPHA)
  
  write.scenario.outputs.final(
    scenario.id = scenario.id,
    scenario.outdir = scenario.outdir,
    trial.rows = trial.rows,
    method.rows = method.rows,
    pointwise.rows = pointwise.rows,
    power.summary = power.summary,
    average.statistics = average.statistics,
    method.comparison = method.comparison,
    pointwise.summary = pointwise.summary
  )
  
  if (isTRUE(SAVE_EXAMPLE_PLOTS)) {
    plot.scenario.power.final(power.summary, scenario.outdir, scenario.id)
    plot.scenario.pvalues.final(method.rows, scenario.outdir, scenario.id)
    plot.scenario.max.fixed.final(method.rows, scenario.outdir, scenario.id)
  }
  
  list(
    trial = trial.rows,
    method = method.rows,
    pointwise = pointwise.rows,
    power = power.summary,
    averages = average.statistics,
    method.comparison = method.comparison,
    pointwise.summary = pointwise.summary
  )
}

# Global plots and runner.
plot.global.power.final <- function(power.table, outdir) {
  if (is.null(power.table) || nrow(power.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  d <- power.table[power.table$method %in% key.methods, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  scenarios <- unique(d$scenario_id)
  for (side in c("one_sided", "two_sided")) {
    mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios),
                  dimnames = list(method_label_final(key.methods), scenario_display_label_final(scenarios)))
    for (i in seq_len(nrow(d))) {
      rr <- match(d$method[i], key.methods)
      cc <- match(d$scenario_id[i], scenarios)
      mat[rr, cc] <- if (side == "one_sided") d$rejection_proportion_one_sided[i] else d$rejection_proportion_two_sided[i]
    }
    fname <- if (side == "one_sided") "GLOBAL_power_comparison_key_methods_one_sided.png" else "GLOBAL_power_comparison_key_methods_two_sided.png"
    ok <- safe_png_open(file.path(outdir, fname), width = 3200, height = 1700)
    if (ok) {
      old.par <- graphics::par(no.readonly = TRUE)
      graphics::par(mar = c(9, 5, 4, 1))
      graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                        ylab = paste0("Rejection proportion, alpha = ", ALPHA),
                        main = paste0("Power / rejection proportion: ", gsub("_", "-", side)),
                        legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.65), cex.names = 0.8)
      graphics::abline(h = ALPHA, lty = 2)
      graphics::par(old.par)
      safe_dev_off(ok)
    }
  }
  # Compatibility filename: use the one-sided plot as the original-style default.
  one.file <- file.path(outdir, "GLOBAL_power_comparison_key_methods_one_sided.png")
  compat.file <- file.path(outdir, "GLOBAL_power_comparison_key_methods.png")
  if (file.exists(one.file)) try(file.copy(one.file, compat.file, overwrite = TRUE), silent = TRUE)
  invisible(NULL)
}

plot.global.pvalue.final <- function(method.table, outdir) {
  if (is.null(method.table) || nrow(method.table) == 0) return(invisible(NULL))
  agg <- summarise.method.power.final(method.table, alpha = ALPHA)
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  d <- agg[agg$method %in% key.methods, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  ok <- safe_png_open(file.path(outdir, "GLOBAL_permutation_pvalue_comparison.png"), width = 2800, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(one_sided = finite_or_na(d$median_p_one_sided), two_sided = finite_or_na(d$median_p_two_sided))
  colnames(mat) <- method_label_final(d$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = "Median p-value", main = "Median p-values across all scenarios",
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

plot.global.fixed.final <- function(method.table, outdir) {
  if (is.null(method.table) || nrow(method.table) == 0) return(invisible(NULL))
  wr <- method.table[method.table$statistic_type == "Win ratio", , drop = FALSE]
  if (nrow(wr) == 0) return(invisible(NULL))
  agg <- summarise.method.power.final(wr, alpha = ALPHA)
  ok <- safe_png_open(file.path(outdir, "GLOBAL_fixed_selected_pvalue_comparison.png"), width = 2800, height = 1500)
  if (!ok) return(invisible(NULL))
  old.par <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old.par); safe_dev_off(ok) }, add = TRUE)
  mat <- rbind(
    max_one_sided = finite_or_na(agg$median_p_one_sided),
    fixed_one_sided = finite_or_na(agg$mean_fixed_p_one_sided),
    max_two_sided = finite_or_na(agg$median_p_two_sided),
    fixed_two_sided = finite_or_na(agg$mean_fixed_p_two_sided)
  )
  colnames(mat) <- method_label_final(agg$method)
  graphics::par(mar = c(10, 5, 4, 1))
  graphics::barplot(mat, beside = TRUE, las = 2, ylim = c(0, 1),
                    ylab = "p-value", main = "Max vs fixed-selected WR p-values",
                    legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.7), cex.names = 0.8)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(NULL)
}

run.all.scenarios.final <- function(scenario.grid = build.final.15.scenario.grid(),
                                    nsim = NSIM,
                                    B = B_PERM,
                                    outdir = OUTDIR,
                                    master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                                    resume = RESUME_IF_EXISTS) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  scenario.grid$scenario_index <- seq_len(nrow(scenario.grid))
  scenario.grid$scenario_label <- scenario_display_label_final(scenario.grid$scenario_id)
  safe_write_csv(scenario.grid, file.path(outdir, "GLOBAL_scenario_config.csv"))
  
  settings.table <- data.frame(
    setting = c(
      "NSIM", "B_PERM", "MASTER_SEED", "ALPHA", "number_of_scenarios",
      "WR_outputs", "WR_one_sided_statistic", "WR_two_sided_statistic",
      "WR_adaptive_one_sided_rule", "WR_adaptive_two_sided_rule", "WR_permutation_scheme",
      "logrank_outputs", "one_sided_logrank_direction",
      "p_grid_low_exploratory", "p_grid_primary", "p_grid_full", "clinical_t_months",
      "checkpoint_every", "resume_if_exists", "save_example_plots", "output_folder"
    ),
    value = c(
      as.character(nsim), as.character(B), as.character(master.seed), as.character(ALPHA), as.character(nrow(scenario.grid)),
      "one-sided and two-sided for every WR pathway",
      "WR; larger values favor treatment",
      "abs(log(WR)); right-tail permutation comparison on distance from the null value 1",
      "max WR, selection repeated inside each permutation",
      "max abs(log(WR)), with the same selection repeated inside each permutation",
      "treatment-label permutation",
      "one-sided and two-sided for death and composite endpoints",
      "treatment benefit, Cox HR < 1 for ARM = 1",
      if (exists("P_GRID_EXPLORATORY")) paste0(sprintf("%.2f", min(P_GRID_EXPLORATORY)), " to ", sprintf("%.2f", max(P_GRID_EXPLORATORY))) else NA_character_,
      if (exists("P_GRID_PRIMARY")) paste0(sprintf("%.2f", min(P_GRID_PRIMARY)), " to ", sprintf("%.2f", max(P_GRID_PRIMARY))) else NA_character_,
      if (exists("P_GRID_ALL")) paste0(sprintf("%.2f", min(P_GRID_ALL)), " to ", sprintf("%.2f", max(P_GRID_ALL)), "; no p = 0") else NA_character_,
      if (exists("CLINICAL_T_MONTHS")) paste(CLINICAL_T_MONTHS, collapse = ", ") else "1, 3, 6, 12, 18, 24",
      as.character(CHECKPOINT_EVERY), as.character(RESUME_IF_EXISTS), as.character(SAVE_EXAMPLE_PLOTS), outdir
    ),
    stringsAsFactors = FALSE
  )
  safe_write_csv(settings.table, file.path(outdir, "GLOBAL_settings.csv"))
  
  cat("\n===== FINAL Simulation settings =====\n")
  cat("Output folder:", outdir, "\n")
  cat("NSIM =", nsim, "; B_PERM =", B, "; scenarios =", nrow(scenario.grid), "\n")
  cat("WR outputs: one-sided and two-sided\n")
  cat("WR two-sided statistic: abs(log(WR)); adaptive rule: max abs(log(WR)) repeated inside each permutation\n")
  cat("Log-rank outputs: one-sided benefit and two-sided\n")
  
  all.trial <- data.frame()
  all.method <- data.frame()
  all.pointwise <- data.frame()
  all.power <- data.frame()
  all.averages <- data.frame()
  all.method.comparison <- data.frame()
  all.pointwise.summary <- data.frame()
  
  for (i in seq_len(nrow(scenario.grid))) {
    scenario.row <- scenario.grid[i, , drop = FALSE]
    out.i <- run.scenario.final(
      scenario.row = scenario.row,
      nsim = nsim,
      B = B,
      master.seed = master.seed,
      outdir = outdir,
      resume = resume,
      checkpoint.every = CHECKPOINT_EVERY,
      verbose = VERBOSE
    )
    
    all.trial <- rbind_fill_base(all.trial, out.i$trial)
    all.method <- rbind_fill_base(all.method, out.i$method)
    all.pointwise <- rbind_fill_base(all.pointwise, out.i$pointwise)
    all.power <- rbind_fill_base(all.power, out.i$power)
    all.averages <- rbind_fill_base(all.averages, out.i$averages)
    all.method.comparison <- rbind_fill_base(all.method.comparison, out.i$method.comparison)
    all.pointwise.summary <- rbind_fill_base(all.pointwise.summary, out.i$pointwise.summary)
    
    # New explicit global names
    safe_write_csv(all.trial, file.path(outdir, "GLOBAL_trial_summary.csv"))
    safe_write_csv(all.method, file.path(outdir, "GLOBAL_method_trial_results.csv"))
    safe_write_csv(all.pointwise, file.path(outdir, "GLOBAL_pointwise_results.csv"))
    safe_write_csv(all.power, file.path(outdir, "GLOBAL_power_summary.csv"))
    safe_write_csv(all.averages, file.path(outdir, "GLOBAL_average_statistics.csv"))
    safe_write_csv(all.method.comparison, file.path(outdir, "GLOBAL_pathway_method_comparison.csv"))
    safe_write_csv(all.pointwise.summary, file.path(outdir, "GLOBAL_pointwise_summary.csv"))
    
    # Old-style compatibility global names
    safe_write_csv(all.method, file.path(outdir, "GLOBAL_all_scenarios_raw_results.csv"))
    safe_write_csv(all.pointwise, file.path(outdir, "GLOBAL_all_scenarios_pointwise_results.csv"))
    safe_write_csv(all.power, file.path(outdir, "GLOBAL_all_scenarios_power_summary.csv"))
    safe_write_csv(all.averages, file.path(outdir, "GLOBAL_all_scenarios_average_statistics.csv"))
    safe_write_csv(all.method.comparison, file.path(outdir, "GLOBAL_all_scenarios_pathway_method_comparison.csv"))
    safe_write_csv(all.pointwise.summary, file.path(outdir, "GLOBAL_all_scenarios_pointwise_average_pvalue_tie_summary.csv"))
  }
  
  if (isTRUE(SAVE_EXAMPLE_PLOTS)) {
    plot.global.power.final(all.power, outdir)
    plot.global.pvalue.final(all.method, outdir)
    plot.global.fixed.final(all.method, outdir)
    # Compatibility placeholder for low-p plot: low-p summaries are in the CSVs; this figure focuses on selected p.
    try({
      wr <- all.method[all.method$statistic_type == "Win ratio", , drop = FALSE]
      if (nrow(wr) > 0) {
        agg <- summarise.method.power.final(wr, alpha = ALPHA)
        ok <- safe_png_open(file.path(outdir, "GLOBAL_low_p_selection_frequency.png"), width = 2600, height = 1500)
        if (ok) {
          old.par <- graphics::par(no.readonly = TRUE)
          graphics::par(mar = c(10, 5, 4, 1))
          vals <- finite_or_na(agg$mean_selected_p_one_sided)
          names(vals) <- method_label_final(agg$method)
          graphics::barplot(vals, las = 2, ylim = c(0, 1), ylab = "Mean selected p, one-sided", main = "Selected p summary")
          graphics::par(old.par)
          safe_dev_off(ok)
        }
      }
    }, silent = TRUE)
  }
  
  manifest <- data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE)
  safe_write_csv(manifest, file.path(outdir, "GLOBAL_output_manifest.csv"))
  
  full.bundle <- list(
    scenario.grid = scenario.grid,
    trial = all.trial,
    raw = all.method,
    method = all.method,
    pointwise = all.pointwise,
    power = all.power,
    averages = all.averages,
    method.comparison = all.method.comparison,
    pointwise.summary = all.pointwise.summary,
    settings = list(
      NSIM = nsim,
      B_PERM = B,
      MASTER_SEED = master.seed,
      ALPHA = ALPHA,
      P_GRID_EXPLORATORY = if (exists("P_GRID_EXPLORATORY")) P_GRID_EXPLORATORY else NULL,
      P_GRID_PRIMARY = if (exists("P_GRID_PRIMARY")) P_GRID_PRIMARY else NULL,
      P_GRID_ALL = if (exists("P_GRID_ALL")) P_GRID_ALL else NULL,
      CLINICAL_T_MONTHS = if (exists("CLINICAL_T_MONTHS")) CLINICAL_T_MONTHS else NULL,
      WR_sidedness = "one-sided and two-sided",
      WR_two_sided_statistic = "abs(log(WR)); right-tail permutation comparison on distance from the null value 1",
      logrank_sidedness = "one-sided benefit and two-sided"
    )
  )
  safe_save_rds(full.bundle, file.path(outdir, "GLOBAL_full_results_bundle.rds"))
  
  cat("\n===== All scenarios complete =====\n")
  cat("Output folder:", outdir, "\n")
  cat("Full-run settings: NSIM =", nsim, ", B_PERM =", B, "\n")
  cat("Scenario count:", nrow(scenario.grid), "\n")
  cat("Main tables:\n")
  cat("  GLOBAL_scenario_config.csv\n")
  cat("  GLOBAL_settings.csv\n")
  cat("  GLOBAL_method_trial_results.csv / GLOBAL_all_scenarios_raw_results.csv\n")
  cat("  GLOBAL_power_summary.csv / GLOBAL_all_scenarios_power_summary.csv\n")
  cat("  GLOBAL_average_statistics.csv / GLOBAL_all_scenarios_average_statistics.csv\n")
  cat("  GLOBAL_pathway_method_comparison.csv / GLOBAL_all_scenarios_pathway_method_comparison.csv\n")
  cat("  GLOBAL_pointwise_results.csv / GLOBAL_all_scenarios_pointwise_results.csv\n")
  cat("  GLOBAL_pointwise_summary.csv / GLOBAL_all_scenarios_pointwise_average_pvalue_tie_summary.csv\n")
  cat("  GLOBAL_full_results_bundle.rds\n")
  cat("Main figures:\n")
  cat("  GLOBAL_power_comparison_key_methods.png\n")
  cat("  GLOBAL_power_comparison_key_methods_one_sided.png\n")
  cat("  GLOBAL_power_comparison_key_methods_two_sided.png\n")
  cat("  GLOBAL_permutation_pvalue_comparison.png\n")
  cat("  GLOBAL_fixed_selected_pvalue_comparison.png\n")
  cat("  GLOBAL_low_p_selection_frequency.png\n")
  cat("Each scenario folder contains old-style Sxx_* files plus new explicit files and example_sim1_* tables.\n")
  
  invisible(c(full.bundle, list(manifest = manifest)))
}

run.all.scenarios <- function(scenario.grid = build.final.15.scenario.grid(),
                              nsim = NSIM,
                              B = B_PERM,
                              outdir = OUTDIR,
                              master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L) {
  run.all.scenarios.final(
    scenario.grid = scenario.grid,
    nsim = nsim,
    B = B,
    outdir = outdir,
    master.seed = master.seed,
    resume = RESUME_IF_EXISTS
  )
}


# Trial-level parallelization only.
# at most 4 worker R sessions, with each worker limited 

get_integer_env_final <- function(name, default) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) return(as.integer(default))
  out <- suppressWarnings(as.integer(val))
  if (is.na(out) || out < 1L) as.integer(default) else out
}

logical_env_final <- function(name, default = TRUE) {
  val <- Sys.getenv(name, unset = NA_character_)
  if (is.na(val) || !nzchar(val)) return(isTRUE(default))
  !(toupper(val) %in% c("FALSE", "F", "0", "NO", "N"))
}

physical_cores_final <- tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_)
if (is.na(physical_cores_final) || physical_cores_final < 1L) {
  physical_cores_final <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 4L)
}
if (is.na(physical_cores_final) || physical_cores_final < 1L) physical_cores_final <- 4L

PARALLEL_TRIALS <- logical_env_final("PARALLEL_TRIALS", TRUE)
SAFE_WORKER_CAP <- get_integer_env_final("SAFE_WORKER_CAP", 4L)
N_WORKERS <- get_integer_env_final(
  "N_WORKERS",
  max(1L, min(SAFE_WORKER_CAP, as.integer(physical_cores_final) - 2L))
)
N_WORKERS <- max(1L, min(N_WORKERS, as.integer(physical_cores_final)))
PAUSE_BETWEEN_BATCHES_SECONDS <- as.numeric(Sys.getenv("PAUSE_BETWEEN_BATCHES_SECONDS", unset = "0"))
if (!is.finite(PAUSE_BETWEEN_BATCHES_SECONDS) || PAUSE_BETWEEN_BATCHES_SECONDS < 0) {
  PAUSE_BETWEEN_BATCHES_SECONDS <- 0
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

CPP_FAST_WR_CODE_FOR_PARALLEL <- r"---(

// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

int count_hosp_until_cpp(const NumericVector& hosp_times,
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

double safe_wr_cpp(double win_score, double loss_score, double total_pairs, double eps) {
  return ((win_score / total_pairs) + eps) / ((loss_score / total_pairs) + eps);
}

// [[Rcpp::export]]
List fast_wr_core_revised_cpp(NumericVector futime,
                              IntegerVector cnsr,
                              IntegerVector arm,
                              NumericVector hosp_times,
                              IntegerVector hosp_start,
                              IntegerVector hosp_len,
                              NumericVector p_grid,
                              NumericVector t_grid,
                              double eps = 1e-8) {
  int n = futime.size();
  std::vector<int> trt;
  std::vector<int> ctrl;
  trt.reserve(n);
  ctrl.reserve(n);

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) trt.push_back(i);
    if (arm[i] == 0) ctrl.push_back(i);
  }

  int n_trt = trt.size();
  int n_ctrl = ctrl.size();
  double total_pairs = static_cast<double>(n_trt) * static_cast<double>(n_ctrl);
  if (total_pairs <= 0) stop("Need at least one treatment patient and one control patient.");

  int K = t_grid.size();
  int P = p_grid.size();

  double D1_win = 0.0, D1_loss = 0.0, H2_win = 0.0, H2_loss = 0.0;
  double H1_win = 0.0, H1_loss = 0.0, D2_win = 0.0, D2_loss = 0.0;
  double true_tie_pairs = 0.0;

  NumericVector t_wins(K), t_losses(K), pr_win_t(K), pr_loss_t(K), pr_tie_t(K), WRt(K);

  for (int ii = 0; ii < n_trt; ++ii) {
    int ti = trt[ii];
    for (int jj = 0; jj < n_ctrl; ++jj) {
      int cj = ctrl[jj];

      int death_sign = 0;
      double ft = futime[ti];
      double fc = futime[cj];
      int dt = cnsr[ti];
      int dc = cnsr[cj];

      if (dt == 1 && dc == 1) {
        if (ft > fc) death_sign = 1;
        else if (ft < fc) death_sign = -1;
      } else if (dt == 0 && dc == 1 && ft >= fc) {
        death_sign = 1;
      } else if (dt == 1 && dc == 0 && fc >= ft) {
        death_sign = -1;
      }

      double common_t = ft < fc ? ft : fc;
      int ht = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, ti, common_t);
      int hc = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, cj, common_t);
      int hosp_sign = 0;
      if (ht < hc) hosp_sign = 1;
      else if (ht > hc) hosp_sign = -1;

      if (death_sign > 0) D1_win += 1.0;
      else if (death_sign < 0) D1_loss += 1.0;
      else if (hosp_sign > 0) H2_win += 1.0;
      else if (hosp_sign < 0) H2_loss += 1.0;
      else true_tie_pairs += 1.0;

      if (hosp_sign > 0) H1_win += 1.0;
      else if (hosp_sign < 0) H1_loss += 1.0;
      else if (death_sign > 0) D2_win += 1.0;
      else if (death_sign < 0) D2_loss += 1.0;

      for (int kk = 0; kk < K; ++kk) {
        int sign_t = 0;
        if (death_sign != 0 && std::fabs(ft - fc) >= t_grid[kk]) {
          sign_t = death_sign;
        } else {
          sign_t = hosp_sign;
        }
        if (sign_t > 0) t_wins[kk] += 1.0;
        else if (sign_t < 0) t_losses[kk] += 1.0;
      }
    }
  }

  NumericVector WR_death_first(P), win_death_first(P), loss_death_first(P);
  NumericVector WR_hosp_first(P), win_hosp_first(P), loss_hosp_first(P);
  for (int pp = 0; pp < P; ++pp) {
    double p = p_grid[pp];
    win_death_first[pp] = p * D1_win + (1.0 - p) * H2_win;
    loss_death_first[pp] = p * D1_loss + (1.0 - p) * H2_loss;
    WR_death_first[pp] = safe_wr_cpp(win_death_first[pp], loss_death_first[pp], total_pairs, eps);

    win_hosp_first[pp] = p * H1_win + (1.0 - p) * D2_win;
    loss_hosp_first[pp] = p * H1_loss + (1.0 - p) * D2_loss;
    WR_hosp_first[pp] = safe_wr_cpp(win_hosp_first[pp], loss_hosp_first[pp], total_pairs, eps);
  }

  double ordinary_win_score = D1_win + H2_win;
  double ordinary_loss_score = D1_loss + H2_loss;
  double ordinaryWR = safe_wr_cpp(ordinary_win_score, ordinary_loss_score, total_pairs, eps);

  for (int kk = 0; kk < K; ++kk) {
    pr_win_t[kk] = t_wins[kk] / total_pairs;
    pr_loss_t[kk] = t_losses[kk] / total_pairs;
    pr_tie_t[kk] = 1.0 - pr_win_t[kk] - pr_loss_t[kk];
    WRt[kk] = safe_wr_cpp(t_wins[kk], t_losses[kk], total_pairs, eps);
  }

  return List::create(
    Named("ordinaryWR") = ordinaryWR,
    Named("ordinary_win_score") = ordinary_win_score,
    Named("ordinary_loss_score") = ordinary_loss_score,
    Named("total_pairs") = total_pairs,
    Named("D1_win") = D1_win,
    Named("D1_loss") = D1_loss,
    Named("H2_win") = H2_win,
    Named("H2_loss") = H2_loss,
    Named("H1_win") = H1_win,
    Named("H1_loss") = H1_loss,
    Named("D2_win") = D2_win,
    Named("D2_loss") = D2_loss,
    Named("true_tie_pairs") = true_tie_pairs,
    Named("p_grid") = p_grid,
    Named("WR_death_first") = WR_death_first,
    Named("win_death_first") = win_death_first,
    Named("loss_death_first") = loss_death_first,
    Named("WR_hosp_first") = WR_hosp_first,
    Named("win_hosp_first") = win_hosp_first,
    Named("loss_hosp_first") = loss_hosp_first,
    Named("t_grid") = t_grid,
    Named("WRt") = WRt,
    Named("t_wins") = t_wins,
    Named("t_losses") = t_losses,
    Named("pr_win_t") = pr_win_t,
    Named("pr_loss_t") = pr_loss_t,
    Named("pr_tie_t") = pr_tie_t
  );
}
)---"

parallel_export_names_final <- function() {
  exclude <- c(
    "fast_wr_core_revised_cpp", "out", "cl", "cluster", ".Random.seed",
    "scenario.grid.current"
  )
  setdiff(ls(envir = .GlobalEnv), exclude)
}

make.parallel.cluster.final <- function(n.workers = N_WORKERS, verbose = TRUE) {
  n.workers <- as.integer(n.workers)
  if (!isTRUE(PARALLEL_TRIALS) || n.workers <= 1L) return(NULL)
  
  cl <- tryCatch(
    parallel::makePSOCKcluster(n.workers),
    error = function(e) {
      warning("Could not start PSOCK cluster; falling back to serial calculation: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(cl)) return(NULL)
  
  ok <- tryCatch({
    parallel::clusterExport(cl, parallel_export_names_final(), envir = .GlobalEnv)
    parallel::clusterEvalQ(cl, {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
      suppressPackageStartupMessages({
        library(survival)
        library(Rcpp)
      })
      Rcpp::sourceCpp(code = CPP_FAST_WR_CODE_FOR_PARALLEL)
      NULL
    })
    TRUE
  }, error = function(e) {
    warning("Could not initialize worker sessions; falling back to serial calculation: ", conditionMessage(e))
    FALSE
  })
  
  if (!isTRUE(ok)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    return(NULL)
  }
  
  if (isTRUE(verbose)) {
    cat("  Parallel trials enabled with", n.workers, "worker R sessions; each worker uses 1 compute thread.
")
  }
  cl
}

safe.run.one.simulation.worker.final <- function(s,
                                                 scenario.row,
                                                 master.seed,
                                                 B,
                                                 scenario.outdir,
                                                 scenario.label) {
  seed <- as.integer(master.seed + scenario.row$scenario_index * 1000000L + s * 104729L)
  
  tryCatch(
    run.one.simulation.final(
      scenario.row = scenario.row,
      sim.index = s,
      seed = seed,
      B = B,
      save.example = (isTRUE(SAVE_EXAMPLE_PLOTS) && s == 1),
      scenario.outdir = scenario.outdir
    ),
    error = function(e) {
      err.method <- data.frame(
        scenario_id = as.character(scenario.row$scenario_id),
        scenario_label = scenario.label,
        scenario_type = scenario.row$scenario_type,
        description = scenario.row$description,
        sim_index = s,
        method = NA_character_,
        method_label = NA_character_,
        statistic_type = NA_character_,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      list(
        trial.row = data.frame(
          scenario_id = as.character(scenario.row$scenario_id),
          scenario_label = scenario.label,
          scenario_type = scenario.row$scenario_type,
          description = scenario.row$description,
          sim_index = s,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        ),
        method.rows = err.method,
        pointwise.rows = data.frame()
      )
    }
  )
}

read_existing_or_checkpoint_final <- function(rds.file, csv.file) {
  if (file.exists(rds.file)) {
    return(tryCatch(readRDS(rds.file), error = function(e) data.frame()))
  }
  if (file.exists(csv.file)) {
    return(tryCatch(utils::read.csv(csv.file, stringsAsFactors = FALSE), error = function(e) data.frame()))
  }
  data.frame()
}

# Override scenario runner: same outputs as before, but pending trials are run
# in parallel batches. Seeds are unchanged, so results are reproducible and do
# not depend on worker scheduling.
run.scenario.final <- function(scenario.row,
                               nsim = NSIM,
                               B = B_PERM,
                               master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                               outdir = OUTDIR,
                               resume = RESUME_IF_EXISTS,
                               checkpoint.every = CHECKPOINT_EVERY,
                               verbose = VERBOSE) {
  scenario.label <- scenario_display_label_final(scenario.row$scenario_id)
  scenario.outdir <- file.path(outdir, paste0(sprintf("%02d", scenario.row$scenario_index), "_", scenario.row$scenario_id))
  dir.create(scenario.outdir, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(verbose)) {
    cat("\n===== Scenario", scenario.row$scenario_index, ":", scenario.label, "=====\n")
    cat("HR =", scenario.row$HR,
        "; hosp scale =", scenario.row$evt.rate.scale.param.trt,
        "; FU =", scenario.row$max.FU,
        "; censor =", scenario.row$censor.rate, "\n")
    cat("NSIM =", nsim, "; B_PERM =", B, "\n")
  }
  
  scenario.id <- as.character(scenario.row$scenario_id)
  raw.rds <- file.path(scenario.outdir, paste0(scenario.id, "_raw_results.rds"))
  pointwise.rds <- file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.rds"))
  trial.rds <- file.path(scenario.outdir, paste0(scenario.id, "_trial_summary.rds"))
  
  existing.trial <- data.frame()
  existing.method <- data.frame()
  existing.pointwise <- data.frame()
  completed <- 0L
  
  if (isTRUE(resume)) {
    existing.method <- read_existing_or_checkpoint_final(raw.rds, file.path(scenario.outdir, "checkpoint_method_trial_results.csv"))
    existing.pointwise <- read_existing_or_checkpoint_final(pointwise.rds, file.path(scenario.outdir, "checkpoint_pointwise_results.csv"))
    existing.trial <- read_existing_or_checkpoint_final(trial.rds, file.path(scenario.outdir, "checkpoint_trial_summary.csv"))
    
    if (nrow(existing.method) > 0 && "sim_index" %in% names(existing.method)) {
      completed <- suppressWarnings(max(as.integer(existing.method$sim_index), na.rm = TRUE))
      if (!is.finite(completed)) completed <- 0L
    }
    completed <- min(as.integer(completed), as.integer(nsim))
    if (isTRUE(verbose) && completed > 0L) {
      cat("Resuming", scenario.id, "from simulated trial", completed + 1L, "\n")
    }
  }
  
  if (completed >= nsim) {
    trial.rows <- existing.trial
    method.rows <- existing.method
    pointwise.rows <- existing.pointwise
  } else {
    pending <- seq.int(completed + 1L, nsim)
    trial.list <- list()
    method.list <- list()
    pointwise.list <- list()
    
    cl <- make.parallel.cluster.final(n.workers = N_WORKERS, verbose = verbose)
    on.exit({ if (!is.null(cl)) try(parallel::stopCluster(cl), silent = TRUE) }, add = TRUE)
    
    batch.starts <- seq(1L, length(pending), by = checkpoint.every)
    for (bb in seq_along(batch.starts)) {
      batch.idx <- pending[batch.starts[bb]:min(length(pending), batch.starts[bb] + checkpoint.every - 1L)]
      if (isTRUE(verbose)) {
        mode.label <- if (!is.null(cl)) paste0("parallel, ", N_WORKERS, " workers") else "serial fallback"
        cat("  running trials", min(batch.idx), "to", max(batch.idx), "of", nsim, "(", mode.label, ")
")
      }
      
      res.batch <- if (!is.null(cl)) {
        parallel::parLapplyLB(
          cl,
          batch.idx,
          safe.run.one.simulation.worker.final,
          scenario.row = scenario.row,
          master.seed = master.seed,
          B = B,
          scenario.outdir = scenario.outdir,
          scenario.label = scenario.label
        )
      } else {
        lapply(
          batch.idx,
          safe.run.one.simulation.worker.final,
          scenario.row = scenario.row,
          master.seed = master.seed,
          B = B,
          scenario.outdir = scenario.outdir,
          scenario.label = scenario.label
        )
      }
      
      trial.list <- c(trial.list, lapply(res.batch, `[[`, "trial.row"))
      method.list <- c(method.list, lapply(res.batch, `[[`, "method.rows"))
      pointwise.list <- c(pointwise.list, lapply(res.batch, `[[`, "pointwise.rows"))
      
      trial.tmp <- rbind_fill_base(existing.trial, trial.list)
      method.tmp <- rbind_fill_base(existing.method, method.list)
      point.tmp <- rbind_fill_base(existing.pointwise, pointwise.list)
      
      safe_write_csv(trial.tmp, file.path(scenario.outdir, "checkpoint_trial_summary.csv"))
      safe_write_csv(method.tmp, file.path(scenario.outdir, "checkpoint_method_trial_results.csv"))
      safe_write_csv(point.tmp, file.path(scenario.outdir, "checkpoint_pointwise_results.csv"))
      safe_save_rds(trial.tmp, trial.rds)
      safe_save_rds(method.tmp, raw.rds)
      safe_save_rds(point.tmp, pointwise.rds)
      
      if (isTRUE(verbose)) {
        cat("  checkpoint written through trial", max(batch.idx), "of", nsim, "\n")
      }
      if (PAUSE_BETWEEN_BATCHES_SECONDS > 0) Sys.sleep(PAUSE_BETWEEN_BATCHES_SECONDS)
    }
    
    trial.rows <- rbind_fill_base(existing.trial, trial.list)
    method.rows <- rbind_fill_base(existing.method, method.list)
    pointwise.rows <- rbind_fill_base(existing.pointwise, pointwise.list)
  }
  
  power.summary <- summarise.method.power.final(method.rows, alpha = ALPHA)
  average.statistics <- make.average.statistics.final(method.rows, alpha = ALPHA)
  method.comparison <- make.pathway.method.comparison.final(power.summary)
  pointwise.summary <- summarise.pointwise.final(pointwise.rows, alpha = ALPHA)
  
  write.scenario.outputs.final(
    scenario.id = scenario.id,
    scenario.outdir = scenario.outdir,
    trial.rows = trial.rows,
    method.rows = method.rows,
    pointwise.rows = pointwise.rows,
    power.summary = power.summary,
    average.statistics = average.statistics,
    method.comparison = method.comparison,
    pointwise.summary = pointwise.summary
  )
  
  if (isTRUE(SAVE_EXAMPLE_PLOTS)) {
    plot.scenario.power.final(power.summary, scenario.outdir, scenario.id)
    plot.scenario.pvalues.final(method.rows, scenario.outdir, scenario.id)
    plot.scenario.max.fixed.final(method.rows, scenario.outdir, scenario.id)
  }
  
  list(
    trial = trial.rows,
    method = method.rows,
    pointwise = pointwise.rows,
    power = power.summary,
    averages = average.statistics,
    method.comparison = method.comparison,
    pointwise.summary = pointwise.summary
  )
}

# Add a small helper file inside the output folder so the parallel safety
write.parallel.settings.final <- function(outdir = OUTDIR) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  x <- data.frame(
    setting = c("PARALLEL_TRIALS", "N_WORKERS", "SAFE_WORKER_CAP", "physical_cores_detected",
                "worker_thread_limit", "parallel_level", "pause_between_batches_seconds"),
    value = c(as.character(PARALLEL_TRIALS), as.character(N_WORKERS), as.character(SAFE_WORKER_CAP),
              as.character(physical_cores_final), "1 BLAS/OpenMP thread per worker",
              "simulated trials parallelized; permutations remain inside each trial",
              as.character(PAUSE_BETWEEN_BATCHES_SECONDS)),
    stringsAsFactors = FALSE
  )
  safe_write_csv(x, file.path(outdir, "GLOBAL_parallel_settings.csv"))
  invisible(x)
}

# Override run.all.scenarios.final only to write the extra parallel-settings CSV
# before using the original full-output global runner body above through the
# existing function name. We preserve all global/scenario-level outputs.
run.all.scenarios.final.original.parallel_safe <- run.all.scenarios.final
run.all.scenarios.final <- function(scenario.grid = build.final.15.scenario.grid(),
                                    nsim = NSIM,
                                    B = B_PERM,
                                    outdir = OUTDIR,
                                    master.seed = if (exists("MASTER_SEED")) MASTER_SEED else 2026L,
                                    resume = RESUME_IF_EXISTS) {
  write.parallel.settings.final(outdir)
  run.all.scenarios.final.original.parallel_safe(
    scenario.grid = scenario.grid,
    nsim = nsim,
    B = B,
    outdir = outdir,
    master.seed = master.seed,
    resume = resume
  )
}

cat("\n===== SAFE parallel acceleration settings =====\n")
cat("PARALLEL_TRIALS =", PARALLEL_TRIALS, "; N_WORKERS =", N_WORKERS,
    "; physical cores detected =", physical_cores_final, "\n")
cat("Each worker is limited to 1 BLAS/OpenMP thread. Set N_WORKERS=1 or PARALLEL_TRIALS=FALSE to run serially.\n")

# run block
if (!exists("RUN_FULL_SIMULATION")) {
  RUN_FULL_SIMULATION <- !identical(Sys.getenv("RUN_FULL_SIMULATION"), "FALSE")
}

if (isTRUE(RUN_FULL_SIMULATION)) {
  scenario.grid.current <- build.final.15.scenario.grid()
  out <- run.all.scenarios.final(
    scenario.grid = scenario.grid.current,
    nsim = NSIM,
    B = B_PERM,
    outdir = OUTDIR,
    master.seed = MASTER_SEED,
    resume = RESUME_IF_EXISTS
  )
} else {
  message("RUN_FULL_SIMULATION is FALSE: functions were loaded, but the full simulation was not started.")
}

if (!exists("OUTDIR")) {
  OUTDIR <- "simulation_WR_15scenarios_WR_and_LR_one_two_sided_abs_log_1000x500_FULL_OUTPUT_PARALLEL_SAFE"
}
ALPHA <- if (exists("ALPHA")) ALPHA else 0.05

key.methods <- c(
  "ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR",
  "maxWRp_primary", "maxOrderWR_primary", "maxWRt",
  "logrank_death", "logrank_composite"
)

short_method_label <- function(method) {
  map <- c(
    ordinaryWR = "WR-D\np=0.50",
    traditionalWR_hosp_first = "WR-H\np=0.50",
    traditionalOrderWR = "Order WR\np=0.50",
    maxWRp_primary = "Max WR(p)\nD-first",
    maxWRp_low = "Max WR(p)\np<0.50",
    maxWRp_full = "Max WR(p)\nfull grid",
    maxOrderWR_primary = "Max order\nWR(p)",
    maxOrderWR_full = "Max order\nfull grid",
    maxWRt = "Max WR(t)",
    logrank_death = "LR death",
    logrank_composite = "LR composite"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out
}

scenario_short_label <- function(x) {
  x <- as.character(x)
  map <- c(
    S01_equivalence_equal_arms = "S01 Equal",
    S02_similarity_small_benefit_both = "S02 Similar +",
    S03_similarity_small_harm_both = "S03 Similar -",
    S04_superiority_moderate_better_both = "S04 Sup mod",
    S05_superiority_strong_better_both = "S05 Sup strong",
    S06_inferiority_moderate_worse_both = "S06 Inf mod",
    S07_inferiority_strong_worse_both = "S07 Inf strong",
    S08_death_benefit_only = "S08 Death +",
    S09_hosp_benefit_only = "S09 Hosp +",
    S10_death_benefit_hosp_harm = "S10 Death + / Hosp -",
    S11_death_harm_hosp_benefit = "S11 Death - / Hosp +",
    S12_high_random_censoring_superiority = "S12 High censor",
    S13_low_censoring_long_followup_superiority = "S13 Long FU",
    S14_strong_hosp_benefit_only = "S14 Strong hosp +",
    S15_weak_death_benefit_strong_hosp_benefit = "S15 Weak death + strong hosp +"
  )
  hit <- x %in% names(map)
  x[hit] <- unname(map[x[hit]])
  x
}

title_short <- function(x) {
  x <- as.character(x)[1]
  sx <- if (grepl("^S[0-9]{2}", x)) sub("^(S[0-9]{2}).*$", "\\1", x) else x
  sim <- if (grepl("example_sim[0-9]+", x)) sub(".*(example_sim[0-9]+).*", "\\1", x) else ""
  if (grepl("^S[0-9]{2}", x) && nzchar(sim)) return(paste(sx, sim))
  if (grepl("^S[0-9]{2}", x)) return(sx)
  x
}

finite_or_na <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- NA_real_
  x
}

read_csv_if_exists <- function(file) {
  if (length(file) == 0 || is.na(file) || !file.exists(file)) return(NULL)
  tryCatch(utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
}

barplot_png <- function(file, mat, ylab, main, legend = TRUE,
                        width = 3200, height = 1700, res = 160,
                        las = 1, mar = c(6, 5, 4, 1), legend.cex = 0.8,
                        cex.names = 0.9) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(file, width = width, height = height, res = res)
  old <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::par(old); grDevices::dev.off() }, add = TRUE)
  graphics::par(mar = mar)
  graphics::barplot(mat, beside = TRUE, las = las, ylim = c(0, 1),
                    ylab = ylab, main = main,
                    legend.text = if (legend) TRUE else NULL,
                    args.legend = if (legend) list(x = "topright", bty = "n", cex = legend.cex) else NULL,
                    cex.names = cex.names)
  graphics::abline(h = ALPHA, lty = 2)
  invisible(file)
}

summarise_method_rows <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(data.frame())
  nums <- c(
    "permutation_p_value_one_sided", "permutation_p_value_two_sided",
    "fixed_parameter_p_value_one_sided", "fixed_parameter_p_value_two_sided",
    "rejected_one_sided", "rejected_two_sided",
    "selected_p_one_sided", "selected_p_two_sided"
  )
  for (nm in nums) if (!nm %in% names(d)) d[[nm]] <- NA_real_
  d$rejected_one_sided <- as.numeric(d$rejected_one_sided)
  d$rejected_two_sided <- as.numeric(d$rejected_two_sided)
  split_keys <- interaction(d$scenario_id, d$method, drop = TRUE, lex.order = TRUE)
  pieces <- lapply(split(d, split_keys), function(x) {
    data.frame(
      scenario_id = x$scenario_id[1],
      scenario_label = if ("scenario_label" %in% names(x)) x$scenario_label[1] else scenario_short_label(x$scenario_id[1]),
      method = x$method[1],
      method_label = short_method_label(x$method[1]),
      statistic_type = if ("statistic_type" %in% names(x)) x$statistic_type[1] else NA_character_,
      median_p_one_sided = suppressWarnings(stats::median(finite_or_na(x$permutation_p_value_one_sided), na.rm = TRUE)),
      median_p_two_sided = suppressWarnings(stats::median(finite_or_na(x$permutation_p_value_two_sided), na.rm = TRUE)),
      mean_fixed_p_one_sided = suppressWarnings(mean(finite_or_na(x$fixed_parameter_p_value_one_sided), na.rm = TRUE)),
      mean_fixed_p_two_sided = suppressWarnings(mean(finite_or_na(x$fixed_parameter_p_value_two_sided), na.rm = TRUE)),
      rejection_proportion_one_sided = suppressWarnings(mean(finite_or_na(x$rejected_one_sided), na.rm = TRUE)),
      rejection_proportion_two_sided = suppressWarnings(mean(finite_or_na(x$rejected_two_sided), na.rm = TRUE)),
      mean_selected_p_one_sided = suppressWarnings(mean(finite_or_na(x$selected_p_one_sided), na.rm = TRUE)),
      mean_selected_p_two_sided = suppressWarnings(mean(finite_or_na(x$selected_p_two_sided), na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

plot_example_pvalues <- function(pvalue.file) {
  d <- read_csv_if_exists(pvalue.file)
  if (is.null(d) || nrow(d) == 0 || !"method" %in% names(d)) return(invisible(FALSE))
  d <- d[d$method %in% key.methods, , drop = FALSE]
  if (nrow(d) == 0) return(invisible(FALSE))
  one <- if ("max_permutation_p_value_one_sided" %in% names(d)) d$max_permutation_p_value_one_sided else d$permutation_p_value_one_sided
  two <- if ("max_permutation_p_value_two_sided" %in% names(d)) d$max_permutation_p_value_two_sided else d$permutation_p_value_two_sided
  mat <- rbind(one_sided = finite_or_na(one), two_sided = finite_or_na(two))
  colnames(mat) <- short_method_label(d$method)
  prefix <- sub("_pvalue_table\\.csv$", "", basename(pvalue.file))
  out <- file.path(dirname(pvalue.file), paste0(prefix, "_pvalue_key_methods.png"))
  barplot_png(out, mat, ylab = "p-value", main = paste0(title_short(prefix), ": one- vs two-sided p-values"))
  invisible(TRUE)
}

plot_scenario_from_tables <- function(scen.dir) {
  pwr <- read_csv_if_exists(file.path(scen.dir, "power_summary.csv"))
  raw <- read_csv_if_exists(file.path(scen.dir, "method_trial_results.csv"))
  if (is.null(raw)) raw <- read_csv_if_exists(list.files(scen.dir, pattern = "_raw_results\\.csv$", full.names = TRUE)[1])
  sid <- sub("^[0-9]+_", "", basename(scen.dir))
  
  if (!is.null(pwr) && nrow(pwr) > 0 && "method" %in% names(pwr)) {
    d <- pwr[pwr$method %in% key.methods, , drop = FALSE]
    if (nrow(d) > 0) {
      mat <- rbind(one_sided = finite_or_na(d$rejection_proportion_one_sided),
                   two_sided = finite_or_na(d$rejection_proportion_two_sided))
      colnames(mat) <- short_method_label(d$method)
      barplot_png(file.path(scen.dir, paste0(sid, "_power_comparison.png")), mat,
                  ylab = paste0("Rejection proportion, alpha = ", ALPHA),
                  main = paste0(title_short(sid), ": one- vs two-sided power"))
    }
  }
  
  if (!is.null(raw) && nrow(raw) > 0 && "method" %in% names(raw)) {
    agg <- summarise_method_rows(raw)
    agg <- agg[agg$method %in% key.methods, , drop = FALSE]
    if (nrow(agg) > 0) {
      mat <- rbind(one_sided = finite_or_na(agg$median_p_one_sided),
                   two_sided = finite_or_na(agg$median_p_two_sided))
      colnames(mat) <- short_method_label(agg$method)
      barplot_png(file.path(scen.dir, paste0(sid, "_permutation_pvalue_comparison.png")), mat,
                  ylab = "Median permutation/log-rank p-value",
                  main = paste0(title_short(sid), ": median p-values"))
      
      wr <- agg[agg$statistic_type == "Win ratio" | is.na(agg$statistic_type), , drop = FALSE]
      if (nrow(wr) > 0) {
        mat <- rbind(max_one_sided = finite_or_na(wr$median_p_one_sided),
                     fixed_one_sided = finite_or_na(wr$mean_fixed_p_one_sided),
                     max_two_sided = finite_or_na(wr$median_p_two_sided),
                     fixed_two_sided = finite_or_na(wr$mean_fixed_p_two_sided))
        colnames(mat) <- short_method_label(wr$method)
        barplot_png(file.path(scen.dir, paste0(sid, "_max_vs_fixed_selected_pvalues.png")), mat,
                    ylab = "p-value", main = paste0(title_short(sid), ": max vs fixed selected p-values"),
                    legend.cex = 0.75)
      }
    }
  }
  
  ex <- list.files(scen.dir, pattern = "_example_sim[0-9]+_pvalue_table\\.csv$", full.names = TRUE)
  if (length(ex) > 0) invisible(vapply(ex, plot_example_pvalues, logical(1)))
  invisible(TRUE)
}

plot_global_from_tables <- function(outdir = OUTDIR) {
  pwr <- read_csv_if_exists(file.path(outdir, "GLOBAL_power_summary.csv"))
  raw <- read_csv_if_exists(file.path(outdir, "GLOBAL_method_trial_results.csv"))
  if (is.null(raw)) raw <- read_csv_if_exists(file.path(outdir, "GLOBAL_all_scenarios_raw_results.csv"))
  
  if (!is.null(pwr) && nrow(pwr) > 0) {
    d <- pwr[pwr$method %in% key.methods, , drop = FALSE]
    scenarios <- unique(d$scenario_id)
    for (side in c("one_sided", "two_sided")) {
      mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios),
                    dimnames = list(short_method_label(key.methods), scenario_short_label(scenarios)))
      for (i in seq_len(nrow(d))) {
        rr <- match(d$method[i], key.methods); cc <- match(d$scenario_id[i], scenarios)
        mat[rr, cc] <- if (side == "one_sided") d$rejection_proportion_one_sided[i] else d$rejection_proportion_two_sided[i]
      }
      fname <- if (side == "one_sided") "GLOBAL_power_comparison_key_methods_one_sided.png" else "GLOBAL_power_comparison_key_methods_two_sided.png"
      barplot_png(file.path(outdir, fname), mat,
                  ylab = paste0("Rejection proportion, alpha = ", ALPHA),
                  main = paste0("Power / rejection proportion: ", gsub("_", "-", side)),
                  width = 3400, height = 1800, las = 2, mar = c(7, 5, 4, 1), legend.cex = 0.65, cex.names = 0.8)
    }
    one.file <- file.path(outdir, "GLOBAL_power_comparison_key_methods_one_sided.png")
    compat.file <- file.path(outdir, "GLOBAL_power_comparison_key_methods.png")
    if (file.exists(one.file)) try(file.copy(one.file, compat.file, overwrite = TRUE), silent = TRUE)
  }
  
  if (!is.null(raw) && nrow(raw) > 0) {
    agg <- summarise_method_rows(raw)
    agg <- agg[agg$method %in% key.methods, , drop = FALSE]
    if (nrow(agg) > 0) {
      mat <- rbind(one_sided = finite_or_na(agg$median_p_one_sided), two_sided = finite_or_na(agg$median_p_two_sided))
      colnames(mat) <- short_method_label(agg$method)
      barplot_png(file.path(outdir, "GLOBAL_permutation_pvalue_comparison.png"), mat,
                  ylab = "Median p-value", main = "Median p-values across all scenarios",
                  width = 3400, height = 1700)
      
      wr <- agg[agg$statistic_type == "Win ratio" | is.na(agg$statistic_type), , drop = FALSE]
      mat <- rbind(max_one_sided = finite_or_na(wr$median_p_one_sided),
                   fixed_one_sided = finite_or_na(wr$mean_fixed_p_one_sided),
                   max_two_sided = finite_or_na(wr$median_p_two_sided),
                   fixed_two_sided = finite_or_na(wr$mean_fixed_p_two_sided))
      colnames(mat) <- short_method_label(wr$method)
      barplot_png(file.path(outdir, "GLOBAL_fixed_selected_pvalue_comparison.png"), mat,
                  ylab = "p-value", main = "Max vs fixed-selected WR p-values",
                  width = 3400, height = 1700, legend.cex = 0.7)
    }
  }
  invisible(TRUE)
}

replot_all_short_labels <- function(outdir = OUTDIR) {
  if (!dir.exists(outdir)) stop("Output folder does not exist: ", outdir)
  scen.dirs <- list.dirs(outdir, recursive = FALSE, full.names = TRUE)
  scen.dirs <- scen.dirs[grepl("/[0-9]{2}_S[0-9]{2}_", gsub("\\\\", "/", scen.dirs))]
  for (d in scen.dirs) {
    message("Replotting scenario figures: ", basename(d))
    try(plot_scenario_from_tables(d), silent = TRUE)
  }
  message("Replotting global figures")
  try(plot_global_from_tables(outdir), silent = TRUE)
  message("Done. Numeric CSV/RDS results were not changed; PNG figures were overwritten/created.")
  invisible(TRUE)
}

# Run automatically when sourced, unless AUTO_REPLOT=FALSE was set.
if (!identical(toupper(Sys.getenv("AUTO_REPLOT", "TRUE")), "FALSE")) {
  replot_all_short_labels(OUTDIR)
}
