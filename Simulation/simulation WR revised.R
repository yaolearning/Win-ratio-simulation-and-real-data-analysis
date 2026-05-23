#Packages and global settings
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_if_missing("survival")
install_if_missing("Rcpp")

library(survival)
library(Rcpp)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

OUTDIR <- "simulation WR revised max outputs"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Main simulation controls for the planned full run.
# The default full setting is NSIM = 100 simulated trials per scenario
# and B_PERM = 500 treatment-label permutations per simulated trial.
# The script is checkpointed by scenario, so it can be restarted.
# The following controls can be overridden before sourcing the file.
# Example quick test before source(): NSIM <- 2; B_PERM <- 20
if (!exists("NSIM")) NSIM <- 100
if (!exists("B_PERM")) B_PERM <- 500
if (!exists("MASTER_SEED")) MASTER_SEED <- 2026
if (!exists("ALPHA")) ALPHA <- 0.05
if (!exists("CHECKPOINT_EVERY")) CHECKPOINT_EVERY <- 5
if (!exists("RESUME_IF_EXISTS")) RESUME_IF_EXISTS <- TRUE
if (!exists("VERBOSE")) VERBOSE <- TRUE
if (!exists("SAVE_EXAMPLE_PLOTS")) SAVE_EXAMPLE_PLOTS <- TRUE

# Weighted p grids.
# p is the weight assigned to the first endpoint in the selected hierarchy.
# For death-first order: p weights death/survival; 1-p weights hospitalization.
# For hospitalization-first order: p weights hospitalization; 1-p weights death.
# No p = 0 value is used. Values below 0.5 are exploratory only;
# the primary weighted analysis remains p in [0.5, 1].
P_GRID_EXPLORATORY <- seq(0.01, 0.49, by = 0.01)
P_GRID_LOW <- P_GRID_EXPLORATORY
P_GRID_PRIMARY <- seq(0.50, 1.00, by = 0.01)
P_GRID_ALL <- sort(unique(c(P_GRID_EXPLORATORY, P_GRID_PRIMARY)))

# Clinically interpretable backup candidates for a time-difference threshold t.
# These are not claimed to be universal clinical margins. They are candidate
# values for sensitivity analysis and should be justified for the disease area.
CLINICAL_T_MONTHS <- c(1, 3, 6, 12, 18, 24)
MAX_T_GRID_SIZE <- 10

EPS_WR <- 1e-8


#Data simulation
simulate.one.dataset <-
  function(N = 50,                           # number of patients per arm;
                                             # can be a vector, e.g., N=c(50,100), which
                                             # means 50 patients in control arm
                                             # and 100 patients in treatment arm
           mort.rate.ctrl = -log(0.6),       # mortality rate in control arm
                                             # (1-yr survival rate=60% in control arm)
           mort.rate.trt = -log(0.6) * 0.6,  # mortality rate in treatment arm
                                             # (HR=0.6)
           evt.rate.shape.param = 5,
           evt.rate.scale.param.ctr = 1,     # avg. freq. of hospitalization in
                                             # control arm = evt.rate.shape.param
                                             # * evt.rate.scale.param.ctr
                                             # (5 event/yr)
           evt.rate.scale.param.trt = 1/2,   # avg. freq. of hospitalization in
                                             # treatment arm =
                                             # evt.rate.shape.param *
                                             # evt.rate.scale.param.trt
                                             # (2.5 event/yr)
           max.followup = 1)                 # maximum follow-up = 1 years
  {
    if (length(N) == 1) N <- rep(N, 2)
    rr <- N[2] / N[1]  # randomization ratio; kept for consistency with old code

    # Subject.ID, Arm, Followup.Time, Censor.Ind, Survival.Time,
    # Censor.Time, Freq.of.Hosp, No.of.Hosp
    table.output <- cbind(
      1:sum(N),
      rep(0:1, N),
      rep(max.followup, sum(N)),
      rep(0, sum(N)),
      c(rexp(N[1], rate = mort.rate.ctrl),
        rexp(N[2], rate = mort.rate.trt)),
      rep(Inf, sum(N)),
      c(rgamma(N[1], shape = evt.rate.shape.param,
               scale = evt.rate.scale.param.ctr),
        rgamma(N[2], shape = evt.rate.shape.param,
               scale = evt.rate.scale.param.trt)),
      rep(0, sum(N))
    )

    table.output <- data.frame(table.output)
    colnames(table.output) <- c(
      "SUBJID", "ARM", "FUTIME", "CNSR", "SURVTIME",
      "CNSRTIME", "FREQHOSP", "NUMHOSP"
    )

    for (i in seq_len(sum(N))) {
      table.output$CNSR[i] <-
        ifelse(table.output$SURVTIME[i] < table.output$FUTIME[i] &&
                 table.output$SURVTIME[i] < table.output$CNSRTIME[i], 1, 0)
      table.output$FUTIME[i] <-
        min(table.output$FUTIME[i],
            table.output$SURVTIME[i],
            table.output$CNSRTIME[i])
    }

    hosp.times.list <- list()
    for (i in seq_len(sum(N))) {
      hosp.times.list[[i]] <- NA
      cum.time <- 0
      followup.time <- table.output$FUTIME[i]
      hosp.rate <- table.output$FREQHOSP[i]

      while (TRUE) {
        new.hosp.time <- rexp(1, rate = hosp.rate)
        cum.time <- cum.time + new.hosp.time
        if (cum.time < followup.time) {
          hosp.times.list[[i]] <- c(hosp.times.list[[i]], new.hosp.time)
        } else {
          break
        }
      }
    }

    for (i in seq_len(sum(N))) {
      if (length(hosp.times.list[[i]]) > 1) {
        hosp.times.list[[i]] <- hosp.times.list[[i]][-1]
        table.output$NUMHOSP[i] <- length(hosp.times.list[[i]])
      } else {
        hosp.times.list[[i]] <- numeric(0)
        table.output$NUMHOSP[i] <- 0
      }
    }

    list(table.output = table.output, hosp.times.list = hosp.times.list)
  }

# Optional wrapper for random censoring scenarios.
# The core simulation function above is intentionally not changed.
apply.random.censoring <- function(ds, censor.rate = 0, seed = NULL) {
  if (is.null(censor.rate) || is.na(censor.rate) || censor.rate <= 0) return(ds)
  if (!is.null(seed)) set.seed(seed)

  tab <- ds$table.output
  n <- nrow(tab)
  ctime <- rexp(n, rate = censor.rate)

  for (i in seq_len(n)) {
    if (ctime[i] < tab$FUTIME[i]) {
      tab$FUTIME[i] <- ctime[i]
      tab$CNSRTIME[i] <- ctime[i]
      tab$CNSR[i] <- 0

      x <- ds$hosp.times.list[[i]]
      if (length(x) == 0 || all(is.na(x))) {
        ds$hosp.times.list[[i]] <- numeric(0)
        tab$NUMHOSP[i] <- 0
      } else {
        abs.x <- cumsum(as.numeric(x[!is.na(x)]))
        keep <- abs.x <= tab$FUTIME[i]
        abs.keep <- abs.x[keep]
        if (length(abs.keep) == 0) {
          ds$hosp.times.list[[i]] <- numeric(0)
          tab$NUMHOSP[i] <- 0
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

#Preparation 
prepare.ds.fast <- function(ds) {
  # Convert hospitalization gap times to absolute hospitalization times.
  # This avoids repeated cumsum() inside pairwise loops.
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
  ds$hosp.start <- starts  # zero-based starts for C++
  ds$hosp.len <- lens
  ds
}

safe_mean <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

safe_wr <- function(win.score, loss.score, total.pairs, eps = EPS_WR) {
  ((win.score / total.pairs) + eps) / ((loss.score / total.pairs) + eps)
}

extract_nearest <- function(x, value) {
  which.min(abs(x - value))
}

# Fast Rcpp for matrix
Rcpp::sourceCpp(code = '
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
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
  if (total_pairs <= 0) {
    stop("Need at least one treatment patient and one control patient.");
  }

  int K = t_grid.size();
  int P = p_grid.size();

  // Death-first hierarchy counts.
  double D1_win = 0.0;
  double D1_loss = 0.0;
  double H2_win = 0.0;
  double H2_loss = 0.0;

  // Hospitalization-first hierarchy counts.
  double H1_win = 0.0;
  double H1_loss = 0.0;
  double D2_win = 0.0;
  double D2_loss = 0.0;

  // Endpoint signs and true unresolved ties.
  double death_win_pairs = 0.0;
  double death_loss_pairs = 0.0;
  double hosp_win_pairs = 0.0;
  double hosp_loss_pairs = 0.0;
  double true_tie_pairs = 0.0;

  std::vector<double> win_diff_t(K + 1, 0.0);
  std::vector<double> loss_diff_t(K + 1, 0.0);

  auto add_range = [&](std::vector<double>& diff, int start, int end, double value) {
    if (start > end) return;
    diff[start] += value;
    if (end + 1 < static_cast<int>(diff.size())) diff[end + 1] -= value;
  };

  for (int ii = 0; ii < n_trt; ++ii) {
    int pt1 = trt[ii];
    double fu1 = futime[pt1];
    int c1_original = cnsr[pt1];

    for (int jj = 0; jj < n_ctrl; ++jj) {
      int pt2 = ctrl[jj];
      double fu2 = futime[pt2];
      int c2_original = cnsr[pt2];
      double min_fu = fu1 < fu2 ? fu1 : fu2;

      int hosp1 = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, pt1, min_fu);
      int hosp2 = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, pt2, min_fu);

      int hosp_sign = 0;
      if (hosp2 > hosp1) hosp_sign = 1;       // treatment has fewer hospitalizations
      else if (hosp2 < hosp1) hosp_sign = -1; // treatment has more hospitalizations

      int c1 = c1_original;
      int c2 = c2_original;
      if (fu1 < fu2) c2 = 0;
      if (fu2 < fu1) c1 = 0;

      int death_sign = 0;
      if (c1 == 0 && c2 == 1) death_sign = 1;       // treatment survives longer
      else if (c1 == 1 && c2 == 0) death_sign = -1; // treatment dies earlier

      if (death_sign > 0) death_win_pairs += 1.0;
      if (death_sign < 0) death_loss_pairs += 1.0;
      if (hosp_sign > 0) hosp_win_pairs += 1.0;
      if (hosp_sign < 0) hosp_loss_pairs += 1.0;
      if (death_sign == 0 && hosp_sign == 0) true_tie_pairs += 1.0;

      // Death-first hierarchy.
      if (death_sign > 0) {
        D1_win += 1.0;
      } else if (death_sign < 0) {
        D1_loss += 1.0;
      } else {
        if (hosp_sign > 0) H2_win += 1.0;
        else if (hosp_sign < 0) H2_loss += 1.0;
      }

      // Hospitalization-first hierarchy.
      if (hosp_sign > 0) {
        H1_win += 1.0;
      } else if (hosp_sign < 0) {
        H1_loss += 1.0;
      } else {
        if (death_sign > 0) D2_win += 1.0;
        else if (death_sign < 0) D2_loss += 1.0;
      }

      // Threshold pathway: endpoint 1 is survival/death time difference.
      // If |survival difference| > t, use survival difference;
      // otherwise use hospitalization count.
      double diff_surv = fu1 - fu2;
      double abs_diff = diff_surv >= 0.0 ? diff_surv : -diff_surv;

      int k = 0;
      while (k < K && t_grid[k] < abs_diff) {
        k++;
      }

      if (k > 0) {
        if (diff_surv > 0.0) {
          add_range(win_diff_t, 0, k - 1, 1.0);
        } else if (diff_surv < 0.0) {
          add_range(loss_diff_t, 0, k - 1, 1.0);
        }
      }

      if (k < K) {
        if (hosp_sign > 0) {
          add_range(win_diff_t, k, K - 1, 1.0);
        } else if (hosp_sign < 0) {
          add_range(loss_diff_t, k, K - 1, 1.0);
        }
      }
    }
  }

  NumericVector WR_death_first(P), win_death_first(P), loss_death_first(P);
  NumericVector WR_hosp_first(P), win_hosp_first(P), loss_hosp_first(P);

  for (int a = 0; a < P; ++a) {
    double p = p_grid[a];

    win_death_first[a] = p * D1_win + (1.0 - p) * H2_win;
    loss_death_first[a] = p * D1_loss + (1.0 - p) * H2_loss;
    WR_death_first[a] = safe_wr_cpp(win_death_first[a], loss_death_first[a], total_pairs, eps);

    win_hosp_first[a] = p * H1_win + (1.0 - p) * D2_win;
    loss_hosp_first[a] = p * H1_loss + (1.0 - p) * D2_loss;
    WR_hosp_first[a] = safe_wr_cpp(win_hosp_first[a], loss_hosp_first[a], total_pairs, eps);
  }

  double ordinary_win_score = 0.5 * D1_win + 0.5 * H2_win;
  double ordinary_loss_score = 0.5 * D1_loss + 0.5 * H2_loss;
  double ordinaryWR = safe_wr_cpp(ordinary_win_score, ordinary_loss_score, total_pairs, eps);

  NumericVector t_wins(K), t_losses(K), WRt(K), pr_win_t(K), pr_loss_t(K), pr_tie_t(K);
  double cur_win = 0.0;
  double cur_loss = 0.0;
  for (int kk = 0; kk < K; ++kk) {
    cur_win += win_diff_t[kk];
    cur_loss += loss_diff_t[kk];
    t_wins[kk] = cur_win;
    t_losses[kk] = cur_loss;
    pr_win_t[kk] = cur_win / total_pairs;
    pr_loss_t[kk] = cur_loss / total_pairs;
    pr_tie_t[kk] = 1.0 - pr_win_t[kk] - pr_loss_t[kk];
    WRt[kk] = safe_wr_cpp(cur_win, cur_loss, total_pairs, eps);
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
    Named("death_win_pairs") = death_win_pairs,
    Named("death_loss_pairs") = death_loss_pairs,
    Named("hosp_win_pairs") = hosp_win_pairs,
    Named("hosp_loss_pairs") = hosp_loss_pairs,
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

max_from_curve <- function(curve, method, lower, upper) {
  ix <- which(curve$p >= lower & curve$p <= upper)
  if (length(ix) == 0) {
    return(data.frame(
      method = method,
      order = NA_character_,
      p = NA_real_,
      WR = NA_real_,
      win.score = NA_real_,
      loss.score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  k <- ix[which.max(curve$WR[ix])]
  data.frame(
    method = method,
    order = as.character(curve$order[k]),
    p = curve$p[k],
    WR = curve$WR[k],
    win.score = curve$win.score[k],
    loss.score = curve$loss.score[k],
    stringsAsFactors = FALSE
  )
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

  weighted.death.first$win.pairs <- ifelse(
    weighted.death.first$p >= 1,
    D1.win,
    D1.win + H2.win
  )
  weighted.death.first$loss.pairs <- ifelse(
    weighted.death.first$p >= 1,
    D1.loss,
    D1.loss + H2.loss
  )
  weighted.death.first$tie.count <- total.pairs -
    weighted.death.first$win.pairs - weighted.death.first$loss.pairs
  weighted.death.first$tie.pr <- weighted.death.first$tie.count / total.pairs

  weighted.hosp.first <- data.frame(
    order = "hospitalization_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_hosp_first),
    win.score = as.numeric(core$win_hosp_first),
    loss.score = as.numeric(core$loss_hosp_first),
    stringsAsFactors = FALSE
  )

  weighted.hosp.first$win.pairs <- ifelse(
    weighted.hosp.first$p >= 1,
    H1.win,
    H1.win + D2.win
  )
  weighted.hosp.first$loss.pairs <- ifelse(
    weighted.hosp.first$p >= 1,
    H1.loss,
    H1.loss + D2.loss
  )
  weighted.hosp.first$tie.count <- total.pairs -
    weighted.hosp.first$win.pairs - weighted.hosp.first$loss.pairs
  weighted.hosp.first$tie.pr <- weighted.hosp.first$tie.count / total.pairs

  weighted.all.orders <- rbind(weighted.death.first, weighted.hosp.first)

  threshold <- data.frame(
    t = as.numeric(core$t_grid),
    t.months = as.numeric(core$t_grid) * 12,
    WR = as.numeric(core$WRt),
    wins = as.numeric(core$t_wins),
    losses = as.numeric(core$t_losses),
    pr.win = as.numeric(core$pr_win_t),
    pr.loss = as.numeric(core$pr_loss_t),
    pr.tie = as.numeric(core$pr_tie_t),
    stringsAsFactors = FALSE
  )
  threshold$tie.count <- total.pairs - threshold$wins - threshold$losses
  best.t.idx <- which.max(threshold$WR)

  p.low.values <- p.grid[p.grid < 0.5]
  p.primary.values <- p.grid[p.grid >= 0.5]
  p.low.lower <- if (length(p.low.values) > 0) min(p.low.values) else NA_real_
  p.low.upper <- if (length(p.low.values) > 0) max(p.low.values) else NA_real_
  p.primary.lower <- if (length(p.primary.values) > 0) min(p.primary.values) else NA_real_
  p.primary.upper <- if (length(p.primary.values) > 0) max(p.primary.values) else NA_real_

  list(
    ordinaryWR = as.numeric(core$ordinaryWR),
    ordinary.win.score = as.numeric(core$ordinary_win_score),
    ordinary.loss.score = as.numeric(core$ordinary_loss_score),

    weighted.death.first = weighted.death.first,
    weighted.hosp.first = weighted.hosp.first,
    weighted.all.orders = weighted.all.orders,

    max.weighted.primary = max_from_curve(weighted.death.first, "maxWRp_primary_death_first", p.primary.lower, p.primary.upper),
    max.weighted.low = max_from_curve(weighted.death.first, "maxWRp_low_death_first", p.low.lower, p.low.upper),
    max.weighted.full = max_from_curve(weighted.death.first, "maxWRp_full_death_first", min(p.grid), max(p.grid)),

    max.order.primary = max_from_curve(weighted.all.orders, "maxOrderWR_primary", p.primary.lower, p.primary.upper),
    max.order.full = max_from_curve(weighted.all.orders, "maxOrderWR_full", min(p.grid), max(p.grid)),

    threshold = threshold,
    max.threshold = threshold[best.t.idx, , drop = FALSE],

    counts = list(
      total.pairs = as.numeric(core$total_pairs),
      D1.win = as.numeric(core$D1_win),
      D1.loss = as.numeric(core$D1_loss),
      H2.win = as.numeric(core$H2_win),
      H2.loss = as.numeric(core$H2_loss),
      H1.win = as.numeric(core$H1_win),
      H1.loss = as.numeric(core$H1_loss),
      D2.win = as.numeric(core$D2_win),
      D2.loss = as.numeric(core$D2_loss),
      death.win.pairs = as.numeric(core$death_win_pairs),
      death.loss.pairs = as.numeric(core$death_loss_pairs),
      hosp.win.pairs = as.numeric(core$hosp_win_pairs),
      hosp.loss.pairs = as.numeric(core$hosp_loss_pairs),
      true.tie.pairs = as.numeric(core$true_tie_pairs),
      true.tie.pr = as.numeric(core$true_tie_pairs) / as.numeric(core$total_pairs)
    )
  )
}

wr_at_order_p <- function(engine.out, order, p) {
  if (is.na(order) || is.na(p)) return(NA_real_)
  curve <- if (order == "death_first") engine.out$weighted.death.first else engine.out$weighted.hosp.first
  idx <- extract_nearest(curve$p, p)
  curve$WR[idx]
}

threshold_at_t <- function(engine.out, t) {
  if (is.na(t)) return(NA_real_)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold$WR[idx]
}

threshold_tie_at_t <- function(engine.out, t) {
  if (is.na(t)) return(NA_real_)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold$pr.tie[idx]
}


#Composite endpoint and threshold-grid
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
      median.num.hosp = median(d$NUMHOSP, na.rm = TRUE),
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

run.logrank.tests <- function(ds) {
  tab <- ds$table.output

  lr.death <- survival::survdiff(
    survival::Surv(FUTIME, CNSR) ~ ARM,
    data = tab
  )
  p.death <- pchisq(lr.death$chisq, df = length(lr.death$n) - 1, lower.tail = FALSE)

  comp <- compute.composite.endpoint.data(ds)
  lr.comp <- survival::survdiff(
    survival::Surv(time, event) ~ ARM,
    data = comp
  )
  p.comp <- pchisq(lr.comp$chisq, df = length(lr.comp$n) - 1, lower.tail = FALSE)

  data.frame(
    method = c("Log-rank death endpoint", "Log-rank composite endpoint"),
    statistic = c(lr.death$chisq, lr.comp$chisq),
    p.value = c(p.death, p.comp),
    note = c(
      "Death / survival endpoint only",
      "Time to earliest of death or first hospitalization"
    ),
    stringsAsFactors = FALSE
  )
}

choose.threshold.grid.primary <- function(ds,
                                          clinical.months = CLINICAL_T_MONTHS,
                                          max.K = MAX_T_GRID_SIZE) {
  # Only first-endpoint information is used here.
  # No hospitalization time, recurrent-event count, or composite endpoint time is used.
  tab <- ds$table.output
  max.months <- max(tab$FUTIME, na.rm = TRUE) * 12

  primary.months <- tab$FUTIME * 12
  death.months <- tab$FUTIME[tab$CNSR == 1] * 12

  pair.abs.diff <- numeric(0)
  if (length(primary.months) >= 2) {
    if (length(primary.months) <= 300) {
      pair.abs.diff <- as.vector(stats::dist(primary.months, method = "manhattan"))
    } else {
      set.seed(999)
      id1 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      id2 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      pair.abs.diff <- abs(primary.months[id1] - primary.months[id2])
    }
  }

  diagnostics <- data.frame(
    quantity = c(
      "first endpoint follow-up mean",
      "first endpoint follow-up median",
      "observed death time mean among deaths",
      "observed death time median among deaths",
      "absolute pairwise first-endpoint difference mean",
      "absolute pairwise first-endpoint difference median",
      "absolute pairwise first-endpoint difference q25",
      "absolute pairwise first-endpoint difference q75"
    ),
    months = c(
      safe_mean(primary.months),
      safe_quantile(primary.months, 0.50),
      safe_mean(death.months),
      safe_quantile(death.months, 0.50),
      safe_mean(pair.abs.diff),
      safe_quantile(pair.abs.diff, 0.50),
      safe_quantile(pair.abs.diff, 0.25),
      safe_quantile(pair.abs.diff, 0.75)
    ),
    source = "first endpoint only",
    stringsAsFactors = FALSE
  )

  clinical <- clinical.months[clinical.months > 0 & clinical.months <= max.months]

  priority <- c(
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference mean"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference median"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q25"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q75"],
    0.25 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    0.50 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    diagnostics$months[diagnostics$quantity == "observed death time median among deaths"]
  )

  priority <- round(priority, digits = 1)
  priority <- priority[is.finite(priority) & !is.na(priority)]
  priority <- priority[priority > 0 & priority <= max.months]

  selected <- unique(sort(round(clinical, digits = 1)))
  for (x in priority) {
    if (length(selected) >= max.K) break
    selected <- unique(sort(c(selected, x)))
  }

  if (length(selected) == 0) {
    selected <- round(seq(max.months / 4, max.months, length.out = min(4, max.K)), 1)
    selected <- selected[selected > 0]
  }

  if (length(selected) > max.K) {
    keep.idx <- unique(round(seq(1, length(selected), length.out = max.K)))
    selected <- selected[keep.idx]
  }

  selected <- unique(sort(selected))

  threshold.table <- data.frame(
    t.months = selected,
    t.years = selected / 12,
    source = ifelse(selected %in% round(clinical, 1), "clinical candidate", "first-endpoint data-informed"),
    stringsAsFactors = FALSE
  )

  list(
    t.grid = selected / 12,
    threshold.table = threshold.table,
    diagnostics = diagnostics
  )
}

#Permutation tests: max-statistic and fixed-selected p-values
perm.test.revised <- function(ds,
                              B = B_PERM,
                              seed = MASTER_SEED,
                              p.grid = P_GRID_ALL,
                              t.grid = CLINICAL_T_MONTHS / 12,
                              verbose = FALSE) {
  ds <- prepare.ds.fast(ds)
  set.seed(seed)

  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR",
    "maxWRp_primary",
    "maxWRp_low",
    "maxWRp_full",
    "maxOrderWR_primary",
    "maxOrderWR_full",
    "maxWRt"
  )

  fixed.names <- c(
    "ordinaryWR",
    "fixedWRp_primary",
    "fixedWRp_low",
    "fixedWRp_full",
    "fixedOrder_primary",
    "fixedOrder_full",
    "fixedWRt"
  )

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )

  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names))
  colnames(perm.max) <- max.names
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names))
  colnames(perm.fixed) <- fixed.names

  # Store pointwise permutation curves so the script can report
  # average p-values and average tie counts at each candidate p and t.
  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.p.primary = NA_real_,
    selected.p.low = NA_real_,
    selected.p.full = NA_real_,
    selected.order.primary = NA_character_,
    selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_,
    selected.order.p.full = NA_real_,
    selected.t = NA_real_,
    stringsAsFactors = FALSE
  )

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

    perm.max[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.max[b, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
    perm.max[b, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
    perm.max[b, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
    perm.max[b, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
    perm.max[b, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
    perm.max[b, "maxWRt"] <- out.b$max.threshold$WR[1]

    perm.fixed[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.fixed[b, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
    perm.fixed[b, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
    perm.fixed[b, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
    perm.fixed[b, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
    perm.fixed[b, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
    perm.fixed[b, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

    selected.perm$selected.p.primary[b] <- out.b$max.weighted.primary$p[1]
    selected.perm$selected.p.low[b] <- out.b$max.weighted.low$p[1]
    selected.perm$selected.p.full[b] <- out.b$max.weighted.full$p[1]
    selected.perm$selected.order.primary[b] <- out.b$max.order.primary$order[1]
    selected.perm$selected.order.p.primary[b] <- out.b$max.order.primary$p[1]
    selected.perm$selected.order.full[b] <- out.b$max.order.full$order[1]
    selected.perm$selected.order.p.full[b] <- out.b$max.order.full$p[1]
    selected.perm$selected.t[b] <- out.b$max.threshold$t[1]
  }

  calc_perm_p <- function(perm.vec, obs.value) {
    (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  }

  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))

  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k])
  })
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k])
  })
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) {
    calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k])
  })

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(
      death_first = perm.pointwise.death.first,
      hospitalization_first = perm.pointwise.hosp.first,
      threshold = perm.pointwise.threshold
    ),
    p.value.pointwise = list(
      death_first = p.value.pointwise.death.first,
      hospitalization_first = p.value.pointwise.hosp.first,
      threshold = p.value.pointwise.threshold
    ),
    selected.fixed = list(
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    p.grid = p.grid,
    t.grid = t.grid
  )
}

#Output tables and plotting helpers
fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 0.001, "<0.001", sprintf("%.4f", x)))
}

save_png <- function(outdir, filename, width = 1800, height = 1100, res = 160) {
  png(file.path(outdir, filename), width = width, height = height, res = res)
}

make.pvalue.table <- function(test.out, logrank.results = NULL) {
  max.tab <- data.frame(
    method = names(test.out$T.obs.max),
    statistic.type = "max-statistic permutation p-value",
    observed.statistic = as.numeric(test.out$T.obs.max),
    p.value = as.numeric(test.out$p.value.max),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  fixed.tab <- data.frame(
    method = names(test.out$T.obs.fixed),
    statistic.type = "fixed selected-parameter permutation p-value",
    observed.statistic = as.numeric(test.out$T.obs.fixed),
    p.value = as.numeric(test.out$p.value.fixed),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  out <- rbind(max.tab, fixed.tab)

  if (!is.null(logrank.results)) {
    lr.tab <- data.frame(
      method = logrank.results$method,
      statistic.type = "log-rank chi-square p-value",
      observed.statistic = logrank.results$statistic,
      p.value = logrank.results$p.value,
      B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.tab)
  }

  out
}

make.max.vs.fixed.table <- function(test.out) {
  data.frame(
    method = c(
      "Traditional WR",
      "Weighted p: death first, p in [0.5, 1]",
      "Weighted p: death first, p < 0.5",
      "Weighted p: death first, full p grid",
      "Maximum order WR: p in [0.5, 1]",
      "Maximum order WR: full p grid",
      "Threshold WR(t)"
    ),
    selected.parameter = c(
      "p = 0.5, death-first",
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), ", death-first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), ", death-first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), ", death-first"),
      paste0("order = ", test.out$selected.fixed$order.primary,
             ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", test.out$selected.fixed$order.full,
             ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    max.statistic = as.numeric(test.out$T.obs.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    fixed.parameter.statistic = as.numeric(test.out$T.obs.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    stringsAsFactors = FALSE
  )
}


get_selected_weighted_row <- function(engine.out, order, p) {
  if (is.na(order) || is.na(p)) return(NULL)
  curve <- if (order == "death_first") engine.out$weighted.death.first else engine.out$weighted.hosp.first
  idx <- extract_nearest(curve$p, p)
  curve[idx, , drop = FALSE]
}

get_selected_threshold_row <- function(engine.out, t) {
  if (is.na(t)) return(NULL)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold[idx, , drop = FALSE]
}

make.selected.max.table <- function(test.out) {
  obs <- test.out$observed

  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)

  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])

  data.frame(
    method = c(
      "Traditional WR",
      "Weighted p: death first, p in [0.5, 1]",
      "Weighted p: death first, p < 0.5",
      "Weighted p: death first, full p grid",
      "Maximum order WR: p in [0.5, 1]",
      "Maximum order WR: full p grid",
      "Threshold WR(t)"
    ),
    selected.order = c(
      "death_first",
      "death_first",
      "death_first",
      "death_first",
      test.out$selected.fixed$order.primary,
      test.out$selected.fixed$order.full,
      NA_character_
    ),
    selected.p = c(
      0.5,
      test.out$selected.fixed$p.primary,
      test.out$selected.fixed$p.low,
      test.out$selected.fixed$p.full,
      test.out$selected.fixed$order.p.primary,
      test.out$selected.fixed$order.p.full,
      NA_real_
    ),
    selected.t.months = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      test.out$selected.fixed$t.months
    ),
    observed.max.statistic = as.numeric(test.out$T.obs.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    selected.tie.count = c(
      get_tie_count(w.traditional),
      get_tie_count(w.primary),
      get_tie_count(w.low),
      get_tie_count(w.full),
      get_tie_count(o.primary),
      get_tie_count(o.full),
      get_tie_count(t.row)
    ),
    selected.tie.pr = c(
      get_tie_pr(w.traditional),
      get_tie_pr(w.primary),
      get_tie_pr(w.low),
      get_tie_pr(w.full),
      get_tie_pr(o.primary),
      get_tie_pr(o.full),
      get_tie_pr(t.row)
    ),
    stringsAsFactors = FALSE
  )
}

make.pointwise.table <- function(test.out) {
  obs <- test.out$observed

  death <- obs$weighted.death.first
  death$family <- "weighted_p"
  death$parameter.name <- "p"
  death$parameter <- death$p
  death$parameter.months <- NA_real_
  death$parameter.label <- paste0("p = ", sprintf("%.2f", death$p))
  death$observed.WR <- death$WR
  death$pointwise.p.value <- as.numeric(test.out$p.value.pointwise$death_first)

  hosp <- obs$weighted.hosp.first
  hosp$family <- "weighted_p"
  hosp$parameter.name <- "p"
  hosp$parameter <- hosp$p
  hosp$parameter.months <- NA_real_
  hosp$parameter.label <- paste0("p = ", sprintf("%.2f", hosp$p))
  hosp$observed.WR <- hosp$WR
  hosp$pointwise.p.value <- as.numeric(test.out$p.value.pointwise$hospitalization_first)

  thr <- obs$threshold
  thr.out <- data.frame(
    order = NA_character_,
    p = NA_real_,
    WR = thr$WR,
    win.score = thr$wins,
    loss.score = thr$losses,
    win.pairs = thr$wins,
    loss.pairs = thr$losses,
    tie.count = thr$tie.count,
    tie.pr = thr$pr.tie,
    family = "threshold_t",
    parameter.name = "t",
    parameter = thr$t,
    parameter.months = thr$t.months,
    parameter.label = paste0("t = ", sprintf("%.1f", thr$t.months), " months"),
    observed.WR = thr$WR,
    pointwise.p.value = as.numeric(test.out$p.value.pointwise$threshold),
    stringsAsFactors = FALSE
  )

  keep.cols <- c(
    "family", "order", "parameter.name", "parameter", "parameter.months",
    "parameter.label", "observed.WR", "pointwise.p.value", "win.score",
    "loss.score", "win.pairs", "loss.pairs", "tie.count", "tie.pr"
  )

  out <- rbind(
    death[, keep.cols, drop = FALSE],
    hosp[, keep.cols, drop = FALSE],
    thr.out[, keep.cols, drop = FALSE]
  )
  out$total.pairs <- test.out$observed$counts$total.pairs
  out
}

plot.weighted.death.first <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_weighted_death_first_WRp.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  tab <- engine.out$weighted.death.first
  plot(tab$p, tab$WR, type = "l", lwd = 2,
       xlab = "p: death/survival weight",
       ylab = "WR(p)",
       main = "Death-first weighted WR(p)")
  abline(v = 0.5, lty = 3, lwd = 2)
  points(engine.out$max.weighted.primary$p, engine.out$max.weighted.primary$WR, pch = 19)
  points(engine.out$max.weighted.low$p, engine.out$max.weighted.low$WR, pch = 17)
  legend("topright",
         legend = c(
           paste0("primary max p = ", sprintf("%.2f", engine.out$max.weighted.primary$p),
                  ", WR = ", sprintf("%.3f", engine.out$max.weighted.primary$WR)),
           paste0("low-p max p = ", sprintf("%.2f", engine.out$max.weighted.low$p),
                  ", WR = ", sprintf("%.3f", engine.out$max.weighted.low$WR))
         ),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.order.curves <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_maximum_order_WRp.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  dtab <- engine.out$weighted.death.first
  htab <- engine.out$weighted.hosp.first
  yr <- range(c(dtab$WR, htab$WR), finite = TRUE)
  plot(dtab$p, dtab$WR, type = "l", lwd = 2, ylim = yr,
       xlab = "p: first-endpoint weight",
       ylab = "WR(p)",
       main = "Effect of hierarchy order on maximum WR")
  lines(htab$p, htab$WR, lwd = 2, lty = 2)
  points(engine.out$max.order.primary$p, engine.out$max.order.primary$WR, pch = 19)
  abline(v = 0.5, lty = 3, lwd = 2)
  legend("topright",
         legend = c(
           "death first",
           "hospitalization first",
           paste0("selected order = ", engine.out$max.order.primary$order),
           paste0("selected p = ", sprintf("%.2f", engine.out$max.order.primary$p)),
           paste0("max WR = ", sprintf("%.3f", engine.out$max.order.primary$WR))
         ),
         lty = c(1, 2, NA, NA, NA),
         pch = c(NA, NA, NA, NA, NA),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.threshold.curve <- function(engine.out, outdir, prefix, p.value = NA_real_) {
  save_png(outdir, paste0(prefix, "_threshold_WRt.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  tab <- engine.out$threshold
  plot(tab$t.months, tab$WR, type = "b", pch = 19, lwd = 2,
       xlab = "threshold t (months)",
       ylab = "WR(t)",
       main = "Threshold WR(t)")
  abline(v = engine.out$max.threshold$t.months, lty = 2, lwd = 2)
  legend("topright",
         legend = c(
           paste0("selected t = ", sprintf("%.1f", engine.out$max.threshold$t.months), " months"),
           paste0("max WR = ", sprintf("%.3f", engine.out$max.threshold$WR)),
           paste0("max-stat p = ", fmt_p(p.value))
         ),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.max.vs.fixed.pvalues <- function(max.vs.fixed.table, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_max_vs_fixed_pvalues.png"), width = 2100, height = 1150)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(7.5, 5.2, 3.2, 1.2))
  mat <- rbind(
    "Max-statistic" = max.vs.fixed.table$max.statistic.p.value,
    "Fixed-selected" = max.vs.fixed.table$fixed.parameter.p.value
  )
  bp <- barplot(mat, beside = TRUE,
                names.arg = c("Trad WR", "p>=0.5", "p<0.5", "full p", "order p>=0.5", "order full", "threshold"),
                ylim = c(0, 1),
                ylab = "Permutation p-value",
                main = "Max-statistic vs fixed-selected p-values",
                legend.text = TRUE,
                args.legend = list(x = "topright", bty = "n"),
                las = 2, cex.names = 0.75)
  abline(h = ALPHA, lty = 2, lwd = 2)
  text(x = bp, y = pmin(mat + 0.035, 0.98), labels = sprintf("%.3f", mat), cex = 0.65)
  par(old.par)
  dev.off()
}

plot.km.with.threshold <- function(ds, engine.out, outdir, prefix) {
  tab <- ds$table.output
  comp <- compute.composite.endpoint.data(ds)
  selected.t <- engine.out$max.threshold$t[1]

  save_png(outdir, paste0(prefix, "_KM_death_with_selected_t.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  fit.death <- survival::survfit(survival::Surv(FUTIME, CNSR) ~ ARM, data = tab)
  plot(fit.death, lwd = 2, mark.time = TRUE,
       xlab = "Time (years)", ylab = "Survival probability",
       main = "Kaplan-Meier curve: death endpoint")
  abline(v = selected.t, lty = 2, lwd = 2)
  legend("bottomleft",
         legend = c("ARM 0: control", "ARM 1: treatment",
                    paste0("selected t = ", sprintf("%.1f", selected.t * 12), " months")),
         lty = c(1, 1, 2), lwd = c(2, 2, 2), bty = "n", cex = 0.85)
  par(old.par)
  dev.off()

  save_png(outdir, paste0(prefix, "_KM_composite_with_selected_t.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  fit.comp <- survival::survfit(survival::Surv(time, event) ~ ARM, data = comp)
  plot(fit.comp, lwd = 2, mark.time = TRUE,
       xlab = "Time (years)", ylab = "Composite-event-free probability",
       main = "Kaplan-Meier curve: composite endpoint")
  abline(v = selected.t, lty = 2, lwd = 2)
  legend("bottomleft",
         legend = c("ARM 0: control", "ARM 1: treatment",
                    paste0("selected t = ", sprintf("%.1f", selected.t * 12), " months")),
         lty = c(1, 1, 2), lwd = c(2, 2, 2), bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}


#Scenario configuration
build.scenario.grid <- function() {
  base.mort <- -log(0.6)

  #Scenario interpretation:
  #equivalence/null: treatment and control have the same generating parameters.
  #similarity: treatment is only slightly better or slightly worse.
  #superiority: treatment is meaningfully better for the configured endpoint(s).
  #inferiority: treatment is meaningfully worse for the configured endpoint(s).
  #mixed effects: one endpoint improves while the other worsens.
  #censoring stress tests: same superiority signal, but censoring/follow-up differs.
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
    "S13_low_censoring_long_followup_superiority"
  )

  # Short labels used in tables and plots. They intentionally avoid underscores.
  scenario_label <- c(
    "S01 Equal",
    "S02 Similar +",
    "S03 Similar -",
    "S04 Sup mod",
    "S05 Sup strong",
    "S06 Inf mod",
    "S07 Inf strong",
    "S08 Death +",
    "S09 Hosp +",
    "S10 Death + / Hosp -",
    "S11 Death - / Hosp +",
    "S12 High censor",
    "S13 Long FU"
  )

  scenario_type <- c(
    "equivalence/null",
    "similarity",
    "similarity",
    "superiority",
    "superiority",
    "inferiority",
    "inferiority",
    "component-specific superiority",
    "component-specific superiority",
    "mixed effects",
    "mixed effects",
    "superiority with high censoring",
    "superiority with long follow-up"
  )

  description <- c(
    "Treatment and control have equal death and hospitalization distributions.",
    "Treatment is slightly better for both death and hospitalization; near-null/similarity case.",
    "Treatment is slightly worse for both death and hospitalization; near-null/similarity case.",
    "Treatment is moderately better for both death and hospitalization.",
    "Treatment is strongly better for both death and hospitalization.",
    "Treatment is moderately worse for both death and hospitalization.",
    "Treatment is strongly worse for both death and hospitalization.",
    "Treatment improves death only; hospitalization distribution is equal.",
    "Treatment improves hospitalization only; death distribution is equal.",
    "Treatment improves death but worsens hospitalization.",
    "Treatment worsens death but improves hospitalization.",
    "Treatment improves both endpoints, but random censoring is high.",
    "Treatment improves both endpoints with longer follow-up and low random censoring."
  )

  n.scen <- length(scenario_id)

  # Notes on parameters:
  #   HR < 1 means lower mortality hazard in treatment.
  #   evt.rate.scale.param.trt < evt.rate.scale.param.ctr means fewer recurrent
  #   hospitalizations on average in treatment because the gamma mean is
  #   shape * scale.
  data.frame(
    scenario_index = seq_along(scenario_id),
    scenario_id = scenario_id,
    scenario_label = scenario_label,
    scenario_type = scenario_type,
    description = description,
    N0 = rep(50, n.scen),
    N1 = rep(50, n.scen),
    mort.rate.ctrl = rep(base.mort, n.scen),

    # Mortality scenario parameter.
    HR = c(
      1.00,  # equivalence/null
      0.90,  # similarity, small benefit
      1.10,  # similarity, small harm
      0.75,  # superiority, moderate benefit
      0.60,  # superiority, strong benefit
      1.25,  # inferiority, moderate harm
      1.40,  # inferiority, strong harm
      0.60,  # death benefit only
      1.00,  # hospitalization benefit only
      0.60,  # death benefit + hospitalization harm
      1.40,  # death harm + hospitalization benefit
      0.60,  # high censoring superiority
      0.60   # low censoring, long follow-up superiority
    ),

    evt.rate.shape.param = rep(5, n.scen),
    evt.rate.scale.param.ctr = rep(1.0, n.scen),

    # Recurrent hospitalization scenario parameter.
    evt.rate.scale.param.trt = c(
      1.00,  # equivalence/null
      0.90,  # similarity, small benefit
      1.10,  # similarity, small harm
      0.75,  # superiority, moderate benefit
      0.50,  # superiority, strong benefit
      1.25,  # inferiority, moderate harm
      1.50,  # inferiority, strong harm
      1.00,  # death benefit only
      0.50,  # hospitalization benefit only
      1.50,  # death benefit + hospitalization harm
      0.50,  # death harm + hospitalization benefit
      0.50,  # high censoring superiority
      0.50   # low censoring, long follow-up superiority
    ),

    max.FU = c(
      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3
    ),

    censor.rate = c(
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.5, 0.05
    ),

    stringsAsFactors = FALSE
  )
}


# Single-trial extraction and scenario summaries
get_comp_arm_value <- function(comp.stats, arm, varname) {
  z <- comp.stats[comp.stats$ARM == arm, varname]
  if (length(z) == 0) return(NA_real_)
  as.numeric(z[1])
}

get_scenario_field <- function(scenario.row, field, default = NA_character_) {
  if (field %in% names(scenario.row)) {
    value <- scenario.row[[field]][1]
    if (length(value) == 0 || is.na(value)) return(default)
    return(as.character(value))
  }
  default
}

extract.trial.row <- function(scenario.row, sim.index, test.out, logrank.results, comp.stats,
                              threshold.info, error.message = NA_character_) {
  obs <- test.out$observed

  death.first.primary <- obs$max.weighted.primary
  hosp.first.primary <- max_from_curve(obs$weighted.hosp.first, "hosp_first_primary", 0.50, 1.00)
  order.diff <- hosp.first.primary$WR[1] - death.first.primary$WR[1]

  w.traditional.row <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.primary$p[1])
  w.low.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.low$p[1])
  w.full.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.full$p[1])
  o.primary.row <- get_selected_weighted_row(obs, obs$max.order.primary$order[1], obs$max.order.primary$p[1])
  o.full.row <- get_selected_weighted_row(obs, obs$max.order.full$order[1], obs$max.order.full$p[1])

  data.frame(
    scenario_id = scenario.row$scenario_id,
    scenario_type = get_scenario_field(scenario.row, "scenario_type"),
    description = get_scenario_field(scenario.row, "description"),
    sim = sim.index,
    error.message = error.message,

    ordinaryWR = test.out$T.obs.max["ordinaryWR"],
    traditional.tie.count = w.traditional.row$tie.count[1],
    traditional.tie.pr = w.traditional.row$tie.pr[1],
    max.pvalue.ordinaryWR = test.out$p.value.max["ordinaryWR"],
    fixed.pvalue.ordinaryWR = test.out$p.value.fixed["ordinaryWR"],

    maxWRp_primary = test.out$T.obs.max["maxWRp_primary"],
    selected.p.primary = obs$max.weighted.primary$p[1],
    weighted.primary.tie.count = w.primary.row$tie.count[1],
    weighted.primary.tie.pr = w.primary.row$tie.pr[1],
    max.pvalue.maxWRp_primary = test.out$p.value.max["maxWRp_primary"],
    fixed.pvalue.maxWRp_primary = test.out$p.value.fixed["fixedWRp_primary"],

    maxWRp_low = test.out$T.obs.max["maxWRp_low"],
    selected.p.low = obs$max.weighted.low$p[1],
    weighted.low.tie.count = w.low.row$tie.count[1],
    weighted.low.tie.pr = w.low.row$tie.pr[1],
    max.pvalue.maxWRp_low = test.out$p.value.max["maxWRp_low"],
    fixed.pvalue.maxWRp_low = test.out$p.value.fixed["fixedWRp_low"],

    maxWRp_full = test.out$T.obs.max["maxWRp_full"],
    selected.p.full = obs$max.weighted.full$p[1],
    weighted.full.tie.count = w.full.row$tie.count[1],
    weighted.full.tie.pr = w.full.row$tie.pr[1],
    selected.full.is.low.p = as.integer(obs$max.weighted.full$p[1] < 0.5),
    max.pvalue.maxWRp_full = test.out$p.value.max["maxWRp_full"],
    fixed.pvalue.maxWRp_full = test.out$p.value.fixed["fixedWRp_full"],

    death.first.primary.WR = death.first.primary$WR[1],
    hosp.first.primary.WR = hosp.first.primary$WR[1],
    order.diff.hosp.minus.death = order.diff,
    order.changed.primary = as.integer(obs$max.order.primary$order[1] == "hospitalization_first"),
    maxOrderWR_primary = test.out$T.obs.max["maxOrderWR_primary"],
    selected.order.primary = obs$max.order.primary$order[1],
    selected.order.p.primary = obs$max.order.primary$p[1],
    order.primary.tie.count = o.primary.row$tie.count[1],
    order.primary.tie.pr = o.primary.row$tie.pr[1],
    max.pvalue.maxOrderWR_primary = test.out$p.value.max["maxOrderWR_primary"],
    fixed.pvalue.maxOrderWR_primary = test.out$p.value.fixed["fixedOrder_primary"],

    maxOrderWR_full = test.out$T.obs.max["maxOrderWR_full"],
    selected.order.full = obs$max.order.full$order[1],
    selected.order.p.full = obs$max.order.full$p[1],
    order.full.tie.count = o.full.row$tie.count[1],
    order.full.tie.pr = o.full.row$tie.pr[1],
    selected.order.full.is.low.p = as.integer(obs$max.order.full$p[1] < 0.5),
    max.pvalue.maxOrderWR_full = test.out$p.value.max["maxOrderWR_full"],
    fixed.pvalue.maxOrderWR_full = test.out$p.value.fixed["fixedOrder_full"],

    maxWRt = test.out$T.obs.max["maxWRt"],
    selected.t.years = obs$max.threshold$t[1],
    selected.t.months = obs$max.threshold$t.months[1],
    threshold.tie.count = obs$max.threshold$tie.count[1],
    threshold.tie.pr = obs$max.threshold$pr.tie[1],
    max.pvalue.maxWRt = test.out$p.value.max["maxWRt"],
    fixed.pvalue.maxWRt = test.out$p.value.fixed["fixedWRt"],

    true.hierarchical.tie.count = obs$counts$true.tie.pairs,
    true.hierarchical.tie.pr = obs$counts$true.tie.pr,
    total.pairs = obs$counts$total.pairs,

    logrank.death.statistic = logrank.results$statistic[logrank.results$method == "Log-rank death endpoint"],
    logrank.death.p = logrank.results$p.value[logrank.results$method == "Log-rank death endpoint"],
    logrank.composite.statistic = logrank.results$statistic[logrank.results$method == "Log-rank composite endpoint"],
    logrank.composite.p = logrank.results$p.value[logrank.results$method == "Log-rank composite endpoint"],

    death.event.rate.control = get_comp_arm_value(comp.stats, 0, "death.event.rate"),
    death.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "death.event.rate"),
    composite.event.rate.control = get_comp_arm_value(comp.stats, 0, "composite.event.rate"),
    composite.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "composite.event.rate"),
    mean.num.hosp.control = get_comp_arm_value(comp.stats, 0, "mean.num.hosp"),
    mean.num.hosp.treatment = get_comp_arm_value(comp.stats, 1, "mean.num.hosp"),

    threshold.grid.months = paste(sprintf("%.1f", threshold.info$threshold.table$t.months), collapse = ", "),
    stringsAsFactors = FALSE
  )
}


# Base-R row bind that tolerates error rows with fewer columns.
rbind_fill_base <- function(a, b) {
  if (nrow(a) == 0) return(b)
  cols <- union(names(a), names(b))
  for (cc in setdiff(cols, names(a))) a[[cc]] <- NA
  for (cc in setdiff(cols, names(b))) b[[cc]] <- NA
  rbind(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}

extract.pointwise.rows <- function(scenario.row, sim.index, test.out) {
  out <- make.pointwise.table(test.out)
  out$scenario_id <- scenario.row$scenario_id
  out$scenario_type <- get_scenario_field(scenario.row, "scenario_type")
  out$description <- get_scenario_field(scenario.row, "description")
  out$sim <- sim.index
  out[, c("scenario_id", "scenario_type", "description", "sim",
           setdiff(names(out), c("scenario_id", "scenario_type", "description", "sim"))),
      drop = FALSE]
}

summarise.pointwise.results <- function(pointwise, alpha = ALPHA) {
  if (nrow(pointwise) == 0) return(data.frame())

  keys <- unique(pointwise[, c(
    "scenario_id", "family", "order", "parameter.name",
    "parameter", "parameter.months", "parameter.label"
  ), drop = FALSE])

  out.list <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    key <- keys[i, , drop = FALSE]
    keep <- pointwise$scenario_id == key$scenario_id &
      pointwise$family == key$family &
      pointwise$parameter.name == key$parameter.name &
      abs(pointwise$parameter - key$parameter) < 1e-10

    if (is.na(key$order)) {
      keep <- keep & is.na(pointwise$order)
    } else {
      keep <- keep & pointwise$order == key$order
    }

    d <- pointwise[keep, , drop = FALSE]

    out.list[[i]] <- data.frame(
      scenario_id = key$scenario_id,
      scenario_type = if ("scenario_type" %in% names(d)) d$scenario_type[1] else NA_character_,
      description = if ("description" %in% names(d)) d$description[1] else NA_character_,
      family = key$family,
      order = key$order,
      parameter.name = key$parameter.name,
      parameter = key$parameter,
      parameter.months = key$parameter.months,
      parameter.label = key$parameter.label,
      nsim.available = nrow(d),
      mean.observed.WR = mean(d$observed.WR, na.rm = TRUE),
      sd.observed.WR = safe_sd(d$observed.WR),
      mean.pointwise.p.value = mean(d$pointwise.p.value, na.rm = TRUE),
      median.pointwise.p.value = median(d$pointwise.p.value, na.rm = TRUE),
      rejection.proportion.pointwise = mean(d$pointwise.p.value < alpha, na.rm = TRUE),
      mean.win.score = mean(d$win.score, na.rm = TRUE),
      mean.loss.score = mean(d$loss.score, na.rm = TRUE),
      mean.win.pairs = mean(d$win.pairs, na.rm = TRUE),
      mean.loss.pairs = mean(d$loss.pairs, na.rm = TRUE),
      mean.tie.count = mean(d$tie.count, na.rm = TRUE),
      mean.tie.pr = mean(d$tie.pr, na.rm = TRUE),
      mean.total.pairs = mean(d$total.pairs, na.rm = TRUE),
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, out.list)
  rownames(out) <- NULL
  out
}


#Scenario
run.one.simulation <- function(scenario.row, sim.index, seed, B = B_PERM, save.example = FALSE,
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

  ds <- apply.random.censoring(ds, censor.rate = scenario.row$censor.rate, seed = seed + 17)
  ds <- prepare.ds.fast(ds)

  threshold.info <- choose.threshold.grid.primary(ds)

  test.out <- perm.test.revised(
    ds = ds,
    B = B,
    seed = seed + 100000,
    p.grid = P_GRID_ALL,
    t.grid = threshold.info$t.grid,
    verbose = FALSE
  )

  logrank.results <- run.logrank.tests(ds)
  comp.stats <- composite.statistics(ds)

  if (save.example) {
    save.one.dataset.outputs(
      ds = ds,
      test.out = test.out,
      threshold.info = threshold.info,
      logrank.results = logrank.results,
      comp.stats = comp.stats,
      outdir = scenario.outdir,
      prefix = paste0("example_sim", sim.index)
    )
    saveRDS(
      list(ds = ds, threshold.info = threshold.info, test.out = test.out,
           logrank.results = logrank.results, comp.stats = comp.stats),
      file.path(scenario.outdir, paste0("example_sim", sim.index, "_full_object.rds"))
    )
  }

  summary.row <- extract.trial.row(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out,
    logrank.results = logrank.results,
    comp.stats = comp.stats,
    threshold.info = threshold.info
  )

  pointwise.rows <- extract.pointwise.rows(
    scenario.row = scenario.row,
    sim.index = sim.index,
    test.out = test.out
  )

  list(summary.row = summary.row, pointwise.rows = pointwise.rows)
}


#Final reporting layer: clean labels, scenario comparisons,
scenario.display.label <- function(x) {
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
    S13_low_censoring_long_followup_superiority = "S13 Long FU"
  )
  out <- as.character(x)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

method.long.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional WR",
    maxWRp_primary = "Weighted max WR, p >= 0.5",
    maxWRp_low = "Weighted max WR, p < 0.5",
    maxWRp_full = "Weighted max WR, full p grid",
    maxOrderWR_primary = "Maximum-order WR, p >= 0.5",
    maxOrderWR_full = "Maximum-order WR, full p grid",
    maxWRt = "Threshold max WR(t)",
    logrank_death = "Log-rank death",
    logrank_composite = "Log-rank composite",
    fixedWRp_primary = "Weighted fixed-selected WR, p >= 0.5",
    fixedWRp_low = "Weighted fixed-selected WR, p < 0.5",
    fixedWRp_full = "Weighted fixed-selected WR, full p grid",
    fixedOrderWR_primary = "Maximum-order fixed-selected WR, p >= 0.5",
    fixedOrderWR_full = "Maximum-order fixed-selected WR, full p grid",
    fixedWRt = "Threshold fixed-selected WR(t)"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

method.display.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional\nWR",
    maxWRp_primary = "Weighted\nmax",
    maxWRp_low = "Low-p\nmax",
    maxWRp_full = "Full-p\nmax",
    maxOrderWR_primary = "Order\nmax",
    maxOrderWR_full = "Order full\nmax",
    maxWRt = "Threshold\nmax",
    logrank_death = "LR\ndeath",
    logrank_composite = "LR\ncomposite",
    fixedWRp_primary = "Weighted\nfixed",
    fixedWRp_low = "Low-p\nfixed",
    fixedWRp_full = "Full-p\nfixed",
    fixedOrderWR_primary = "Order\nfixed",
    fixedOrderWR_full = "Order full\nfixed",
    fixedWRt = "Threshold\nfixed"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

clean.order.label <- function(x) {
  out <- as.character(x)
  out <- gsub("death_first", "death first", out)
  out <- gsub("hospitalization_first", "hospitalization first", out)
  out <- gsub("_", " ", out)
  out
}

safe.mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe.median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe.col <- function(df, col) {
  if (is.na(col) || !(col %in% names(df))) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col]])
}

mode.string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

get.power.value <- function(power.table, method) {
  if (is.na(method) || is.null(power.table) || nrow(power.table) == 0) return(NA_real_)
  z <- power.table$rejection.proportion[power.table$method == method]
  if (length(z) == 0) return(NA_real_)
  as.numeric(z[1])
}


# tables and figures
make.pathway.results.table <- function(test.out, logrank.results = NULL) {
  obs <- test.out$observed

  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)

  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])

  method.id <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )

  out <- data.frame(
    method = method.id,
    method_label = method.long.label(method.id),
    statistic_type = "Win ratio",
    statistic_name = c(
      "WR(p = 0.50), death first",
      "max WR(p), death first, p in [0.50, 1.00]",
      "max WR(p), death first, p in [0.01, 0.49]",
      "max WR(p), death first, p in [0.01, 1.00]",
      "max WR(p), selected order, p in [0.50, 1.00]",
      "max WR(p), selected order, p in [0.01, 1.00]",
      "max WR(t)"
    ),
    observed_statistic = as.numeric(test.out$T.obs.max[method.id]),
    selected_order = c(
      "death_first", "death_first", "death_first", "death_first",
      test.out$selected.fixed$order.primary,
      test.out$selected.fixed$order.full,
      NA_character_
    ),
    selected_p = c(
      0.5,
      test.out$selected.fixed$p.primary,
      test.out$selected.fixed$p.low,
      test.out$selected.fixed$p.full,
      test.out$selected.fixed$order.p.primary,
      test.out$selected.fixed$order.p.full,
      NA_real_
    ),
    selected_t_months = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      test.out$selected.fixed$t.months
    ),
    selected_value = c(
      "p = 0.50; death first",
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), "; death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary),
             "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full),
             "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    permutation_p_value = as.numeric(test.out$p.value.max[method.id]),
    selected_parameter_p_value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
      "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
    )]),
    tie_count = c(
      get_tie_count(w.traditional), get_tie_count(w.primary), get_tie_count(w.low),
      get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full),
      get_tie_count(t.row)
    ),
    tie_proportion = c(
      get_tie_pr(w.traditional), get_tie_pr(w.primary), get_tie_pr(w.low),
      get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full),
      get_tie_pr(t.row)
    ),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  if (!is.null(logrank.results)) {
    lr.death <- logrank.results[logrank.results$method == "Log-rank death endpoint", , drop = FALSE]
    lr.comp <- logrank.results[logrank.results$method == "Log-rank composite endpoint", , drop = FALSE]
    lr.rows <- data.frame(
      method = c("logrank_death", "logrank_composite"),
      method_label = method.long.label(c("logrank_death", "logrank_composite")),
      statistic_type = "Log-rank chi-square",
      statistic_name = c("Death endpoint", "Composite endpoint"),
      observed_statistic = c(lr.death$statistic[1], lr.comp$statistic[1]),
      selected_order = NA_character_,
      selected_p = NA_real_,
      selected_t_months = NA_real_,
      selected_value = NA_character_,
      permutation_p_value = c(lr.death$p.value[1], lr.comp$p.value[1]),
      selected_parameter_p_value = NA_real_,
      tie_count = NA_real_,
      tie_proportion = NA_real_,
      B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.rows)
  }

  out
}

plot.pvalue.bar <- function(tab, outdir, filename, title, p.col = "permutation_p_value",
                            method.order = NULL, show.alpha = TRUE) {
  if (is.null(tab) || nrow(tab) == 0 || !(p.col %in% names(tab))) return(invisible(NULL))
  d <- tab
  if (!is.null(method.order)) d <- d[d$method %in% method.order, , drop = FALSE]
  d <- d[!is.na(d[[p.col]]), , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  if (!is.null(method.order)) {
    d$method <- factor(d$method, levels = method.order)
    d <- d[order(d$method), , drop = FALSE]
  }
  d <- d[order(d[[p.col]], decreasing = TRUE), , drop = FALSE]

  save_png(outdir, filename, width = 2100, height = 1150)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 13.0, 3.2, 1.2))
  bp <- barplot(d[[p.col]], horiz = TRUE,
                names.arg = method.long.label(as.character(d$method)),
                xlim = c(0, 1), las = 1,
                xlab = "P-value",
                main = title,
                cex.names = 0.86)
  if (show.alpha) abline(v = ALPHA, lty = 2, lwd = 2)
  text(x = pmin(d[[p.col]] + 0.035, 0.96), y = bp,
       labels = fmt_p(d[[p.col]]), cex = 0.82)
  if (show.alpha) legend("bottomright", legend = paste0("alpha = ", ALPHA), bty = "n")
  par(old.par)
  dev.off()
}

plot.example.pvalue.figures <- function(pathway.tab, outdir, prefix) {
  all.methods <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt",
    "logrank_death", "logrank_composite"
  )
  wr.methods <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  weighted.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full")
  logrank.methods <- c("logrank_death", "logrank_composite")
  threshold.methods <- c("maxWRt")

  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_all_methods_with_logrank.png"),
                  "All method p-values", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_WR_pathways_only_no_logrank.png"),
                  "WR pathway p-values", "permutation_p_value", wr.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_traditional_and_weighted_only.png"),
                  "Traditional and weighted WR p-values", "permutation_p_value", weighted.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_logrank_pvalues_separate.png"),
                  "Log-rank p-values", "permutation_p_value", logrank.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_threshold_pvalue_separate.png"),
                  "Threshold pathway p-value", "permutation_p_value", threshold.methods)
}

plot.permutation.null.revised <- function(test.out, method, outdir, prefix) {
  if (!(method %in% colnames(test.out$T.perm.max))) return(invisible(NULL))
  x <- as.numeric(test.out$T.perm.max[, method])
  obs.value <- as.numeric(test.out$T.obs.max[method])
  p.value <- as.numeric(test.out$p.value.max[method])

  save_png(outdir, paste0(prefix, "_perm_null_", method, ".png"), width = 1600, height = 1000)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.2, 1.2))
  hist(x, breaks = 35,
       main = paste0("Permutation null: ", method.long.label(method)),
       xlab = "Permutation statistic",
       ylab = "Frequency")
  abline(v = obs.value, lty = 2, lwd = 2)
  legend("topright",
         legend = c(
           paste0("observed = ", sprintf("%.3f", obs.value)),
           paste0("p = ", fmt_p(p.value)),
           paste0("B = ", test.out$B)
         ), bty = "n")
  par(old.par)
  dev.off()
}

plot.example.permutation.nulls <- function(test.out, outdir, prefix) {
  for (m in c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")) {
    plot.permutation.null.revised(test.out, m, outdir, prefix)
  }
}

# Override one-example output to include the DIG-style plots plus the new pathway table.
save.one.dataset.outputs <- function(ds, test.out, threshold.info, logrank.results,
                                     comp.stats, outdir, prefix = "example") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  obs <- test.out$observed

  write.csv(obs$weighted.death.first, file.path(outdir, paste0(prefix, "_weighted_death_first_curve.csv")), row.names = FALSE)
  write.csv(obs$weighted.hosp.first, file.path(outdir, paste0(prefix, "_weighted_hospitalization_first_curve.csv")), row.names = FALSE)
  write.csv(obs$threshold, file.path(outdir, paste0(prefix, "_threshold_curve.csv")), row.names = FALSE)
  write.csv(threshold.info$diagnostics, file.path(outdir, paste0(prefix, "_threshold_diagnostics_first_endpoint_only.csv")), row.names = FALSE)
  write.csv(threshold.info$threshold.table, file.path(outdir, paste0(prefix, "_threshold_grid.csv")), row.names = FALSE)
  write.csv(logrank.results, file.path(outdir, paste0(prefix, "_logrank_results.csv")), row.names = FALSE)
  write.csv(comp.stats, file.path(outdir, paste0(prefix, "_composite_statistics.csv")), row.names = FALSE)
  write.csv(make.pvalue.table(test.out, logrank.results), file.path(outdir, paste0(prefix, "_pvalue_table.csv")), row.names = FALSE)
  write.csv(make.selected.max.table(test.out), file.path(outdir, paste0(prefix, "_selected_max_parameter_table.csv")), row.names = FALSE)
  write.csv(make.pointwise.table(test.out), file.path(outdir, paste0(prefix, "_pointwise_pvalue_and_tie_table.csv")), row.names = FALSE)

  pathway.tab <- make.pathway.results.table(test.out, logrank.results)
  write.csv(pathway.tab, file.path(outdir, paste0(prefix, "_pathway_results_table.csv")), row.names = FALSE)

  max.vs.fixed.table <- make.max.vs.fixed.table(test.out)
  write.csv(max.vs.fixed.table, file.path(outdir, paste0(prefix, "_max_vs_fixed_parameter_pvalues.csv")), row.names = FALSE)

  plot.weighted.death.first(obs, outdir, prefix)
  plot.order.curves(obs, outdir, prefix)
  plot.threshold.curve(obs, outdir, prefix, p.value = test.out$p.value.max["maxWRt"])
  plot.max.vs.fixed.pvalues(max.vs.fixed.table, outdir, prefix)
  plot.km.with.threshold(ds, obs, outdir, prefix)
  plot.example.pvalue.figures(pathway.tab, outdir, prefix)
  plot.example.permutation.nulls(test.out, outdir, prefix)
}

############################################################
# Scenario summary tables and scenario-level figures
############################################################

summarise.scenario.results <- function(res, alpha = ALPHA) {
  ok <- res[is.na(res$error.message) | res$error.message == "", , drop = FALSE]
  if (nrow(ok) == 0) return(list(power = data.frame(), averages = data.frame()))

  method <- c(
    "ordinaryWR",
    "maxWRp_primary", "fixedWRp_primary",
    "maxWRp_low", "fixedWRp_low",
    "maxWRp_full", "fixedWRp_full",
    "maxOrderWR_primary", "fixedOrderWR_primary",
    "maxOrderWR_full", "fixedOrderWR_full",
    "maxWRt", "fixedWRt",
    "logrank_death", "logrank_composite"
  )
  pcols <- c(
    "max.pvalue.ordinaryWR",
    "max.pvalue.maxWRp_primary", "fixed.pvalue.maxWRp_primary",
    "max.pvalue.maxWRp_low", "fixed.pvalue.maxWRp_low",
    "max.pvalue.maxWRp_full", "fixed.pvalue.maxWRp_full",
    "max.pvalue.maxOrderWR_primary", "fixed.pvalue.maxOrderWR_primary",
    "max.pvalue.maxOrderWR_full", "fixed.pvalue.maxOrderWR_full",
    "max.pvalue.maxWRt", "fixed.pvalue.maxWRt",
    "logrank.death.p", "logrank.composite.p"
  )

  power <- data.frame(
    method = method,
    method_label = method.long.label(method),
    rejection.proportion = sapply(pcols, function(cc) mean(safe.col(ok, cc) < alpha, na.rm = TRUE)),
    alpha = alpha,
    nsim.available = nrow(ok),
    stringsAsFactors = FALSE
  )

  averages <- data.frame(
    nsim.available = nrow(ok),
    mean.ordinaryWR = safe.mean(ok$ordinaryWR),
    mean.max.pvalue.ordinaryWR = safe.mean(ok$max.pvalue.ordinaryWR),
    mean.fixed.pvalue.ordinaryWR = safe.mean(ok$fixed.pvalue.ordinaryWR),
    mean.traditional.tie.count = safe.mean(ok$traditional.tie.count),
    mean.traditional.tie.pr = safe.mean(ok$traditional.tie.pr),

    mean.maxWRp.primary = safe.mean(ok$maxWRp_primary),
    mean.selected.p.primary = safe.mean(ok$selected.p.primary),
    mean.max.pvalue.maxWRp.primary = safe.mean(ok$max.pvalue.maxWRp_primary),
    mean.fixed.pvalue.maxWRp.primary = safe.mean(ok$fixed.pvalue.maxWRp_primary),
    mean.weighted.primary.tie.count = safe.mean(ok$weighted.primary.tie.count),
    mean.weighted.primary.tie.pr = safe.mean(ok$weighted.primary.tie.pr),

    mean.maxWRp.low = safe.mean(ok$maxWRp_low),
    mean.selected.p.low = safe.mean(ok$selected.p.low),
    mean.max.pvalue.maxWRp.low = safe.mean(ok$max.pvalue.maxWRp_low),
    mean.fixed.pvalue.maxWRp.low = safe.mean(ok$fixed.pvalue.maxWRp_low),
    mean.weighted.low.tie.count = safe.mean(ok$weighted.low.tie.count),
    mean.weighted.low.tie.pr = safe.mean(ok$weighted.low.tie.pr),

    mean.maxWRp.full = safe.mean(ok$maxWRp_full),
    mean.selected.p.full = safe.mean(ok$selected.p.full),
    prop.full.selected.p.below.0.5 = mean(ok$selected.full.is.low.p == 1, na.rm = TRUE),
    mean.max.pvalue.maxWRp.full = safe.mean(ok$max.pvalue.maxWRp_full),
    mean.fixed.pvalue.maxWRp.full = safe.mean(ok$fixed.pvalue.maxWRp_full),
    mean.weighted.full.tie.count = safe.mean(ok$weighted.full.tie.count),
    mean.weighted.full.tie.pr = safe.mean(ok$weighted.full.tie.pr),

    mean.maxOrderWR.primary = safe.mean(ok$maxOrderWR_primary),
    prop.order.changed.primary = mean(ok$order.changed.primary == 1, na.rm = TRUE),
    mean.order.diff.hosp.minus.death = safe.mean(ok$order.diff.hosp.minus.death),
    mean.selected.order.p.primary = safe.mean(ok$selected.order.p.primary),
    mode.selected.order.primary = clean.order.label(mode.string(ok$selected.order.primary)),
    mean.max.pvalue.maxOrderWR.primary = safe.mean(ok$max.pvalue.maxOrderWR_primary),
    mean.fixed.pvalue.maxOrderWR.primary = safe.mean(ok$fixed.pvalue.maxOrderWR_primary),
    mean.order.primary.tie.count = safe.mean(ok$order.primary.tie.count),
    mean.order.primary.tie.pr = safe.mean(ok$order.primary.tie.pr),

    mean.maxOrderWR.full = safe.mean(ok$maxOrderWR_full),
    mean.selected.order.p.full = safe.mean(ok$selected.order.p.full),
    mode.selected.order.full = clean.order.label(mode.string(ok$selected.order.full)),
    prop.order.full.selected.p.below.0.5 = mean(ok$selected.order.full.is.low.p == 1, na.rm = TRUE),
    mean.max.pvalue.maxOrderWR.full = safe.mean(ok$max.pvalue.maxOrderWR_full),
    mean.fixed.pvalue.maxOrderWR.full = safe.mean(ok$fixed.pvalue.maxOrderWR_full),
    mean.order.full.tie.count = safe.mean(ok$order.full.tie.count),
    mean.order.full.tie.pr = safe.mean(ok$order.full.tie.pr),

    mean.maxWRt = safe.mean(ok$maxWRt),
    mean.selected.t.months = safe.mean(ok$selected.t.months),
    sd.selected.t.months = safe_sd(ok$selected.t.months),
    mean.max.pvalue.maxWRt = safe.mean(ok$max.pvalue.maxWRt),
    mean.fixed.pvalue.maxWRt = safe.mean(ok$fixed.pvalue.maxWRt),
    mean.threshold.tie.count = safe.mean(ok$threshold.tie.count),
    mean.threshold.tie.pr = safe.mean(ok$threshold.tie.pr),

    mean.true.hierarchical.tie.count = safe.mean(ok$true.hierarchical.tie.count),
    mean.true.hierarchical.tie.pr = safe.mean(ok$true.hierarchical.tie.pr),
    mean.logrank.death.statistic = safe.mean(ok$logrank.death.statistic),
    mean.logrank.death.p = safe.mean(ok$logrank.death.p),
    mean.logrank.composite.statistic = safe.mean(ok$logrank.composite.statistic),
    mean.logrank.composite.p = safe.mean(ok$logrank.composite.p),
    mean.death.event.rate.control = safe.mean(ok$death.event.rate.control),
    mean.death.event.rate.treatment = safe.mean(ok$death.event.rate.treatment),
    mean.composite.event.rate.control = safe.mean(ok$composite.event.rate.control),
    mean.composite.event.rate.treatment = safe.mean(ok$composite.event.rate.treatment),
    mean.num.hosp.control = safe.mean(ok$mean.num.hosp.control),
    mean.num.hosp.treatment = safe.mean(ok$mean.num.hosp.treatment),
    stringsAsFactors = FALSE
  )

  list(power = power, averages = averages)
}

make.scenario.method.comparison <- function(res, power.table, alpha = ALPHA) {
  ok <- res[is.na(res$error.message) | res$error.message == "", , drop = FALSE]
  if (nrow(ok) == 0) return(data.frame())

  specs <- data.frame(
    method = c(
      "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
      "maxOrderWR_primary", "maxOrderWR_full", "maxWRt",
      "logrank_death", "logrank_composite"
    ),
    pathway = c(
      "Traditional WR", "Weighted WR", "Weighted WR", "Weighted WR",
      "Maximum-order WR", "Maximum-order WR", "Threshold WR",
      "Log-rank", "Log-rank"
    ),
    statistic.column = c(
      "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
      "maxOrderWR_primary", "maxOrderWR_full", "maxWRt",
      "logrank.death.statistic", "logrank.composite.statistic"
    ),
    permutation.p.column = c(
      "max.pvalue.ordinaryWR", "max.pvalue.maxWRp_primary", "max.pvalue.maxWRp_low", "max.pvalue.maxWRp_full",
      "max.pvalue.maxOrderWR_primary", "max.pvalue.maxOrderWR_full", "max.pvalue.maxWRt",
      "logrank.death.p", "logrank.composite.p"
    ),
    selected.pvalue.column = c(
      "fixed.pvalue.ordinaryWR", "fixed.pvalue.maxWRp_primary", "fixed.pvalue.maxWRp_low", "fixed.pvalue.maxWRp_full",
      "fixed.pvalue.maxOrderWR_primary", "fixed.pvalue.maxOrderWR_full", "fixed.pvalue.maxWRt",
      NA_character_, NA_character_
    ),
    selected.p.column = c(
      NA_character_, "selected.p.primary", "selected.p.low", "selected.p.full",
      "selected.order.p.primary", "selected.order.p.full", NA_character_,
      NA_character_, NA_character_
    ),
    selected.t.column = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "selected.t.months", NA_character_, NA_character_),
    selected.order.column = c(NA_character_, NA_character_, NA_character_, NA_character_, "selected.order.primary", "selected.order.full", NA_character_, NA_character_, NA_character_),
    tie.count.column = c(
      "traditional.tie.count", "weighted.primary.tie.count", "weighted.low.tie.count", "weighted.full.tie.count",
      "order.primary.tie.count", "order.full.tie.count", "threshold.tie.count",
      NA_character_, NA_character_
    ),
    tie.pr.column = c(
      "traditional.tie.pr", "weighted.primary.tie.pr", "weighted.low.tie.pr", "weighted.full.tie.pr",
      "order.primary.tie.pr", "order.full.tie.pr", "threshold.tie.pr",
      NA_character_, NA_character_
    ),
    fixed.power.method = c(
      "ordinaryWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
      "fixedOrderWR_primary", "fixedOrderWR_full", "fixedWRt",
      NA_character_, NA_character_
    ),
    stringsAsFactors = FALSE
  )

  rows <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i, , drop = FALSE]
    stat <- safe.col(ok, sp$statistic.column)
    pp <- safe.col(ok, sp$permutation.p.column)
    fixed.pp <- if (!is.na(sp$selected.pvalue.column)) safe.col(ok, sp$selected.pvalue.column) else rep(NA_real_, nrow(ok))
    sel.p <- if (!is.na(sp$selected.p.column)) safe.col(ok, sp$selected.p.column) else rep(NA_real_, nrow(ok))
    sel.t <- if (!is.na(sp$selected.t.column)) safe.col(ok, sp$selected.t.column) else rep(NA_real_, nrow(ok))
    tie.count <- if (!is.na(sp$tie.count.column)) safe.col(ok, sp$tie.count.column) else rep(NA_real_, nrow(ok))
    tie.pr <- if (!is.na(sp$tie.pr.column)) safe.col(ok, sp$tie.pr.column) else rep(NA_real_, nrow(ok))
    sel.order <- if (!is.na(sp$selected.order.column) && sp$selected.order.column %in% names(ok)) as.character(ok[[sp$selected.order.column]]) else rep(NA_character_, nrow(ok))

    selected.summary <- NA_character_
    if (sp$method == "ordinaryWR") selected.summary <- "p = 0.50; death first"
    if (sp$method %in% c("maxWRp_primary", "maxWRp_low", "maxWRp_full")) {
      selected.summary <- paste0("mean p = ", sprintf("%.2f", safe.mean(sel.p)), "; death first")
    }
    if (sp$method %in% c("maxOrderWR_primary", "maxOrderWR_full")) {
      selected.summary <- paste0("mode order = ", clean.order.label(mode.string(sel.order)),
                                 "; mean p = ", sprintf("%.2f", safe.mean(sel.p)))
    }
    if (sp$method == "maxWRt") {
      selected.summary <- paste0("mean t = ", sprintf("%.1f", safe.mean(sel.t)), " months")
    }

    rows[[i]] <- data.frame(
      scenario_id = ok$scenario_id[1],
      scenario_label = scenario.display.label(ok$scenario_id[1]),
      scenario_type = ok$scenario_type[1],
      description = ok$description[1],
      method = sp$method,
      method_label = method.long.label(sp$method),
      pathway = sp$pathway,
      statistic_type = ifelse(grepl("logrank", sp$method), "Log-rank chi-square", "Win ratio"),
      nsim.available = nrow(ok),
      mean_statistic = safe.mean(stat),
      sd_statistic = safe_sd(stat),
      median_statistic = safe.median(stat),
      selected_value_summary = selected.summary,
      mean_selected_p = safe.mean(sel.p),
      median_selected_p = safe.median(sel.p),
      mean_selected_t_months = safe.mean(sel.t),
      median_selected_t_months = safe.median(sel.t),
      selected_order_mode = clean.order.label(mode.string(sel.order)),
      prop_hospitalization_first = ifelse(all(is.na(sel.order)), NA_real_, mean(sel.order == "hospitalization_first", na.rm = TRUE)),
      mean_permutation_p_value = safe.mean(pp),
      median_permutation_p_value = safe.median(pp),
      power_permutation_p_value = mean(pp < alpha, na.rm = TRUE),
      mean_selected_parameter_p_value = safe.mean(fixed.pp),
      median_selected_parameter_p_value = safe.median(fixed.pp),
      power_selected_parameter_p_value = mean(fixed.pp < alpha, na.rm = TRUE),
      power_selected_parameter_from_power_table = get.power.value(power.table, sp$fixed.power.method),
      mean_tie_count = safe.mean(tie.count),
      mean_tie_proportion = safe.mean(tie.pr),
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

plot.scenario.method.comparison <- function(method.table, outdir, prefix) {
  if (nrow(method.table) == 0) return(invisible(NULL))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  title.label <- scenario.display.label(prefix)

  # Scenario-level permutation p-value comparison for all pathways and log-rank tests.
  pvalue.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite")
  pvtab <- method.table[method.table$method %in% pvalue.methods, , drop = FALSE]
  if (nrow(pvtab) > 0) {
    pvtab$method <- factor(pvtab$method, levels = pvalue.methods)
    pvtab <- pvtab[order(pvtab$method), , drop = FALSE]
    save_png(outdir, paste0(prefix, "_permutation_pvalue_comparison.png"), width = 2300, height = 1200)
    old.par <- par(no.readonly = TRUE)
    par(mar = c(8.8, 5.2, 3.2, 1.2))
    bp <- barplot(pvtab$mean_permutation_p_value,
                  names.arg = method.display.label(pvtab$method),
                  ylim = c(0, 1), las = 1,
                  ylab = "Mean permutation / log-rank p-value",
                  main = paste0(title.label, ": pathway p-value comparison"),
                  cex.names = 0.82)
    abline(h = ALPHA, lty = 2, lwd = 2)
    text(bp, pmin(pvtab$mean_permutation_p_value + 0.04, 0.98),
         labels = sprintf("%.3f", pvtab$mean_permutation_p_value), cex = 0.68)
    par(old.par)
    dev.off()
  }

  # Scenario-level max-statistic vs fixed-selected p-value comparison for WR pathways.
  fixed.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  ftab <- method.table[method.table$method %in% fixed.methods, , drop = FALSE]
  if (nrow(ftab) > 0) {
    ftab$method <- factor(ftab$method, levels = fixed.methods)
    ftab <- ftab[order(ftab$method), , drop = FALSE]
    mat <- rbind(
      "Max-statistic" = ftab$mean_permutation_p_value,
      "Fixed-selected" = ftab$mean_selected_parameter_p_value
    )
    colnames(mat) <- method.display.label(ftab$method)
    save_png(outdir, paste0(prefix, "_max_vs_fixed_selected_pvalues.png"), width = 2300, height = 1200)
    old.par <- par(no.readonly = TRUE)
    par(mar = c(8.8, 5.2, 3.2, 1.2))
    bp <- barplot(mat, beside = TRUE, ylim = c(0, 1), las = 1,
                  ylab = "Mean p-value",
                  main = paste0(title.label, ": max-statistic vs fixed-selected p-values"),
                  legend.text = TRUE,
                  args.legend = list(x = "topright", bty = "n"),
                  cex.names = 0.82)
    abline(h = ALPHA, lty = 2, lwd = 2)
    text(bp, pmin(mat + 0.04, 0.98), labels = sprintf("%.3f", mat), cex = 0.60)
    par(old.par)
    dev.off()
  }

  # Scenario-level power comparison for all pathways and log-rank tests.
  powtab <- pvtab
  if (nrow(powtab) > 0) {
    save_png(outdir, paste0(prefix, "_power_comparison.png"), width = 2300, height = 1200)
    old.par <- par(no.readonly = TRUE)
    par(mar = c(8.8, 5.2, 3.2, 1.2))
    bp <- barplot(powtab$power_permutation_p_value,
                  names.arg = method.display.label(powtab$method),
                  ylim = c(0, 1), las = 1,
                  ylab = paste0("Rejection proportion at alpha = ", ALPHA),
                  main = paste0(title.label, ": power / rejection proportion"),
                  cex.names = 0.82)
    abline(h = ALPHA, lty = 2, lwd = 2)
    text(bp, pmin(powtab$power_permutation_p_value + 0.04, 0.98),
         labels = sprintf("%.2f", powtab$power_permutation_p_value), cex = 0.68)
    par(old.par)
    dev.off()
  }
}

# Disable pointwise p-value/tie graphics. The tables remain available.
plot.scenario.pointwise.summary <- function(pointwise.summary, outdir, prefix) {
  invisible(NULL)
}

run.scenario <- function(scenario.row,
                         nsim = NSIM,
                         B = B_PERM,
                         master.seed = MASTER_SEED,
                         outdir = OUTDIR,
                         resume = RESUME_IF_EXISTS,
                         checkpoint.every = CHECKPOINT_EVERY,
                         verbose = VERBOSE,
                         scenario.index = NULL) {
  scenario.id <- scenario.row$scenario_id
  scenario.outdir <- file.path(outdir, scenario.id)
  if (!dir.exists(scenario.outdir)) dir.create(scenario.outdir, recursive = TRUE)

  raw.rds <- file.path(scenario.outdir, paste0(scenario.id, "_raw_results.rds"))
  raw.csv <- file.path(scenario.outdir, paste0(scenario.id, "_raw_results.csv"))
  pointwise.rds <- file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.rds"))
  pointwise.csv <- file.path(scenario.outdir, paste0(scenario.id, "_pointwise_results.csv"))

  if (resume && file.exists(raw.rds)) {
    res <- readRDS(raw.rds)
    completed <- if (nrow(res) == 0) 0 else max(res$sim, na.rm = TRUE)
    pointwise.res <- if (file.exists(pointwise.rds)) readRDS(pointwise.rds) else data.frame()
    if (verbose) cat("Resuming", scenario.id, "from simulation", completed + 1, "\n")
  } else {
    res <- data.frame()
    pointwise.res <- data.frame()
    completed <- 0
  }

  if (!is.null(scenario.index) && !is.na(scenario.index)) {
    scenario.seed.index <- scenario.index
  } else {
    default.ids <- build.scenario.grid()$scenario_id
    scenario.seed.index <- match(scenario.id, default.ids)
    if (is.na(scenario.seed.index)) {
      scenario.seed.index <- (sum(utf8ToInt(as.character(scenario.id))) %% 10000) + 1
    }
  }

  if (completed >= nsim) {
    if (verbose) cat("Scenario", scenario.id, "already complete.\n")
  } else {
    for (s in seq.int(completed + 1, nsim)) {
      if (verbose) cat("Scenario", scenario.id, ": simulation", s, "of", nsim, "\n")
      sim.seed <- master.seed + scenario.seed.index * 1000000 + s * 10007

      sim.out <- tryCatch(
        run.one.simulation(
          scenario.row = scenario.row,
          sim.index = s,
          seed = sim.seed,
          B = B,
          save.example = SAVE_EXAMPLE_PLOTS && s == 1,
          scenario.outdir = scenario.outdir
        ),
        error = function(e) {
          warning("Simulation failed: ", scenario.id, " sim ", s, ": ", conditionMessage(e))
          list(
            summary.row = data.frame(
              scenario_id = scenario.id,
              scenario_type = get_scenario_field(scenario.row, "scenario_type"),
              description = get_scenario_field(scenario.row, "description"),
              sim = s,
              error.message = conditionMessage(e),
              stringsAsFactors = FALSE
            ),
            pointwise.rows = data.frame()
          )
        }
      )

      res <- rbind_fill_base(res, sim.out$summary.row)
      if (nrow(sim.out$pointwise.rows) > 0) {
        pointwise.res <- rbind_fill_base(pointwise.res, sim.out$pointwise.rows)
      }

      if (s %% checkpoint.every == 0 || s == nsim) {
        saveRDS(res, raw.rds)
        write.csv(res, raw.csv, row.names = FALSE)
        saveRDS(pointwise.res, pointwise.rds)
        write.csv(pointwise.res, pointwise.csv, row.names = FALSE)
      }
    }
  }

  summary.out <- summarise.scenario.results(res, alpha = ALPHA)
  pointwise.summary <- summarise.pointwise.results(pointwise.res, alpha = ALPHA)
  method.comparison <- make.scenario.method.comparison(res, summary.out$power, alpha = ALPHA)

  if (nrow(summary.out$power) > 0) {
    summary.out$power$scenario_id <- scenario.id
    summary.out$power$scenario_label <- scenario.display.label(scenario.id)
    summary.out$power$scenario_type <- get_scenario_field(scenario.row, "scenario_type")
    summary.out$power$description <- get_scenario_field(scenario.row, "description")
    summary.out$averages$scenario_id <- scenario.id
    summary.out$averages$scenario_label <- scenario.display.label(scenario.id)
    summary.out$averages$scenario_type <- get_scenario_field(scenario.row, "scenario_type")
    summary.out$averages$description <- get_scenario_field(scenario.row, "description")
    write.csv(summary.out$power, file.path(scenario.outdir, paste0(scenario.id, "_power_summary.csv")), row.names = FALSE)
    write.csv(summary.out$averages, file.path(scenario.outdir, paste0(scenario.id, "_average_statistics.csv")), row.names = FALSE)
  }

  if (nrow(method.comparison) > 0) {
    write.csv(method.comparison, file.path(scenario.outdir, paste0(scenario.id, "_pathway_method_comparison.csv")), row.names = FALSE)
    plot.scenario.method.comparison(method.comparison, scenario.outdir, scenario.id)
  }

  if (nrow(pointwise.summary) > 0) {
    pointwise.summary$scenario_id <- scenario.id
    pointwise.summary$scenario_label <- scenario.display.label(scenario.id)
    pointwise.summary$scenario_type <- get_scenario_field(scenario.row, "scenario_type")
    pointwise.summary$description <- get_scenario_field(scenario.row, "description")
    write.csv(pointwise.summary, file.path(scenario.outdir, paste0(scenario.id, "_pointwise_average_pvalue_tie_summary.csv")), row.names = FALSE)
  }

  list(
    raw = res,
    pointwise = pointwise.res,
    power = summary.out$power,
    averages = summary.out$averages,
    method.comparison = method.comparison,
    pointwise.summary = pointwise.summary
  )
}


# Global figures
plot.global.power.clean <- function(power.table, outdir) {
  if (nrow(power.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  tab <- power.table[power.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  scenarios <- unique(tab$scenario_id)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios),
                dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  for (i in seq_len(nrow(tab))) {
    mat[method.display.label(tab$method[i]), scenario.display.label(tab$scenario_id[i])] <- tab$rejection.proportion[i]
  }
  save_png(outdir, "GLOBAL_power_comparison_key_methods.png", width = 3100, height = 1500)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2,
          ylab = paste0("Rejection proportion at alpha = ", ALPHA),
          main = "Power / rejection proportion across scenarios",
          legend.text = TRUE,
          args.legend = list(x = "topright", bty = "n", cex = 0.72),
          cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2)
  par(old.par)
  dev.off()
}

plot.global.permutation.pvalues.clean <- function(method.table, outdir) {
  if (nrow(method.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  tab <- method.table[method.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  scenarios <- unique(tab$scenario_id)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios),
                dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  for (i in seq_len(nrow(tab))) {
    mat[method.display.label(tab$method[i]), scenario.display.label(tab$scenario_id[i])] <- tab$mean_permutation_p_value[i]
  }
  save_png(outdir, "GLOBAL_permutation_pvalue_comparison.png", width = 3100, height = 1500)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2,
          ylab = "Mean permutation / log-rank p-value",
          main = "Permutation p-value comparison across scenarios",
          legend.text = TRUE,
          args.legend = list(x = "topright", bty = "n", cex = 0.72),
          cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2)
  par(old.par)
  dev.off()
}

plot.global.max.fixed.pvalues.clean <- function(method.table, outdir) {
  if (nrow(method.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")
  tab <- method.table[method.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  scenarios <- unique(tab$scenario_id)
  mat.max <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios),
                    dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  mat.fix <- mat.max
  for (i in seq_len(nrow(tab))) {
    rr <- method.display.label(tab$method[i])
    cc <- scenario.display.label(tab$scenario_id[i])
    mat.max[rr, cc] <- tab$mean_permutation_p_value[i]
    mat.fix[rr, cc] <- tab$mean_selected_parameter_p_value[i]
  }
  # Two panels are intentionally saved as separate files to keep labels readable.
  save_png(outdir, "GLOBAL_fixed_selected_pvalue_comparison.png", width = 3000, height = 1450)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat.fix, beside = TRUE, ylim = c(0, 1), las = 2,
          ylab = "Mean fixed-selected permutation p-value",
          main = "Fixed-selected p-value comparison across scenarios",
          legend.text = TRUE,
          args.legend = list(x = "topright", bty = "n", cex = 0.75),
          cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2)
  par(old.par)
  dev.off()
}

plot.global.low.p.clean <- function(avg.table, outdir) {
  if (nrow(avg.table) == 0 || !("prop.full.selected.p.below.0.5" %in% names(avg.table))) return(invisible(NULL))
  labs <- if ("scenario_label" %in% names(avg.table)) avg.table$scenario_label else scenario.display.label(avg.table$scenario_id)
  save_png(outdir, "GLOBAL_low_p_selection_frequency.png", width = 2500, height = 1300)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  bp <- barplot(avg.table$prop.full.selected.p.below.0.5,
                names.arg = labs,
                ylim = c(0, 1), las = 2,
                ylab = "Proportion of simulations",
                main = "How often full-grid maximum selects p < 0.5",
                cex.names = 0.84)
  text(bp, pmin(avg.table$prop.full.selected.p.below.0.5 + 0.04, 0.98),
       labels = sprintf("%.2f", avg.table$prop.full.selected.p.below.0.5), cex = 0.75)
  par(old.par)
  dev.off()
}

run.all.scenarios <- function(scenario.grid = build.scenario.grid(),
                              nsim = NSIM,
                              B = B_PERM,
                              outdir = OUTDIR,
                              master.seed = MASTER_SEED) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  if (!("scenario_label" %in% names(scenario.grid))) {
    scenario.grid$scenario_label <- scenario.display.label(scenario.grid$scenario_id)
  }
  write.csv(scenario.grid, file.path(outdir, "GLOBAL_scenario_config.csv"), row.names = FALSE)

  settings.table <- data.frame(
    setting = c(
      "NSIM", "B_PERM", "MASTER_SEED", "ALPHA", "number_of_scenarios",
      "p_grid_low_exploratory", "p_grid_primary", "p_grid_full",
      "clinical_t_months", "checkpoint_every", "resume_if_exists", "save_example_plots"
    ),
    value = c(
      as.character(nsim),
      as.character(B),
      as.character(master.seed),
      as.character(ALPHA),
      as.character(nrow(scenario.grid)),
      paste0(sprintf("%.2f", min(P_GRID_EXPLORATORY)), " to ", sprintf("%.2f", max(P_GRID_EXPLORATORY))),
      paste0(sprintf("%.2f", min(P_GRID_PRIMARY)), " to ", sprintf("%.2f", max(P_GRID_PRIMARY))),
      paste0(sprintf("%.2f", min(P_GRID_ALL)), " to ", sprintf("%.2f", max(P_GRID_ALL)), "; no p = 0"),
      paste(CLINICAL_T_MONTHS, collapse = ", "),
      as.character(CHECKPOINT_EVERY),
      as.character(RESUME_IF_EXISTS),
      as.character(SAVE_EXAMPLE_PLOTS)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(settings.table, file.path(outdir, "GLOBAL_settings.csv"), row.names = FALSE)

  cat("\n===== Simulation settings =====\n")
  cat("Output folder:", outdir, "\n")
  cat("NSIM =", nsim, "; B_PERM =", B, "; scenarios =", nrow(scenario.grid), "\n")
  cat("Primary weighted p grid: [", min(P_GRID_PRIMARY), ", ", max(P_GRID_PRIMARY), "]\n", sep = "")
  cat("Exploratory low-p grid: [", min(P_GRID_EXPLORATORY), ", ", max(P_GRID_EXPLORATORY), "]\n", sep = "")

  all.raw <- data.frame()
  all.pointwise <- data.frame()
  all.power <- data.frame()
  all.averages <- data.frame()
  all.method.comparison <- data.frame()
  all.pointwise.summary <- data.frame()

  for (i in seq_len(nrow(scenario.grid))) {
    scenario.row <- scenario.grid[i, , drop = FALSE]
    out.i <- run.scenario(
      scenario.row = scenario.row,
      nsim = nsim,
      B = B,
      master.seed = master.seed,
      outdir = outdir,
      resume = RESUME_IF_EXISTS,
      checkpoint.every = CHECKPOINT_EVERY,
      verbose = VERBOSE,
      scenario.index = i
    )

    all.raw <- rbind_fill_base(all.raw, out.i$raw)
    all.pointwise <- rbind_fill_base(all.pointwise, out.i$pointwise)
    all.power <- rbind_fill_base(all.power, out.i$power)
    all.averages <- rbind_fill_base(all.averages, out.i$averages)
    all.method.comparison <- rbind_fill_base(all.method.comparison, out.i$method.comparison)
    all.pointwise.summary <- rbind_fill_base(all.pointwise.summary, out.i$pointwise.summary)

    write.csv(all.raw, file.path(outdir, "GLOBAL_all_scenarios_raw_results.csv"), row.names = FALSE)
    write.csv(all.pointwise, file.path(outdir, "GLOBAL_all_scenarios_pointwise_results.csv"), row.names = FALSE)
    write.csv(all.power, file.path(outdir, "GLOBAL_all_scenarios_power_summary.csv"), row.names = FALSE)
    write.csv(all.averages, file.path(outdir, "GLOBAL_all_scenarios_average_statistics.csv"), row.names = FALSE)
    write.csv(all.method.comparison, file.path(outdir, "GLOBAL_all_scenarios_pathway_method_comparison.csv"), row.names = FALSE)
    write.csv(all.pointwise.summary, file.path(outdir, "GLOBAL_all_scenarios_pointwise_average_pvalue_tie_summary.csv"), row.names = FALSE)
  }

  plot.global.power.clean(all.power, outdir)
  plot.global.permutation.pvalues.clean(all.method.comparison, outdir)
  plot.global.max.fixed.pvalues.clean(all.method.comparison, outdir)
  plot.global.low.p.clean(all.averages, outdir)

  manifest <- data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE)
  write.csv(manifest, file.path(outdir, "GLOBAL_output_manifest.csv"), row.names = FALSE)

  saveRDS(
    list(
      scenario.grid = scenario.grid,
      raw = all.raw,
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
        P_GRID_EXPLORATORY = P_GRID_EXPLORATORY,
        P_GRID_PRIMARY = P_GRID_PRIMARY,
        P_GRID_ALL = P_GRID_ALL,
        CLINICAL_T_MONTHS = CLINICAL_T_MONTHS
      )
    ),
    file.path(outdir, "GLOBAL_full_results_bundle.rds")
  )

  cat("\n===== All scenarios complete =====\n")
  cat("Output folder:", outdir, "\n")
  cat("Full-run settings: NSIM =", nsim, ", B_PERM =", B, "\n")
  cat("Scenario count:", nrow(scenario.grid), "\n")
  cat("Main tables:\n")
  cat("  GLOBAL_scenario_config.csv\n")
  cat("  GLOBAL_settings.csv\n")
  cat("  GLOBAL_all_scenarios_raw_results.csv\n")
  cat("  GLOBAL_all_scenarios_power_summary.csv\n")
  cat("  GLOBAL_all_scenarios_average_statistics.csv\n")
  cat("  GLOBAL_all_scenarios_pathway_method_comparison.csv\n")
  cat("  GLOBAL_all_scenarios_pointwise_results.csv\n")
  cat("  GLOBAL_all_scenarios_pointwise_average_pvalue_tie_summary.csv\n")
  cat("Main figures:\n")
  cat("  GLOBAL_power_comparison_key_methods.png\n")
  cat("  GLOBAL_permutation_pvalue_comparison.png\n")
  cat("  GLOBAL_fixed_selected_pvalue_comparison.png\n")
  cat("  GLOBAL_low_p_selection_frequency.png\n")
  cat("Each scenario folder contains:\n")
  cat("  Sxx_pathway_method_comparison.csv\n")
  cat("  Sxx_average_statistics.csv  # table only, no average-statistic plots\n")
  cat("  Sxx_permutation_pvalue_comparison.png\n")
  cat("  Sxx_max_vs_fixed_selected_pvalues.png\n")
  cat("  Sxx_power_comparison.png\n")
  cat("  example_sim1_* DIG-style figures and tables\n")

  invisible(list(
    raw = all.raw,
    pointwise = all.pointwise,
    power = all.power,
    averages = all.averages,
    method.comparison = all.method.comparison,
    pointwise.summary = all.pointwise.summary
  ))
}

#Traditional order-sensitivity add-on
#Fixed p = 0.50, compare death-first vs hospitalization-first.
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
    wins = as.numeric(core$t_wins),
    losses = as.numeric(core$t_losses),
    pr.win = as.numeric(core$pr_win_t),
    pr.loss = as.numeric(core$pr_loss_t),
    pr.tie = as.numeric(core$pr_tie_t),
    stringsAsFactors = FALSE
  )
  threshold$tie.count <- total.pairs - threshold$wins - threshold$losses
  best.t.idx <- which.max(threshold$WR)

  p.low.values <- p.grid[p.grid < 0.5]
  p.primary.values <- p.grid[p.grid >= 0.5]
  p.low.lower <- if (length(p.low.values) > 0) min(p.low.values) else NA_real_
  p.low.upper <- if (length(p.low.values) > 0) max(p.low.values) else NA_real_
  p.primary.lower <- if (length(p.primary.values) > 0) min(p.primary.values) else NA_real_
  p.primary.upper <- if (length(p.primary.values) > 0) max(p.primary.values) else NA_real_

  traditional.death.first <- max_from_curve(weighted.death.first, "ordinaryWR_death_first", 0.50, 0.50)
  traditional.hosp.first <- max_from_curve(weighted.hosp.first, "ordinaryWR_hosp_first", 0.50, 0.50)
  max.traditional.order <- max_from_curve(weighted.all.orders, "traditionalOrderWR", 0.50, 0.50)

  list(
    ordinaryWR = as.numeric(core$ordinaryWR),
    ordinary.win.score = as.numeric(core$ordinary_win_score),
    ordinary.loss.score = as.numeric(core$ordinary_loss_score),
    traditional.death.first = traditional.death.first,
    traditional.hosp.first = traditional.hosp.first,
    max.traditional.order = max.traditional.order,
    weighted.death.first = weighted.death.first,
    weighted.hosp.first = weighted.hosp.first,
    weighted.all.orders = weighted.all.orders,
    max.weighted.primary = max_from_curve(weighted.death.first, "maxWRp_primary_death_first", p.primary.lower, p.primary.upper),
    max.weighted.low = max_from_curve(weighted.death.first, "maxWRp_low_death_first", p.low.lower, p.low.upper),
    max.weighted.full = max_from_curve(weighted.death.first, "maxWRp_full_death_first", min(p.grid), max(p.grid)),
    max.order.primary = max_from_curve(weighted.all.orders, "maxOrderWR_primary", p.primary.lower, p.primary.upper),
    max.order.full = max_from_curve(weighted.all.orders, "maxOrderWR_full", min(p.grid), max(p.grid)),
    threshold = threshold,
    max.threshold = threshold[best.t.idx, , drop = FALSE],
    counts = list(
      total.pairs = as.numeric(core$total_pairs),
      D1.win = as.numeric(core$D1_win),
      D1.loss = as.numeric(core$D1_loss),
      H2.win = as.numeric(core$H2_win),
      H2.loss = as.numeric(core$H2_loss),
      H1.win = as.numeric(core$H1_win),
      H1.loss = as.numeric(core$H1_loss),
      D2.win = as.numeric(core$D2_win),
      D2.loss = as.numeric(core$D2_loss),
      death.win.pairs = as.numeric(core$death_win_pairs),
      death.loss.pairs = as.numeric(core$death_loss_pairs),
      hosp.win.pairs = as.numeric(core$hosp_win_pairs),
      hosp.loss.pairs = as.numeric(core$hosp_loss_pairs),
      true.tie.pairs = as.numeric(core$true_tie_pairs),
      true.tie.pr = as.numeric(core$true_tie_pairs) / as.numeric(core$total_pairs)
    )
  )
}

traditional_order_stat <- function(engine.out) {
  death.wr <- wr_at_order_p(engine.out, "death_first", 0.50)
  hosp.wr <- wr_at_order_p(engine.out, "hospitalization_first", 0.50)
  selected.order <- if (!is.na(hosp.wr) && (is.na(death.wr) || hosp.wr > death.wr)) {
    "hospitalization_first"
  } else {
    "death_first"
  }
  data.frame(
    WR = max(c(death.wr, hosp.wr), na.rm = TRUE),
    order = selected.order,
    p = 0.50,
    death.first.WR = death.wr,
    hosp.first.WR = hosp.wr,
    stringsAsFactors = FALSE
  )
}

perm.test.revised <- function(ds,
                              B = B_PERM,
                              seed = MASTER_SEED,
                              p.grid = P_GRID_ALL,
                              t.grid = CLINICAL_T_MONTHS / 12,
                              verbose = FALSE) {
  ds <- prepare.ds.fast(ds)
  set.seed(seed)
  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR",
    "traditionalWR_hosp_first",
    "traditionalOrderWR",
    "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.names <- c(
    "ordinaryWR",
    "fixedTraditionalHospFirst",
    "fixedTraditionalOrderWR",
    "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )

  fixed.traditional.order <- obs$max.traditional.order$order[1]
  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    traditionalWR_hosp_first = wr_at_order_p(obs, "hospitalization_first", 0.50),
    traditionalOrderWR = obs$max.traditional.order$WR[1],
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )

  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedTraditionalHospFirst = wr_at_order_p(obs, "hospitalization_first", 0.50),
    fixedTraditionalOrderWR = wr_at_order_p(obs, fixed.traditional.order, 0.50),
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names)); colnames(perm.max) <- max.names
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names)); colnames(perm.fixed) <- fixed.names

  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.traditional.order = NA_character_,
    selected.p.primary = NA_real_, selected.p.low = NA_real_, selected.p.full = NA_real_,
    selected.order.primary = NA_character_, selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_, selected.order.p.full = NA_real_,
    selected.t = NA_real_, stringsAsFactors = FALSE
  )

  perm.seeds <- seed + seq_len(B) * 1009L
  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b == B || b %% 50 == 0)) cat("  Permutation", b, "of", B, "\n")
    set.seed(perm.seeds[b])
    ds.b <- ds
    ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)
    out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)

    perm.pointwise.death.first[b, ] <- out.b$weighted.death.first$WR
    perm.pointwise.hosp.first[b, ] <- out.b$weighted.hosp.first$WR
    perm.pointwise.threshold[b, ] <- out.b$threshold$WR

    perm.max[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.max[b, "traditionalWR_hosp_first"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
    perm.max[b, "traditionalOrderWR"] <- out.b$max.traditional.order$WR[1]
    perm.max[b, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
    perm.max[b, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
    perm.max[b, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
    perm.max[b, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
    perm.max[b, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
    perm.max[b, "maxWRt"] <- out.b$max.threshold$WR[1]

    perm.fixed[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.fixed[b, "fixedTraditionalHospFirst"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
    perm.fixed[b, "fixedTraditionalOrderWR"] <- wr_at_order_p(out.b, fixed.traditional.order, 0.50)
    perm.fixed[b, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
    perm.fixed[b, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
    perm.fixed[b, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
    perm.fixed[b, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
    perm.fixed[b, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
    perm.fixed[b, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

    selected.perm$selected.traditional.order[b] <- out.b$max.traditional.order$order[1]
    selected.perm$selected.p.primary[b] <- out.b$max.weighted.primary$p[1]
    selected.perm$selected.p.low[b] <- out.b$max.weighted.low$p[1]
    selected.perm$selected.p.full[b] <- out.b$max.weighted.full$p[1]
    selected.perm$selected.order.primary[b] <- out.b$max.order.primary$order[1]
    selected.perm$selected.order.p.primary[b] <- out.b$max.order.primary$p[1]
    selected.perm$selected.order.full[b] <- out.b$max.order.full$order[1]
    selected.perm$selected.order.p.full[b] <- out.b$max.order.full$p[1]
    selected.perm$selected.t[b] <- out.b$max.threshold$t[1]
  }

  calc_perm_p <- function(perm.vec, obs.value) (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))
  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k]))
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k]))
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k]))

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(death_first = perm.pointwise.death.first, hospitalization_first = perm.pointwise.hosp.first, threshold = perm.pointwise.threshold),
    p.value.pointwise = list(death_first = p.value.pointwise.death.first, hospitalization_first = p.value.pointwise.hosp.first, threshold = p.value.pointwise.threshold),
    selected.fixed = list(
      traditional.order = fixed.traditional.order,
      traditional.p = 0.50,
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    p.grid = p.grid,
    t.grid = t.grid
  )
}

method.long.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional WR, death first",
    traditionalWR_hosp_first = "Traditional WR, hospitalization first",
    traditionalOrderWR = "Traditional order-selected WR",
    maxWRp_primary = "Weighted max WR, p >= 0.5",
    maxWRp_low = "Weighted max WR, p < 0.5",
    maxWRp_full = "Weighted max WR, full p grid",
    maxOrderWR_primary = "Maximum-order weighted WR, p >= 0.5",
    maxOrderWR_full = "Maximum-order weighted WR, full p grid",
    maxWRt = "Threshold max WR(t)",
    logrank_death = "Log-rank death",
    logrank_composite = "Log-rank composite",
    fixedTraditionalHospFirst = "Traditional fixed WR, hospitalization first",
    fixedTraditionalOrderWR = "Traditional fixed-selected order WR",
    fixedWRp_primary = "Weighted fixed-selected WR, p >= 0.5",
    fixedWRp_low = "Weighted fixed-selected WR, p < 0.5",
    fixedWRp_full = "Weighted fixed-selected WR, full p grid",
    fixedOrderWR_primary = "Maximum-order fixed-selected WR, p >= 0.5",
    fixedOrderWR_full = "Maximum-order fixed-selected WR, full p grid",
    fixedWRt = "Threshold fixed-selected WR(t)"
  )
  out <- as.character(method); hit <- out %in% names(map); out[hit] <- unname(map[out[hit]]); out[!hit] <- gsub("_", " ", out[!hit]); out
}

method.display.label <- function(method) {
  map <- c(
    ordinaryWR = "Trad death\nWR",
    traditionalWR_hosp_first = "Trad hosp\nWR",
    traditionalOrderWR = "Trad order\nmax",
    maxWRp_primary = "Weighted\nmax",
    maxWRp_low = "Low-p\nmax",
    maxWRp_full = "Full-p\nmax",
    maxOrderWR_primary = "Weighted order\nmax",
    maxOrderWR_full = "Order full\nmax",
    maxWRt = "Threshold\nmax",
    logrank_death = "LR\ndeath",
    logrank_composite = "LR\ncomposite",
    fixedTraditionalHospFirst = "Trad hosp\nfixed",
    fixedTraditionalOrderWR = "Trad order\nfixed",
    fixedWRp_primary = "Weighted\nfixed",
    fixedWRp_low = "Low-p\nfixed",
    fixedWRp_full = "Full-p\nfixed",
    fixedOrderWR_primary = "Weighted order\nfixed",
    fixedOrderWR_full = "Order full\nfixed",
    fixedWRt = "Threshold\nfixed"
  )
  out <- as.character(method); hit <- out %in% names(map); out[hit] <- unname(map[out[hit]]); out[!hit] <- gsub("_", " ", out[!hit]); out
}

make.max.vs.fixed.table <- function(test.out) {
  method.ids <- c(
    "ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR",
    "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.ids <- c(
    "ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR",
    "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )
  data.frame(
    method.id = method.ids,
    method = method.long.label(method.ids),
    selected.parameter = c(
      "p = 0.50, death first",
      "p = 0.50, hospitalization first",
      paste0("selected order = ", clean.order.label(test.out$selected.fixed$traditional.order), "; p = 0.50"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), ", death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), ", death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), ", death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary), ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full), ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    max.statistic = as.numeric(test.out$T.obs.max[method.ids]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[method.ids]),
    fixed.parameter.statistic = as.numeric(test.out$T.obs.fixed[fixed.ids]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[fixed.ids]),
    stringsAsFactors = FALSE
  )
}

make.selected.max.table <- function(test.out) {
  obs <- test.out$observed
  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional <- get_selected_weighted_row(obs, test.out$selected.fixed$traditional.order, 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)
  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])
  method.ids <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  fixed.ids <- c("ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full", "fixedOrder_primary", "fixedOrder_full", "fixedWRt")
  data.frame(
    method = method.long.label(method.ids),
    method.id = method.ids,
    selected.order = c("death_first", "hospitalization_first", test.out$selected.fixed$traditional.order, "death_first", "death_first", "death_first", test.out$selected.fixed$order.primary, test.out$selected.fixed$order.full, NA_character_),
    selected.p = c(0.5, 0.5, 0.5, test.out$selected.fixed$p.primary, test.out$selected.fixed$p.low, test.out$selected.fixed$p.full, test.out$selected.fixed$order.p.primary, test.out$selected.fixed$order.p.full, NA_real_),
    selected.t.months = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, test.out$selected.fixed$t.months),
    observed.max.statistic = as.numeric(test.out$T.obs.max[method.ids]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[method.ids]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[fixed.ids]),
    selected.tie.count = c(get_tie_count(w.traditional), get_tie_count(h.traditional), get_tie_count(to.traditional), get_tie_count(w.primary), get_tie_count(w.low), get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full), get_tie_count(t.row)),
    selected.tie.pr = c(get_tie_pr(w.traditional), get_tie_pr(h.traditional), get_tie_pr(to.traditional), get_tie_pr(w.primary), get_tie_pr(w.low), get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full), get_tie_pr(t.row)),
    stringsAsFactors = FALSE
  )
}

make.pathway.results.table <- function(test.out, logrank.results = NULL) {
  obs <- test.out$observed
  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional <- get_selected_weighted_row(obs, test.out$selected.fixed$traditional.order, 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)
  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])
  method.id <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  fixed.id <- c("ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full", "fixedOrder_primary", "fixedOrder_full", "fixedWRt")
  out <- data.frame(
    method = method.id,
    method_label = method.long.label(method.id),
    statistic_type = "Win ratio",
    statistic_name = c(
      "WR(p = 0.50), death first",
      "WR(p = 0.50), hospitalization first",
      "max WR(order), p = 0.50",
      "max WR(p), death first, p in [0.50, 1.00]",
      "max WR(p), death first, p in [0.01, 0.49]",
      "max WR(p), death first, p in [0.01, 1.00]",
      "max WR(order, p), p in [0.50, 1.00]",
      "max WR(order, p), p in [0.01, 1.00]",
      "max WR(t)"
    ),
    observed_statistic = as.numeric(test.out$T.obs.max[method.id]),
    selected_order = c("death_first", "hospitalization_first", test.out$selected.fixed$traditional.order, "death_first", "death_first", "death_first", test.out$selected.fixed$order.primary, test.out$selected.fixed$order.full, NA_character_),
    selected_p = c(0.5, 0.5, 0.5, test.out$selected.fixed$p.primary, test.out$selected.fixed$p.low, test.out$selected.fixed$p.full, test.out$selected.fixed$order.p.primary, test.out$selected.fixed$order.p.full, NA_real_),
    selected_t_months = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, test.out$selected.fixed$t.months),
    selected_value = c(
      "p = 0.50; death first",
      "p = 0.50; hospitalization first",
      paste0("selected order = ", clean.order.label(test.out$selected.fixed$traditional.order), "; p = 0.50"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), "; death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary), "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full), "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    permutation_p_value = as.numeric(test.out$p.value.max[method.id]),
    selected_parameter_p_value = as.numeric(test.out$p.value.fixed[fixed.id]),
    tie_count = c(get_tie_count(w.traditional), get_tie_count(h.traditional), get_tie_count(to.traditional), get_tie_count(w.primary), get_tie_count(w.low), get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full), get_tie_count(t.row)),
    tie_proportion = c(get_tie_pr(w.traditional), get_tie_pr(h.traditional), get_tie_pr(to.traditional), get_tie_pr(w.primary), get_tie_pr(w.low), get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full), get_tie_pr(t.row)),
    B = test.out$B,
    stringsAsFactors = FALSE
  )
  if (!is.null(logrank.results)) {
    lr.death <- logrank.results[logrank.results$method == "Log-rank death endpoint", , drop = FALSE]
    lr.comp <- logrank.results[logrank.results$method == "Log-rank composite endpoint", , drop = FALSE]
    lr.rows <- data.frame(
      method = c("logrank_death", "logrank_composite"),
      method_label = method.long.label(c("logrank_death", "logrank_composite")),
      statistic_type = "Log-rank chi-square",
      statistic_name = c("Death endpoint", "Composite endpoint"),
      observed_statistic = c(lr.death$statistic[1], lr.comp$statistic[1]),
      selected_order = NA_character_, selected_p = NA_real_, selected_t_months = NA_real_, selected_value = NA_character_,
      permutation_p_value = c(lr.death$p.value[1], lr.comp$p.value[1]),
      selected_parameter_p_value = NA_real_, tie_count = NA_real_, tie_proportion = NA_real_, B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.rows)
  }
  out
}

plot.order.curves <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_maximum_order_WRp.png"), width = 1900, height = 1150)
  old.par <- par(no.readonly = TRUE); par(mar = c(5.2, 5.2, 3.4, 1.2))
  dtab <- engine.out$weighted.death.first; htab <- engine.out$weighted.hosp.first
  yr <- range(c(dtab$WR, htab$WR), finite = TRUE)
  plot(dtab$p, dtab$WR, type = "l", lwd = 2, ylim = yr,
       xlab = "p: first-endpoint weight", ylab = "WR(p)", main = "Hierarchy order comparison")
  lines(htab$p, htab$WR, lwd = 2, lty = 2)
  points(0.50, engine.out$traditional.death.first$WR[1], pch = 1, cex = 1.5)
  points(0.50, engine.out$traditional.hosp.first$WR[1], pch = 2, cex = 1.5)
  points(engine.out$max.traditional.order$p[1], engine.out$max.traditional.order$WR[1], pch = 17, cex = 1.3)
  points(engine.out$max.order.primary$p, engine.out$max.order.primary$WR, pch = 19)
  abline(v = 0.5, lty = 3, lwd = 2)
  legend("topright",
         legend = c("death first", "hospitalization first",
                    paste0("traditional selected = ", clean.order.label(engine.out$max.traditional.order$order[1])),
                    paste0("traditional order WR = ", sprintf("%.3f", engine.out$max.traditional.order$WR[1])),
                    paste0("weighted selected = ", clean.order.label(engine.out$max.order.primary$order[1])),
                    paste0("weighted p = ", sprintf("%.2f", engine.out$max.order.primary$p[1]))),
         lty = c(1, 2, NA, NA, NA, NA), pch = c(NA, NA, 17, NA, 19, NA), bty = "n", cex = 0.82)
  par(old.par); dev.off()
}

plot.max.vs.fixed.pvalues <- function(max.vs.fixed.table, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_max_vs_fixed_pvalues.png"), width = 2500, height = 1250)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  mat <- rbind("Max-statistic" = max.vs.fixed.table$max.statistic.p.value, "Fixed-selected" = max.vs.fixed.table$fixed.parameter.p.value)
  labels <- if ("method.id" %in% names(max.vs.fixed.table)) method.display.label(max.vs.fixed.table$method.id) else max.vs.fixed.table$method
  bp <- barplot(mat, beside = TRUE, names.arg = labels, ylim = c(0, 1), ylab = "Permutation p-value", main = "Max-statistic vs fixed-selected p-values", legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), las = 2, cex.names = 0.74)
  abline(h = ALPHA, lty = 2, lwd = 2)
  text(x = bp, y = pmin(mat + 0.035, 0.98), labels = sprintf("%.3f", mat), cex = 0.58)
  par(old.par); dev.off()
}

plot.example.pvalue.figures <- function(pathway.tab, outdir, prefix) {
  all.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite")
  wr.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  weighted.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full")
  logrank.methods <- c("logrank_death", "logrank_composite")
  threshold.methods <- c("maxWRt")
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_all_methods_with_logrank.png"), "All method p-values", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_WR_pathways_only_no_logrank.png"), "WR pathway p-values", "permutation_p_value", wr.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_traditional_and_weighted_only.png"), "Traditional/order and weighted WR p-values", "permutation_p_value", weighted.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_logrank_pvalues_separate.png"), "Log-rank p-values", "permutation_p_value", logrank.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_threshold_pvalue_separate.png"), "Threshold pathway p-value", "permutation_p_value", threshold.methods)
}

plot.example.permutation.nulls <- function(test.out, outdir, prefix) {
  for (m in c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")) {
    plot.permutation.null.revised(test.out, m, outdir, prefix)
  }
}

extract.trial.row <- function(scenario.row, sim.index, test.out, logrank.results, comp.stats,
                              threshold.info, error.message = NA_character_) {
  obs <- test.out$observed
  death.first.primary <- obs$max.weighted.primary
  hosp.first.primary <- max_from_curve(obs$weighted.hosp.first, "hosp_first_primary", 0.50, 1.00)
  order.diff <- hosp.first.primary$WR[1] - death.first.primary$WR[1]

  w.traditional.row <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional.row <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional.row <- get_selected_weighted_row(obs, obs$max.traditional.order$order[1], 0.5)
  w.primary.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.primary$p[1])
  w.low.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.low$p[1])
  w.full.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.full$p[1])
  o.primary.row <- get_selected_weighted_row(obs, obs$max.order.primary$order[1], obs$max.order.primary$p[1])
  o.full.row <- get_selected_weighted_row(obs, obs$max.order.full$order[1], obs$max.order.full$p[1])

  data.frame(
    scenario_id = scenario.row$scenario_id,
    scenario_type = get_scenario_field(scenario.row, "scenario_type"),
    description = get_scenario_field(scenario.row, "description"),
    sim = sim.index,
    error.message = error.message,

    ordinaryWR = test.out$T.obs.max["ordinaryWR"],
    traditional.tie.count = w.traditional.row$tie.count[1],
    traditional.tie.pr = w.traditional.row$tie.pr[1],
    max.pvalue.ordinaryWR = test.out$p.value.max["ordinaryWR"],
    fixed.pvalue.ordinaryWR = test.out$p.value.fixed["ordinaryWR"],

    traditionalWR_hosp_first = test.out$T.obs.max["traditionalWR_hosp_first"],
    traditional.hosp.first.tie.count = h.traditional.row$tie.count[1],
    traditional.hosp.first.tie.pr = h.traditional.row$tie.pr[1],
    max.pvalue.traditionalWR_hosp_first = test.out$p.value.max["traditionalWR_hosp_first"],
    fixed.pvalue.traditionalWR_hosp_first = test.out$p.value.fixed["fixedTraditionalHospFirst"],

    traditional.death.first.WR = w.traditional.row$WR[1],
    traditional.hosp.first.WR = h.traditional.row$WR[1],
    traditional.order.diff.hosp.minus.death = h.traditional.row$WR[1] - w.traditional.row$WR[1],
    traditional.order.changed = as.integer(obs$max.traditional.order$order[1] == "hospitalization_first"),
    traditionalOrderWR = test.out$T.obs.max["traditionalOrderWR"],
    selected.traditional.order = obs$max.traditional.order$order[1],
    traditional.order.tie.count = to.traditional.row$tie.count[1],
    traditional.order.tie.pr = to.traditional.row$tie.pr[1],
    max.pvalue.traditionalOrderWR = test.out$p.value.max["traditionalOrderWR"],
    fixed.pvalue.traditionalOrderWR = test.out$p.value.fixed["fixedTraditionalOrderWR"],

    maxWRp_primary = test.out$T.obs.max["maxWRp_primary"],
    selected.p.primary = obs$max.weighted.primary$p[1],
    weighted.primary.tie.count = w.primary.row$tie.count[1],
    weighted.primary.tie.pr = w.primary.row$tie.pr[1],
    max.pvalue.maxWRp_primary = test.out$p.value.max["maxWRp_primary"],
    fixed.pvalue.maxWRp_primary = test.out$p.value.fixed["fixedWRp_primary"],

    maxWRp_low = test.out$T.obs.max["maxWRp_low"],
    selected.p.low = obs$max.weighted.low$p[1],
    weighted.low.tie.count = w.low.row$tie.count[1],
    weighted.low.tie.pr = w.low.row$tie.pr[1],
    max.pvalue.maxWRp_low = test.out$p.value.max["maxWRp_low"],
    fixed.pvalue.maxWRp_low = test.out$p.value.fixed["fixedWRp_low"],

    maxWRp_full = test.out$T.obs.max["maxWRp_full"],
    selected.p.full = obs$max.weighted.full$p[1],
    weighted.full.tie.count = w.full.row$tie.count[1],
    weighted.full.tie.pr = w.full.row$tie.pr[1],
    selected.full.is.low.p = as.integer(obs$max.weighted.full$p[1] < 0.5),
    max.pvalue.maxWRp_full = test.out$p.value.max["maxWRp_full"],
    fixed.pvalue.maxWRp_full = test.out$p.value.fixed["fixedWRp_full"],

    death.first.primary.WR = death.first.primary$WR[1],
    hosp.first.primary.WR = hosp.first.primary$WR[1],
    order.diff.hosp.minus.death = order.diff,
    order.changed.primary = as.integer(obs$max.order.primary$order[1] == "hospitalization_first"),
    maxOrderWR_primary = test.out$T.obs.max["maxOrderWR_primary"],
    selected.order.primary = obs$max.order.primary$order[1],
    selected.order.p.primary = obs$max.order.primary$p[1],
    order.primary.tie.count = o.primary.row$tie.count[1],
    order.primary.tie.pr = o.primary.row$tie.pr[1],
    max.pvalue.maxOrderWR_primary = test.out$p.value.max["maxOrderWR_primary"],
    fixed.pvalue.maxOrderWR_primary = test.out$p.value.fixed["fixedOrder_primary"],

    maxOrderWR_full = test.out$T.obs.max["maxOrderWR_full"],
    selected.order.full = obs$max.order.full$order[1],
    selected.order.p.full = obs$max.order.full$p[1],
    order.full.tie.count = o.full.row$tie.count[1],
    order.full.tie.pr = o.full.row$tie.pr[1],
    selected.order.full.is.low.p = as.integer(obs$max.order.full$p[1] < 0.5),
    max.pvalue.maxOrderWR_full = test.out$p.value.max["maxOrderWR_full"],
    fixed.pvalue.maxOrderWR_full = test.out$p.value.fixed["fixedOrder_full"],

    maxWRt = test.out$T.obs.max["maxWRt"],
    selected.t.years = obs$max.threshold$t[1],
    selected.t.months = obs$max.threshold$t.months[1],
    threshold.tie.count = obs$max.threshold$tie.count[1],
    threshold.tie.pr = obs$max.threshold$pr.tie[1],
    max.pvalue.maxWRt = test.out$p.value.max["maxWRt"],
    fixed.pvalue.maxWRt = test.out$p.value.fixed["fixedWRt"],

    true.hierarchical.tie.count = obs$counts$true.tie.pairs,
    true.hierarchical.tie.pr = obs$counts$true.tie.pr,
    total.pairs = obs$counts$total.pairs,
    logrank.death.statistic = logrank.results$statistic[logrank.results$method == "Log-rank death endpoint"],
    logrank.death.p = logrank.results$p.value[logrank.results$method == "Log-rank death endpoint"],
    logrank.composite.statistic = logrank.results$statistic[logrank.results$method == "Log-rank composite endpoint"],
    logrank.composite.p = logrank.results$p.value[logrank.results$method == "Log-rank composite endpoint"],
    death.event.rate.control = get_comp_arm_value(comp.stats, 0, "death.event.rate"),
    death.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "death.event.rate"),
    composite.event.rate.control = get_comp_arm_value(comp.stats, 0, "composite.event.rate"),
    composite.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "composite.event.rate"),
    mean.num.hosp.control = get_comp_arm_value(comp.stats, 0, "mean.num.hosp"),
    mean.num.hosp.treatment = get_comp_arm_value(comp.stats, 1, "mean.num.hosp"),
    threshold.grid.months = paste(sprintf("%.1f", threshold.info$threshold.table$t.months), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

summarise.scenario.results <- function(res, alpha = ALPHA) {
  ok <- res[is.na(res$error.message) | res$error.message == "", , drop = FALSE]
  if (nrow(ok) == 0) return(list(power = data.frame(), averages = data.frame()))
  method <- c(
    "ordinaryWR", "traditionalWR_hosp_first", "fixedTraditionalHospFirst", "traditionalOrderWR", "fixedTraditionalOrderWR",
    "maxWRp_primary", "fixedWRp_primary", "maxWRp_low", "fixedWRp_low", "maxWRp_full", "fixedWRp_full",
    "maxOrderWR_primary", "fixedOrderWR_primary", "maxOrderWR_full", "fixedOrderWR_full", "maxWRt", "fixedWRt", "logrank_death", "logrank_composite"
  )
  pcols <- c(
    "max.pvalue.ordinaryWR", "max.pvalue.traditionalWR_hosp_first", "fixed.pvalue.traditionalWR_hosp_first", "max.pvalue.traditionalOrderWR", "fixed.pvalue.traditionalOrderWR",
    "max.pvalue.maxWRp_primary", "fixed.pvalue.maxWRp_primary", "max.pvalue.maxWRp_low", "fixed.pvalue.maxWRp_low", "max.pvalue.maxWRp_full", "fixed.pvalue.maxWRp_full",
    "max.pvalue.maxOrderWR_primary", "fixed.pvalue.maxOrderWR_primary", "max.pvalue.maxOrderWR_full", "fixed.pvalue.maxOrderWR_full", "max.pvalue.maxWRt", "fixed.pvalue.maxWRt", "logrank.death.p", "logrank.composite.p"
  )
  power <- data.frame(method = method, method_label = method.long.label(method), rejection.proportion = sapply(pcols, function(cc) mean(safe.col(ok, cc) < alpha, na.rm = TRUE)), alpha = alpha, nsim.available = nrow(ok), stringsAsFactors = FALSE)
  averages <- data.frame(
    nsim.available = nrow(ok),
    mean.ordinaryWR = safe.mean(ok$ordinaryWR), mean.max.pvalue.ordinaryWR = safe.mean(ok$max.pvalue.ordinaryWR), mean.fixed.pvalue.ordinaryWR = safe.mean(ok$fixed.pvalue.ordinaryWR), mean.traditional.tie.count = safe.mean(ok$traditional.tie.count), mean.traditional.tie.pr = safe.mean(ok$traditional.tie.pr),
    mean.traditionalWR.hosp.first = safe.mean(ok$traditionalWR_hosp_first), mean.max.pvalue.traditionalWR.hosp.first = safe.mean(ok$max.pvalue.traditionalWR_hosp_first), mean.fixed.pvalue.traditionalWR.hosp.first = safe.mean(ok$fixed.pvalue.traditionalWR_hosp_first), mean.traditional.hosp.first.tie.count = safe.mean(ok$traditional.hosp.first.tie.count), mean.traditional.hosp.first.tie.pr = safe.mean(ok$traditional.hosp.first.tie.pr),
    mean.traditional.death.first.WR = safe.mean(ok$traditional.death.first.WR), mean.traditional.hosp.first.WR = safe.mean(ok$traditional.hosp.first.WR), mean.traditional.order.diff.hosp.minus.death = safe.mean(ok$traditional.order.diff.hosp.minus.death), prop.traditional.order.changed = mean(ok$traditional.order.changed == 1, na.rm = TRUE), mode.selected.traditional.order = clean.order.label(mode.string(ok$selected.traditional.order)), mean.traditionalOrderWR = safe.mean(ok$traditionalOrderWR), mean.max.pvalue.traditionalOrderWR = safe.mean(ok$max.pvalue.traditionalOrderWR), mean.fixed.pvalue.traditionalOrderWR = safe.mean(ok$fixed.pvalue.traditionalOrderWR), mean.traditional.order.tie.count = safe.mean(ok$traditional.order.tie.count), mean.traditional.order.tie.pr = safe.mean(ok$traditional.order.tie.pr),
    mean.maxWRp.primary = safe.mean(ok$maxWRp_primary), mean.selected.p.primary = safe.mean(ok$selected.p.primary), mean.max.pvalue.maxWRp.primary = safe.mean(ok$max.pvalue.maxWRp_primary), mean.fixed.pvalue.maxWRp.primary = safe.mean(ok$fixed.pvalue.maxWRp_primary), mean.weighted.primary.tie.count = safe.mean(ok$weighted.primary.tie.count), mean.weighted.primary.tie.pr = safe.mean(ok$weighted.primary.tie.pr),
    mean.maxWRp.low = safe.mean(ok$maxWRp_low), mean.selected.p.low = safe.mean(ok$selected.p.low), mean.max.pvalue.maxWRp.low = safe.mean(ok$max.pvalue.maxWRp_low), mean.fixed.pvalue.maxWRp.low = safe.mean(ok$fixed.pvalue.maxWRp_low), mean.weighted.low.tie.count = safe.mean(ok$weighted.low.tie.count), mean.weighted.low.tie.pr = safe.mean(ok$weighted.low.tie.pr),
    mean.maxWRp.full = safe.mean(ok$maxWRp_full), mean.selected.p.full = safe.mean(ok$selected.p.full), prop.full.selected.p.below.0.5 = mean(ok$selected.full.is.low.p == 1, na.rm = TRUE), mean.max.pvalue.maxWRp.full = safe.mean(ok$max.pvalue.maxWRp_full), mean.fixed.pvalue.maxWRp.full = safe.mean(ok$fixed.pvalue.maxWRp_full), mean.weighted.full.tie.count = safe.mean(ok$weighted.full.tie.count), mean.weighted.full.tie.pr = safe.mean(ok$weighted.full.tie.pr),
    mean.maxOrderWR.primary = safe.mean(ok$maxOrderWR_primary), prop.order.changed.primary = mean(ok$order.changed.primary == 1, na.rm = TRUE), mean.order.diff.hosp.minus.death = safe.mean(ok$order.diff.hosp.minus.death), mean.selected.order.p.primary = safe.mean(ok$selected.order.p.primary), mode.selected.order.primary = clean.order.label(mode.string(ok$selected.order.primary)), mean.max.pvalue.maxOrderWR.primary = safe.mean(ok$max.pvalue.maxOrderWR_primary), mean.fixed.pvalue.maxOrderWR.primary = safe.mean(ok$fixed.pvalue.maxOrderWR_primary), mean.order.primary.tie.count = safe.mean(ok$order.primary.tie.count), mean.order.primary.tie.pr = safe.mean(ok$order.primary.tie.pr),
    mean.maxOrderWR.full = safe.mean(ok$maxOrderWR_full), mean.selected.order.p.full = safe.mean(ok$selected.order.p.full), mode.selected.order.full = clean.order.label(mode.string(ok$selected.order.full)), prop.order.full.selected.p.below.0.5 = mean(ok$selected.order.full.is.low.p == 1, na.rm = TRUE), mean.max.pvalue.maxOrderWR.full = safe.mean(ok$max.pvalue.maxOrderWR_full), mean.fixed.pvalue.maxOrderWR.full = safe.mean(ok$fixed.pvalue.maxOrderWR_full), mean.order.full.tie.count = safe.mean(ok$order.full.tie.count), mean.order.full.tie.pr = safe.mean(ok$order.full.tie.pr),
    mean.maxWRt = safe.mean(ok$maxWRt), mean.selected.t.months = safe.mean(ok$selected.t.months), sd.selected.t.months = safe_sd(ok$selected.t.months), mean.max.pvalue.maxWRt = safe.mean(ok$max.pvalue.maxWRt), mean.fixed.pvalue.maxWRt = safe.mean(ok$fixed.pvalue.maxWRt), mean.threshold.tie.count = safe.mean(ok$threshold.tie.count), mean.threshold.tie.pr = safe.mean(ok$threshold.tie.pr),
    mean.true.hierarchical.tie.count = safe.mean(ok$true.hierarchical.tie.count), mean.true.hierarchical.tie.pr = safe.mean(ok$true.hierarchical.tie.pr), mean.logrank.death.statistic = safe.mean(ok$logrank.death.statistic), mean.logrank.death.p = safe.mean(ok$logrank.death.p), mean.logrank.composite.statistic = safe.mean(ok$logrank.composite.statistic), mean.logrank.composite.p = safe.mean(ok$logrank.composite.p),
    mean.death.event.rate.control = safe.mean(ok$death.event.rate.control), mean.death.event.rate.treatment = safe.mean(ok$death.event.rate.treatment), mean.composite.event.rate.control = safe.mean(ok$composite.event.rate.control), mean.composite.event.rate.treatment = safe.mean(ok$composite.event.rate.treatment), mean.num.hosp.control = safe.mean(ok$mean.num.hosp.control), mean.num.hosp.treatment = safe.mean(ok$mean.num.hosp.treatment), stringsAsFactors = FALSE)
  list(power = power, averages = averages)
}

make.scenario.method.comparison <- function(res, power.table, alpha = ALPHA) {
  ok <- res[is.na(res$error.message) | res$error.message == "", , drop = FALSE]
  if (nrow(ok) == 0) return(data.frame())
  specs <- data.frame(
    method = c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite"),
    pathway = c("Traditional WR", "Traditional WR", "Traditional order WR", "Weighted WR", "Weighted WR", "Weighted WR", "Maximum-order weighted WR", "Maximum-order weighted WR", "Threshold WR", "Log-rank", "Log-rank"),
    statistic.column = c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank.death.statistic", "logrank.composite.statistic"),
    permutation.p.column = c("max.pvalue.ordinaryWR", "max.pvalue.traditionalWR_hosp_first", "max.pvalue.traditionalOrderWR", "max.pvalue.maxWRp_primary", "max.pvalue.maxWRp_low", "max.pvalue.maxWRp_full", "max.pvalue.maxOrderWR_primary", "max.pvalue.maxOrderWR_full", "max.pvalue.maxWRt", "logrank.death.p", "logrank.composite.p"),
    selected.pvalue.column = c("fixed.pvalue.ordinaryWR", "fixed.pvalue.traditionalWR_hosp_first", "fixed.pvalue.traditionalOrderWR", "fixed.pvalue.maxWRp_primary", "fixed.pvalue.maxWRp_low", "fixed.pvalue.maxWRp_full", "fixed.pvalue.maxOrderWR_primary", "fixed.pvalue.maxOrderWR_full", "fixed.pvalue.maxWRt", NA_character_, NA_character_),
    selected.p.column = c(NA_character_, NA_character_, NA_character_, "selected.p.primary", "selected.p.low", "selected.p.full", "selected.order.p.primary", "selected.order.p.full", NA_character_, NA_character_, NA_character_),
    selected.t.column = c(NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "selected.t.months", NA_character_, NA_character_),
    selected.order.column = c(NA_character_, NA_character_, "selected.traditional.order", NA_character_, NA_character_, NA_character_, "selected.order.primary", "selected.order.full", NA_character_, NA_character_, NA_character_),
    tie.count.column = c("traditional.tie.count", "traditional.hosp.first.tie.count", "traditional.order.tie.count", "weighted.primary.tie.count", "weighted.low.tie.count", "weighted.full.tie.count", "order.primary.tie.count", "order.full.tie.count", "threshold.tie.count", NA_character_, NA_character_),
    tie.pr.column = c("traditional.tie.pr", "traditional.hosp.first.tie.pr", "traditional.order.tie.pr", "weighted.primary.tie.pr", "weighted.low.tie.pr", "weighted.full.tie.pr", "order.primary.tie.pr", "order.full.tie.pr", "threshold.tie.pr", NA_character_, NA_character_),
    fixed.power.method = c("ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full", "fixedOrderWR_primary", "fixedOrderWR_full", "fixedWRt", NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )
  rows <- vector("list", nrow(specs))
  for (i in seq_len(nrow(specs))) {
    sp <- specs[i, , drop = FALSE]
    stat <- safe.col(ok, sp$statistic.column); pp <- safe.col(ok, sp$permutation.p.column)
    fixed.pp <- if (!is.na(sp$selected.pvalue.column)) safe.col(ok, sp$selected.pvalue.column) else rep(NA_real_, nrow(ok))
    sel.p <- if (!is.na(sp$selected.p.column)) safe.col(ok, sp$selected.p.column) else rep(NA_real_, nrow(ok))
    sel.t <- if (!is.na(sp$selected.t.column)) safe.col(ok, sp$selected.t.column) else rep(NA_real_, nrow(ok))
    tie.count <- if (!is.na(sp$tie.count.column)) safe.col(ok, sp$tie.count.column) else rep(NA_real_, nrow(ok))
    tie.pr <- if (!is.na(sp$tie.pr.column)) safe.col(ok, sp$tie.pr.column) else rep(NA_real_, nrow(ok))
    sel.order <- if (!is.na(sp$selected.order.column) && sp$selected.order.column %in% names(ok)) as.character(ok[[sp$selected.order.column]]) else rep(NA_character_, nrow(ok))
    selected.summary <- NA_character_
    if (sp$method == "ordinaryWR") selected.summary <- "p = 0.50; death first"
    if (sp$method == "traditionalWR_hosp_first") selected.summary <- "p = 0.50; hospitalization first"
    if (sp$method == "traditionalOrderWR") selected.summary <- paste0("mode order = ", clean.order.label(mode.string(sel.order)), "; p = 0.50")
    if (sp$method %in% c("maxWRp_primary", "maxWRp_low", "maxWRp_full")) selected.summary <- paste0("mean p = ", sprintf("%.2f", safe.mean(sel.p)), "; death first")
    if (sp$method %in% c("maxOrderWR_primary", "maxOrderWR_full")) selected.summary <- paste0("mode order = ", clean.order.label(mode.string(sel.order)), "; mean p = ", sprintf("%.2f", safe.mean(sel.p)))
    if (sp$method == "maxWRt") selected.summary <- paste0("mean t = ", sprintf("%.1f", safe.mean(sel.t)), " months")
    rows[[i]] <- data.frame(scenario_id = ok$scenario_id[1], scenario_label = scenario.display.label(ok$scenario_id[1]), scenario_type = ok$scenario_type[1], description = ok$description[1], method = sp$method, method_label = method.long.label(sp$method), pathway = sp$pathway, statistic_type = ifelse(grepl("logrank", sp$method), "Log-rank chi-square", "Win ratio"), nsim.available = nrow(ok), mean_statistic = safe.mean(stat), sd_statistic = safe_sd(stat), median_statistic = safe.median(stat), selected_value_summary = selected.summary, mean_selected_p = safe.mean(sel.p), median_selected_p = safe.median(sel.p), mean_selected_t_months = safe.mean(sel.t), median_selected_t_months = safe.median(sel.t), selected_order_mode = clean.order.label(mode.string(sel.order)), prop_hospitalization_first = ifelse(all(is.na(sel.order)), NA_real_, mean(sel.order == "hospitalization_first", na.rm = TRUE)), mean_permutation_p_value = safe.mean(pp), median_permutation_p_value = safe.median(pp), power_permutation_p_value = mean(pp < alpha, na.rm = TRUE), mean_selected_parameter_p_value = safe.mean(fixed.pp), median_selected_parameter_p_value = safe.median(fixed.pp), power_selected_parameter_p_value = mean(fixed.pp < alpha, na.rm = TRUE), power_selected_parameter_from_power_table = get.power.value(power.table, sp$fixed.power.method), mean_tie_count = safe.mean(tie.count), mean_tie_proportion = safe.mean(tie.pr), alpha = alpha, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL; out
}

plot.global.power.clean <- function(power.table, outdir) {
  if (nrow(power.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  tab <- power.table[power.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  scenarios <- unique(tab$scenario_id)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios), dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  for (i in seq_len(nrow(tab))) mat[method.display.label(tab$method[i]), scenario.display.label(tab$scenario_id[i])] <- tab$rejection.proportion[i]
  save_png(outdir, "GLOBAL_power_comparison_key_methods.png", width = 3400, height = 1600)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2, ylab = paste0("Rejection proportion at alpha = ", ALPHA), main = "Power / rejection proportion across scenarios", legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.65), cex.names = 0.82)
  abline(h = ALPHA, lty = 2, lwd = 2); par(old.par); dev.off()
}

plot.global.mean.pvalues.clean <- function(method.table, outdir) {
  if (nrow(method.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  tab <- method.table[method.table$method %in% key.methods, , drop = FALSE]
  scenarios <- unique(tab$scenario_id)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios), dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  for (i in seq_len(nrow(tab))) mat[method.display.label(tab$method[i]), scenario.display.label(tab$scenario_id[i])] <- tab$mean_permutation_p_value[i]
  save_png(outdir, "GLOBAL_permutation_pvalue_comparison.png", width = 3400, height = 1600)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2, ylab = "Mean permutation p-value", main = "Mean permutation p-value across scenarios", legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.65), cex.names = 0.82)
  abline(h = ALPHA, lty = 2, lwd = 2); par(old.par); dev.off()
}

plot.global.fixed.pvalues.clean <- function(method.table, outdir) {
  if (nrow(method.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")
  tab <- method.table[method.table$method %in% key.methods, , drop = FALSE]
  scenarios <- unique(tab$scenario_id)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(scenarios), dimnames = list(method.display.label(key.methods), scenario.display.label(scenarios)))
  for (i in seq_len(nrow(tab))) mat[method.display.label(tab$method[i]), scenario.display.label(tab$scenario_id[i])] <- tab$mean_selected_parameter_p_value[i]
  save_png(outdir, "GLOBAL_fixed_selected_pvalue_comparison.png", width = 3100, height = 1500)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2, ylab = "Mean fixed-selected p-value", main = "Mean fixed-selected p-value across scenarios", legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.72), cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2); par(old.par); dev.off()
}

# Use the updated global comparison functions in the original runner.
plot.global.permutation.pvalues.clean <- plot.global.mean.pvalues.clean
plot.global.max.fixed.pvalues.clean <- plot.global.fixed.pvalues.clean

#Run block
out <- run.all.scenarios(
  nsim = NSIM,
  B = B_PERM,
  scenario.grid = build.scenario.grid()
)
