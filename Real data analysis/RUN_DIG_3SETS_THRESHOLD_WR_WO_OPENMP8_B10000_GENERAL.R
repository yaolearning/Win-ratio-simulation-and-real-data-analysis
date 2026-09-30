# ===========================================================================
# DIG THREE-ENDPOINT SETS: threshold-maximized WR/WO with permutation inference
#
# Current analysis settings:
#   - Three endpoint sets, all including death
#   - Death is the ONLY thresholded endpoint
#   - Death thresholds: 0, 7, 14, 21 days
#   - WR and WO
#   - Optimization over ordering, weighting, thresholds, and joint combinations
#   - Adaptive statistics are re-maximized independently in EVERY permutation
#   - Default B = 10,000 permutations
#   - Default OpenMP threads = 8
#
# This GitHub/general version uses RELATIVE paths and does not assume a local
# Windows drive. Put this R file and dig3_3sets_threshold_openmp_engine.cpp
# in the same directory, or update CPP_FILE below.
# ===========================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# ------------------------------- SETTINGS ----------------------------------
THRESHOLDS_DAYS <- c(0, 7, 14, 21)
B_PERM <- 10000L
N_THREADS <- 8L
MASTER_SEED <- 20260923L
CHECKPOINT_EVERY <- 25L

# Set TRUE for a quick debugging run.
PILOT_ONLY <- FALSE

# Output folder is created relative to the current working directory.
OUTPUT_ROOT <- file.path(
  getwd(),
  "DIG_3endpoint_THREE_SETS_threshold_max"
)

SETS_TO_RUN <- c(
  "Death_FirstHosp_RecurrentHosp",
  "Death_FirstHosp_SVA",
  "Death_RecurrentHosp_SVA"
)

# No FirstHosp_RecurrentHosp_SVA set in this analysis.

# C++ engine expected in the same working directory as this R script.
CPP_FILE <- file.path(getwd(), "dig3_3sets_threshold_openmp_engine.cpp")
# ---------------------------------------------------------------------------

if (PILOT_ONLY) {
  B_PERM <- 20L
  OUTPUT_ROOT <- paste0(OUTPUT_ROOT, "_PILOT")
} else {
  OUTPUT_ROOT <- paste0(OUTPUT_ROOT, "_B", B_PERM)
}

for (pkg in c("asympTest", "Rcpp", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Install package: ", pkg)
  }
}

library(Rcpp)
library(data.table)

if (!file.exists(CPP_FILE)) {
  stop(
    "C++ file not found: ", CPP_FILE,
    "\nPlace dig3_3sets_threshold_openmp_engine.cpp in the working directory ",
    "or update CPP_FILE."
  )
}

dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OMP_THREAD_LIMIT = as.character(N_THREADS)
)

Rcpp::sourceCpp(CPP_FILE, rebuild = TRUE, verbose = FALSE)

cat("OpenMP threads available:", dig_openmp_threads_cpp(), "\n")
if (dig_openmp_threads_cpp() == 1L && N_THREADS > 1L) {
  warning("OpenMP is unavailable in this R compiler setup; runtime may be long.")
}

data("DIGdata", package = "asympTest")
dig <- DIGdata

required_vars <- c(
  "TRTMT", "DEATH", "DEATHDAY",
  "HOSP", "HOSPDAYS", "NHOSP", "SVA"
)
stopifnot(all(required_vars %in% names(dig)))

arm <- suppressWarnings(as.integer(as.character(dig$TRTMT)))

if (!identical(sort(unique(arm)), c(0L, 1L))) {
  stop("Check DIG TRTMT coding; expected 0/1.")
}

cat("Treatment coding: 0=Placebo, 1=Digoxin; counts:\n")
print(table(arm))

if (anyNA(dig[, required_vars])) {
  stop("Missing endpoint/arm values: review missing-data handling before running.")
}

if (!all(dig$DEATH %in% c(0, 1)) || !all(dig$HOSP %in% c(0, 1))) {
  stop("Event indicators must be 0/1.")
}

if (!all(dig$SVA %in% c(0, 1))) {
  stop("Check SVA binary coding.")
}

if (any(dig$NHOSP < 0 | dig$NHOSP != floor(dig$NHOSP))) {
  stop("NHOSP must be a nonnegative integer hospitalization count.")
}

cat("NHOSP recurrent-hospitalization count summary (verify definition):\n")
print(summary(dig$NHOSP))

safe_ratio <- function(a, b) {
  if (a == 0 && b == 0) return(1)
  if (b == 0) return(Inf)
  a / b
}

safe_abslog <- function(x) {
  if (x == 0 || is.infinite(x)) return(Inf)
  abs(log(x))
}

# Two-sided comparison uses |log(WR)| or |log(WO)|.
comparison_score <- function(v, side) {
  if (side == "one") v else vapply(v, safe_abslog, numeric(1))
}

run_mode <- function(mode) {
  stopifnot(mode %in% SETS_TO_RUN)

  m <- 3L
  odir <- file.path(OUTPUT_ROOT, mode)
  dir.create(odir, recursive = TRUE, showWarnings = FALSE)

  ckpt_file <- file.path(odir, "CHECKPOINT.rds")
  complete_file <- file.path(odir, "COMPLETE.rds")
  n <- nrow(dig)

  tm <- matrix(NA_real_, n, m)
  em <- matrix(0L, n, m)
  vm <- matrix(NA_real_, n, m)

  # Endpoint 1 is death in all three sets and is the ONLY thresholded endpoint.
  tm[, 1] <- as.numeric(dig$DEATHDAY)
  em[, 1] <- as.integer(dig$DEATH)

  if (mode == "Death_FirstHosp_RecurrentHosp") {
    tm[, 2] <- as.numeric(dig$HOSPDAYS)
    em[, 2] <- as.integer(dig$HOSP)
    vm[, 3] <- as.numeric(dig$NHOSP)

    endpoint_names <- c(
      "Death",
      "First all-cause hospitalization",
      "Recurrent hospitalization (NHOSP count)"
    )
    type_code <- c(1L, 1L, 2L)

  } else if (mode == "Death_FirstHosp_SVA") {
    tm[, 2] <- as.numeric(dig$HOSPDAYS)
    em[, 2] <- as.integer(dig$HOSP)
    vm[, 3] <- as.numeric(dig$SVA)

    endpoint_names <- c(
      "Death",
      "First all-cause hospitalization",
      "SVA hospitalization"
    )
    type_code <- c(1L, 1L, 2L)

  } else if (mode == "Death_RecurrentHosp_SVA") {
    vm[, 2] <- as.numeric(dig$NHOSP)
    vm[, 3] <- as.numeric(dig$SVA)

    endpoint_names <- c(
      "Death",
      "Recurrent hospitalization (NHOSP count)",
      "SVA hospitalization"
    )
    type_code <- c(1L, 2L, 2L)

  } else {
    stop("Unrecognized endpoint set: ", mode)
  }

  cat(
    "\nEndpoint set: ", mode, "\n",
    paste(seq_along(endpoint_names), endpoint_names, sep = ": ", collapse = "\n"),
    "\n",
    sep = ""
  )

  # Six endpoint orders; natural order listed first.
  orders <- rbind(
    c(1L, 2L, 3L),
    c(1L, 3L, 2L),
    c(2L, 1L, 3L),
    c(2L, 3L, 1L),
    c(3L, 1L, 2L),
    c(3L, 2L, 1L)
  )
  no <- nrow(orders)

  # Prespecified weight vertices.
  weights <- rbind(
    c(1, 0, 0),
    c(0.5, 0.5, 0),
    c(1/3, 1/3, 1/3)
  )
  nw <- nrow(weights)

  # Equal weighting row = "original" weighting.
  ew <- nw

  # Threshold vectors follow endpoint identity, regardless of endpoint order.
  # Only death receives thresholds; other endpoints remain at 0.
  tg <- data.table(
    t_death = THRESHOLDS_DAYS,
    t_hosp = 0
  )
  setorder(tg, t_death, t_hosp)

  th <- cbind(
    as.matrix(tg),
    t_sva = rep(0, nrow(tg))
  )
  storage.mode(th) <- "double"

  nt <- nrow(th)
  base_t <- which(tg$t_death == 0 & tg$t_hosp == 0)

  if (length(base_t) != 1L) {
    stop("Exactly one zero/zero threshold vector required.")
  }

  method_names <- c(
    "Original",
    "M_ord",
    "M_wt",
    "M_ord_wt",
    "M_thresh",
    "M_ord_thresh",
    "M_wt_thresh",
    "M_ord_wt_thresh"
  )

  # Candidate rows in the same order as C++ output: threshold, order, weight.
  meta <- CJ(
    tid = seq_len(nt),
    oid = seq_len(no),
    wid = seq_len(nw),
    sorted = FALSE
  )
  setorder(meta, tid, oid, wid)

  meta[, candidate_id := .I]

  meta[, selected_order := vapply(
    oid,
    function(o) paste(orders[o, ], collapse = "->"),
    character(1)
  )]

  meta[, selected_weight := vapply(
    wid,
    function(w) {
      paste0("(", paste(sprintf("%.3f", weights[w, ]), collapse = ","), ")")
    },
    character(1)
  )]

  meta[, t_death_days := tg$t_death[tid]]
  meta[, t_hosp_days := tg$t_hosp[tid]]

  for (j in seq_len(m)) {
    meta[[paste0("p", j)]] <- weights[meta$wid, j]
  }

  # Prespecified candidate space for each adaptive method.
  method_ids <- list(
    Original = meta[tid == base_t & oid == 1 & wid == ew, candidate_id],
    M_ord = meta[tid == base_t & wid == ew, candidate_id],
    M_wt = meta[tid == base_t & oid == 1, candidate_id],
    M_ord_wt = meta[tid == base_t, candidate_id],
    M_thresh = meta[oid == 1 & wid == ew, candidate_id],
    M_ord_thresh = meta[wid == ew, candidate_id],
    M_wt_thresh = meta[oid == 1, candidate_id],
    M_ord_wt_thresh = meta$candidate_id
  )

  stopifnot(
    all(th[, 2] == 0),
    all(th[, 3] == 0),
    identical(
      sort(unique(th[, 1])),
      sort(THRESHOLDS_DAYS)
    )
  )

  evaluate <- function(trt) {
    z <- dig_evaluate_candidates_cpp(
      as.integer(trt),
      as.integer(type_code),
      tm,
      em,
      vm,
      orders,
      th,
      as.integer(N_THREADS)
    )

    out <- copy(meta)
    np <- z$n_pairs

    W <- matrix(NA_real_, nrow(out), m)
    L <- W
    U <- W

    for (t in seq_len(nt)) {
      for (o in seq_len(no)) {
        r <- (t - 1L) * no + o
        ix <- which(out$tid == t & out$oid == o)

        for (j in seq_len(m)) {
          W[ix, j] <- z$wins[r, j]
          L[ix, j] <- z$losses[r, j]
          U[ix, j] <- z$unresolved_after[r, j]
        }
      }
    }

    out[, c("weighted_wins", "weighted_losses", "WR", "WO") :=
          list(0, 0, 0, 0)]

    for (j in seq_len(m)) {
      out$weighted_wins <-
        out$weighted_wins + out[[paste0("p", j)]] * W[, j]

      out$weighted_losses <-
        out$weighted_losses + out[[paste0("p", j)]] * L[, j]

      out[[paste0("tie_after_rank", j, "_prop")]] <- U[, j] / np
    }

    final_ties <- U[, m]

    out$WR <- mapply(
      safe_ratio,
      out$weighted_wins,
      out$weighted_losses
    )

    # Current WO convention retained from the latest DIG analysis:
    # final unresolved pairs contribute 0.5 to both treatment and control.
    out$WO <- mapply(
      safe_ratio,
      out$weighted_wins + 0.5 * final_ties,
      out$weighted_losses + 0.5 * final_ties
    )

    out
  }

  cat(
    "\n", mode, " endpoints: ", n, " subjects; ",
    nt, " threshold vectors; ",
    no, " orders; ",
    nw, " weights; ",
    nrow(meta), " candidates.\n",
    sep = ""
  )
  cat("Threshold units: DAYS. Threshold is applied ONLY to DEATH.\n")

  fwrite(tg, file.path(odir, "threshold_grid.csv"))
  fwrite(meta, file.path(odir, "candidate_metadata.csv"))

  obs <- evaluate(arm)

  fwrite(
    obs,
    file.path(odir, "observed_all_candidates.csv")
  )

  fwrite(
    obs[tid == base_t & wid == ew],
    file.path(odir, "observed_fixed_orders.csv")
  )

  # One-sided output is retained for internal checking.
  # Manuscript-facing 3-endpoint reporting should use TWO-SIDED results only.
  specs <- CJ(
    measure = c("WR", "WO"),
    method = method_names,
    side = c("one", "two"),
    sorted = FALSE
  )

  # Deterministic tie-breaking follows candidate order:
  # threshold, then endpoint order, then weight.
  obs_rows <- vector("list", nrow(specs))

  for (i in seq_len(nrow(specs))) {
    a <- specs[i]
    pool <- obs[method_ids[[a$method]]]

    scores <- comparison_score(
      pool[[a$measure]],
      a$side
    )

    if (anyNA(scores)) {
      stop("NA statistic in observed candidates.")
    }

    k <- which.max(scores)
    best <- pool[k]

    z <- data.table(
      measure = a$measure,
      method = a$method,
      side = a$side,
      n_candidates = nrow(pool),
      max_statistic = best[[a$measure]],
      comparison_statistic = scores[k],
      selected_candidate_id = best$candidate_id,
      selected_order = best$selected_order,
      selected_weight = best$selected_weight,
      t_death_days = best$t_death_days,
      t_hosp_days = best$t_hosp_days
    )

    for (j in seq_len(m)) {
      z[[paste0("p", j)]] <- best[[paste0("p", j)]]
      z[[paste0("tie_after_rank", j, "_prop")]] <-
        best[[paste0("tie_after_rank", j, "_prop")]]
    }

    obs_rows[[i]] <- z
  }

  observed <- rbindlist(obs_rows)

  fwrite(
    observed,
    file.path(odir, "observed_method_selections.csv")
  )

  # Store config so checkpoint/resume cannot silently mix analyses.
  config <- list(
    mode = mode,
    endpoint_names = endpoint_names,
    type_code = type_code,
    B = B_PERM,
    seed = MASTER_SEED,
    threads = N_THREADS,
    threshold_days = THRESHOLDS_DAYS,
    weights = weights,
    orders = orders,
    n = n,
    group_sizes = as.integer(table(arm)),
    m = m
  )

  if (file.exists(complete_file)) {
    completed <- readRDS(complete_file)

    if (!identical(completed$config, config) ||
        !identical(completed$observed, observed)) {
      stop(
        "COMPLETE.rds belongs to a different analysis/configuration. ",
        "Move it before running."
      )
    }

    cat("Already completed: ", complete_file, "\n", sep = "")
    return(invisible(completed$results))
  }

  if (file.exists(ckpt_file)) {
    ck <- readRDS(ckpt_file)

    if (!identical(ck$config, config)) {
      stop(
        "Existing checkpoint has a DIFFERENT configuration. ",
        "Move/delete it intentionally."
      )
    }

    if (!identical(ck$observed, observed)) {
      stop(
        "Existing checkpoint has different observed data/results; ",
        "cannot resume."
      )
    }

    b_start <- ck$last_b + 1L
    exceeds <- ck$exceeds
    tie_sum <- ck$tie_sum

    cat(
      "Resuming from permutation ",
      b_start,
      " / ",
      B_PERM,
      "\n",
      sep = ""
    )

  } else {
    b_start <- 1L
    exceeds <- integer(nrow(observed))
    tie_sum <- matrix(0, nrow(observed), m)
  }

  save_checkpoint <- function(b) {
    tmp <- paste0(ckpt_file, ".tmp")

    saveRDS(
      list(
        config = config,
        observed = observed,
        last_b = b,
        exceeds = exceeds,
        tie_sum = tie_sum
      ),
      tmp
    )

    if (file.exists(ckpt_file)) {
      file.remove(ckpt_file)
    }

    if (!file.rename(tmp, ckpt_file)) {
      stop("Checkpoint rename failed: ", tmp)
    }
  }

  start <- Sys.time()

  if (b_start <= B_PERM) {
    for (b in seq.int(b_start, B_PERM)) {

      # Seed-per-permutation makes checkpoint/resume exactly reproducible.
      set.seed(MASTER_SEED + b)

      pc <- evaluate(
        sample(arm, length(arm), replace = FALSE)
      )

      for (i in seq_len(nrow(observed))) {
        a <- observed[i]
        pool <- pc[method_ids[[a$method]]]

        scores <- comparison_score(
          pool[[a$measure]],
          a$side
        )

        k <- which.max(scores)

        exceeds[i] <- exceeds[i] +
          as.integer(scores[k] >= a$comparison_statistic)

        for (j in seq_len(m)) {
          tie_sum[i, j] <- tie_sum[i, j] +
            pool[[paste0("tie_after_rank", j, "_prop")]][k]
        }
      }

      if (
        b == 1L ||
        b %% CHECKPOINT_EVERY == 0L ||
        b == B_PERM
      ) {
        save_checkpoint(b)

        mins <- as.numeric(
          difftime(Sys.time(), start, units = "mins")
        )

        cat(
          mode, ": ",
          b, "/", B_PERM,
          "; elapsed since resume ",
          round(mins, 1),
          " min\n",
          sep = ""
        )
      }
    }
  }

  final <- copy(observed)

  final[, permutation_exceed_count := exceeds]
  final[, permutation_p_value :=
          (permutation_exceed_count + 1) / (B_PERM + 1)]
  final[, significant_0_05 := permutation_p_value < 0.05]
  final[, B_permutations := B_PERM]
  final[, endpoint_set := mode]
  final[, endpoint_1 := endpoint_names[1]]
  final[, endpoint_2 := endpoint_names[2]]
  final[, endpoint_3 := endpoint_names[3]]
  final[, n_treatment := sum(arm == 1L)]
  final[, n_control := sum(arm == 0L)]

  for (j in seq_len(m)) {
    final[[paste0(
      "mean_permuted_tie_after_rank",
      j,
      "_prop"
    )]] <- tie_sum[, j] / B_PERM
  }

  setorder(final, measure, side, method)

  fwrite(
    final,
    file.path(
      odir,
      sprintf("DIG_%dENDPOINT_MAX_WR_WO_B%d.csv", m, B_PERM)
    )
  )

  # Compact summary table for review/meetings.
  keep <- c(
    "endpoint_set",
    "endpoint_1",
    "endpoint_2",
    "endpoint_3",
    "measure",
    "method",
    "side",
    "n_candidates",
    "max_statistic",
    "permutation_p_value",
    "selected_order",
    "selected_weight",
    "t_death_days",
    "t_hosp_days",
    paste0("p", seq_len(m)),
    paste0("tie_after_rank", seq_len(m), "_prop"),
    "B_permutations"
  )

  fwrite(
    final[, ..keep],
    file.path(
      odir,
      sprintf("DIG_%dENDPOINT_SUMMARY_B%d.csv", m, B_PERM)
    )
  )

  # Manuscript-facing two-sided-only summary.
  fwrite(
    final[side == "two", ..keep],
    file.path(
      odir,
      sprintf("DIG_%dENDPOINT_SUMMARY_TWOSIDED_B%d.csv", m, B_PERM)
    )
  )

  saveRDS(
    list(
      config = config,
      observed = observed,
      results = final,
      threshold_grid = tg
    ),
    complete_file
  )

  cat("DONE: ", odir, "\n", sep = "")
  invisible(final)
}

# Each endpoint set has its own checkpoint and output directory.
ALL_RESULTS <- lapply(SETS_TO_RUN, run_mode)
names(ALL_RESULTS) <- SETS_TO_RUN

combined <- rbindlist(
  ALL_RESULTS,
  fill = TRUE,
  use.names = TRUE
)

combined_keep <- c(
  "endpoint_set",
  "endpoint_1",
  "endpoint_2",
  "endpoint_3",
  "measure",
  "method",
  "side",
  "n_candidates",
  "max_statistic",
  "permutation_p_value",
  "selected_order",
  "selected_weight",
  "t_death_days",
  "t_hosp_days",
  "p1",
  "p2",
  "p3",
  "tie_after_rank1_prop",
  "tie_after_rank2_prop",
  "tie_after_rank3_prop",
  "B_permutations"
)

fwrite(
  combined[, ..combined_keep],
  file.path(
    OUTPUT_ROOT,
    sprintf("DIG_THREE_SETS_COMBINED_SUMMARY_B%d.csv", B_PERM)
  )
)

fwrite(
  combined,
  file.path(
    OUTPUT_ROOT,
    sprintf("DIG_THREE_SETS_COMBINED_FULL_B%d.csv", B_PERM)
  )
)

# Current manuscript-facing 3-endpoint results are TWO-SIDED only.
fwrite(
  combined[side == "two", ..combined_keep],
  file.path(
    OUTPUT_ROOT,
    sprintf("DIG_THREE_SETS_COMBINED_TWOSIDED_B%d.csv", B_PERM)
  )
)

cat("\nAll THREE requested endpoint sets complete. Fourth set excluded.\n")
