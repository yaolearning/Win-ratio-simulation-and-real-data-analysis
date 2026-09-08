# ============================================================================
# RUN_DIG_3ENDPOINT_STANDARD_MAX_WR_WO_OPENMP8.R
# DIG 3-endpoint standard maximized WR/WO analysis
# Correct permutation re-optimization + OpenMP 8-thread pairwise engine
#
# Endpoints:
#   1) Death                            : DEATHDAY + DEATH
#   2) First all-cause hospitalization : HOSPDAYS + HOSP
#   3) SVA hospitalization             : SVA (adverse binary; 1 worse)
#
# Methods:
#   Original
#   Six fixed orders (descriptive only)
#   M_ord
#   M_wt
#   M_ord_wt
#
# Threshold methods are intentionally not run in this script.
# ============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

# ----------------------------- User settings -----------------------------
B_PERM <- 500L
MASTER_SEED <- 20260908L
N_THREADS <- 8L
CHECKPOINT_EVERY <- 25L
FORCE_RERUN <- TRUE

# Put BOTH this R file and dig3_openmp_engine.cpp in D:/Win Ratio Research
CODE_DIR <- "D:/Win Ratio Research"
CPP_FILE <- file.path(CODE_DIR, "dig3_openmp_engine.cpp")

# Keep using your previous output folder
OUTPUT_DIR <- "D:/Win Ratio Research/DIG_3endpoint_37tests_results_openmp8"
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "optimized_standard_checkpoints")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Packages ----------------------------------
need <- c("asympTest", "Rcpp", "data.table")
for (p in need) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, type = "binary")
  }
}
library(asympTest)
library(Rcpp)
library(data.table)

# ----------------------------- Load DIG ----------------------------------
data("DIGdata", package = "asympTest")
dig <- DIGdata
cat("\nDIG loaded. N =", nrow(dig), "\n")

# ---------------------- Detect treatment column --------------------------
# Try common names first. If none are found, stop and print likely candidates.
arm_candidates <- c("TRTMT", "TRT", "TREAT", "TREATMENT", "ARM", "TRT_AB", "DIG")
ARM_SOURCE <- arm_candidates[arm_candidates %in% names(dig)][1]

if (is.na(ARM_SOURCE)) {
  cat("\nCould not auto-detect treatment column. Candidate columns:\n")
  print(grep("TRT|TREAT|ARM|DIG|PLACEBO", names(dig), value = TRUE, ignore.case = TRUE))
  stop("Set ARM_SOURCE manually near the top of the script.")
}
cat("Treatment column detected:", ARM_SOURCE, "\n")
print(table(dig[[ARM_SOURCE]], useNA = "ifany"))

# -------------------------- Check endpoints -------------------------------
required_vars <- c("DEATH", "DEATHDAY", "HOSP", "HOSPDAYS", "SVA")
missing_vars <- setdiff(required_vars, names(dig))
if (length(missing_vars) > 0) {
  cat("\nRelevant available columns:\n")
  print(grep("DEATH|HOSP|SVA|WHF|DAY", names(dig), value = TRUE, ignore.case = TRUE))
  stop("Missing required endpoint variable(s): ", paste(missing_vars, collapse = ", "))
}

# -------------------------- Recode treatment ------------------------------
arm_raw <- suppressWarnings(as.numeric(as.character(dig[[ARM_SOURCE]])))
uniq <- sort(unique(arm_raw[!is.na(arm_raw)]))

if (identical(uniq, c(0, 1))) {
  ARM <- as.integer(arm_raw)
} else if (identical(uniq, c(1, 2))) {
  # Assumption: 1=treatment, 2=control. If your labels indicate otherwise, change here.
  ARM <- ifelse(arm_raw == 1, 1L, 0L)
} else {
  stop("Unexpected treatment coding in ", ARM_SOURCE, ": ", paste(uniq, collapse = ", "))
}
cat("Recoded arm counts (1=treatment, 0=control):\n")
print(table(ARM))

# -------------------------- Endpoint matrices -----------------------------
N <- nrow(dig)
M_ENDPOINTS <- 3L
TYPE_CODE <- c(1L, 1L, 2L)   # time, time, adverse lower-is-better

TIME_MAT <- matrix(NA_real_, nrow = N, ncol = 3)
EVENT_MAT <- matrix(0L, nrow = N, ncol = 3)
VALUE_MAT <- matrix(NA_real_, nrow = N, ncol = 3)

TIME_MAT[, 1] <- as.numeric(dig$DEATHDAY)
EVENT_MAT[, 1] <- as.integer(as.numeric(dig$DEATH) == 1)

TIME_MAT[, 2] <- as.numeric(dig$HOSPDAYS)
EVENT_MAT[, 2] <- as.integer(as.numeric(dig$HOSP) == 1)

VALUE_MAT[, 3] <- as.numeric(dig$SVA)

ENDPOINT_NAMES <- c("Death", "First all-cause hospitalization", "SVA hospitalization")
cat("\nEndpoint summaries:\n")
cat("Death:\n"); print(table(EVENT_MAT[,1]))
cat("First all-cause hospitalization:\n"); print(table(EVENT_MAT[,2]))
cat("SVA:\n"); print(table(VALUE_MAT[,3], useNA = "ifany"))

# -------------------------- Orders and weights -----------------------------
ORDERS <- rbind(
  c(1L,2L,3L), c(1L,3L,2L), c(2L,1L,3L),
  c(2L,3L,1L), c(3L,1L,2L), c(3L,2L,1L)
)
ORDER_NAMES <- apply(ORDERS, 1, paste, collapse = "->")

WEIGHTS <- rbind(
  c(1,0,0),
  c(0.5,0.5,0),
  c(1/3,1/3,1/3)
)
WEIGHT_NAMES <- apply(WEIGHTS, 1, function(x) {
  paste0("(", paste(sprintf("%.3f", x), collapse = ","), ")")
})
NATURAL_ORDER_ID <- 1L
EQUAL_WEIGHT_ID <- 3L
MAIN_METHODS <- c("Original", "M_ord", "M_wt", "M_ord_wt")

# ---------------------------- Compile C++ ---------------------------------
if (!file.exists(CPP_FILE)) stop("C++ file not found: ", CPP_FILE)
Sys.setenv(OMP_NUM_THREADS = as.character(N_THREADS))
Sys.setenv(OMP_THREAD_LIMIT = as.character(N_THREADS))
cat("\nCompiling:", CPP_FILE, "\n")
Rcpp::sourceCpp(CPP_FILE, rebuild = TRUE, verbose = FALSE)
cat("OpenMP max threads reported by engine:", openmp_max_threads_cpp(), "\n")

# ----------------------------- Helpers ------------------------------------
ratio_safe <- function(num, den) {
  if (is.na(num) || is.na(den)) return(NA_real_)
  if (num == 0 && den == 0) return(1)
  if (den == 0 && num > 0) return(Inf)
  if (num == 0 && den > 0) return(0)
  num / den
}

abslog_safe <- function(x) {
  if (is.na(x)) return(NA_real_)
  if (!is.finite(x) || x <= 0) return(Inf)
  abs(log(x))
}

calculate_candidate_statistics <- function(engine_result) {
  wins <- engine_result$wins
  losses <- engine_result$losses
  ties <- engine_result$ties
  n_pairs <- engine_result$n_pairs

  rows <- vector("list", nrow(ORDERS) * nrow(WEIGHTS))
  rr <- 1L
  for (o in seq_len(nrow(ORDERS))) {
    for (w in seq_len(nrow(WEIGHTS))) {
      ww <- WEIGHTS[w,]
      weighted_win <- sum(ww * wins[o,])
      weighted_loss <- sum(ww * losses[o,])
      WR <- ratio_safe(weighted_win, weighted_loss)
      WO <- ratio_safe(weighted_win + 0.5 * ties[o], weighted_loss + 0.5 * ties[o])

      rows[[rr]] <- data.table(
        order_id = o,
        order = ORDER_NAMES[o],
        weight_id = w,
        weight = WEIGHT_NAMES[w],
        p1 = ww[1], p2 = ww[2], p3 = ww[3],
        WR = WR,
        WO = WO,
        WR_abslog = abslog_safe(WR),
        WO_abslog = abslog_safe(WO),
        weighted_win = weighted_win,
        weighted_loss = weighted_loss,
        ties = ties[o],
        tie_proportion = ties[o] / n_pairs,
        win_rank1 = wins[o,1], win_rank2 = wins[o,2], win_rank3 = wins[o,3],
        loss_rank1 = losses[o,1], loss_rank2 = losses[o,2], loss_rank3 = losses[o,3],
        n_pairs = n_pairs
      )
      rr <- rr + 1L
    }
  }
  rbindlist(rows)
}

evaluate_dataset <- function(arm_vector) {
  eng <- evaluate_six_orders_openmp(
    arm = as.integer(arm_vector),
    type_code = as.integer(TYPE_CODE),
    time_mat = TIME_MAT,
    event_mat = EVENT_MAT,
    value_mat = VALUE_MAT,
    orders = ORDERS,
    n_threads = N_THREADS
  )
  calculate_candidate_statistics(eng)
}

get_method_candidates <- function(candidates, method_name) {
  if (method_name == "Original") {
    return(candidates[order_id == NATURAL_ORDER_ID & weight_id == EQUAL_WEIGHT_ID])
  }
  if (method_name == "M_ord") {
    return(candidates[weight_id == EQUAL_WEIGHT_ID])
  }
  if (method_name == "M_wt") {
    return(candidates[order_id == NATURAL_ORDER_ID])
  }
  if (method_name == "M_ord_wt") {
    return(candidates)
  }
  stop("Unknown method: ", method_name)
}

select_best <- function(candidates, measure = c("WR","WO"), side = c("one","two")) {
  measure <- match.arg(measure)
  side <- match.arg(side)
  stat <- if (side == "one") candidates[[measure]] else candidates[[paste0(measure, "_abslog")]]
  if (all(is.na(stat))) stop("All candidate statistics are NA.")
  candidates[which.max(stat)]
}

# --------------------------- Observed data --------------------------------
cat("\nEvaluating observed data...\n")
obs_start <- Sys.time()
OBS_CANDIDATES <- evaluate_dataset(ARM)
cat("Observed evaluation seconds:", round(as.numeric(difftime(Sys.time(), obs_start, units="secs")), 2), "\n")

fwrite(OBS_CANDIDATES, file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_observed_18_candidates.csv"))
FIXED_ORDER_RESULTS <- OBS_CANDIDATES[weight_id == EQUAL_WEIGHT_ID]
fwrite(FIXED_ORDER_RESULTS, file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_3ENDPOINT_FIXED_ORDERS.csv"))

obs_rows <- list(); rr <- 1L
for (measure_name in c("WR","WO")) {
  for (method_name in MAIN_METHODS) {
    d <- get_method_candidates(OBS_CANDIDATES, method_name)
    for (side_name in c("one","two")) {
      best <- select_best(d, measure_name, side_name)
      comp_stat <- if (side_name == "one") best[[measure_name]] else best[[paste0(measure_name, "_abslog")]]
      obs_rows[[rr]] <- data.table(
        measure = measure_name,
        method = method_name,
        side = side_name,
        statistic = best[[measure_name]],
        comparison_statistic = comp_stat,
        order = best$order,
        weight = best$weight,
        p1 = best$p1, p2 = best$p2, p3 = best$p3,
        tie_proportion = best$tie_proportion
      )
      rr <- rr + 1L
    }
  }
}
OBS_SELECTIONS <- rbindlist(obs_rows)
fwrite(OBS_SELECTIONS, file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_observed_method_selections.csv"))
print(OBS_SELECTIONS)

# --------------------------- Permutations ---------------------------------
set.seed(MASTER_SEED)
PERM_ARM_MATRIX <- replicate(B_PERM, sample(ARM, replace = FALSE), simplify = "matrix")
checkpoint_file <- file.path(CHECKPOINT_DIR, paste0("standard_max_WR_WO_B", B_PERM, "_checkpoint.rds"))

keys <- paste(OBS_SELECTIONS$measure, OBS_SELECTIONS$method, OBS_SELECTIONS$side, sep = "__")
EXCEED <- setNames(rep(0L, length(keys)), keys)
PERM_SELECTIONS <- data.table()
start_b <- 1L

if (file.exists(checkpoint_file) && !FORCE_RERUN) {
  chk <- readRDS(checkpoint_file)
  start_b <- chk$next_b
  EXCEED <- chk$EXCEED
  PERM_SELECTIONS <- chk$PERM_SELECTIONS
  cat("Resuming from permutation", start_b, "\n")
}

cat("\nStarting", B_PERM, "permutations with", N_THREADS, "OpenMP threads...\n")
perm_start <- Sys.time()

if (start_b <= B_PERM) {
  for (b in start_b:B_PERM) {
    PERM_CANDIDATES <- evaluate_dataset(PERM_ARM_MATRIX[,b])
    perm_rows <- list(); pr <- 1L

    for (measure_name in c("WR","WO")) {
      for (method_name in MAIN_METHODS) {
        d <- get_method_candidates(PERM_CANDIDATES, method_name)

        for (side_name in c("one","two")) {
          best <- select_best(d, measure_name, side_name)
          comp_stat <- if (side_name == "one") best[[measure_name]] else best[[paste0(measure_name, "_abslog")]]

          obs_row <- OBS_SELECTIONS[
            measure == measure_name & method == method_name & side == side_name
          ]

          key <- paste(measure_name, method_name, side_name, sep = "__")
          if (!is.na(comp_stat) && !is.na(obs_row$comparison_statistic) && comp_stat >= obs_row$comparison_statistic) {
            EXCEED[key] <- EXCEED[key] + 1L
          }

          perm_rows[[pr]] <- data.table(
            b = b,
            measure = measure_name,
            method = method_name,
            side = side_name,
            statistic = best[[measure_name]],
            comparison_statistic = comp_stat,
            selected_order = best$order,
            selected_weight = best$weight,
            selected_p1 = best$p1,
            selected_p2 = best$p2,
            selected_p3 = best$p3,
            tie_proportion = best$tie_proportion
          )
          pr <- pr + 1L
        }
      }
    }

    PERM_SELECTIONS <- rbind(PERM_SELECTIONS, rbindlist(perm_rows), fill = TRUE)

    if (b == 1L || b %% 5L == 0L) {
      elapsed_min <- as.numeric(difftime(Sys.time(), perm_start, units = "mins"))
      done <- b - start_b + 1L
      avg_sec <- elapsed_min * 60 / done
      eta_min <- avg_sec * (B_PERM - b) / 60
      cat(sprintf("Permutation %d/%d | elapsed %.2f min | ETA %.2f min\n", b, B_PERM, elapsed_min, eta_min))
    }

    if (b %% CHECKPOINT_EVERY == 0L || b == B_PERM) {
      saveRDS(list(
        next_b = b + 1L,
        EXCEED = EXCEED,
        PERM_SELECTIONS = PERM_SELECTIONS
      ), checkpoint_file)
      cat("Checkpoint saved at", b, "\n")
    }
  }
}

# --------------------------- Final results --------------------------------
FINAL_RESULTS <- copy(OBS_SELECTIONS)
FINAL_RESULTS[, key := paste(measure, method, side, sep = "__")]
FINAL_RESULTS[, exceed := as.integer(EXCEED[key])]
FINAL_RESULTS[, permutation_p := (exceed + 1) / (B_PERM + 1)]
FINAL_RESULTS[, significant_0.05 := permutation_p < 0.05]

perm_ties <- PERM_SELECTIONS[, .(
  mean_perm_tie_proportion = mean(tie_proportion, na.rm = TRUE)
), by = .(measure, method, side)]
FINAL_RESULTS <- merge(FINAL_RESULTS, perm_ties, by = c("measure","method","side"), all.x = TRUE, sort = FALSE)

FINAL_RESULTS[, n_candidates := fifelse(
  method == "Original", 1L,
  fifelse(method == "M_ord", 6L,
          fifelse(method == "M_wt", 3L, 18L))
)]

final_file <- file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_3ENDPOINT_WR_WO_RESULTS_B500.csv")
fwrite(FINAL_RESULTS, final_file)
fwrite(PERM_SELECTIONS, file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_3ENDPOINT_PERMUTATION_SELECTIONS_B500.csv"))

saveRDS(list(
  final_results = FINAL_RESULTS,
  fixed_orders = FIXED_ORDER_RESULTS,
  observed_candidates = OBS_CANDIDATES,
  permutation_selections = PERM_SELECTIONS,
  B = B_PERM,
  seed = MASTER_SEED,
  threads = N_THREADS,
  endpoint_names = ENDPOINT_NAMES,
  orders = ORDERS,
  weights = WEIGHTS
), file.path(OUTPUT_DIR, "OPTIMIZED_STANDARD_3ENDPOINT_FULL_RESULT_B500.rds"))

cat("\n================ FINAL RESULTS ================\n")
print(FINAL_RESULTS[, .(
  measure, method, side, n_candidates, statistic, permutation_p,
  order, weight, tie_proportion, mean_perm_tie_proportion
)])
cat("\nSaved to:\n", final_file, "\n")
cat("Total runtime (minutes):", round(as.numeric(difftime(Sys.time(), perm_start, units = "mins")), 2), "\n")
