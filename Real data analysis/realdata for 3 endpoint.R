# Clean environment
rm(list = ls())
options(stringsAsFactors = FALSE)


#USER CONFIGURATION
#Input and output paths
DATA_PATH <- "C:/Users/Lenovo/Desktop/your_real_data.csv"
DATA_OBJECT_NAME <- NULL

OUTPUT_DIR <- "C:/Users/Lenovo/Desktop/realdata_3endpoint_WR_WO_single_permutation_multi_threshold"

# Data format
#Use "wide" if each subject has one row.
#Use "event_long" if the data has one row per event, e.g. ID/time/status/trt.
DATA_FORMAT <- "wide"  # "wide" or "event_long"

#wide format:
#ID_COL, ARM_COL, and endpoint columns in ENDPOINTS_WIDE.
#ARM_COL must be coded as 1=treatment A and 0=control B.
ID_COL <- "SUBJID"
ARM_COL <- "ARM"

#Endpoint types supported:
#type = "time"       adverse time-to-event endpoint; larger event-free time is better
#                     needs time_col and event_col, event_col coded 1=event, 0=censoring
#type = "count"      recurrent/count adverse endpoint; lower count is better; needs count_col
#type = "binary"     adverse binary endpoint; lower value is better; needs value_col
#type = "continuous" continuous endpoint; set higher_better TRUE/FALSE; needs value_col
#For 3 endpoints, put three entries here.
#For 2 endpoints, put two entries here; the code still works.
ENDPOINTS_WIDE <- list(
  list(name = "Outcome 1", type = "time",  time_col = "FUTIME", event_col = "CNSR"),
  list(name = "Outcome 2", type = "count", count_col = "NUMHOSP"),
  # Replace this with your real third endpoint. Delete it for two-endpoint data.
  list(name = "Outcome 3", type = "count", count_col = "ENDPOINT3_COUNT")
)

# Event-long format
# Required columns for event_long format:
# EVENT_ID_COL, EVENT_ARM_COL, EVENT_TIME_COL, EVENT_STATUS_COL.
# status code 0 is treated as censoring/follow-up row.
EVENT_ID_COL <- "patid"
EVENT_ARM_COL <- "trt_ab"
EVENT_TIME_COL <- "time"
EVENT_STATUS_COL <- "status"

# For event_long format, define endpoints by status_code.
# type = "time"  : first event of that status code; otherwise censored at max follow-up
# type = "count" : number of events with that status code
EVENT_ENDPOINTS <- list(
  list(name = "Death", type = "time",  status_code = 1),
  list(name = "Hospitalization", type = "count", status_code = 2),
  ## Replace or delete depending on your data
  list(name = "Endpoint 3", type = "count", status_code = 3)
)



#Thresholds are applied only to time-to-event endpoints.
#They follow the endpoint number, not the rank in a selected ordering.
#If THRESHOLD_TIME_ENDPOINTS = "auto", every endpoint with type = "time"
#receives its own threshold grid. Therefore, if there are two time-to-event
#endpoints, the threshold search uses two thresholds, one for each time endpoint.
#You may also specify a subset manually, e.g. c(1, 3). Every listed endpoint
THRESHOLD_TIME_ENDPOINTS <- "auto"  # "auto" or an integer vector such as c(1, 3)


#Units must match the data time unit.
#Non-time endpoints are ignored and always use threshold 0.
#Example:
#THRESHOLD_GRID_BY_ENDPOINT <- list(
#`1` = c(0, 3, 6, 12, 18, 24),
#`2` = NULL,
#`3` = c(0, 6, 12, 24)
THRESHOLD_GRID_BY_ENDPOINT <- list(
  `1` = NULL,
  `2` = NULL,
  `3` = NULL
)

TIME_UNIT_LABEL <- "months"
MAX_AUTO_THRESHOLDS_PER_ENDPOINT <- 10L

#Weight search configuration
#WEIGHT_MODE = "vertices" uses vertices of the ordered weight simplex.
#For 3 endpoints with p1 >= p2 >= p3 >= 0 and sum=1, vertices are:
#(1,0,0), (0.5,0.5,0), (1/3,1/3,1/3).
WEIGHT_MODE <- "vertices"  # "vertices" or "grid"
WEIGHT_STEP <- 0.05

#Default is ordered weights: the first endpoint in the selected order gets
#the largest weight, the second endpoint gets the second largest weight, and
#the last endpoint gets the smallest weight.
WEIGHT_CONSTRAINT <- "ordered"  # "ordered" or "simplex"


#permutation configuration
B_PERM <- 500L
SEED <- 20260818
CHECKPOINT_EVERY <- 25L
FORCE_RERUN <- FALSE


#Method/output configuration
REPORT_ALL_FIXED_ORDERS <- TRUE
RUN_WR <- TRUE
RUN_WO <- TRUE
RUN_LOGRANK <- TRUE
RUN_COMPOSITE_LOGRANK_IF_POSSIBLE <- TRUE


#Packages
need_pkgs <- c("data.table", "survival", "ggplot2")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(data.table)
library(survival)
library(ggplot2)

for (d in c(OUTPUT_DIR, file.path(OUTPUT_DIR, "checkpoints"), file.path(OUTPUT_DIR, "figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

set.seed(SEED)

#Helper functions
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

ratio_safe <- function(num, den) {
  if (is.na(num) || is.na(den)) return(NA_real_)
  if (den == 0 && num == 0) return(1)
  if (den == 0 && num > 0) return(Inf)
  if (num == 0 && den > 0) return(0)
  num / den
}

abslog_safe <- function(x) {
  if (is.na(x)) return(NA_real_)
  if (is.infinite(x)) return(Inf)
  if (x <= 0) return(Inf)
  abs(log(x))
}

format_order <- function(order_vec) paste(order_vec, collapse = "->")
parse_order <- function(order_key) as.integer(strsplit(as.character(order_key), "->", fixed = TRUE)[[1]])

all_permutations <- function(x) {
  x <- as.integer(x)
  if (length(x) == 1) return(list(x))
  out <- list()
  idx <- 1L
  for (i in seq_along(x)) {
    rest <- x[-i]
    sub <- all_permutations(rest)
    for (s in sub) {
      out[[idx]] <- c(x[i], s)
      idx <- idx + 1L
    }
  }
  out
}

clean_method_file_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

make_weight_string <- function(w) paste0("(", paste(sprintf("%.3f", w), collapse = ","), ")")

make_threshold_string <- function(thr_vec) {
  paste0("(", paste(sprintf("%.4f", as.numeric(thr_vec)), collapse = ","), ")")
}


#Read and prepare data

read_any_data <- function(path, object_name = NULL) {
  if (!file.exists(path)) stop("DATA_PATH not found: ", path)
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    return(as.data.frame(data.table::fread(path, data.table = FALSE)))
  }
  
  if (ext == "rds") {
    obj <- readRDS(path)
    if (!is.data.frame(obj)) stop("RDS file is not a data frame.")
    return(as.data.frame(obj))
  }
  
  if (ext %in% c("rdata", "rda")) {
    env <- new.env()
    nm <- load(path, envir = env)
    if (!is.null(object_name)) {
      if (!exists(object_name, envir = env)) stop("DATA_OBJECT_NAME not found in RData: ", object_name)
      obj <- get(object_name, envir = env)
      if (!is.data.frame(obj)) stop("DATA_OBJECT_NAME is not a data frame.")
      return(as.data.frame(obj))
    }
    df_names <- nm[sapply(nm, function(z) is.data.frame(get(z, envir = env)))]
    if (length(df_names) == 0) stop("No data frame found in RData.")
    if (length(df_names) > 1) {
      message("Multiple data frames found in RData. Using the first one: ", df_names[1])
      message("Set DATA_OBJECT_NAME if this is not intended.")
    }
    return(as.data.frame(get(df_names[1], envir = env)))
  }
  
  stop("Unsupported file extension: ", ext)
}

prepare_wide_data <- function(raw, id_col, arm_col, endpoints) {
  if (!id_col %in% names(raw)) stop("ID_COL not found: ", id_col)
  if (!arm_col %in% names(raw)) stop("ARM_COL not found: ", arm_col)
  
  d <- data.frame(
    id = raw[[id_col]],
    arm = safe_num(raw[[arm_col]]),
    stringsAsFactors = FALSE
  )
  
  if (!all(d$arm %in% c(0, 1))) {
    stop("ARM_COL must be coded as 1=treatment A and 0=control B.")
  }
  
  endpoint_specs <- list()
  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]
    ep_type <- tolower(ep$type)
    ep_name <- if (!is.null(ep$name)) ep$name else paste0("Outcome ", j)
    
    if (ep_type == "time") {
      if (!ep$time_col %in% names(raw)) stop("Missing time column for endpoint ", j, ": ", ep$time_col)
      if (!ep$event_col %in% names(raw)) stop("Missing event column for endpoint ", j, ": ", ep$event_col)
      tcol <- paste0("e", j, "_time")
      ecol <- paste0("e", j, "_event")
      d[[tcol]] <- safe_num(raw[[ep$time_col]])
      d[[ecol]] <- as.integer(safe_num(raw[[ep$event_col]]) == 1)
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "time", time_col = tcol, event_col = ecol)
      
    } else if (ep_type == "count") {
      if (!ep$count_col %in% names(raw)) stop("Missing count column for endpoint ", j, ": ", ep$count_col)
      ccol <- paste0("e", j, "_count")
      d[[ccol]] <- safe_num(raw[[ep$count_col]])
      d[[ccol]][is.na(d[[ccol]])] <- 0
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "count", count_col = ccol)
      
    } else if (ep_type == "binary") {
      if (!ep$value_col %in% names(raw)) stop("Missing binary value column for endpoint ", j, ": ", ep$value_col)
      vcol <- paste0("e", j, "_value")
      d[[vcol]] <- safe_num(raw[[ep$value_col]])
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "binary", value_col = vcol)
      
    } else if (ep_type == "continuous") {
      if (!ep$value_col %in% names(raw)) stop("Missing continuous value column for endpoint ", j, ": ", ep$value_col)
      vcol <- paste0("e", j, "_value")
      d[[vcol]] <- safe_num(raw[[ep$value_col]])
      higher_better <- if (!is.null(ep$higher_better)) isTRUE(ep$higher_better) else TRUE
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "continuous", value_col = vcol, higher_better = higher_better)
      
    } else {
      stop("Unsupported endpoint type for endpoint ", j, ": ", ep$type)
    }
  }
  
  d <- d[!is.na(d$id) & !is.na(d$arm), , drop = FALSE]
  if (any(duplicated(d$id))) stop("Wide data must have one row per subject. Duplicated IDs found.")
  
  list(data = d, endpoints = endpoint_specs)
}

prepare_event_long_data <- function(raw, id_col, arm_col, time_col, status_col, endpoints) {
  required <- c(id_col, arm_col, time_col, status_col)
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) stop("Missing event-long columns: ", paste(missing, collapse = ", "))
  
  dt <- as.data.table(raw)
  dt[, .id_tmp := get(id_col)]
  dt[, .arm_tmp := safe_num(get(arm_col))]
  dt[, .time_tmp := safe_num(get(time_col))]
  dt[, .status_tmp := safe_num(get(status_col))]
  
  subj <- dt[, .(
    id = first_non_missing(.id_tmp),
    arm = first_non_missing(.arm_tmp),
    followup_time = max(.time_tmp, na.rm = TRUE)
  ), by = .id_tmp]
  subj[, .id_tmp := NULL]
  
  if (!all(subj$arm %in% c(0, 1))) stop("EVENT_ARM_COL must be coded as 1=treatment A and 0=control B.")
  
  endpoint_specs <- list()
  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]
    ep_type <- tolower(ep$type)
    ep_name <- if (!is.null(ep$name)) ep$name else paste0("Outcome ", j)
    code <- ep$status_code
    
    if (ep_type == "time") {
      tmp <- dt[.status_tmp == code, .(event_time = min(.time_tmp, na.rm = TRUE)), by = .id_tmp]
      names(tmp)[names(tmp) == ".id_tmp"] <- "id"
      subj <- merge(subj, tmp, by = "id", all.x = TRUE)
      tcol <- paste0("e", j, "_time")
      ecol <- paste0("e", j, "_event")
      subj[[ecol]] <- ifelse(is.na(subj$event_time), 0L, 1L)
      subj[[tcol]] <- ifelse(is.na(subj$event_time), subj$followup_time, subj$event_time)
      subj$event_time <- NULL
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "time", time_col = tcol, event_col = ecol)
      
    } else if (ep_type == "count") {
      tmp <- dt[, .(event_count = sum(.status_tmp == code, na.rm = TRUE)), by = .id_tmp]
      names(tmp)[names(tmp) == ".id_tmp"] <- "id"
      subj <- merge(subj, tmp, by = "id", all.x = TRUE)
      ccol <- paste0("e", j, "_count")
      subj[[ccol]] <- subj$event_count
      subj[[ccol]][is.na(subj[[ccol]])] <- 0
      subj$event_count <- NULL
      endpoint_specs[[j]] <- list(id = j, name = ep_name, type = "count", count_col = ccol)
      
    } else {
      stop("For event_long data, endpoint type must be 'time' or 'count'.")
    }
  }
  
  subj <- as.data.frame(subj)
  list(data = subj, endpoints = endpoint_specs)
}

raw_data <- read_any_data(DATA_PATH, DATA_OBJECT_NAME)

if (DATA_FORMAT == "wide") {
  prepared <- prepare_wide_data(raw_data, ID_COL, ARM_COL, ENDPOINTS_WIDE)
} else if (DATA_FORMAT == "event_long") {
  prepared <- prepare_event_long_data(raw_data, EVENT_ID_COL, EVENT_ARM_COL, EVENT_TIME_COL, EVENT_STATUS_COL, EVENT_ENDPOINTS)
} else {
  stop("DATA_FORMAT must be 'wide' or 'event_long'.")
}

subject_data <- prepared$data
endpoint_specs <- prepared$endpoints
M_ENDPOINTS <- length(endpoint_specs)
if (M_ENDPOINTS < 2 || M_ENDPOINTS > 3) stop("This script supports 2 or 3 endpoints.")

cat("Prepared subject-level data:\n")
cat("  n subjects:", nrow(subject_data), "\n")
cat("  treatment A:", sum(subject_data$arm == 1), "\n")
cat("  control B:", sum(subject_data$arm == 0), "\n")
cat("  endpoints:", M_ENDPOINTS, "\n")


#endpoint comparison
compare_endpoint_pairs <- function(data, idx_A, idx_B, endpoint_spec, threshold = 0) {
  n <- length(idx_A)
  out <- integer(n)  # +1 treatment wins, -1 control wins, 0 tie/unresolved
  threshold <- ifelse(is.na(threshold), 0, threshold)
  
  if (endpoint_spec$type == "time") {
    ta <- safe_num(data[[endpoint_spec$time_col]][idx_A])
    tb <- safe_num(data[[endpoint_spec$time_col]][idx_B])
    ea <- as.integer(data[[endpoint_spec$event_col]][idx_A] == 1)
    eb <- as.integer(data[[endpoint_spec$event_col]][idx_B] == 1)
    
    both_event <- ea == 1 & eb == 1 & !is.na(ta) & !is.na(tb)
    out[both_event & (ta - tb > threshold)] <- 1L
    out[both_event & (tb - ta > threshold)] <- -1L
    
    ## A censored after B has event: A has known longer event-free time.
    a_cens_b_event <- ea == 0 & eb == 1 & !is.na(ta) & !is.na(tb)
    out[a_cens_b_event & (ta - tb > threshold)] <- 1L
    
    ## B censored after A has event: B has known longer event-free time.
    a_event_b_cens <- ea == 1 & eb == 0 & !is.na(ta) & !is.na(tb)
    out[a_event_b_cens & (tb - ta > threshold)] <- -1L
    
    return(out)
  }
  
  if (endpoint_spec$type == "count") {
    ca <- safe_num(data[[endpoint_spec$count_col]][idx_A])
    cb <- safe_num(data[[endpoint_spec$count_col]][idx_B])
    ok <- !is.na(ca) & !is.na(cb)
    #adverse count: lower is better. Count endpoints do not receive thresholding here.
    out[ok & (cb - ca > 0)] <- 1L
    out[ok & (ca - cb > 0)] <- -1L
    return(out)
  }
  
  if (endpoint_spec$type == "binary") {
    va <- safe_num(data[[endpoint_spec$value_col]][idx_A])
    vb <- safe_num(data[[endpoint_spec$value_col]][idx_B])
    ok <- !is.na(va) & !is.na(vb)
    #adverse binary: lower is better. Binary endpoints do not receive thresholding here.
    out[ok & (vb - va > 0)] <- 1L
    out[ok & (va - vb > 0)] <- -1L
    return(out)
  }
  
  if (endpoint_spec$type == "continuous") {
    va <- safe_num(data[[endpoint_spec$value_col]][idx_A])
    vb <- safe_num(data[[endpoint_spec$value_col]][idx_B])
    ok <- !is.na(va) & !is.na(vb)
    higher_better <- isTRUE(endpoint_spec$higher_better)
    #Continuous endpoints do not receive thresholding here unless separately modified.
    if (higher_better) {
      out[ok & (va - vb > 0)] <- 1L
      out[ok & (vb - va > 0)] <- -1L
    } else {
      out[ok & (vb - va > 0)] <- 1L
      out[ok & (va - vb > 0)] <- -1L
    }
    return(out)
  }
  
  stop("Unsupported endpoint type: ", endpoint_spec$type)
}

get_pair_indices <- function(arm_vec) {
  idx_A <- which(arm_vec == 1)
  idx_B <- which(arm_vec == 0)
  if (length(idx_A) == 0 || length(idx_B) == 0) stop("Both arms must have at least one subject.")
  list(
    A = rep(idx_A, each = length(idx_B)),
    B = rep(idx_B, times = length(idx_A)),
    n_pairs = length(idx_A) * length(idx_B)
  )
}

counts_for_candidate <- function(data, arm_vec, order_vec, weights, thresholds_by_endpoint, pairs = NULL) {
  if (is.null(pairs)) pairs <- get_pair_indices(arm_vec)
  idx_A <- pairs$A
  idx_B <- pairs$B
  n_pairs <- pairs$n_pairs
  
  wins_by_rank <- rep(0, M_ENDPOINTS)
  losses_by_rank <- rep(0, M_ENDPOINTS)
  unresolved <- rep(TRUE, n_pairs)
  
  for (r in seq_along(order_vec)) {
    endpoint_id <- order_vec[r]
    ep <- endpoint_specs[[endpoint_id]]
    
    #Threshold follows the endpoint itself and is applied only to time-to-event endpoints.
    thr <- if (identical(tolower(ep$type), "time")) thresholds_by_endpoint[endpoint_id] else 0
    if (is.na(thr)) thr <- 0
    
    cmp <- compare_endpoint_pairs(data, idx_A, idx_B, ep, threshold = thr)
    unresolved_idx <- which(unresolved)
    if (length(unresolved_idx) == 0) break
    cmp_u <- cmp[unresolved_idx]
    
    wins_by_rank[r] <- sum(cmp_u == 1L, na.rm = TRUE)
    losses_by_rank[r] <- sum(cmp_u == -1L, na.rm = TRUE)
    
    resolved_local <- which(cmp_u != 0L)
    if (length(resolved_local) > 0) {
      unresolved[unresolved_idx[resolved_local]] <- FALSE
    }
  }
  
  tie_count <- sum(unresolved)
  weights <- as.numeric(weights)
  if (length(weights) < M_ENDPOINTS) weights <- c(weights, rep(0, M_ENDPOINTS - length(weights)))
  weights <- weights[seq_len(M_ENDPOINTS)]
  
  weighted_win <- sum(weights * wins_by_rank)
  weighted_loss <- sum(weights * losses_by_rank)
  
  WR_statistic <- ratio_safe(weighted_win, weighted_loss)
  WO_statistic <- ratio_safe(weighted_win + 0.5 * tie_count,
                             weighted_loss + 0.5 * tie_count)
  
  list(
    WR_statistic = WR_statistic,
    WO_statistic = WO_statistic,
    WR_abslog = abslog_safe(WR_statistic),
    WO_abslog = abslog_safe(WO_statistic),
    weighted_win = weighted_win,
    weighted_loss = weighted_loss,
    tie_count = tie_count,
    tie_proportion = tie_count / n_pairs,
    n_pairs = n_pairs,
    wins_by_rank = wins_by_rank,
    losses_by_rank = losses_by_rank
  )
}


#Threshold and weight grids
auto_threshold_grid <- function(data, endpoint_spec, max_thresholds = 10L) {
  ## Threshold grid is defined only for time-to-event endpoints.
  if (is.null(endpoint_spec) || !identical(tolower(endpoint_spec$type), "time")) return(0)
  
  vals <- c(0)
  t <- safe_num(data[[endpoint_spec$time_col]])
  e <- as.integer(data[[endpoint_spec$event_col]] == 1)
  t_event <- t[e == 1 & !is.na(t)]
  t_all <- t[!is.na(t)]
  
  if (length(t_all) == 0) return(0)
  max_t <- max(t_all, na.rm = TRUE)
  
  clinical <- c(1, 3, 6, 12, 18, 24, 36, 48, 60)
  clinical <- clinical[clinical > 0 & clinical < max_t]
  
  if (length(t_event) >= 2) {
    if (length(t_event) > 800) t_event <- sample(t_event, 800)
    diffs <- abs(as.vector(stats::dist(t_event)))
    diffs <- diffs[is.finite(diffs) & diffs > 0]
    qs <- as.numeric(stats::quantile(diffs, probs = c(0.10, 0.25, 0.50, 0.75, 0.90), na.rm = TRUE))
    vals <- c(vals, clinical, qs)
  } else {
    vals <- c(vals, clinical)
  }
  
  vals <- sort(unique(round(vals[is.finite(vals) & vals >= 0], 4)))
  if (length(vals) > max_thresholds) {
    keep_idx <- unique(round(seq(1, length(vals), length.out = max_thresholds)))
    vals <- vals[keep_idx]
    vals <- sort(unique(c(0, vals)))
  }
  vals
}

resolve_threshold_time_endpoints <- function(endpoint_specs, user_setting) {
  time_ids <- which(vapply(endpoint_specs, function(ep) identical(tolower(ep$type), "time"), logical(1)))
  
  if (length(time_ids) == 0) return(integer(0))
  
  if (is.character(user_setting) && length(user_setting) == 1 && tolower(user_setting) == "auto") {
    return(time_ids)
  }
  
  ids <- as.integer(user_setting)
  if (any(is.na(ids))) stop("THRESHOLD_TIME_ENDPOINTS must be 'auto' or an integer vector of endpoint IDs.")
  bad <- setdiff(ids, time_ids)
  if (length(bad) > 0) {
    stop("Thresholds can only be assigned to time-to-event endpoints. Invalid endpoint(s): ", paste(bad, collapse = ", "))
  }
  ids
}

get_manual_threshold_grid <- function(endpoint_id) {
  val <- THRESHOLD_GRID_BY_ENDPOINT[[as.character(endpoint_id)]]
  if (is.null(val)) return(NULL)
  val <- sort(unique(as.numeric(val)))
  val <- val[is.finite(val) & val >= 0]
  if (!0 %in% val) val <- sort(unique(c(0, val)))
  val
}

make_threshold_grid_df <- function(data, endpoint_specs) {
  active_ids <- resolve_threshold_time_endpoints(endpoint_specs, THRESHOLD_TIME_ENDPOINTS)
  m <- length(endpoint_specs)
  
  values_list <- vector("list", m)
  source_vec <- rep("none", m)
  
  for (j in seq_len(m)) {
    colname <- paste0("threshold_e", j)
    if (j %in% active_ids) {
      manual <- get_manual_threshold_grid(j)
      if (!is.null(manual)) {
        values_list[[j]] <- manual
        source_vec[j] <- "manual"
      } else {
        values_list[[j]] <- auto_threshold_grid(data, endpoint_specs[[j]], MAX_AUTO_THRESHOLDS_PER_ENDPOINT)
        source_vec[j] <- "auto"
      }
    } else {
      values_list[[j]] <- 0
      source_vec[j] <- if (identical(tolower(endpoint_specs[[j]]$type), "time")) "inactive_time_endpoint" else "not_time_endpoint"
    }
    names(values_list)[j] <- colname
  }
  
  dt <- do.call(data.table::CJ, c(values_list, list(sorted = FALSE)))
  thr_cols <- paste0("threshold_e", seq_len(m))
  dt[, threshold_key := apply(.SD, 1, function(z) paste(sprintf("%.4f", as.numeric(z)), collapse = "|")), .SDcols = thr_cols]
  
  info <- data.table(
    endpoint = seq_len(m),
    endpoint_name = vapply(endpoint_specs, function(ep) ep$name, character(1)),
    endpoint_type = vapply(endpoint_specs, function(ep) ep$type, character(1)),
    threshold_active = seq_len(m) %in% active_ids,
    threshold_source = source_vec,
    threshold_values = vapply(seq_len(m), function(j) paste(values_list[[j]], collapse = ","), character(1))
  )
  
  list(grid = dt, info = info, active_ids = active_ids)
}

make_weight_grid <- function(m, mode = "vertices", step = 0.05, constraint = "ordered") {
  if (m == 2) {
    if (mode == "vertices") {
      if (constraint == "ordered") {
        w <- rbind(c(1, 0), c(0.5, 0.5))
      } else {
        w <- rbind(c(1, 0), c(0, 1), c(0.5, 0.5))
      }
    } else {
      p1 <- seq(0, 1, by = step)
      w <- cbind(p1, 1 - p1)
      if (constraint == "ordered") w <- w[w[, 1] >= w[, 2] - 1e-9, , drop = FALSE]
    }
  } else if (m == 3) {
    if (mode == "vertices") {
      if (constraint == "ordered") {
        w <- rbind(c(1, 0, 0), c(0.5, 0.5, 0), c(1/3, 1/3, 1/3))
      } else {
        w <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1/3, 1/3, 1/3))
      }
    } else {
      vals <- seq(0, 1, by = step)
      rows <- list()
      idx <- 1L
      for (p1 in vals) {
        for (p2 in vals) {
          p3 <- 1 - p1 - p2
          if (p3 < -1e-9) next
          if (p3 < 0) p3 <- 0
          ww <- c(p1, p2, p3)
          if (abs(sum(ww) - 1) > 1e-8) next
          if (constraint == "ordered" && !(ww[1] >= ww[2] - 1e-9 && ww[2] >= ww[3] - 1e-9)) next
          rows[[idx]] <- ww
          idx <- idx + 1L
        }
      }
      w <- do.call(rbind, rows)
    }
  } else {
    stop("Only 2 or 3 endpoints are supported.")
  }
  
  w <- unique(round(w, 8))
  colnames(w) <- paste0("p", seq_len(ncol(w)))
  as.data.frame(w)
}

threshold_objects <- make_threshold_grid_df(subject_data, endpoint_specs)
THRESHOLD_GRID_DF <- threshold_objects$grid
THRESHOLD_INFO <- threshold_objects$info
THRESHOLD_ACTIVE_IDS <- threshold_objects$active_ids

WEIGHT_GRID <- make_weight_grid(M_ENDPOINTS, WEIGHT_MODE, WEIGHT_STEP, WEIGHT_CONSTRAINT)
EQUAL_WEIGHTS <- rep(1 / M_ENDPOINTS, M_ENDPOINTS)
NATURAL_ORDER <- seq_len(M_ENDPOINTS)
ALL_ORDERS <- all_permutations(seq_len(M_ENDPOINTS))

fwrite(THRESHOLD_INFO, file.path(OUTPUT_DIR, "REALDATA_threshold_grid_by_endpoint.csv"))
fwrite(THRESHOLD_GRID_DF, file.path(OUTPUT_DIR, "REALDATA_threshold_candidate_combinations.csv"))
fwrite(as.data.table(WEIGHT_GRID), file.path(OUTPUT_DIR, "REALDATA_weight_grid.csv"))

cat("Threshold-active time endpoint(s):", ifelse(length(THRESHOLD_ACTIVE_IDS) == 0, "none", paste(THRESHOLD_ACTIVE_IDS, collapse = ", ")), "\n")
cat("Threshold candidate combinations:", nrow(THRESHOLD_GRID_DF), "\n")
cat("Weight grid:\n")
print(WEIGHT_GRID)

#Candidate grids and method definitions
zero_threshold_grid_df <- function(m) {
  dt <- data.table(dummy = 1L)
  for (j in seq_len(m)) dt[[paste0("threshold_e", j)]] <- 0
  dt[, dummy := NULL]
  dt[, threshold_key := paste(rep("0.0000", m), collapse = "|")]
  dt
}

candidate_key_dt <- function(dt) {
  p_cols <- paste0("p", seq_len(M_ENDPOINTS))
  t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
  p_part <- apply(as.data.frame(dt[, ..p_cols]), 1, function(z) paste(sprintf("%.8f", as.numeric(z)), collapse = "|"))
  t_part <- apply(as.data.frame(dt[, ..t_cols]), 1, function(z) paste(sprintf("%.4f", as.numeric(z)), collapse = "|"))
  paste(dt$order_key, p_part, t_part, sep = "__")
}

make_candidate_grid <- function(orders, weights_df, threshold_df) {
  rows <- list()
  idx <- 1L
  t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
  
  for (ord in orders) {
    ord_key <- format_order(ord)
    for (i in seq_len(nrow(weights_df))) {
      ww <- as.numeric(unlist(weights_df[i, seq_len(M_ENDPOINTS), drop = FALSE], use.names = FALSE))
      for (j in seq_len(nrow(threshold_df))) {
        row <- data.table(order_key = ord_key)
        for (k in seq_len(M_ENDPOINTS)) row[[paste0("p", k)]] <- ww[k]
        for (k in seq_len(M_ENDPOINTS)) row[[paste0("threshold_e", k)]] <- threshold_df[[paste0("threshold_e", k)]][j]
        rows[[idx]] <- row
        idx <- idx + 1L
      }
    }
  }
  
  out <- rbindlist(rows, fill = TRUE)
  out[, candidate_key := candidate_key_dt(out)]
  unique(out, by = "candidate_key")
}

make_methods <- function() {
  methods <- list()
  
  eq_df <- as.data.frame(as.list(EQUAL_WEIGHTS))
  names(eq_df) <- paste0("p", seq_len(M_ENDPOINTS))
  
  zero_thr <- zero_threshold_grid_df(M_ENDPOINTS)
  active_thr <- THRESHOLD_GRID_DF
  
  natural_grid <- make_candidate_grid(list(NATURAL_ORDER), eq_df, zero_thr)
  methods[[length(methods) + 1L]] <- list(
    id = "original",
    notation = ifelse(M_ENDPOINTS == 2, "WR/WO", "WR/WO^{(3)}"),
    description = "Fixed natural order, equal weights, no threshold",
    grid = natural_grid
  )
  
  if (REPORT_ALL_FIXED_ORDERS) {
    for (ord in ALL_ORDERS) {
      if (identical(ord, NATURAL_ORDER)) next
      ord_key <- format_order(ord)
      methods[[length(methods) + 1L]] <- list(
        id = paste0("fixed_order_", gsub("->", "_", ord_key, fixed = TRUE)),
        notation = paste0("Fixed order ", ord_key),
        description = paste0("Fixed order ", ord_key, ", equal weights, no threshold"),
        grid = make_candidate_grid(list(ord), eq_df, zero_thr)
      )
    }
  }
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord", "M_ord^{(3)}"),
    description = "Maximize over endpoint order; equal weights; no threshold",
    grid = make_candidate_grid(ALL_ORDERS, eq_df, zero_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_wt",
    notation = ifelse(M_ENDPOINTS == 2, "M_wt^0.5", "M_wt^{(3)}"),
    description = "Fixed natural order; maximize over ordered weights; no threshold",
    grid = make_candidate_grid(list(NATURAL_ORDER), WEIGHT_GRID, zero_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord,wt^0.5", "M_ord,wt^{(3)}"),
    description = "Maximize over endpoint order and ordered weights; no threshold",
    grid = make_candidate_grid(ALL_ORDERS, WEIGHT_GRID, zero_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_thresh", "M_thresh^{(3)}"),
    description = "Fixed natural order; equal weights; maximize threshold(s) on time-to-event endpoint(s)",
    grid = make_candidate_grid(list(NATURAL_ORDER), eq_df, active_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_wt_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_wt,thresh^0.5", "M_wt,thresh^{(3)}"),
    description = "Fixed natural order; maximize ordered weights and threshold(s) on time-to-event endpoint(s)",
    grid = make_candidate_grid(list(NATURAL_ORDER), WEIGHT_GRID, active_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord,thresh", "M_ord,thresh^{(3)}"),
    description = "Maximize over endpoint order and threshold(s) on time-to-event endpoint(s); equal weights",
    grid = make_candidate_grid(ALL_ORDERS, eq_df, active_thr)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord,wt,thresh^0.5", "M_ord,wt,thresh^{(3)}"),
    description = "Full adaptive: maximize over order, ordered weights, and time-endpoint threshold(s)",
    grid = make_candidate_grid(ALL_ORDERS, WEIGHT_GRID, active_thr)
  )
  
  methods
}

METHODS <- make_methods()

method_overview <- rbindlist(lapply(METHODS, function(m) {
  data.table(
    method_id = m$id,
    notation = m$notation,
    description = m$description,
    n_candidates = nrow(m$grid)
  )
}))
fwrite(method_overview, file.path(OUTPUT_DIR, "REALDATA_method_overview.csv"))

ALL_CANDIDATES <- rbindlist(lapply(METHODS, function(m) m$grid), fill = TRUE)
ALL_CANDIDATES <- unique(ALL_CANDIDATES, by = "candidate_key")
fwrite(ALL_CANDIDATES, file.path(OUTPUT_DIR, "REALDATA_all_unique_candidates.csv"))

cat("Unique candidates evaluated per observed/permutation dataset:", nrow(ALL_CANDIDATES), "\n")


#Candidate evaluation and selection

evaluate_all_candidates <- function(data, arm_vec, candidate_grid) {
  pairs <- get_pair_indices(arm_vec)
  rows <- vector("list", nrow(candidate_grid))
  t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
  
  for (i in seq_len(nrow(candidate_grid))) {
    cg <- candidate_grid[i]
    ord <- parse_order(cg$order_key)
    weights <- as.numeric(unlist(cg[, paste0("p", seq_len(M_ENDPOINTS)), with = FALSE], use.names = FALSE))
    thresholds_by_endpoint <- as.numeric(unlist(cg[, ..t_cols], use.names = FALSE))
    
    cc <- counts_for_candidate(
      data = data,
      arm_vec = arm_vec,
      order_vec = ord,
      weights = weights,
      thresholds_by_endpoint = thresholds_by_endpoint,
      pairs = pairs
    )
    
    win_cols <- as.list(rep(NA_real_, 3))
    loss_cols <- as.list(rep(NA_real_, 3))
    names(win_cols) <- paste0("win_rank", 1:3)
    names(loss_cols) <- paste0("loss_rank", 1:3)
    for (r in seq_len(M_ENDPOINTS)) {
      win_cols[[paste0("win_rank", r)]] <- cc$wins_by_rank[r]
      loss_cols[[paste0("loss_rank", r)]] <- cc$losses_by_rank[r]
    }
    
    rows[[i]] <- c(
      as.list(cg),
      list(
        WR_statistic = cc$WR_statistic,
        WO_statistic = cc$WO_statistic,
        WR_abslog = cc$WR_abslog,
        WO_abslog = cc$WO_abslog,
        weighted_win = cc$weighted_win,
        weighted_loss = cc$weighted_loss,
        tie_count = cc$tie_count,
        tie_proportion = cc$tie_proportion,
        n_pairs = cc$n_pairs
      ),
      win_cols,
      loss_cols
    )
  }
  
  rbindlist(rows, fill = TRUE)
}

select_best_candidate <- function(counts_dt, candidate_keys, measure = c("WR", "WO"), side = c("one", "two")) {
  measure <- match.arg(measure)
  side <- match.arg(side)
  
  d <- counts_dt[candidate_key %in% candidate_keys]
  if (nrow(d) == 0) stop("No candidate rows found for selection.")
  
  stat_col <- if (measure == "WR") "WR_statistic" else "WO_statistic"
  abs_col <- if (measure == "WR") "WR_abslog" else "WO_abslog"
  
  if (side == "one") {
    ord <- order(d[[stat_col]], decreasing = TRUE, na.last = TRUE)
  } else {
    ord <- order(d[[abs_col]], decreasing = TRUE, na.last = TRUE)
  }
  
  d[ord[1]]
}

get_stat_value <- function(row, measure = c("WR", "WO"), side = c("one", "two")) {
  measure <- match.arg(measure)
  side <- match.arg(side)
  if (measure == "WR" && side == "one") return(as.numeric(row$WR_statistic[1]))
  if (measure == "WO" && side == "one") return(as.numeric(row$WO_statistic[1]))
  if (measure == "WR" && side == "two") return(as.numeric(row$WR_abslog[1]))
  if (measure == "WO" && side == "two") return(as.numeric(row$WO_abslog[1]))
}

threshold_summary_string <- function(row) {
  t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
  vals <- as.numeric(unlist(row[, ..t_cols], use.names = FALSE))
  paste(paste0("e", seq_len(M_ENDPOINTS), "=", sprintf("%.4f", vals)), collapse = ";")
}

flatten_selection_for_perm <- function(row, measure, method_id, side) {
  t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
  vals <- as.numeric(unlist(row[, ..t_cols], use.names = FALSE))
  out <- data.table(
    measure = measure,
    method_id = method_id,
    side = side,
    statistic = if (measure == "WR") row$WR_statistic else row$WO_statistic,
    abslog_statistic = if (measure == "WR") row$WR_abslog else row$WO_abslog,
    selected_order = row$order_key,
    selected_p1 = row$p1,
    selected_p2 = row$p2,
    selected_p3 = if (M_ENDPOINTS == 3) row$p3 else NA_real_,
    selected_thresholds = paste(vals, collapse = ","),
    tie_count = row$tie_count,
    tie_proportion = row$tie_proportion
  )
  for (j in seq_len(M_ENDPOINTS)) out[[paste0("selected_t_e", j)]] <- vals[j]
  out
}


#One permutation process for all WR/WO methods
run_all_methods_single_permutation_process <- function() {
  measures_to_run <- c(if (RUN_WR) "WR", if (RUN_WO) "WO")
  if (length(measures_to_run) == 0) stop("At least one of RUN_WR or RUN_WO must be TRUE.")
  
  final_path <- file.path(OUTPUT_DIR, "checkpoints", "ALL_METHODS_SINGLE_PERMUTATION_final.rds")
  checkpoint_path <- file.path(OUTPUT_DIR, "checkpoints", "ALL_METHODS_SINGLE_PERMUTATION_checkpoint.rds")
  
  if (file.exists(final_path) && !FORCE_RERUN) {
    return(readRDS(final_path))
  }
  
  cat("\nEvaluating observed data for all unique candidates...\n")
  obs_counts <- evaluate_all_candidates(subject_data, subject_data$arm, ALL_CANDIDATES)
  fwrite(obs_counts, file.path(OUTPUT_DIR, "REALDATA_observed_candidates_all_unique.csv"))
  
  combo_rows <- list()
  combo_idx <- 1L
  for (measure in measures_to_run) {
    for (m in METHODS) {
      keys <- m$grid$candidate_key
      obs_one <- select_best_candidate(obs_counts, keys, measure, "one")
      obs_two <- select_best_candidate(obs_counts, keys, measure, "two")
      
      combo_rows[[combo_idx]] <- data.table(
        combo_id = paste(measure, m$id, sep = "__"),
        measure = measure,
        method_id = m$id,
        notation = m$notation,
        description = m$description,
        n_candidates = nrow(m$grid),
        obs_one_stat = get_stat_value(obs_one, measure, "one"),
        obs_two_abslog = get_stat_value(obs_two, measure, "two"),
        obs_one_candidate = obs_one$candidate_key,
        obs_two_candidate = obs_two$candidate_key
      )
      combo_idx <- combo_idx + 1L
    }
  }
  combo_dt <- rbindlist(combo_rows, fill = TRUE)
  
  exceed_one <- setNames(rep(0L, nrow(combo_dt)), combo_dt$combo_id)
  exceed_two <- setNames(rep(0L, nrow(combo_dt)), combo_dt$combo_id)
  perm_selected <- data.table()
  start_b <- 1L
  
  if (file.exists(checkpoint_path) && !FORCE_RERUN) {
    chk <- readRDS(checkpoint_path)
    start_b <- chk$next_b
    exceed_one <- chk$exceed_one
    exceed_two <- chk$exceed_two
    perm_selected <- chk$perm_selected
    cat("  Resuming from permutation", start_b, "\n")
  }
  
  PERM_ARM_MATRIX <- replicate(B_PERM, sample(subject_data$arm), simplify = "matrix")
  
  if (start_b <= B_PERM) {
    for (b in start_b:B_PERM) {
      perm_arm <- PERM_ARM_MATRIX[, b]
      perm_counts <- evaluate_all_candidates(subject_data, perm_arm, ALL_CANDIDATES)
      
      perm_rows_b <- list()
      rr <- 1L
      for (measure in measures_to_run) {
        for (m in METHODS) {
          keys <- m$grid$candidate_key
          combo_id_value <- paste(measure, m$id, sep = "__")
          obs_row <- combo_dt[combo_id == combo_id_value]
          
          perm_one <- select_best_candidate(perm_counts, keys, measure, "one")
          perm_two <- select_best_candidate(perm_counts, keys, measure, "two")
          
          perm_one_stat <- get_stat_value(perm_one, measure, "one")
          perm_two_abslog <- get_stat_value(perm_two, measure, "two")
          
          if (!is.na(perm_one_stat) && !is.na(obs_row$obs_one_stat) && perm_one_stat >= obs_row$obs_one_stat) {
            exceed_one[combo_id_value] <- exceed_one[combo_id_value] + 1L
          }
          if (!is.na(perm_two_abslog) && !is.na(obs_row$obs_two_abslog) && perm_two_abslog >= obs_row$obs_two_abslog) {
            exceed_two[combo_id_value] <- exceed_two[combo_id_value] + 1L
          }
          
          one_flat <- flatten_selection_for_perm(perm_one, measure, m$id, "one")
          two_flat <- flatten_selection_for_perm(perm_two, measure, m$id, "two")
          one_flat[, b := b]
          two_flat[, b := b]
          perm_rows_b[[rr]] <- one_flat; rr <- rr + 1L
          perm_rows_b[[rr]] <- two_flat; rr <- rr + 1L
        }
      }
      
      perm_selected <- rbind(perm_selected, rbindlist(perm_rows_b, fill = TRUE), fill = TRUE)
      
      if (b %% CHECKPOINT_EVERY == 0 || b == B_PERM) {
        saveRDS(
          list(
            next_b = b + 1L,
            exceed_one = exceed_one,
            exceed_two = exceed_two,
            perm_selected = perm_selected,
            combo_dt = combo_dt
          ),
          checkpoint_path
        )
        cat("  completed permutation", b, "of", B_PERM, "\n")
      }
    }
  }
  
  result <- list(
    obs_counts = obs_counts,
    combo_dt = combo_dt,
    exceed_one = exceed_one,
    exceed_two = exceed_two,
    perm_selected = perm_selected,
    B = B_PERM,
    methods = METHODS,
    all_candidates = ALL_CANDIDATES
  )
  
  saveRDS(result, final_path)
  result
}

analysis_out <- run_all_methods_single_permutation_process()


#final result tables
get_obs_selection <- function(obs_counts, key_value) {
  obs_counts[candidate_key == key_value][1]
}

flatten_final_result <- function(analysis_out) {
  combo_dt <- copy(analysis_out$combo_dt)
  obs_counts <- analysis_out$obs_counts
  perm_selected <- analysis_out$perm_selected
  B <- analysis_out$B
  
  rows <- list()
  for (i in seq_len(nrow(combo_dt))) {
    cd <- combo_dt[i]
    combo_id <- cd$combo_id
    
    obs_one <- get_obs_selection(obs_counts, cd$obs_one_candidate)
    obs_two <- get_obs_selection(obs_counts, cd$obs_two_candidate)
    
    perm_one <- perm_selected[measure == cd$measure & method_id == cd$method_id & side == "one"]
    perm_two <- perm_selected[measure == cd$measure & method_id == cd$method_id & side == "two"]
    
    p_one <- (as.integer(analysis_out$exceed_one[combo_id]) + 1) / (B + 1)
    p_two <- (as.integer(analysis_out$exceed_two[combo_id]) + 1) / (B + 1)
    
    t_cols <- paste0("threshold_e", seq_len(M_ENDPOINTS))
    t_one_vals <- as.numeric(unlist(obs_one[, ..t_cols], use.names = FALSE))
    t_two_vals <- as.numeric(unlist(obs_two[, ..t_cols], use.names = FALSE))
    
    row <- data.table(
      measure = cd$measure,
      method_id = cd$method_id,
      notation = cd$notation,
      description = cd$description,
      n_candidates = cd$n_candidates,
      B_perm = B,
      
      statistic_one = if (cd$measure == "WR") obs_one$WR_statistic else obs_one$WO_statistic,
      p_one = p_one,
      reject_one_0.05 = p_one < 0.05,
      exceed_one = as.integer(analysis_out$exceed_one[combo_id]),
      selected_order_one = obs_one$order_key,
      selected_p1_one = obs_one$p1,
      selected_p2_one = obs_one$p2,
      selected_p3_one = if (M_ENDPOINTS == 3) obs_one$p3 else NA_real_,
      selected_thresholds_one = paste(t_one_vals, collapse = ","),
      weighted_win_one = obs_one$weighted_win,
      weighted_loss_one = obs_one$weighted_loss,
      tie_count_one = obs_one$tie_count,
      tie_proportion_one = obs_one$tie_proportion,
      win_rank1_one = obs_one$win_rank1,
      win_rank2_one = obs_one$win_rank2,
      win_rank3_one = obs_one$win_rank3,
      loss_rank1_one = obs_one$loss_rank1,
      loss_rank2_one = obs_one$loss_rank2,
      loss_rank3_one = obs_one$loss_rank3,
      mean_perm_tie_count_one = mean(perm_one$tie_count, na.rm = TRUE),
      mean_perm_tie_proportion_one = mean(perm_one$tie_proportion, na.rm = TRUE),
      
      statistic_two = if (cd$measure == "WR") obs_two$WR_statistic else obs_two$WO_statistic,
      abslog_statistic_two = if (cd$measure == "WR") obs_two$WR_abslog else obs_two$WO_abslog,
      p_two = p_two,
      reject_two_0.05 = p_two < 0.05,
      exceed_two = as.integer(analysis_out$exceed_two[combo_id]),
      selected_order_two = obs_two$order_key,
      selected_p1_two = obs_two$p1,
      selected_p2_two = obs_two$p2,
      selected_p3_two = if (M_ENDPOINTS == 3) obs_two$p3 else NA_real_,
      selected_thresholds_two = paste(t_two_vals, collapse = ","),
      weighted_win_two = obs_two$weighted_win,
      weighted_loss_two = obs_two$weighted_loss,
      tie_count_two = obs_two$tie_count,
      tie_proportion_two = obs_two$tie_proportion,
      win_rank1_two = obs_two$win_rank1,
      win_rank2_two = obs_two$win_rank2,
      win_rank3_two = obs_two$win_rank3,
      loss_rank1_two = obs_two$loss_rank1,
      loss_rank2_two = obs_two$loss_rank2,
      loss_rank3_two = obs_two$loss_rank3,
      mean_perm_tie_count_two = mean(perm_two$tie_count, na.rm = TRUE),
      mean_perm_tie_proportion_two = mean(perm_two$tie_proportion, na.rm = TRUE)
    )
    
    for (j in seq_len(M_ENDPOINTS)) {
      row[[paste0("selected_t_e", j, "_one")]] <- t_one_vals[j]
      row[[paste0("selected_t_e", j, "_two")]] <- t_two_vals[j]
    }
    
    rows[[i]] <- row
  }
  
  rbindlist(rows, fill = TRUE)
}

method_results <- flatten_final_result(analysis_out)
fwrite(method_results, file.path(OUTPUT_DIR, "REALDATA_all_WR_WO_method_results.csv"))
fwrite(analysis_out$perm_selected, file.path(OUTPUT_DIR, "REALDATA_permutation_selected_results_all_methods.csv"))

#Endpoint and data summaries
endpoint_summary_rows <- list()
for (j in seq_along(endpoint_specs)) {
  ep <- endpoint_specs[[j]]
  if (ep$type == "time") {
    endpoint_summary_rows[[j]] <- data.table(
      endpoint = j,
      name = ep$name,
      type = ep$type,
      n_events_A = sum(subject_data$arm == 1 & subject_data[[ep$event_col]] == 1, na.rm = TRUE),
      n_events_B = sum(subject_data$arm == 0 & subject_data[[ep$event_col]] == 1, na.rm = TRUE),
      mean_A = mean(subject_data[[ep$time_col]][subject_data$arm == 1], na.rm = TRUE),
      mean_B = mean(subject_data[[ep$time_col]][subject_data$arm == 0], na.rm = TRUE),
      median_A = median(subject_data[[ep$time_col]][subject_data$arm == 1], na.rm = TRUE),
      median_B = median(subject_data[[ep$time_col]][subject_data$arm == 0], na.rm = TRUE)
    )
  } else if (ep$type == "count") {
    endpoint_summary_rows[[j]] <- data.table(
      endpoint = j,
      name = ep$name,
      type = ep$type,
      n_events_A = sum(subject_data[[ep$count_col]][subject_data$arm == 1] > 0, na.rm = TRUE),
      n_events_B = sum(subject_data[[ep$count_col]][subject_data$arm == 0] > 0, na.rm = TRUE),
      mean_A = mean(subject_data[[ep$count_col]][subject_data$arm == 1], na.rm = TRUE),
      mean_B = mean(subject_data[[ep$count_col]][subject_data$arm == 0], na.rm = TRUE),
      median_A = median(subject_data[[ep$count_col]][subject_data$arm == 1], na.rm = TRUE),
      median_B = median(subject_data[[ep$count_col]][subject_data$arm == 0], na.rm = TRUE)
    )
  } else {
    vcol <- ep$value_col
    endpoint_summary_rows[[j]] <- data.table(
      endpoint = j,
      name = ep$name,
      type = ep$type,
      n_events_A = NA_integer_,
      n_events_B = NA_integer_,
      mean_A = mean(subject_data[[vcol]][subject_data$arm == 1], na.rm = TRUE),
      mean_B = mean(subject_data[[vcol]][subject_data$arm == 0], na.rm = TRUE),
      median_A = median(subject_data[[vcol]][subject_data$arm == 1], na.rm = TRUE),
      median_B = median(subject_data[[vcol]][subject_data$arm == 0], na.rm = TRUE)
    )
  }
}
endpoint_summary <- rbindlist(endpoint_summary_rows, fill = TRUE)
fwrite(endpoint_summary, file.path(OUTPUT_DIR, "REALDATA_endpoint_summary.csv"))

subject_summary <- data.table(
  n_subjects = nrow(subject_data),
  n_treatment_A = sum(subject_data$arm == 1),
  n_control_B = sum(subject_data$arm == 0),
  n_endpoints = M_ENDPOINTS,
  threshold_time_endpoint_ids = paste(THRESHOLD_ACTIVE_IDS, collapse = ","),
  n_threshold_candidate_combinations = nrow(THRESHOLD_GRID_DF),
  B_perm = B_PERM,
  seed = SEED,
  one_permutation_process_for_all_methods = TRUE,
  weight_constraint = WEIGHT_CONSTRAINT,
  weight_mode = WEIGHT_MODE
)
fwrite(subject_summary, file.path(OUTPUT_DIR, "REALDATA_subject_summary.csv"))


#Log-rank tests
run_logrank_one_endpoint <- function(data, ep, label) {
  if (ep$type != "time") {
    return(data.table(test = label, available = FALSE, reason = "endpoint is not time-to-event"))
  }
  
  surv_obj <- survival::Surv(time = data[[ep$time_col]], event = data[[ep$event_col]] == 1)
  fit <- survival::survdiff(surv_obj ~ data$arm)
  p_two <- 1 - stats::pchisq(fit$chisq, df = 1)
  
  cox <- tryCatch(survival::coxph(surv_obj ~ data$arm), error = function(e) NULL)
  if (!is.null(cox)) {
    sm <- summary(cox)
    beta <- sm$coef[1, "coef"]
    se <- sm$coef[1, "se(coef)"]
    z <- beta / se
    hr <- exp(beta)
    #Alternative is treatment benefit for adverse event: HR < 1, so z < 0.
    p_one_benefit <- stats::pnorm(z)
  } else {
    beta <- se <- z <- hr <- p_one_benefit <- NA_real_
  }
  
  data.table(
    test = label,
    available = TRUE,
    reason = NA_character_,
    chisq = fit$chisq,
    p_two_sided = p_two,
    HR_treatment_vs_control = hr,
    beta = beta,
    se = se,
    z = z,
    p_one_sided_benefit = p_one_benefit
  )
}

make_composite_time_data <- function(data, endpoint_specs) {
  if (!all(sapply(endpoint_specs, function(ep) ep$type == "time"))) return(NULL)
  
  times <- sapply(endpoint_specs, function(ep) safe_num(data[[ep$time_col]]))
  events <- sapply(endpoint_specs, function(ep) as.integer(data[[ep$event_col]] == 1))
  
  comp_time <- rep(NA_real_, nrow(data))
  comp_event <- rep(0L, nrow(data))
  
  for (i in seq_len(nrow(data))) {
    ev_times <- times[i, events[i, ] == 1]
    if (length(ev_times) > 0) {
      comp_time[i] <- min(ev_times, na.rm = TRUE)
      comp_event[i] <- 1L
    } else {
      comp_time[i] <- max(times[i, ], na.rm = TRUE)
      comp_event[i] <- 0L
    }
  }
  
  data.frame(id = data$id, arm = data$arm, comp_time = comp_time, comp_event = comp_event)
}

logrank_results <- data.table()
if (RUN_LOGRANK) {
  logrank_results <- rbind(
    logrank_results,
    run_logrank_one_endpoint(subject_data, endpoint_specs[[1]], "Log-rank outcome 1"),
    fill = TRUE
  )
  
  ## Also report log-rank tests for other time-to-event endpoints, if present.
  time_ids <- which(vapply(endpoint_specs, function(ep) ep$type == "time", logical(1)))
  if (length(time_ids) > 0) {
    for (j in setdiff(time_ids, 1L)) {
      logrank_results <- rbind(
        logrank_results,
        run_logrank_one_endpoint(subject_data, endpoint_specs[[j]], paste0("Log-rank outcome ", j)),
        fill = TRUE
      )
    }
  }
  
  if (RUN_COMPOSITE_LOGRANK_IF_POSSIBLE) {
    comp <- make_composite_time_data(subject_data, endpoint_specs)
    if (!is.null(comp)) {
      comp_ep <- list(type = "time", time_col = "comp_time", event_col = "comp_event")
      logrank_results <- rbind(
        logrank_results,
        run_logrank_one_endpoint(comp, comp_ep, "Log-rank composite time-to-first event"),
        fill = TRUE
      )
    } else {
      logrank_results <- rbind(
        logrank_results,
        data.table(test = "Log-rank composite time-to-first event", available = FALSE,
                   reason = "not all endpoints are time-to-event"),
        fill = TRUE
      )
    }
  }
}
fwrite(logrank_results, file.path(OUTPUT_DIR, "REALDATA_logrank_results.csv"))


#Simple figures
plot_dt_one <- copy(method_results)
plot_dt_one[, p_one_label := ifelse(is.na(p_one), NA_character_, sprintf("%.4f", p_one))]
plot_dt_one[, p_two_label := ifelse(is.na(p_two), NA_character_, sprintf("%.4f", p_two))]
plot_dt_one[, method_display := paste(measure, notation)]

p_plot <- ggplot(plot_dt_one, aes(x = reorder(method_display, p_one), y = p_one, fill = measure)) +
  geom_col(position = position_dodge(width = 0.7)) +
  geom_hline(yintercept = 0.05, linetype = 2) +
  coord_flip() +
  labs(
    title = "Real-data permutation p-values: one-sided benefit test",
    x = NULL,
    y = "Permutation p-value"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(OUTPUT_DIR, "figures", "REALDATA_pvalue_barplot_one_sided.png"), p_plot, width = 10, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "figures", "REALDATA_pvalue_barplot_one_sided.pdf"), p_plot, width = 10, height = 7)

if (RUN_LOGRANK && endpoint_specs[[1]]$type == "time") {
  km_df <- subject_data
  ep1 <- endpoint_specs[[1]]
  km_fit <- survival::survfit(survival::Surv(km_df[[ep1$time_col]], km_df[[ep1$event_col]] == 1) ~ arm, data = km_df)
  png(file.path(OUTPUT_DIR, "figures", "REALDATA_KM_outcome1.png"), width = 1800, height = 1400, res = 200)
  plot(km_fit, col = c("blue", "red"), lwd = 2,
       xlab = paste0("Time (", TIME_UNIT_LABEL, ")"), ylab = "Event-free probability",
       main = "Kaplan-Meier curve for outcome 1")
  legend("bottomleft", legend = c("Control B", "Treatment A"), col = c("blue", "red"), lwd = 2, bty = "n")
  dev.off()
}


#README
readme <- c(
  "Real-data WR / WO full adaptive analysis",
  paste0("Generated: ", Sys.time()),
  "",
  "Important interpretation note:",
  "This is a real-data analysis. It reports observed WR/WO statistics and permutation p-values, not statistical power.",
  "Adaptive procedures repeat the full selection inside each treatment-label permutation.",
  "",
  "Major update in this version:",
  "All WR/WO methods are evaluated inside one shared permutation process. The script does NOT rerun a separate permutation loop for each maxing method.",
  "The number of permutations is still B_PERM; the change is computational structure, not the permutation definition.",
  "",
  "Threshold rule:",
  "Thresholds are applied only to time-to-event endpoints. A threshold follows the endpoint number, not its rank in a selected ordering.",
  "If two endpoints are time-to-event endpoints, the threshold search uses two threshold values, one for each time-to-event endpoint.",
  "Count/binary/continuous endpoints are not thresholded in this script.",
  "",
  "Pairwise hierarchy rule:",
  "Pairs unresolved or tied at one endpoint are passed to the next endpoint in the selected ordering.",
  "Pairs unresolved after all endpoints are counted as ties.",
  "",
  "Weight rule:",
  "Default weight constraint is ordered: the first endpoint in the selected order receives the largest weight, the second receives the second largest weight, and the last receives the smallest weight.",
  "For three endpoints with WEIGHT_MODE = vertices and WEIGHT_CONSTRAINT = ordered, the weights are (1,0,0), (0.5,0.5,0), and (1/3,1/3,1/3).",
  "",
  "Main output files:",
  "REALDATA_all_WR_WO_method_results.csv",
  "REALDATA_permutation_selected_results_all_methods.csv",
  "REALDATA_endpoint_summary.csv",
  "REALDATA_subject_summary.csv",
  "REALDATA_logrank_results.csv",
  "REALDATA_method_overview.csv",
  "REALDATA_threshold_grid_by_endpoint.csv",
  "REALDATA_threshold_candidate_combinations.csv",
  "REALDATA_weight_grid.csv",
  "REALDATA_observed_candidates_all_unique.csv",
  "",
  "WR statistic:",
  "sum_r p_r W_A,r / sum_r p_r W_B,r",
  "",
  "WO statistic:",
  "(sum_r p_r W_A,r + 0.5*T) / (sum_r p_r W_B,r + 0.5*T)",
  "",
  paste0("B_PERM = ", B_PERM),
  paste0("SEED = ", SEED),
  paste0("Threshold-active endpoint IDs = ", paste(THRESHOLD_ACTIVE_IDS, collapse = ",")),
  paste0("Number of unique candidate rules evaluated per permutation = ", nrow(ALL_CANDIDATES))
)
writeLines(readme, file.path(OUTPUT_DIR, "README_REALDATA_full_adaptive_WR_WO.txt"))

#Output
cat("\nDone. Outputs saved in:\n")
cat(OUTPUT_DIR, "\n")
cat("\nMain result table:\n")
cat(file.path(OUTPUT_DIR, "REALDATA_all_WR_WO_method_results.csv"), "\n")
cat("\nPermutation selections:\n")
cat(file.path(OUTPUT_DIR, "REALDATA_permutation_selected_results_all_methods.csv"), "\n")
cat("\nLog-rank table:\n")
cat(file.path(OUTPUT_DIR, "REALDATA_logrank_results.csv"), "\n")

