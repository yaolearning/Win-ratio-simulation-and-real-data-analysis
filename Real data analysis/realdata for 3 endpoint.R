#clean environment
rm(list = ls())
options(stringsAsFactors = FALSE)


#Input and output paths
DATA_PATH <- ""
DATA_OBJECT_NAME <- NULL

OUTPUT_DIR <- ""


#Data format
#Use "wide" if each subject has one row.
#Use "event_long" if the data has one row per event, e.g. ID/time/status/trt.
DATA_FORMAT <- "wide"  # "wide" or "event_long"


#Required columns for wide format:
#   ID_COL, ARM_COL, and endpoint columns in ENDPOINTS_WIDE.
# ARM_COL must be coded as 1=treatment A and 0=control B.
ID_COL <- "SUBJID"
ARM_COL <- "ARM"

#Endpoint types supported:
#type = "time"       adverse time-to-event endpoint; larger event-free time is better
#needs time_col and event_col, event_col coded 1=event, 0=censoring
#type = "count"      recurrent/count adverse endpoint; lower count is better
#needs count_col
#type = "binary"     adverse binary endpoint; lower value is better
#needs value_col, e.g. 1=event, 0=no event
#type = "continuous" continuous endpoint; set higher_better TRUE/FALSE
#needs value_col

#For 3 endpoints, put three entries here.
#For 2 endpoints, put two entries here; the code still works.
ENDPOINTS_WIDE <- list(
  list(name = "Outcome 1", type = "time",  time_col = "FUTIME", event_col = "CNSR"),
  list(name = "Outcome 2", type = "count", count_col = "NUMHOSP"),
  ## Replace this with your real third endpoint. Delete it for two-endpoint data.
  list(name = "Outcome 3", type = "count", count_col = "ENDPOINT3_COUNT")
)


#event long format 
#Required columns for event_long format:
#EVENT_ID_COL, EVENT_ARM_COL, EVENT_TIME_COL, EVENT_STATUS_COL.
#status code 0 is treated as censoring/follow-up row.
EVENT_ID_COL <- "patid"
EVENT_ARM_COL <- "trt_ab"
EVENT_TIME_COL <- "time"
EVENT_STATUS_COL <- "status"

#For event_long format, define endpoints by status_code.
#kind = "time"  : first event of that status code; otherwise censored at max follow-up
#kind = "count" : number of events with that status code
EVENT_ENDPOINTS <- list(
  list(name = "Death", type = "time",  status_code = 1),
  list(name = "Hospitalization", type = "count", status_code = 2),
  #Replace or delete depending on your data
  list(name = "Endpoint 3", type = "count", status_code = 3)
)


#Threshold configuration
THRESHOLD_TIME_ENDPOINT <- NULL

#Set THRESHOLD_GRID manually if desired, e.g. c(0, 3, 6, 12, 18, 24).
#Units must match the data time unit. If NULL, the script builds a compact data-driven grid
#from the selected time-to-event endpoint only.
THRESHOLD_GRID <- NULL
TIME_UNIT_LABEL <- "months"
MAX_AUTO_THRESHOLDS <- 10L


#Weight search configuration

#WEIGHT_MODE = "vertices" uses theoretical vertices of the ordered weight simplex.
#For 3 endpoints with p1 >= p2 >= p3 >= 0 and sum=1, vertices are:
#(1,0,0), (0.5,0.5,0), (1/3,1/3,1/3)
#WEIGHT_MODE = "grid" uses a grid over p1,p2,p3.
WEIGHT_MODE <- "vertices"  # "vertices" or "grid"
WEIGHT_STEP <- 0.05
WEIGHT_CONSTRAINT <- "ordered"  # "ordered" or "simplex"

#Permutation configuration
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


#packages
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


# helper function

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

clean_method_file_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

make_p_string <- function(w) {
  paste0("(", paste(sprintf("%.3f", w), collapse = ","), ")")
}


# read data and prepare data
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
    ep_type <- ep$type
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
      stop("Unsupported endpoint type for endpoint ", j, ": ", ep_type)
    }
  }
  
  d <- d[!is.na(d$id) & !is.na(d$arm), , drop = FALSE]
  if (any(duplicated(d$id))) {
    stop("Wide data must have one row per subject. Duplicated IDs found.")
  }
  
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
  
  if (!all(subj$arm %in% c(0, 1))) {
    stop("EVENT_ARM_COL must be coded as 1=treatment A and 0=control B.")
  }
  
  endpoint_specs <- list()
  
  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]
    ep_type <- ep$type
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

resolve_threshold_endpoint <- function(endpoint_specs, requested_endpoint = NULL) {
  time_endpoint_ids <- which(sapply(endpoint_specs, function(ep) identical(tolower(ep$type), "time")))
  
  if (length(time_endpoint_ids) == 0) {
    warning("No time-to-event endpoint was found. Threshold methods will use threshold 0 only.")
    return(NA_integer_)
  }
  
  if (!is.null(requested_endpoint) && !is.na(requested_endpoint)) {
    requested_endpoint <- as.integer(requested_endpoint)
    if (requested_endpoint < 1 || requested_endpoint > length(endpoint_specs)) {
      stop("THRESHOLD_TIME_ENDPOINT is outside the endpoint range.")
    }
    if (!identical(tolower(endpoint_specs[[requested_endpoint]]$type), "time")) {
      stop(
        "THRESHOLD_TIME_ENDPOINT must point to a time-to-event endpoint. ",
        "Endpoint ", requested_endpoint, " has type '", endpoint_specs[[requested_endpoint]]$type, "'. ",
        "Thresholds are not applied to count/binary/continuous endpoints."
      )
    }
    return(requested_endpoint)
  }
  
  if (length(time_endpoint_ids) > 1) {
    warning(
      "Multiple time-to-event endpoints were found. Using endpoint ", time_endpoint_ids[1],
      " by default. Set THRESHOLD_TIME_ENDPOINT to choose a different time endpoint."
    )
  }
  
  time_endpoint_ids[1]
}

THRESHOLD_ENDPOINT <- resolve_threshold_endpoint(endpoint_specs, THRESHOLD_TIME_ENDPOINT)
THRESHOLD_ENDPOINT_NAME <- if (is.na(THRESHOLD_ENDPOINT)) NA_character_ else endpoint_specs[[THRESHOLD_ENDPOINT]]$name
THRESHOLD_ENDPOINT_TYPE <- if (is.na(THRESHOLD_ENDPOINT)) NA_character_ else endpoint_specs[[THRESHOLD_ENDPOINT]]$type

cat("Prepared subject-level data:\n")
cat("  n subjects:", nrow(subject_data), "\n")
cat("  treatment A:", sum(subject_data$arm == 1), "\n")
cat("  control B:", sum(subject_data$arm == 0), "\n")
cat("  endpoints:", M_ENDPOINTS, "\n")


# endpoint comparison rules
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
    ## adverse count: lower is better. In this script, adaptive thresholds are not
    ## assigned to count endpoints, so this branch normally receives threshold = 0.
    out[ok & (cb - ca > threshold)] <- 1L
    out[ok & (ca - cb > threshold)] <- -1L
    return(out)
  }
  
  if (endpoint_spec$type == "binary") {
    va <- safe_num(data[[endpoint_spec$value_col]][idx_A])
    vb <- safe_num(data[[endpoint_spec$value_col]][idx_B])
    ok <- !is.na(va) & !is.na(vb)
    ## adverse binary: lower is better
    out[ok & (vb - va > threshold)] <- 1L
    out[ok & (va - vb > threshold)] <- -1L
    return(out)
  }
  
  if (endpoint_spec$type == "continuous") {
    va <- safe_num(data[[endpoint_spec$value_col]][idx_A])
    vb <- safe_num(data[[endpoint_spec$value_col]][idx_B])
    ok <- !is.na(va) & !is.na(vb)
    higher_better <- isTRUE(endpoint_spec$higher_better)
    if (higher_better) {
      out[ok & (va - vb > threshold)] <- 1L
      out[ok & (vb - va > threshold)] <- -1L
    } else {
      out[ok & (vb - va > threshold)] <- 1L
      out[ok & (va - vb > threshold)] <- -1L
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

counts_for_candidate <- function(data, arm_vec, order_vec, weights, threshold_value,
                                 threshold_endpoint, measure = c("WR", "WO")) {
  measure <- match.arg(measure)
  pairs <- get_pair_indices(arm_vec)
  idx_A <- pairs$A
  idx_B <- pairs$B
  n_pairs <- pairs$n_pairs
  
  wins_by_rank <- rep(0, M_ENDPOINTS)
  losses_by_rank <- rep(0, M_ENDPOINTS)
  unresolved <- rep(TRUE, n_pairs)
  
  for (r in seq_along(order_vec)) {
    endpoint_id <- order_vec[r]
    ep <- endpoint_specs[[endpoint_id]]
    ## The threshold follows the time-to-event endpoint number, not the rank.
    ## It is never applied to count/binary/continuous endpoints.
    thr <- if (!is.na(threshold_endpoint) &&
               endpoint_id == threshold_endpoint &&
               identical(tolower(ep$type), "time")) threshold_value else 0
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
  
  if (measure == "WR") {
    statistic <- ratio_safe(weighted_win, weighted_loss)
  } else {
    statistic <- ratio_safe(weighted_win + 0.5 * tie_count,
                            weighted_loss + 0.5 * tie_count)
  }
  
  out <- list(
    statistic = statistic,
    abslog_statistic = abslog_safe(statistic),
    weighted_win = weighted_win,
    weighted_loss = weighted_loss,
    tie_count = tie_count,
    tie_proportion = tie_count / n_pairs,
    n_pairs = n_pairs,
    wins_by_rank = wins_by_rank,
    losses_by_rank = losses_by_rank
  )
  out
}

#threshold and weight grids
auto_threshold_grid <- function(data, endpoint_spec, max_thresholds = 10L) {
  ## Thresholds are defined only for time-to-event endpoints.
  ## Count, binary, and continuous endpoints are not thresholded.
  if (is.null(endpoint_spec) || !identical(tolower(endpoint_spec$type), "time")) {
    return(0)
  }
  
  vals <- c(0)
  t <- safe_num(data[[endpoint_spec$time_col]])
  e <- as.integer(data[[endpoint_spec$event_col]] == 1)
  t_event <- t[e == 1 & !is.na(t)]
  t_all <- t[!is.na(t)]
  
  if (length(t_all) == 0 || all(!is.finite(t_all))) {
    return(0)
  }
  
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
      if (constraint == "ordered") w <- w[w[, 1] >= w[, 2], , drop = FALSE]
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

if (is.na(THRESHOLD_ENDPOINT)) {
  THRESHOLD_GRID <- 0
} else if (is.null(THRESHOLD_GRID)) {
  THRESHOLD_GRID <- auto_threshold_grid(subject_data, endpoint_specs[[THRESHOLD_ENDPOINT]], MAX_AUTO_THRESHOLDS)
}
THRESHOLD_GRID <- sort(unique(as.numeric(THRESHOLD_GRID)))
if (!0 %in% THRESHOLD_GRID) THRESHOLD_GRID <- sort(unique(c(0, THRESHOLD_GRID)))
if (length(THRESHOLD_GRID) == 0 || all(!is.finite(THRESHOLD_GRID))) THRESHOLD_GRID <- 0

WEIGHT_GRID <- make_weight_grid(M_ENDPOINTS, WEIGHT_MODE, WEIGHT_STEP, WEIGHT_CONSTRAINT)
EQUAL_WEIGHTS <- rep(1 / M_ENDPOINTS, M_ENDPOINTS)
NATURAL_ORDER <- seq_len(M_ENDPOINTS)
ALL_ORDERS <- all_permutations(seq_len(M_ENDPOINTS))

cat("Threshold endpoint:", THRESHOLD_ENDPOINT, "(", THRESHOLD_ENDPOINT_NAME, ",", THRESHOLD_ENDPOINT_TYPE, ")\n")
cat("Threshold grid:", paste(THRESHOLD_GRID, collapse = ", "), "\n")
cat("Weight grid:\n")
print(WEIGHT_GRID)

# method definition
make_candidate_grid <- function(orders, weights_df, thresholds) {
  rows <- list()
  idx <- 1L
  for (ord in orders) {
    ord_key <- format_order(ord)
    for (i in seq_len(nrow(weights_df))) {
      ww <- as.numeric(weights_df[i, seq_len(M_ENDPOINTS), drop = TRUE])
      for (tt in thresholds) {
        rows[[idx]] <- data.frame(
          order_key = ord_key,
          p1 = ww[1],
          p2 = ifelse(M_ENDPOINTS >= 2, ww[2], NA_real_),
          p3 = ifelse(M_ENDPOINTS >= 3, ww[3], NA_real_),
          threshold = tt,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  rbindlist(rows, fill = TRUE)
}

make_methods <- function() {
  methods <- list()
  
  eq_df <- as.data.frame(as.list(EQUAL_WEIGHTS))
  names(eq_df) <- paste0("p", seq_len(M_ENDPOINTS))
  
  natural_grid <- make_candidate_grid(list(NATURAL_ORDER), eq_df, 0)
  methods[[length(methods) + 1L]] <- list(
    id = "original",
    notation = ifelse(M_ENDPOINTS == 2, "WR/WO", "WR/WO_3"),
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
        grid = make_candidate_grid(list(ord), eq_df, 0)
      )
    }
  }
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord",
    notation = "M_ord",
    description = "Maximize over endpoint order; equal weights; no threshold",
    grid = make_candidate_grid(ALL_ORDERS, eq_df, 0)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_wt",
    notation = ifelse(M_ENDPOINTS == 2, "M_wt^0.5", "M_wt^(3)"),
    description = "Fixed natural order; maximize over weights; no threshold",
    grid = make_candidate_grid(list(NATURAL_ORDER), WEIGHT_GRID, 0)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord,wt^0.5", "M_ord,wt^(3)"),
    description = "Maximize over endpoint order and weights; no threshold",
    grid = make_candidate_grid(ALL_ORDERS, WEIGHT_GRID, 0)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_thresh",
    notation = "M_thresh",
    description = "Fixed natural order; equal weights; maximize threshold on the time-to-event endpoint",
    grid = make_candidate_grid(list(NATURAL_ORDER), eq_df, THRESHOLD_GRID)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_wt_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_wt,thresh^0.5", "M_wt,thresh^(3)"),
    description = "Fixed natural order; maximize over weights and threshold on the time-to-event endpoint",
    grid = make_candidate_grid(list(NATURAL_ORDER), WEIGHT_GRID, THRESHOLD_GRID)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_thresh",
    notation = "M_ord,thresh",
    description = "Maximize over endpoint order and threshold on the time-to-event endpoint; equal weights",
    grid = make_candidate_grid(ALL_ORDERS, eq_df, THRESHOLD_GRID)
  )
  
  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt_thresh",
    notation = ifelse(M_ENDPOINTS == 2, "M_ord,wt,thresh^0.5", "M_ord,wt,thresh^(3)"),
    description = "Full adaptive: maximize over order, weights, and threshold on the time-to-event endpoint",
    grid = make_candidate_grid(ALL_ORDERS, WEIGHT_GRID, THRESHOLD_GRID)
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


#candidate evulation
evaluate_candidate_table <- function(data, arm_vec, candidate_grid, measure) {
  rows <- vector("list", nrow(candidate_grid))
  
  for (i in seq_len(nrow(candidate_grid))) {
    cg <- candidate_grid[i]
    ord <- parse_order(cg$order_key)
    weights <- c(cg$p1, cg$p2)
    if (M_ENDPOINTS == 3) weights <- c(cg$p1, cg$p2, cg$p3)
    
    cc <- counts_for_candidate(
      data = data,
      arm_vec = arm_vec,
      order_vec = ord,
      weights = weights,
      threshold_value = cg$threshold,
      threshold_endpoint = THRESHOLD_ENDPOINT,
      measure = measure
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
      list(
        order_key = cg$order_key,
        p1 = cg$p1,
        p2 = cg$p2,
        p3 = ifelse(M_ENDPOINTS == 3, cg$p3, NA_real_),
        threshold_endpoint = THRESHOLD_ENDPOINT,
        threshold = cg$threshold,
        statistic = cc$statistic,
        abslog_statistic = cc$abslog_statistic,
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
  
  as.data.table(rows)
}

select_best_candidate <- function(candidate_table, side = c("one", "two")) {
  side <- match.arg(side)
  d <- copy(candidate_table)
  if (side == "one") {
    ord <- order(d$statistic, decreasing = TRUE, na.last = TRUE)
  } else {
    ord <- order(d$abslog_statistic, decreasing = TRUE, na.last = TRUE)
  }
  d[ord[1]]
}

#permutation
PERM_ARM_MATRIX <- replicate(B_PERM, sample(subject_data$arm), simplify = "matrix")

run_one_method_measure <- function(method_obj, measure) {
  method_id <- method_obj$id
  clean_id <- clean_method_file_name(paste(measure, method_id, sep = "_"))
  checkpoint_path <- file.path(OUTPUT_DIR, "checkpoints", paste0(clean_id, "_checkpoint.rds"))
  final_path <- file.path(OUTPUT_DIR, "checkpoints", paste0(clean_id, "_final.rds"))
  
  if (file.exists(final_path) && !FORCE_RERUN) {
    return(readRDS(final_path))
  }
  
  cat("\nRunning", measure, method_id, "with", nrow(method_obj$grid), "candidates\n")
  
  obs_table <- evaluate_candidate_table(subject_data, subject_data$arm, method_obj$grid, measure)
  obs_one <- select_best_candidate(obs_table, "one")
  obs_two <- select_best_candidate(obs_table, "two")
  
  start_b <- 1L
  perm_results <- data.table()
  exceed_one <- 0L
  exceed_two <- 0L
  
  if (file.exists(checkpoint_path) && !FORCE_RERUN) {
    chk <- readRDS(checkpoint_path)
    start_b <- chk$next_b
    perm_results <- chk$perm_results
    exceed_one <- chk$exceed_one
    exceed_two <- chk$exceed_two
    cat("  Resuming from permutation", start_b, "\n")
  }
  
  if (start_b <= B_PERM) {
    for (b in start_b:B_PERM) {
      perm_arm <- PERM_ARM_MATRIX[, b]
      perm_table <- evaluate_candidate_table(subject_data, perm_arm, method_obj$grid, measure)
      perm_one <- select_best_candidate(perm_table, "one")
      perm_two <- select_best_candidate(perm_table, "two")
      
      if (!is.na(perm_one$statistic) && !is.na(obs_one$statistic) && perm_one$statistic >= obs_one$statistic) {
        exceed_one <- exceed_one + 1L
      }
      if (!is.na(perm_two$abslog_statistic) && !is.na(obs_two$abslog_statistic) && perm_two$abslog_statistic >= obs_two$abslog_statistic) {
        exceed_two <- exceed_two + 1L
      }
      
      perm_results <- rbind(
        perm_results,
        data.table(
          b = b,
          statistic_one = perm_one$statistic,
          statistic_two = perm_two$statistic,
          abslog_two = perm_two$abslog_statistic,
          tie_count_one = perm_one$tie_count,
          tie_count_two = perm_two$tie_count,
          tie_proportion_one = perm_one$tie_proportion,
          tie_proportion_two = perm_two$tie_proportion,
          selected_order_one = perm_one$order_key,
          selected_order_two = perm_two$order_key,
          selected_p1_one = perm_one$p1,
          selected_p2_one = perm_one$p2,
          selected_p3_one = perm_one$p3,
          selected_p1_two = perm_two$p1,
          selected_p2_two = perm_two$p2,
          selected_p3_two = perm_two$p3,
          selected_t_one = perm_one$threshold,
          selected_t_two = perm_two$threshold
        ),
        fill = TRUE
      )
      
      if (b %% CHECKPOINT_EVERY == 0 || b == B_PERM) {
        saveRDS(
          list(
            next_b = b + 1L,
            perm_results = perm_results,
            exceed_one = exceed_one,
            exceed_two = exceed_two
          ),
          checkpoint_path
        )
        cat("  completed permutation", b, "of", B_PERM, "\n")
      }
    }
  }
  
  p_one <- (exceed_one + 1) / (B_PERM + 1)
  p_two <- (exceed_two + 1) / (B_PERM + 1)
  
  out <- list(
    measure = measure,
    method_id = method_id,
    notation = method_obj$notation,
    description = method_obj$description,
    n_candidates = nrow(method_obj$grid),
    obs_candidates = obs_table,
    obs_one = obs_one,
    obs_two = obs_two,
    perm_results = perm_results,
    p_one = p_one,
    p_two = p_two,
    exceed_one = exceed_one,
    exceed_two = exceed_two,
    B = B_PERM
  )
  
  saveRDS(out, final_path)
  out
}

flatten_result <- function(res) {
  one <- res$obs_one
  two <- res$obs_two
  pr <- res$perm_results
  
  data.table(
    measure = res$measure,
    method_id = res$method_id,
    notation = res$notation,
    description = res$description,
    n_candidates = res$n_candidates,
    B_perm = res$B,
    
    statistic_one = one$statistic,
    p_one = res$p_one,
    reject_one_0.05 = res$p_one < 0.05,
    selected_order_one = one$order_key,
    selected_p1_one = one$p1,
    selected_p2_one = one$p2,
    selected_p3_one = one$p3,
    selected_t_one = one$threshold,
    weighted_win_one = one$weighted_win,
    weighted_loss_one = one$weighted_loss,
    tie_count_one = one$tie_count,
    tie_proportion_one = one$tie_proportion,
    win_rank1_one = one$win_rank1,
    win_rank2_one = one$win_rank2,
    win_rank3_one = one$win_rank3,
    loss_rank1_one = one$loss_rank1,
    loss_rank2_one = one$loss_rank2,
    loss_rank3_one = one$loss_rank3,
    mean_perm_tie_count_one = mean(pr$tie_count_one, na.rm = TRUE),
    mean_perm_tie_proportion_one = mean(pr$tie_proportion_one, na.rm = TRUE),
    
    statistic_two = two$statistic,
    abslog_statistic_two = two$abslog_statistic,
    p_two = res$p_two,
    reject_two_0.05 = res$p_two < 0.05,
    selected_order_two = two$order_key,
    selected_p1_two = two$p1,
    selected_p2_two = two$p2,
    selected_p3_two = two$p3,
    selected_t_two = two$threshold,
    weighted_win_two = two$weighted_win,
    weighted_loss_two = two$weighted_loss,
    tie_count_two = two$tie_count,
    tie_proportion_two = two$tie_proportion,
    win_rank1_two = two$win_rank1,
    win_rank2_two = two$win_rank2,
    win_rank3_two = two$win_rank3,
    loss_rank1_two = two$loss_rank1,
    loss_rank2_two = two$loss_rank2,
    loss_rank3_two = two$loss_rank3,
    mean_perm_tie_count_two = mean(pr$tie_count_two, na.rm = TRUE),
    mean_perm_tie_proportion_two = mean(pr$tie_proportion_two, na.rm = TRUE)
  )
}


#WR and WO methods
all_results <- list()
idx <- 1L
measures_to_run <- c(if (RUN_WR) "WR", if (RUN_WO) "WO")

for (measure in measures_to_run) {
  for (m in METHODS) {
    rr <- run_one_method_measure(m, measure)
    all_results[[idx]] <- flatten_result(rr)
    ## Save observed candidates for traceability.
    obs_cand <- copy(rr$obs_candidates)
    obs_cand[, measure := measure]
    obs_cand[, method_id := m$id]
    obs_cand[, notation := m$notation]
    fwrite(obs_cand, file.path(OUTPUT_DIR, paste0("REALDATA_observed_candidates_", measure, "_", clean_method_file_name(m$id), ".csv")))
    idx <- idx + 1L
  }
}

method_results <- rbindlist(all_results, fill = TRUE)
fwrite(method_results, file.path(OUTPUT_DIR, "REALDATA_all_WR_WO_method_results.csv"))

# endpoint and data summaries 
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
  threshold_endpoint = THRESHOLD_ENDPOINT,
  threshold_endpoint_name = THRESHOLD_ENDPOINT_NAME,
  threshold_endpoint_type = THRESHOLD_ENDPOINT_TYPE,
  threshold_grid = paste(THRESHOLD_GRID, collapse = ","),
  B_perm = B_PERM,
  seed = SEED
)
fwrite(subject_summary, file.path(OUTPUT_DIR, "REALDATA_subject_summary.csv"))


#log rank test
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
    ## Alternative is treatment benefit for adverse event: HR < 1, so z < 0.
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


#simple figure
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
  "Main output files:",
  "REALDATA_all_WR_WO_method_results.csv",
  "REALDATA_endpoint_summary.csv",
  "REALDATA_subject_summary.csv",
  "REALDATA_logrank_results.csv",
  "REALDATA_method_overview.csv",
  "REALDATA_observed_candidates_[measure]_[method].csv",
  "",
  "Main selected columns in REALDATA_all_WR_WO_method_results.csv:",
  "statistic_one, p_one, selected_order_one, selected_p1_one, selected_p2_one, selected_p3_one, selected_t_one, tie_count_one, tie_proportion_one",
  "statistic_two, p_two, selected_order_two, selected_p1_two, selected_p2_two, selected_p3_two, selected_t_two, tie_count_two, tie_proportion_two",
  "",
  "Endpoint comparison rules:",
  "For adverse time-to-event endpoints, larger event-free time is better; censored pairs are only ordered when one subject is known to have longer event-free follow-up beyond the threshold.",
  "For adverse count/binary endpoints, smaller values are better.",
  "Pairs unresolved after all endpoints are counted as ties.",
  "",
  "WR statistic:",
  "sum_r p_r W_A,r / sum_r p_r W_B,r",
  "",
  "WO statistic:",
  "(sum_r p_r W_A,r + 0.5*T) / (sum_r p_r W_B,r + 0.5*T)",
  "",
  "Threshold rule:",
  "Thresholds are always applied to the selected time-to-event endpoint number, not to the first rank in the selected order.",
  "Thresholds are never applied to count endpoints; count endpoints are compared directly.",
  paste0("THRESHOLD_ENDPOINT = ", THRESHOLD_ENDPOINT),
  paste0("THRESHOLD_ENDPOINT_NAME = ", THRESHOLD_ENDPOINT_NAME),
  paste0("THRESHOLD_ENDPOINT_TYPE = ", THRESHOLD_ENDPOINT_TYPE),
  paste0("THRESHOLD_GRID = ", paste(THRESHOLD_GRID, collapse = ", ")),
  "",
  "Weight grid:",
  paste(capture.output(print(WEIGHT_GRID)), collapse = "\n")
)
writeLines(readme, file.path(OUTPUT_DIR, "README_REALDATA_full_adaptive_WR_WO.txt"))


#output

cat("\nDone. Outputs saved in:\n")
cat(OUTPUT_DIR, "\n")
cat("\nMain result table:\n")
cat(file.path(OUTPUT_DIR, "REALDATA_all_WR_WO_method_results.csv"), "\n")
cat("\nLog-rank table:\n")
cat(file.path(OUTPUT_DIR, "REALDATA_logrank_results.csv"), "\n")


