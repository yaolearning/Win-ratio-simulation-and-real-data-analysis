threshold_control <- function(time_endpoints = "auto",
                              candidates = NULL,
                              unit = "months",
                              data_probs = c(0.10, 0.25, 0.50, 0.75, 0.90),
                              max_event_times = 800L,
                              max_thresholds = 10L,
                              seed = 999L) {
  if (!(is.character(time_endpoints) && length(time_endpoints) == 1 && tolower(time_endpoints) == "auto")) {
    time_endpoints <- as.integer(time_endpoints)
    if (any(is.na(time_endpoints))) stop("time_endpoints must be 'auto' or endpoint IDs.")
  }

  if (is.list(unit)) {
    units <- unit
  } else if (length(unit) == 1L) {
    units <- unit
  } else {
    units <- as.list(unit)
  }

  structure(
    list(
      time_endpoints = time_endpoints,
      candidates = candidates,
      unit = units,
      data_probs = as.numeric(data_probs),
      max_event_times = as.integer(max_event_times),
      max_thresholds = as.integer(max_thresholds),
      seed = as.integer(seed)
    ),
    class = c("win_threshold_control", "list")
  )
}

.wr_threshold_unit_for_endpoint <- function(control, endpoint_id) {
  u <- control$unit
  if (is.list(u)) {
    val <- u[[as.character(endpoint_id)]]
    if (is.null(val) && length(u) >= endpoint_id) val <- u[[endpoint_id]]
  } else if (length(u) == 1L) {
    val <- u
  } else {
    val <- u[endpoint_id]
  }

  if (is.null(val) || length(val) == 0 || is.na(val[1]) || !nzchar(as.character(val[1]))) {
    return("years")
  }
  .wr_normalize_time_unit(val[1])
}

.wr_threshold_candidates_for_endpoint <- function(control, endpoint_id) {
  x <- control$candidates
  if (is.null(x)) return(numeric(0))

  if (is.list(x)) {
    val <- x[[as.character(endpoint_id)]]
    if (is.null(val) && length(x) >= endpoint_id) val <- x[[endpoint_id]]
  } else {
    if (endpoint_id != 1L) return(numeric(0))
    val <- x
  }

  if (is.null(val)) return(numeric(0))
  val <- as.numeric(val)
  val <- val[is.finite(val) & val >= 0]
  if (length(val) == 0) return(numeric(0))
  sort(unique(.wr_time_to_years(val, .wr_threshold_unit_for_endpoint(control, endpoint_id))))
}

.wr_resolve_threshold_time_endpoints <- function(endpoint_specs, control) {
  time_ids <- which(vapply(endpoint_specs, function(ep) identical(tolower(ep$type), "time"), logical(1)))
  if (length(time_ids) == 0) return(integer(0))

  user_setting <- control$time_endpoints
  if (is.character(user_setting) && length(user_setting) == 1 && tolower(user_setting) == "auto") {
    if (length(endpoint_specs) == 2L) return(intersect(1L, time_ids))
    return(time_ids)
  }

  ids <- as.integer(user_setting)
  bad <- setdiff(ids, time_ids)
  if (length(bad) > 0) {
    stop("Thresholds can only be assigned to time endpoints: ", paste(bad, collapse = ", "))
  }
  unique(ids)
}

.wr_auto_threshold_grid <- function(data,
                                    endpoint_spec,
                                    control,
                                    endpoint_id) {
  if (is.null(endpoint_spec) || !identical(tolower(endpoint_spec$type), "time")) {
    return(numeric(0))
  }

  t <- .wr_safe_num(data[[endpoint_spec$time_col]])
  e <- .wr_safe_num(data[[endpoint_spec$event_col]])
  t_event <- t[e == 1 & is.finite(t)]

  if (length(t_event) < 2) return(numeric(0))

  if (length(t_event) > control$max_event_times) {
    set.seed(as.integer(control$seed + endpoint_id))
    t_event <- sample(t_event, control$max_event_times, replace = FALSE)
  }

  diffs <- abs(as.vector(stats::dist(t_event)))
  diffs <- diffs[is.finite(diffs) & diffs > 0]
  if (length(diffs) == 0) return(numeric(0))

  probs <- sort(unique(as.numeric(control$data_probs)))
  probs <- probs[is.finite(probs) & probs > 0 & probs < 1]
  if (length(probs) == 0) return(numeric(0))

  vals <- as.numeric(stats::quantile(diffs, probs = probs, na.rm = TRUE, names = FALSE))
  vals <- sort(unique(round(vals[is.finite(vals) & vals > 0], 8)))

  if (length(vals) > control$max_thresholds) {
    idx <- unique(round(seq(1, length(vals), length.out = control$max_thresholds)))
    vals <- vals[idx]
  }

  vals
}

.wr_expand_grid <- function(values_list) {
  if (length(values_list) == 0) return(data.frame())
  args <- c(values_list, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  do.call(expand.grid, args)
}

.wr_make_threshold_grid <- function(data,
                                    endpoint_specs,
                                    control = threshold_control()) {
  active_ids <- .wr_resolve_threshold_time_endpoints(endpoint_specs, control)
  m <- length(endpoint_specs)
  values_list <- vector("list", m)
  source_vec <- rep("none", m)
  display_values <- rep("", m)
  display_units <- rep("", m)

  for (j in seq_len(m)) {
    if (j %in% active_ids) {
      manual <- .wr_threshold_candidates_for_endpoint(control, j)
      automatic <- .wr_auto_threshold_grid(data, endpoint_specs[[j]], control, j)
      vals <- sort(unique(c(manual, automatic)))
      if (length(vals) == 0) vals <- 0

      values_list[[j]] <- vals
      if (length(manual) > 0 && length(automatic) > 0) {
        source_vec[j] <- "candidate+data-driven"
      } else if (length(manual) > 0) {
        source_vec[j] <- "candidate"
      } else if (length(automatic) > 0) {
        source_vec[j] <- "data-driven"
      } else {
        source_vec[j] <- "zero_fallback"
      }

      unit_j <- .wr_threshold_unit_for_endpoint(control, j)
      display_units[j] <- unit_j
      display_values[j] <- paste(round(.wr_years_to_unit(vals, unit_j), 8), collapse = ",")
    } else {
      values_list[[j]] <- 0
      source_vec[j] <- if (identical(tolower(endpoint_specs[[j]]$type), "time")) {
        "inactive_time_endpoint"
      } else {
        "not_time_endpoint"
      }
      display_units[j] <- if (identical(tolower(endpoint_specs[[j]]$type), "time")) {
        .wr_threshold_unit_for_endpoint(control, j)
      } else {
        ""
      }
      display_values[j] <- "0"
    }
    names(values_list)[j] <- paste0("threshold_e", j)
  }

  grid <- .wr_expand_grid(values_list)
  threshold_cols <- paste0("threshold_e", seq_len(m))
  if (nrow(grid) > 0) {
    grid$threshold_key <- apply(
      grid[, threshold_cols, drop = FALSE],
      1,
      function(z) paste(sprintf("%.8f", as.numeric(z)), collapse = "|")
    )
  } else {
    grid$threshold_key <- character(0)
  }

  info <- data.frame(
    endpoint = seq_len(m),
    endpoint_name = vapply(endpoint_specs, function(ep) as.character(ep$name), character(1)),
    endpoint_type = vapply(endpoint_specs, function(ep) as.character(ep$type), character(1)),
    threshold_active = seq_len(m) %in% active_ids,
    threshold_source = source_vec,
    threshold_internal_unit = "years",
    threshold_candidate_unit = display_units,
    threshold_values_internal = vapply(seq_len(m), function(j) paste(values_list[[j]], collapse = ","), character(1)),
    threshold_values_candidate_unit = display_values,
    stringsAsFactors = FALSE
  )

  list(
    grid = grid,
    info = info,
    active_ids = active_ids
  )
}
