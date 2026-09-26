endpoint_time <- function(name,
                          time = NULL,
                          event = NULL,
                          unit = "years",
                          status_code = NULL) {
  if (missing(name) || !nzchar(as.character(name)[1])) stop("name is required.")
  if (!is.null(unit)) unit <- .wr_normalize_time_unit(unit)
  structure(
    list(
      name = as.character(name)[1],
      type = "time",
      time_col = if (is.null(time)) NULL else as.character(time)[1],
      event_col = if (is.null(event)) NULL else as.character(event)[1],
      unit = unit,
      status_code = status_code
    ),
    class = c("win_endpoint", "list")
  )
}

endpoint_count <- function(name,
                           count = NULL,
                           comparison = c("observed_total", "pairwise_common_followup"),
                           recurrent_time = NULL,
                           recurrent_id = NULL,
                           followup = NULL,
                           unit = "years",
                           followup_unit = NULL,
                           recurrent_time_type = c("absolute", "gap"),
                           status_code = NULL) {
  if (missing(name) || !nzchar(as.character(name)[1])) stop("name is required.")
  comparison <- match.arg(comparison)
  recurrent_time_type <- match.arg(recurrent_time_type)
  if (!is.null(unit)) unit <- .wr_normalize_time_unit(unit)
  if (!is.null(followup_unit)) followup_unit <- .wr_normalize_time_unit(followup_unit)

  structure(
    list(
      name = as.character(name)[1],
      type = "count",
      count_col = if (is.null(count)) NULL else as.character(count)[1],
      comparison = comparison,
      recurrent_time_col = if (is.null(recurrent_time)) NULL else as.character(recurrent_time)[1],
      recurrent_id_col = if (is.null(recurrent_id)) NULL else as.character(recurrent_id)[1],
      followup_col = if (is.null(followup)) NULL else as.character(followup)[1],
      unit = unit,
      followup_unit = followup_unit,
      recurrent_time_type = recurrent_time_type,
      status_code = status_code
    ),
    class = c("win_endpoint", "list")
  )
}

endpoint_binary <- function(name,
                            value,
                            adverse = TRUE) {
  if (missing(name) || !nzchar(as.character(name)[1])) stop("name is required.")
  if (missing(value) || !nzchar(as.character(value)[1])) stop("value is required.")
  structure(
    list(
      name = as.character(name)[1],
      type = "binary",
      value_col = as.character(value)[1],
      adverse = isTRUE(adverse)
    ),
    class = c("win_endpoint", "list")
  )
}

endpoint_continuous <- function(name,
                                value,
                                higher_better = TRUE) {
  if (missing(name) || !nzchar(as.character(name)[1])) stop("name is required.")
  if (missing(value) || !nzchar(as.character(value)[1])) stop("value is required.")
  structure(
    list(
      name = as.character(name)[1],
      type = "continuous",
      value_col = as.character(value)[1],
      higher_better = isTRUE(higher_better)
    ),
    class = c("win_endpoint", "list")
  )
}

.wr_validate_endpoints <- function(endpoints,
                                   data_format = c("wide", "event_long"),
                                   min_endpoints = 2L,
                                   max_endpoints = 3L) {
  data_format <- match.arg(data_format)
  if (!is.list(endpoints)) stop("endpoints must be a list.")

  m <- length(endpoints)
  if (m < min_endpoints || m > max_endpoints) {
    stop("The analysis currently supports ", min_endpoints, " to ", max_endpoints, " endpoints.")
  }

  allowed <- c("time", "count", "binary", "continuous")

  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]

    if (is.null(ep$type) || !(tolower(ep$type) %in% allowed)) {
      stop("Unsupported endpoint type for endpoint ", j, ".")
    }

    ep$type <- tolower(ep$type)

    if (data_format == "wide") {
      if (ep$type == "time") {
        if (is.null(ep$time_col) || is.null(ep$event_col)) {
          stop("Wide time endpoint ", j, " requires time and event columns.")
        }
      } else if (ep$type == "count") {
        comparison <- if (is.null(ep$comparison)) "observed_total" else ep$comparison

        if (comparison == "observed_total" && is.null(ep$count_col)) {
          stop("Wide observed-total count endpoint ", j, " requires a count column.")
        }

        if (comparison == "pairwise_common_followup" &&
            is.null(ep$count_col) &&
            is.null(ep$recurrent_time_col)) {
          stop(
            "Wide pairwise-common-follow-up count endpoint ",
            j,
            " requires recurrent_time or a count column plus recurrent event times."
          )
        }
      } else {
        if (is.null(ep$value_col)) {
          stop("Wide ", ep$type, " endpoint ", j, " requires a value column.")
        }
      }
    } else {
      if (!(ep$type %in% c("time", "count"))) {
        stop("event_long data currently supports only time and count endpoints.")
      }

      if (is.null(ep$status_code) || length(ep$status_code) != 1 || is.na(ep$status_code)) {
        stop("event_long endpoint ", j, " requires status_code.")
      }
    }
  }

  endpoints
}

.wr_endpoint_first_time_index <- function(endpoints) {
  idx <- which(vapply(
    endpoints,
    function(ep) identical(tolower(ep$type), "time"),
    logical(1)
  ))
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

.wr_recurrent_abs_times <- function(values,
                                    unit = "years",
                                    type = c("absolute", "gap"),
                                    followup = Inf) {
  type <- match.arg(type)
  values <- .wr_parse_numeric_times(values)
  values <- .wr_time_to_years(values, unit)

  if (length(values) == 0) return(numeric(0))

  if (type == "gap") {
    values <- cumsum(values)
  } else {
    values <- sort(values)
  }

  values <- values[
    is.finite(values) &
      values > 0 &
      values <= followup
  ]

  as.numeric(values)
}

.wr_build_wide_recurrent_cache <- function(raw,
                                           prepared,
                                           ep,
                                           endpoint_id,
                                           recurrent_data,
                                           id_col,
                                           default_followup_raw_col,
                                           default_followup_unit) {
  recurrent_time_col <- ep$recurrent_time_col
  recurrent_id_col <- if (is.null(ep$recurrent_id_col)) id_col else ep$recurrent_id_col
  recurrent_unit <- if (is.null(ep$unit)) "years" else ep$unit
  recurrent_type <- if (is.null(ep$recurrent_time_type)) "absolute" else ep$recurrent_time_type

  followup_raw_col <- if (!is.null(ep$followup_col)) {
    ep$followup_col
  } else {
    default_followup_raw_col
  }

  if (is.null(followup_raw_col) || !(followup_raw_col %in% names(raw))) {
    stop(
      "pairwise_common_followup count endpoint ",
      endpoint_id,
      " requires followup= or at least one time endpoint whose time column is the subject follow-up time."
    )
  }

  followup_unit <- if (!is.null(ep$followup_unit)) {
    ep$followup_unit
  } else if (!is.null(ep$followup_col)) {
    recurrent_unit
  } else {
    default_followup_unit
  }

  followup_internal_col <- paste0("e", endpoint_id, "_followup")
  prepared[[followup_internal_col]] <- .wr_time_to_years(
    raw[[followup_raw_col]],
    followup_unit
  )

  ids <- as.character(prepared$id)
  abs_times_list <- vector("list", nrow(prepared))

  if (!is.null(recurrent_data)) {
    if (!is.data.frame(recurrent_data)) recurrent_data <- as.data.frame(recurrent_data)

    if (is.null(recurrent_time_col) || !(recurrent_time_col %in% names(recurrent_data))) {
      stop(
        "Recurrent event-time column not found for endpoint ",
        endpoint_id,
        ": ",
        if (is.null(recurrent_time_col)) "<not specified>" else recurrent_time_col
      )
    }

    if (!(recurrent_id_col %in% names(recurrent_data))) {
      stop(
        "Recurrent subject-ID column not found for endpoint ",
        endpoint_id,
        ": ",
        recurrent_id_col
      )
    }

    recurrent_ids <- as.character(recurrent_data[[recurrent_id_col]])

    for (i in seq_len(nrow(prepared))) {
      z <- recurrent_data[[recurrent_time_col]][recurrent_ids == ids[i]]
      abs_times_list[[i]] <- .wr_recurrent_abs_times(
        z,
        unit = recurrent_unit,
        type = recurrent_type,
        followup = prepared[[followup_internal_col]][i]
      )
    }
  } else {
    if (is.null(recurrent_time_col) || !(recurrent_time_col %in% names(raw))) {
      stop(
        "No recurrent-event table was supplied and subject-level recurrent_time column was not found for endpoint ",
        endpoint_id,
        "."
      )
    }

    raw_id <- as.character(raw[[id_col]])
    match_index <- match(ids, raw_id)

    for (i in seq_len(nrow(prepared))) {
      z <- raw[[recurrent_time_col]][match_index[i]]
      abs_times_list[[i]] <- .wr_recurrent_abs_times(
        z,
        unit = recurrent_unit,
        type = recurrent_type,
        followup = prepared[[followup_internal_col]][i]
      )
    }
  }

  list(
    prepared = prepared,
    cache = list(
      endpoint_id = endpoint_id,
      abs_times_list = abs_times_list,
      followup_col = followup_internal_col,
      recurrent_time_unit = "years"
    ),
    followup_col = followup_internal_col
  )
}

.wr_prepare_wide_data <- function(raw,
                                  id_col,
                                  arm_col,
                                  endpoints,
                                  recurrent_data = NULL) {
  endpoints <- .wr_validate_endpoints(endpoints, "wide")

  if (!is.data.frame(raw)) raw <- as.data.frame(raw)
  if (!(id_col %in% names(raw))) stop("ID column not found: ", id_col)
  if (!(arm_col %in% names(raw))) stop("Treatment column not found: ", arm_col)

  d <- data.frame(
    id = raw[[id_col]],
    arm = .wr_safe_num(raw[[arm_col]]),
    stringsAsFactors = FALSE
  )

  if (any(!is.na(d$arm) & !(d$arm %in% c(0, 1)))) {
    stop("Treatment column must be coded 1 = treatment and 0 = control.")
  }

  keep <- !is.na(d$id) & !is.na(d$arm)
  raw <- raw[keep, , drop = FALSE]
  d <- d[keep, , drop = FALSE]

  if (any(duplicated(d$id))) stop("Wide data must have one row per subject.")

  endpoint_specs <- vector("list", length(endpoints))
  recurrent_cache <- vector("list", length(endpoints))

  first_time_id <- .wr_endpoint_first_time_index(endpoints)
  default_followup_raw_col <- if (!is.na(first_time_id)) {
    endpoints[[first_time_id]]$time_col
  } else {
    NULL
  }
  default_followup_unit <- if (!is.na(first_time_id)) {
    endpoints[[first_time_id]]$unit
  } else {
    "years"
  }

  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]
    ep_type <- tolower(ep$type)
    ep_name <- ep$name

    if (ep_type == "time") {
      if (!(ep$time_col %in% names(raw))) {
        stop("Missing time column for endpoint ", j, ": ", ep$time_col)
      }

      if (!(ep$event_col %in% names(raw))) {
        stop("Missing event column for endpoint ", j, ": ", ep$event_col)
      }

      tcol <- paste0("e", j, "_time")
      ecol <- paste0("e", j, "_event")
      ep_unit <- if (is.null(ep$unit)) "years" else ep$unit

      d[[tcol]] <- .wr_time_to_years(raw[[ep$time_col]], ep_unit)

      event_raw <- .wr_safe_num(raw[[ep$event_col]])
      d[[ecol]] <- ifelse(
        is.na(event_raw),
        NA_integer_,
        as.integer(event_raw == 1)
      )

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "time",
        time_col = tcol,
        event_col = ecol,
        input_unit = .wr_normalize_time_unit(ep_unit),
        analysis_unit = "years"
      )
    } else if (ep_type == "count") {
      comparison <- if (is.null(ep$comparison)) "observed_total" else ep$comparison
      ccol <- paste0("e", j, "_count")

      if (!is.null(ep$count_col) && ep$count_col %in% names(raw)) {
        d[[ccol]] <- .wr_safe_num(raw[[ep$count_col]])
      } else {
        d[[ccol]] <- NA_real_
      }

      followup_internal_col <- NULL

      if (comparison == "pairwise_common_followup") {
        rec <- .wr_build_wide_recurrent_cache(
          raw = raw,
          prepared = d,
          ep = ep,
          endpoint_id = j,
          recurrent_data = recurrent_data,
          id_col = id_col,
          default_followup_raw_col = default_followup_raw_col,
          default_followup_unit = default_followup_unit
        )

        d <- rec$prepared
        recurrent_cache[[j]] <- rec$cache
        followup_internal_col <- rec$followup_col
        d[[ccol]] <- as.numeric(lengths(rec$cache$abs_times_list))
      } else {
        if (is.null(ep$count_col) || !(ep$count_col %in% names(raw))) {
          stop("Missing count column for endpoint ", j, ": ", ep$count_col)
        }
      }

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "count",
        count_col = ccol,
        comparison = comparison,
        recurrent_time_col = ep$recurrent_time_col,
        recurrent_id_col = ep$recurrent_id_col,
        followup_col = followup_internal_col,
        recurrent_time_type = if (is.null(ep$recurrent_time_type)) "absolute" else ep$recurrent_time_type,
        input_unit = if (is.null(ep$unit)) "years" else .wr_normalize_time_unit(ep$unit)
      )
    } else if (ep_type == "binary") {
      if (!(ep$value_col %in% names(raw))) {
        stop("Missing binary value column for endpoint ", j, ": ", ep$value_col)
      }

      vcol <- paste0("e", j, "_value")
      d[[vcol]] <- .wr_safe_num(raw[[ep$value_col]])

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "binary",
        value_col = vcol,
        adverse = if (is.null(ep$adverse)) TRUE else isTRUE(ep$adverse)
      )
    } else if (ep_type == "continuous") {
      if (!(ep$value_col %in% names(raw))) {
        stop("Missing continuous value column for endpoint ", j, ": ", ep$value_col)
      }

      vcol <- paste0("e", j, "_value")
      d[[vcol]] <- .wr_safe_num(raw[[ep$value_col]])

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "continuous",
        value_col = vcol,
        higher_better = if (is.null(ep$higher_better)) TRUE else isTRUE(ep$higher_better)
      )
    }
  }

  list(
    data = d,
    endpoints = endpoint_specs,
    recurrent_cache = recurrent_cache
  )
}

.wr_prepare_event_long_data <- function(raw,
                                        id_col,
                                        arm_col,
                                        time_col,
                                        status_col,
                                        endpoints,
                                        time_unit = "years") {
  endpoints <- .wr_validate_endpoints(endpoints, "event_long")

  if (!is.data.frame(raw)) raw <- as.data.frame(raw)

  required <- c(id_col, arm_col, time_col, status_col)
  missing_cols <- setdiff(required, names(raw))

  if (length(missing_cols) > 0) {
    stop("Missing event-long columns: ", paste(missing_cols, collapse = ", "))
  }

  d0 <- data.frame(
    id = raw[[id_col]],
    arm = .wr_safe_num(raw[[arm_col]]),
    time = .wr_time_to_years(raw[[time_col]], time_unit),
    status = .wr_safe_num(raw[[status_col]]),
    stringsAsFactors = FALSE
  )

  if (any(!is.na(d0$arm) & !(d0$arm %in% c(0, 1)))) {
    stop("Treatment column must be coded 1 = treatment and 0 = control.")
  }

  ids <- unique(d0$id[!is.na(d0$id)])

  subj <- data.frame(
    id = ids,
    arm = vapply(
      ids,
      function(z) .wr_first_non_missing(d0$arm[d0$id == z]),
      numeric(1)
    ),
    followup_time = vapply(
      ids,
      function(z) {
        x <- d0$time[d0$id == z]
        if (length(x) == 0 || all(!is.finite(x))) return(NA_real_)
        max(x[is.finite(x)], na.rm = TRUE)
      },
      numeric(1)
    ),
    stringsAsFactors = FALSE
  )

  endpoint_specs <- vector("list", length(endpoints))
  recurrent_cache <- vector("list", length(endpoints))

  for (j in seq_along(endpoints)) {
    ep <- endpoints[[j]]
    code <- ep$status_code
    ep_name <- ep$name

    if (ep$type == "time") {
      tcol <- paste0("e", j, "_time")
      ecol <- paste0("e", j, "_event")

      event_time <- vapply(
        ids,
        function(z) {
          x <- d0$time[d0$id == z & d0$status == code]
          x <- x[is.finite(x)]
          if (length(x) == 0) return(NA_real_)
          min(x)
        },
        numeric(1)
      )

      subj[[ecol]] <- ifelse(is.na(event_time), 0L, 1L)
      subj[[tcol]] <- ifelse(
        is.na(event_time),
        subj$followup_time,
        event_time
      )

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "time",
        time_col = tcol,
        event_col = ecol,
        input_unit = .wr_normalize_time_unit(time_unit),
        analysis_unit = "years"
      )
    } else if (ep$type == "count") {
      comparison <- if (is.null(ep$comparison)) "observed_total" else ep$comparison
      ccol <- paste0("e", j, "_count")

      abs_times_list <- lapply(
        ids,
        function(z) {
          x <- d0$time[d0$id == z & d0$status == code]
          sort(x[is.finite(x) & x > 0])
        }
      )

      subj[[ccol]] <- as.numeric(lengths(abs_times_list))

      followup_col <- NULL

      if (comparison == "pairwise_common_followup") {
        followup_col <- "followup_time"

        for (i in seq_along(abs_times_list)) {
          abs_times_list[[i]] <- abs_times_list[[i]][
            abs_times_list[[i]] <= subj$followup_time[i]
          ]
        }

        recurrent_cache[[j]] <- list(
          endpoint_id = j,
          abs_times_list = abs_times_list,
          followup_col = followup_col,
          recurrent_time_unit = "years"
        )
      }

      endpoint_specs[[j]] <- list(
        id = j,
        name = ep_name,
        type = "count",
        count_col = ccol,
        comparison = comparison,
        recurrent_time_col = time_col,
        recurrent_id_col = id_col,
        followup_col = followup_col,
        recurrent_time_type = "absolute",
        input_unit = .wr_normalize_time_unit(time_unit)
      )
    }
  }

  list(
    data = subj,
    endpoints = endpoint_specs,
    recurrent_cache = recurrent_cache
  )
}

.wr_prepare_endpoint_cache <- function(data,
                                       endpoint_specs,
                                       recurrent_cache = NULL) {
  n <- nrow(data)
  m <- length(endpoint_specs)

  type_code <- integer(m)
  time_mat <- matrix(NA_real_, nrow = n, ncol = m)
  event_mat <- matrix(NA_integer_, nrow = n, ncol = m)
  value_mat <- matrix(NA_real_, nrow = n, ncol = m)
  followup_mat <- matrix(NA_real_, nrow = n, ncol = m)

  recurrent_times <- vector("list", m)
  recurrent_start <- vector("list", m)
  recurrent_len <- vector("list", m)

  for (j in seq_len(m)) {
    ep <- endpoint_specs[[j]]
    ep_type <- tolower(ep$type)

    recurrent_times[[j]] <- numeric(0)
    recurrent_start[[j]] <- integer(n)
    recurrent_len[[j]] <- integer(n)

    if (ep_type == "time") {
      type_code[j] <- 1L
      time_mat[, j] <- .wr_safe_num(data[[ep$time_col]])

      ev <- .wr_safe_num(data[[ep$event_col]])
      event_mat[, j] <- ifelse(
        is.na(ev),
        NA_integer_,
        as.integer(ev == 1)
      )
    } else if (ep_type == "count") {
      comparison <- if (is.null(ep$comparison)) "observed_total" else ep$comparison
      value_mat[, j] <- .wr_safe_num(data[[ep$count_col]])

      if (comparison == "pairwise_common_followup") {
        type_code[j] <- 7L

        rec <- if (!is.null(recurrent_cache) && length(recurrent_cache) >= j) {
          recurrent_cache[[j]]
        } else {
          NULL
        }

        if (is.null(rec) || is.null(rec$abs_times_list) || is.null(rec$followup_col)) {
          stop(
            "Missing recurrent-event cache for pairwise_common_followup endpoint ",
            j,
            "."
          )
        }

        if (!(rec$followup_col %in% names(data))) {
          stop(
            "Missing follow-up column for pairwise_common_followup endpoint ",
            j,
            ": ",
            rec$followup_col
          )
        }

        followup_mat[, j] <- .wr_safe_num(data[[rec$followup_col]])

        abs_times_list <- rec$abs_times_list
        if (length(abs_times_list) != n) {
          stop("Recurrent-event cache length does not match subject count for endpoint ", j, ".")
        }

        lens <- as.integer(lengths(abs_times_list))
        starts <- if (n == 0L) {
          integer(0)
        } else {
          as.integer(cumsum(c(0L, lens[-n])))
        }

        flat <- as.numeric(unlist(abs_times_list, use.names = FALSE))
        if (length(flat) == 0L) flat <- numeric(0)

        recurrent_times[[j]] <- flat
        recurrent_start[[j]] <- starts
        recurrent_len[[j]] <- lens
      } else {
        type_code[j] <- 2L
      }
    } else if (ep_type == "binary") {
      type_code[j] <- if (isTRUE(ep$adverse)) 3L else 6L
      value_mat[, j] <- .wr_safe_num(data[[ep$value_col]])
    } else if (ep_type == "continuous") {
      type_code[j] <- if (isTRUE(ep$higher_better)) 4L else 5L
      value_mat[, j] <- .wr_safe_num(data[[ep$value_col]])
    } else {
      stop("Unsupported endpoint type: ", ep$type)
    }
  }

  list(
    type_code = type_code,
    time_mat = time_mat,
    event_mat = event_mat,
    value_mat = value_mat,
    followup_mat = followup_mat,
    recurrent_times = recurrent_times,
    recurrent_start = recurrent_start,
    recurrent_len = recurrent_len
  )
}
