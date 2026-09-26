.wr_logrank_one_endpoint <- function(data,
                                     endpoint_spec,
                                     label) {
  if (!identical(endpoint_spec$type, "time")) {
    return(data.frame(
      test = label,
      available = FALSE,
      reason = "endpoint is not time-to-event",
      chisq = NA_real_,
      HR = NA_real_,
      z = NA_real_,
      p_one_sided_benefit = NA_real_,
      p_two_sided = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  time <- data[[endpoint_spec$time_col]]
  event <- data[[endpoint_spec$event_col]]
  arm <- data$arm

  keep <- is.finite(time) & !is.na(event) & !is.na(arm)
  d <- data.frame(
    time = time[keep],
    event = event[keep],
    arm = arm[keep]
  )

  if (nrow(d) == 0 || length(unique(d$arm)) < 2) {
    return(data.frame(
      test = label,
      available = FALSE,
      reason = "insufficient data",
      chisq = NA_real_,
      HR = NA_real_,
      z = NA_real_,
      p_one_sided_benefit = NA_real_,
      p_two_sided = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  surv_object <- survival::Surv(d$time, d$event == 1)
  lr <- tryCatch(
    survival::survdiff(surv_object ~ d$arm),
    error = function(e) NULL
  )

  fit <- tryCatch(
    survival::coxph(surv_object ~ d$arm),
    error = function(e) NULL
  )

  chisq <- if (is.null(lr)) NA_real_ else as.numeric(lr$chisq)
  hr <- NA_real_
  z <- NA_real_

  if (!is.null(fit)) {
    sm <- tryCatch(summary(fit), error = function(e) NULL)
    if (!is.null(sm) && nrow(sm$coef) >= 1L) {
      hr <- exp(as.numeric(sm$coef[1, "coef"]))
      z <- as.numeric(sm$coef[1, "z"])
    }
  }

  p_one <- if (is.finite(z)) stats::pnorm(z) else NA_real_
  p_two <- if (is.finite(z)) {
    2 * stats::pnorm(-abs(z))
  } else if (is.finite(chisq)) {
    stats::pchisq(chisq, df = 1, lower.tail = FALSE)
  } else {
    NA_real_
  }

  data.frame(
    test = label,
    available = TRUE,
    reason = NA_character_,
    chisq = chisq,
    HR = hr,
    z = z,
    p_one_sided_benefit = p_one,
    p_two_sided = p_two,
    stringsAsFactors = FALSE
  )
}

.wr_composite_time_data <- function(data,
                                    endpoint_specs) {
  if (!all(vapply(endpoint_specs, function(ep) identical(ep$type, "time"), logical(1)))) {
    return(NULL)
  }

  times <- do.call(cbind, lapply(endpoint_specs, function(ep) {
    .wr_safe_num(data[[ep$time_col]])
  }))

  events <- do.call(cbind, lapply(endpoint_specs, function(ep) {
    as.integer(.wr_safe_num(data[[ep$event_col]]) == 1)
  }))

  comp_time <- rep(NA_real_, nrow(data))
  comp_event <- rep(0L, nrow(data))

  for (i in seq_len(nrow(data))) {
    event_times <- times[i, events[i, ] == 1, drop = TRUE]
    event_times <- event_times[is.finite(event_times)]

    if (length(event_times) > 0) {
      comp_time[i] <- min(event_times)
      comp_event[i] <- 1L
    } else {
      available_times <- times[i, is.finite(times[i, ]), drop = TRUE]
      comp_time[i] <- if (length(available_times) > 0) max(available_times) else NA_real_
      comp_event[i] <- 0L
    }
  }

  list(
    data = data.frame(
      id = data$id,
      arm = data$arm,
      comp_time = comp_time,
      comp_event = comp_event,
      stringsAsFactors = FALSE
    ),
    endpoint = list(
      type = "time",
      time_col = "comp_time",
      event_col = "comp_event"
    )
  )
}

run_win_logrank <- function(prepared,
                            run = TRUE,
                            composite_if_possible = TRUE) {
  if (!isTRUE(run)) return(data.frame())
  if (!is.list(prepared) || !all(c("data", "endpoints") %in% names(prepared))) {
    return(data.frame())
  }

  rows <- list()
  idx <- 1L

  time_ids <- which(vapply(
    prepared$endpoints,
    function(ep) identical(ep$type, "time"),
    logical(1)
  ))

  for (j in time_ids) {
    rows[[idx]] <- .wr_logrank_one_endpoint(
      data = prepared$data,
      endpoint_spec = prepared$endpoints[[j]],
      label = paste0("Log-rank outcome ", j)
    )
    idx <- idx + 1L
  }

  if (isTRUE(composite_if_possible)) {
    comp <- .wr_composite_time_data(
      data = prepared$data,
      endpoint_specs = prepared$endpoints
    )

    if (!is.null(comp)) {
      rows[[idx]] <- .wr_logrank_one_endpoint(
        data = comp$data,
        endpoint_spec = comp$endpoint,
        label = "Log-rank composite time-to-first event"
      )
    }
  }

  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
