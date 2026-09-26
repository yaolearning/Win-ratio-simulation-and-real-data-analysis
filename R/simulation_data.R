.wr_simulate_one_dataset <- function(n_control = 50L,
                                     n_treatment = 50L,
                                     mort_rate_control = -log(0.6),
                                     hazard_ratio = 1,
                                     hosp_shape = 5,
                                     hosp_scale_control = 1,
                                     hosp_scale_treatment = 1,
                                     followup = 1,
                                     censor_rate = 0,
                                     seed = NULL) {
  n_control <- as.integer(n_control)
  n_treatment <- as.integer(n_treatment)

  if (n_control < 1L || n_treatment < 1L) {
    stop("Both groups must contain at least one subject.")
  }

  if (!is.finite(mort_rate_control) || mort_rate_control <= 0) {
    stop("mort_rate_control must be positive.")
  }

  if (!is.finite(hazard_ratio) || hazard_ratio <= 0) {
    stop("hazard_ratio must be positive.")
  }

  if (!is.finite(hosp_shape) || hosp_shape <= 0) {
    stop("hosp_shape must be positive.")
  }

  if (!is.finite(hosp_scale_control) || hosp_scale_control <= 0 ||
      !is.finite(hosp_scale_treatment) || hosp_scale_treatment <= 0) {
    stop("Hospitalization scale parameters must be positive.")
  }

  if (!is.finite(followup) || followup <= 0) {
    stop("followup must be positive.")
  }

  if (!is.finite(censor_rate) || censor_rate < 0) {
    stop("censor_rate must be non-negative.")
  }

  if (!is.null(seed)) set.seed(as.integer(seed))

  n <- n_control + n_treatment

  table_output <- data.frame(
    SUBJID = seq_len(n),
    ARM = c(rep(0L, n_control), rep(1L, n_treatment)),
    FUTIME = rep(followup, n),
    CNSR = rep(0L, n),
    SURVTIME = c(
      stats::rexp(n_control, rate = mort_rate_control),
      stats::rexp(n_treatment, rate = mort_rate_control * hazard_ratio)
    ),
    CNSRTIME = rep(Inf, n),
    FREQHOSP = c(
      stats::rgamma(
        n_control,
        shape = hosp_shape,
        scale = hosp_scale_control
      ),
      stats::rgamma(
        n_treatment,
        shape = hosp_shape,
        scale = hosp_scale_treatment
      )
    ),
    NUMHOSP = rep(0L, n),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    table_output$CNSR[i] <- as.integer(
      table_output$SURVTIME[i] < table_output$FUTIME[i] &&
        table_output$SURVTIME[i] < table_output$CNSRTIME[i]
    )

    table_output$FUTIME[i] <- min(
      table_output$FUTIME[i],
      table_output$SURVTIME[i],
      table_output$CNSRTIME[i]
    )
  }

  hosp_times_list <- vector("list", n)

  for (i in seq_len(n)) {
    gaps <- numeric(0)
    cumulative_time <- 0
    followup_time <- table_output$FUTIME[i]
    hosp_rate <- table_output$FREQHOSP[i]

    while (is.finite(hosp_rate) && hosp_rate > 0) {
      gap <- stats::rexp(1, rate = hosp_rate)
      cumulative_time <- cumulative_time + gap

      if (cumulative_time < followup_time) {
        gaps <- c(gaps, gap)
      } else {
        break
      }
    }

    hosp_times_list[[i]] <- gaps
    table_output$NUMHOSP[i] <- length(gaps)
  }

  if (censor_rate > 0) {
    if (!is.null(seed)) set.seed(as.integer(seed) + 17L)
    censor_time <- stats::rexp(n, rate = censor_rate)

    for (i in seq_len(n)) {
      if (censor_time[i] < table_output$FUTIME[i]) {
        table_output$FUTIME[i] <- censor_time[i]
        table_output$CNSRTIME[i] <- censor_time[i]
        table_output$CNSR[i] <- 0L

        x <- hosp_times_list[[i]]

        if (length(x) == 0 || all(is.na(x))) {
          hosp_times_list[[i]] <- numeric(0)
          table_output$NUMHOSP[i] <- 0L
        } else {
          absolute_times <- cumsum(as.numeric(x[!is.na(x)]))
          keep <- absolute_times <= table_output$FUTIME[i]
          absolute_keep <- absolute_times[keep]

          if (length(absolute_keep) == 0) {
            hosp_times_list[[i]] <- numeric(0)
            table_output$NUMHOSP[i] <- 0L
          } else {
            hosp_times_list[[i]] <- diff(c(0, absolute_keep))
            table_output$NUMHOSP[i] <- length(absolute_keep)
          }
        }
      }
    }
  }

  list(
    table.output = table_output,
    hosp.times.list = hosp_times_list
  )
}

.wr_simulated_recurrent_table <- function(ds) {
  rows <- vector("list", nrow(ds$table.output))

  for (i in seq_len(nrow(ds$table.output))) {
    gaps <- ds$hosp.times.list[[i]]

    if (length(gaps) == 0) {
      rows[[i]] <- NULL
    } else {
      rows[[i]] <- data.frame(
        SUBJID = rep(ds$table.output$SUBJID[i], length(gaps)),
        HOSPTIME = cumsum(as.numeric(gaps)),
        stringsAsFactors = FALSE
      )
    }
  }

  rows <- rows[!vapply(rows, is.null, logical(1))]

  if (length(rows) == 0) {
    return(data.frame(
      SUBJID = numeric(0),
      HOSPTIME = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

generate_win_dataset <- function(n_control = 50L,
                                 n_treatment = 50L,
                                 mort_rate_control = -log(0.6),
                                 hazard_ratio = 1,
                                 hosp_shape = 5,
                                 hosp_scale_control = 1,
                                 hosp_scale_treatment = 1,
                                 followup = 1,
                                 censor_rate = 0,
                                 seed = 2026L) {
  ds <- .wr_simulate_one_dataset(
    n_control = n_control,
    n_treatment = n_treatment,
    mort_rate_control = mort_rate_control,
    hazard_ratio = hazard_ratio,
    hosp_shape = hosp_shape,
    hosp_scale_control = hosp_scale_control,
    hosp_scale_treatment = hosp_scale_treatment,
    followup = followup,
    censor_rate = censor_rate,
    seed = seed
  )

  out <- list(
    subjects = ds$table.output,
    recurrent_events = .wr_simulated_recurrent_table(ds),
    analysis_data = ds,
    settings = list(
      n_control = as.integer(n_control),
      n_treatment = as.integer(n_treatment),
      mort_rate_control = mort_rate_control,
      hazard_ratio = hazard_ratio,
      hosp_shape = hosp_shape,
      hosp_scale_control = hosp_scale_control,
      hosp_scale_treatment = hosp_scale_treatment,
      followup = followup,
      censor_rate = censor_rate,
      seed = seed
    )
  )

  class(out) <- c("win_simulated_dataset", "list")
  out
}

print.win_simulated_dataset <- function(x, ...) {
  cat("Simulated win-ratio dataset\n")
  cat("Subjects:", nrow(x$subjects), "\n")
  cat("Control:", sum(x$subjects$ARM == 0), "\n")
  cat("Treatment:", sum(x$subjects$ARM == 1), "\n")
  cat("Recurrent events:", nrow(x$recurrent_events), "\n")
  invisible(x)
}
