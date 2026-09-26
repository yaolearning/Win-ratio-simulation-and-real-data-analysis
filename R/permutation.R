.wr_analysis_arm <- function(prepared, backend) {
  if (backend == "two_endpoint_recurrent") {
    return(as.integer(prepared$table.output$ARM))
  }
  as.integer(prepared$data$arm)
}

.wr_set_analysis_arm <- function(prepared, backend, arm) {
  if (backend == "two_endpoint_recurrent") {
    prepared$table.output$ARM <- as.integer(arm)
    return(prepared)
  }
  prepared$data$arm <- as.integer(arm)
  prepared
}

.wr_method_by_id <- function(methods, method_id) {
  idx <- which(vapply(methods, function(x) identical(x$id, method_id), logical(1)))
  if (length(idx) != 1L) stop("Method not found: ", method_id)
  methods[[idx]]
}

.wr_selected_candidate_row <- function(counts,
                                       method,
                                       measure,
                                       side) {
  .wr_select_best_candidate(
    counts = counts,
    candidate_keys = method$grid$candidate_key,
    measure = measure,
    side = side
  )
}

.wr_selection_record <- function(row,
                                 measure,
                                 method,
                                 side,
                                 m) {
  out <- .wr_flatten_selection(
    row = row,
    measure = measure,
    method_id = method$id,
    side = side,
    m = m
  )

  out$notation <- method$notation
  out$description <- method$description
  out$n_candidates <- nrow(method$grid)
  out$candidate_key <- if (nrow(row) > 0) as.character(row$candidate_key[1]) else NA_character_
  out$selected_value <- if (nrow(row) == 0) {
    NA_real_
  } else if (measure == "WR") {
    as.numeric(row$WR_statistic[1])
  } else {
    as.numeric(row$WO_statistic[1])
  }

  log_value <- .wr_log_ratio(out$selected_value)
  out$selected_log_value <- as.numeric(log_value[1])
  out$selected_abs_log_value <- abs(out$selected_log_value)
  out$statistic <- if (side == "one") out$selected_value else out$selected_abs_log_value
  out$selected_direction <- if (is.na(out$selected_log_value)) {
    NA_character_
  } else if (out$selected_log_value >= 0) {
    "upper_benefit"
  } else {
    "lower_harm"
  }

  out
}

.wr_observed_selection_rows <- function(counts,
                                        methods,
                                        method_control,
                                        m) {
  measures <- .wr_measures_to_run(method_control)
  rows <- list()
  idx <- 1L

  for (measure in measures) {
    for (method in methods) {
      for (side in c("one", "two")) {
        selected <- .wr_selected_candidate_row(
          counts = counts,
          method = method,
          measure = measure,
          side = side
        )

        rows[[idx]] <- .wr_selection_record(
          row = selected,
          measure = measure,
          method = method,
          side = side,
          m = m
        )
        idx <- idx + 1L
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$key <- paste(out$measure, out$method_id, out$side, sep = "__")
  out
}

.wr_stat_from_candidate <- function(row,
                                    measure,
                                    side) {
  if (nrow(row) == 0) return(NA_real_)

  if (measure == "WR") {
    if (side == "one") return(as.numeric(row$WR_statistic[1]))
    return(as.numeric(row$WR_abslog[1]))
  }

  if (side == "one") return(as.numeric(row$WO_statistic[1]))
  as.numeric(row$WO_abslog[1])
}

.wr_perm_selection_row <- function(selected,
                                   measure,
                                   method_id,
                                   side,
                                   b,
                                   m) {
  out <- data.frame(
    b = as.integer(b),
    measure = measure,
    method_id = method_id,
    side = side,
    candidate_key = if (nrow(selected) > 0) as.character(selected$candidate_key[1]) else NA_character_,
    statistic = .wr_stat_from_candidate(selected, measure, side),
    selected_value = if (nrow(selected) == 0) {
      NA_real_
    } else if (measure == "WR") {
      as.numeric(selected$WR_statistic[1])
    } else {
      as.numeric(selected$WO_statistic[1])
    },
    selected_order = if (nrow(selected) > 0) as.character(selected$order_key[1]) else NA_character_,
    tie_count = if (nrow(selected) > 0) as.numeric(selected$tie_count[1]) else NA_real_,
    tie_proportion = if (nrow(selected) > 0) as.numeric(selected$tie_proportion[1]) else NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(m)) {
    out[[paste0("selected_p", j)]] <- if (nrow(selected) > 0) {
      as.numeric(selected[[paste0("p", j)]][1])
    } else {
      NA_real_
    }

    out[[paste0("selected_t_e", j)]] <- if (nrow(selected) > 0) {
      as.numeric(selected[[paste0("threshold_e", j)]][1])
    } else {
      NA_real_
    }
  }

  out
}

run_win_permutation <- function(prepared,
                                candidate_grid,
                                methods,
                                method_control,
                                permutation = permutation_control(),
                                engine = engine_control()) {
  backend <- engine$backend
  if (backend == "auto") {
    if (is.list(prepared) &&
        all(c("table.output", "hosp.times.list") %in% names(prepared)) &&
        !is.null(prepared$hosp.times.list)) {
      backend <- "two_endpoint_recurrent"
    } else {
      backend <- "general"
    }
  }

  engine_resolved <- engine
  engine_resolved$backend <- backend

  m <- if (backend == "two_endpoint_recurrent") 2L else length(prepared$endpoints)

  observed_counts <- evaluate_win_candidates(
    prepared = prepared,
    candidate_grid = candidate_grid,
    engine = engine_resolved
  )

  observed <- .wr_observed_selection_rows(
    counts = observed_counts,
    methods = methods,
    method_control = method_control,
    m = m
  )

  observed$adaptive_p_value <- NA_real_
  observed$fixed_selected_p_value <- NA_real_
  observed$mean_perm_tie_count <- NA_real_
  observed$mean_perm_tie_proportion <- NA_real_
  observed$mean_fixed_perm_tie_count <- NA_real_
  observed$mean_fixed_perm_tie_proportion <- NA_real_

  if (!isTRUE(permutation$enabled)) {
    return(list(
      observed_counts = observed_counts,
      observed = observed,
      perm_selected = data.frame(),
      B = 0L,
      backend = backend
    ))
  }

  B <- permutation$B
  keys <- observed$key
  n_keys <- length(keys)

  adaptive_stat <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  fixed_stat <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  adaptive_tie <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  adaptive_tie_pr <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  fixed_tie <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  fixed_tie_pr <- matrix(
    NA_real_,
    nrow = B,
    ncol = n_keys,
    dimnames = list(NULL, keys)
  )

  perm_selected_list <- vector("list", B)
  original_arm <- .wr_analysis_arm(prepared, backend)
  perm_seeds <- as.integer(permutation$seed + seq_len(B) * 1009L)

  for (b in seq_len(B)) {
    set.seed(perm_seeds[b])
    permuted_arm <- sample(original_arm, replace = FALSE)

    prepared_b <- .wr_set_analysis_arm(
      prepared = prepared,
      backend = backend,
      arm = permuted_arm
    )

    counts_b <- evaluate_win_candidates(
      prepared = prepared_b,
      candidate_grid = candidate_grid,
      engine = engine_resolved
    )

    rows_b <- vector("list", n_keys)

    for (j in seq_len(n_keys)) {
      obs_row <- observed[j, , drop = FALSE]
      method <- .wr_method_by_id(methods, obs_row$method_id)
      measure <- as.character(obs_row$measure)
      side <- as.character(obs_row$side)

      selected_b <- .wr_selected_candidate_row(
        counts = counts_b,
        method = method,
        measure = measure,
        side = side
      )

      fixed_b <- counts_b[
        counts_b$candidate_key == obs_row$candidate_key,
        ,
        drop = FALSE
      ]

      adaptive_stat[b, j] <- .wr_stat_from_candidate(selected_b, measure, side)
      fixed_stat[b, j] <- .wr_stat_from_candidate(fixed_b, measure, side)

      adaptive_tie[b, j] <- if (nrow(selected_b) > 0) selected_b$tie_count[1] else NA_real_
      adaptive_tie_pr[b, j] <- if (nrow(selected_b) > 0) selected_b$tie_proportion[1] else NA_real_
      fixed_tie[b, j] <- if (nrow(fixed_b) > 0) fixed_b$tie_count[1] else NA_real_
      fixed_tie_pr[b, j] <- if (nrow(fixed_b) > 0) fixed_b$tie_proportion[1] else NA_real_

      rows_b[[j]] <- .wr_perm_selection_row(
        selected = selected_b,
        measure = measure,
        method_id = obs_row$method_id,
        side = side,
        b = b,
        m = m
      )
    }

    perm_selected_list[[b]] <- do.call(rbind, rows_b)
  }

  perm_selected <- do.call(rbind, perm_selected_list)
  rownames(perm_selected) <- NULL

  for (j in seq_len(n_keys)) {
    obs_stat <- as.numeric(observed$statistic[j])

    observed$adaptive_p_value[j] <- .wr_right_tail_perm_p(
      adaptive_stat[, j],
      obs_stat
    )

    observed$fixed_selected_p_value[j] <- .wr_right_tail_perm_p(
      fixed_stat[, j],
      obs_stat
    )

    observed$mean_perm_tie_count[j] <- .wr_safe_mean(adaptive_tie[, j])
    observed$mean_perm_tie_proportion[j] <- .wr_safe_mean(adaptive_tie_pr[, j])
    observed$mean_fixed_perm_tie_count[j] <- .wr_safe_mean(fixed_tie[, j])
    observed$mean_fixed_perm_tie_proportion[j] <- .wr_safe_mean(fixed_tie_pr[, j])
  }

  list(
    observed_counts = observed_counts,
    observed = observed,
    perm_selected = perm_selected,
    adaptive_stat = adaptive_stat,
    fixed_stat = fixed_stat,
    B = B,
    backend = backend
  )
}
