.wr_get_pair_indices <- function(arm) {
  idx_treatment <- which(arm == 1)
  idx_control <- which(arm == 0)

  if (length(idx_treatment) == 0 || length(idx_control) == 0) {
    stop("Both treatment and control groups must contain at least one subject.")
  }

  list(
    treatment = rep(idx_treatment, each = length(idx_control)),
    control = rep(idx_control, times = length(idx_treatment)),
    n_pairs = length(idx_treatment) * length(idx_control)
  )
}

.wr_count_sorted_times_until <- function(times, t) {
  if (length(times) == 0 || !is.finite(t)) return(0L)
  findInterval(t, times)
}

.wr_compare_endpoint_pairs_R <- function(data,
                                         idx_treatment,
                                         idx_control,
                                         endpoint_spec,
                                         threshold = 0,
                                         recurrent_endpoint = NULL) {
  n <- length(idx_treatment)
  out <- integer(n)
  threshold <- ifelse(is.na(threshold), 0, threshold)

  if (endpoint_spec$type == "time") {
    ta <- .wr_safe_num(data[[endpoint_spec$time_col]][idx_treatment])
    tb <- .wr_safe_num(data[[endpoint_spec$time_col]][idx_control])
    ea <- .wr_safe_num(data[[endpoint_spec$event_col]][idx_treatment])
    eb <- .wr_safe_num(data[[endpoint_spec$event_col]][idx_control])

    both_event <- ea == 1 & eb == 1 & !is.na(ta) & !is.na(tb)
    out[both_event & (ta - tb > threshold)] <- 1L
    out[both_event & (tb - ta > threshold)] <- -1L

    a_censored_b_event <- ea == 0 & eb == 1 & !is.na(ta) & !is.na(tb)
    out[a_censored_b_event & (ta - tb > threshold)] <- 1L

    a_event_b_censored <- ea == 1 & eb == 0 & !is.na(ta) & !is.na(tb)
    out[a_event_b_censored & (tb - ta > threshold)] <- -1L

    return(out)
  }

  if (endpoint_spec$type == "count") {
    comparison <- if (is.null(endpoint_spec$comparison)) {
      "observed_total"
    } else {
      endpoint_spec$comparison
    }

    if (comparison == "pairwise_common_followup") {
      if (is.null(recurrent_endpoint) ||
          is.null(recurrent_endpoint$abs_times_list) ||
          is.null(recurrent_endpoint$followup_col)) {
        stop("Recurrent event times are required for pairwise_common_followup.")
      }

      followup <- .wr_safe_num(data[[recurrent_endpoint$followup_col]])

      for (k in seq_len(n)) {
        ia <- idx_treatment[k]
        ib <- idx_control[k]
        common_followup <- min(followup[ia], followup[ib])

        if (!is.finite(common_followup)) next

        a_count <- .wr_count_sorted_times_until(
          recurrent_endpoint$abs_times_list[[ia]],
          common_followup
        )

        b_count <- .wr_count_sorted_times_until(
          recurrent_endpoint$abs_times_list[[ib]],
          common_followup
        )

        if (b_count > a_count) out[k] <- 1L
        if (a_count > b_count) out[k] <- -1L
      }

      return(out)
    }

    a <- .wr_safe_num(data[[endpoint_spec$count_col]][idx_treatment])
    b <- .wr_safe_num(data[[endpoint_spec$count_col]][idx_control])
    ok <- !is.na(a) & !is.na(b)

    out[ok & (b - a > 0)] <- 1L
    out[ok & (a - b > 0)] <- -1L
    return(out)
  }

  if (endpoint_spec$type == "binary") {
    a <- .wr_safe_num(data[[endpoint_spec$value_col]][idx_treatment])
    b <- .wr_safe_num(data[[endpoint_spec$value_col]][idx_control])
    ok <- !is.na(a) & !is.na(b)

    if (isTRUE(endpoint_spec$adverse)) {
      out[ok & (b - a > 0)] <- 1L
      out[ok & (a - b > 0)] <- -1L
    } else {
      out[ok & (a - b > 0)] <- 1L
      out[ok & (b - a > 0)] <- -1L
    }

    return(out)
  }

  if (endpoint_spec$type == "continuous") {
    a <- .wr_safe_num(data[[endpoint_spec$value_col]][idx_treatment])
    b <- .wr_safe_num(data[[endpoint_spec$value_col]][idx_control])
    ok <- !is.na(a) & !is.na(b)

    if (isTRUE(endpoint_spec$higher_better)) {
      out[ok & (a - b > 0)] <- 1L
      out[ok & (b - a > 0)] <- -1L
    } else {
      out[ok & (b - a > 0)] <- 1L
      out[ok & (a - b > 0)] <- -1L
    }

    return(out)
  }

  stop("Unsupported endpoint type: ", endpoint_spec$type)
}

.wr_counts_for_candidate_R <- function(data,
                                       arm,
                                       endpoint_specs,
                                       recurrent_cache,
                                       order_vec,
                                       weights,
                                       thresholds_by_endpoint,
                                       pairs = NULL) {
  m <- length(endpoint_specs)

  if (is.null(pairs)) {
    pairs <- .wr_get_pair_indices(arm)
  }

  idx_treatment <- pairs$treatment
  idx_control <- pairs$control
  n_pairs <- pairs$n_pairs

  wins_by_rank <- rep(0, m)
  losses_by_rank <- rep(0, m)
  unresolved <- rep(TRUE, n_pairs)

  for (rank in seq_along(order_vec)) {
    endpoint_id <- order_vec[rank]
    endpoint_spec <- endpoint_specs[[endpoint_id]]

    threshold <- if (identical(tolower(endpoint_spec$type), "time")) {
      thresholds_by_endpoint[endpoint_id]
    } else {
      0
    }

    if (is.na(threshold)) threshold <- 0

    recurrent_endpoint <- if (!is.null(recurrent_cache) &&
                              length(recurrent_cache) >= endpoint_id) {
      recurrent_cache[[endpoint_id]]
    } else {
      NULL
    }

    cmp <- .wr_compare_endpoint_pairs_R(
      data = data,
      idx_treatment = idx_treatment,
      idx_control = idx_control,
      endpoint_spec = endpoint_spec,
      threshold = threshold,
      recurrent_endpoint = recurrent_endpoint
    )

    unresolved_idx <- which(unresolved)

    if (length(unresolved_idx) == 0) break

    cmp_unresolved <- cmp[unresolved_idx]

    wins_by_rank[rank] <- sum(
      cmp_unresolved == 1L,
      na.rm = TRUE
    )

    losses_by_rank[rank] <- sum(
      cmp_unresolved == -1L,
      na.rm = TRUE
    )

    resolved_local <- which(cmp_unresolved != 0L)

    if (length(resolved_local) > 0) {
      unresolved[unresolved_idx[resolved_local]] <- FALSE
    }
  }

  tie_count <- sum(unresolved)
  weights <- as.numeric(weights)

  weighted_win <- sum(weights * wins_by_rank)
  weighted_loss <- sum(weights * losses_by_rank)

  wr <- .wr_ratio_safe(
    weighted_win,
    weighted_loss
  )

  wo <- .wr_ratio_safe(
    weighted_win + 0.5 * tie_count,
    weighted_loss + 0.5 * tie_count
  )

  list(
    WR_statistic = wr,
    WO_statistic = wo,
    WR_abslog = .wr_abslog_safe(wr),
    WO_abslog = .wr_abslog_safe(wo),
    weighted_win = weighted_win,
    weighted_loss = weighted_loss,
    tie_count = tie_count,
    tie_proportion = tie_count / n_pairs,
    n_pairs = n_pairs,
    wins_by_rank = wins_by_rank,
    losses_by_rank = losses_by_rank
  )
}

.wr_counts_for_candidate <- function(data,
                                     arm,
                                     endpoint_specs,
                                     recurrent_cache,
                                     cache,
                                     order_vec,
                                     weights,
                                     thresholds_by_endpoint,
                                     pairs) {
  if (exists("cpp_counts_for_candidate_core", mode = "function")) {
    return(cpp_counts_for_candidate_core(
      idx_treatment = as.integer(pairs$treatment),
      idx_control = as.integer(pairs$control),
      order_vec = as.integer(order_vec),
      weights = as.numeric(weights),
      thresholds_by_endpoint = as.numeric(thresholds_by_endpoint),
      type_code = as.integer(cache$type_code),
      time_mat = cache$time_mat,
      event_mat = cache$event_mat,
      value_mat = cache$value_mat,
      followup_mat = cache$followup_mat,
      recurrent_times = cache$recurrent_times,
      recurrent_start = cache$recurrent_start,
      recurrent_len = cache$recurrent_len
    ))
  }

  .wr_counts_for_candidate_R(
    data = data,
    arm = arm,
    endpoint_specs = endpoint_specs,
    recurrent_cache = recurrent_cache,
    order_vec = order_vec,
    weights = weights,
    thresholds_by_endpoint = thresholds_by_endpoint,
    pairs = pairs
  )
}

evaluate_general_candidates <- function(data,
                                        arm,
                                        endpoint_specs,
                                        candidate_grid,
                                        recurrent_cache = NULL) {
  m <- length(endpoint_specs)

  if (!(m %in% c(2L, 3L))) {
    stop("Only 2 or 3 endpoints are supported.")
  }

  if (!is.data.frame(candidate_grid) || nrow(candidate_grid) == 0) {
    stop("candidate_grid must be a non-empty data frame.")
  }

  p_cols <- paste0("p", seq_len(m))
  t_cols <- paste0("threshold_e", seq_len(m))
  required <- c(
    "order_key",
    "candidate_key",
    p_cols,
    t_cols
  )

  missing_cols <- setdiff(
    required,
    names(candidate_grid)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Candidate grid is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  arm <- as.integer(arm)
  pairs <- .wr_get_pair_indices(arm)

  cache <- .wr_prepare_endpoint_cache(
    data = data,
    endpoint_specs = endpoint_specs,
    recurrent_cache = recurrent_cache
  )

  rows <- vector("list", nrow(candidate_grid))

  for (i in seq_len(nrow(candidate_grid))) {
    cg <- candidate_grid[i, , drop = FALSE]

    order_vec <- .wr_parse_order(
      cg$order_key[1]
    )

    weights <- as.numeric(
      cg[1, p_cols, drop = TRUE]
    )

    thresholds <- as.numeric(
      cg[1, t_cols, drop = TRUE]
    )

    counts <- .wr_counts_for_candidate(
      data = data,
      arm = arm,
      endpoint_specs = endpoint_specs,
      recurrent_cache = recurrent_cache,
      cache = cache,
      order_vec = order_vec,
      weights = weights,
      thresholds_by_endpoint = thresholds,
      pairs = pairs
    )

    row <- cg
    row$WR_statistic <- as.numeric(counts$WR_statistic)
    row$WO_statistic <- as.numeric(counts$WO_statistic)
    row$WR_abslog <- as.numeric(counts$WR_abslog)
    row$WO_abslog <- as.numeric(counts$WO_abslog)
    row$weighted_win <- as.numeric(counts$weighted_win)
    row$weighted_loss <- as.numeric(counts$weighted_loss)
    row$tie_count <- as.numeric(counts$tie_count)
    row$tie_proportion <- as.numeric(counts$tie_proportion)
    row$n_pairs <- as.numeric(counts$n_pairs)

    for (rank in seq_len(3L)) {
      row[[paste0("win_rank", rank)]] <- if (rank <= m) {
        as.numeric(counts$wins_by_rank[rank])
      } else {
        NA_real_
      }

      row[[paste0("loss_rank", rank)]] <- if (rank <= m) {
        as.numeric(counts$losses_by_rank[rank])
      } else {
        NA_real_
      }
    }

    rows[[i]] <- row
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
