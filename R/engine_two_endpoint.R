.wr_safe_wr_eps <- function(win_score,
                            loss_score,
                            total_pairs,
                            eps = 1e-8) {
  ((win_score / total_pairs) + eps) / ((loss_score / total_pairs) + eps)
}

.wr_count_hosp_until_R <- function(times, t) {
  if (length(times) == 0) return(0L)
  sum(times <= t)
}

.wr_two_endpoint_base_counts_R <- function(ds) {
  ds <- .wr_prepare_recurrent_cache(ds)
  tab <- ds$table.output
  treatment <- which(tab$ARM == 1)
  control <- which(tab$ARM == 0)

  if (length(treatment) == 0 || length(control) == 0) {
    stop("Both treatment and control groups must contain at least one subject.")
  }

  D1_win <- D1_loss <- H2_win <- H2_loss <- 0
  H1_win <- H1_loss <- D2_win <- D2_loss <- 0

  for (i in treatment) {
    fu1 <- tab$FUTIME[i]
    c1_original <- tab$CNSR[i]

    for (j in control) {
      fu2 <- tab$FUTIME[j]
      c2_original <- tab$CNSR[j]
      common_followup <- min(fu1, fu2)

      hosp1 <- .wr_count_hosp_until_R(ds$hosp.abs.times.list[[i]], common_followup)
      hosp2 <- .wr_count_hosp_until_R(ds$hosp.abs.times.list[[j]], common_followup)

      hosp_sign <- 0L
      if (hosp2 > hosp1) hosp_sign <- 1L
      if (hosp2 < hosp1) hosp_sign <- -1L

      c1 <- c1_original
      c2 <- c2_original
      if (fu1 < fu2) c2 <- 0L
      if (fu2 < fu1) c1 <- 0L

      death_sign <- 0L
      if (c1 == 0 && c2 == 1) death_sign <- 1L
      if (c1 == 1 && c2 == 0) death_sign <- -1L

      if (death_sign > 0) {
        D1_win <- D1_win + 1
      } else if (death_sign < 0) {
        D1_loss <- D1_loss + 1
      } else if (hosp_sign > 0) {
        H2_win <- H2_win + 1
      } else if (hosp_sign < 0) {
        H2_loss <- H2_loss + 1
      }

      if (hosp_sign > 0) {
        H1_win <- H1_win + 1
      } else if (hosp_sign < 0) {
        H1_loss <- H1_loss + 1
      } else if (death_sign > 0) {
        D2_win <- D2_win + 1
      } else if (death_sign < 0) {
        D2_loss <- D2_loss + 1
      }
    }
  }

  list(
    total_pairs = length(treatment) * length(control),
    D1_win = D1_win,
    D1_loss = D1_loss,
    H2_win = H2_win,
    H2_loss = H2_loss,
    H1_win = H1_win,
    H1_loss = H1_loss,
    D2_win = D2_win,
    D2_loss = D2_loss
  )
}

.wr_two_endpoint_threshold_counts_R <- function(ds, t_grid) {
  ds <- .wr_prepare_recurrent_cache(ds)
  tab <- ds$table.output
  treatment <- which(tab$ARM == 1)
  control <- which(tab$ARM == 0)

  if (length(treatment) == 0 || length(control) == 0) {
    stop("Both treatment and control groups must contain at least one subject.")
  }

  t_grid <- as.numeric(t_grid)
  k <- length(t_grid)

  df_first_win <- df_first_loss <- df_second_win <- df_second_loss <- rep(0, k)
  hf_first_win <- hf_first_loss <- hf_second_win <- hf_second_loss <- rep(0, k)

  for (i in treatment) {
    fu1 <- tab$FUTIME[i]
    c1_original <- tab$CNSR[i]

    for (j in control) {
      fu2 <- tab$FUTIME[j]
      c2_original <- tab$CNSR[j]
      common_followup <- min(fu1, fu2)

      hosp1 <- .wr_count_hosp_until_R(ds$hosp.abs.times.list[[i]], common_followup)
      hosp2 <- .wr_count_hosp_until_R(ds$hosp.abs.times.list[[j]], common_followup)

      hosp_sign <- 0L
      if (hosp2 > hosp1) hosp_sign <- 1L
      if (hosp2 < hosp1) hosp_sign <- -1L

      c1 <- c1_original
      c2 <- c2_original
      if (fu1 < fu2) c2 <- 0L
      if (fu2 < fu1) c1 <- 0L

      death_sign <- 0L
      if (c1 == 0 && c2 == 1) death_sign <- 1L
      if (c1 == 1 && c2 == 0) death_sign <- -1L

      abs_diff <- abs(fu1 - fu2)

      for (kk in seq_along(t_grid)) {
        death_threshold_sign <- 0L
        if (death_sign != 0L && abs_diff > t_grid[kk]) {
          death_threshold_sign <- death_sign
        }

        if (death_threshold_sign > 0) {
          df_first_win[kk] <- df_first_win[kk] + 1
        } else if (death_threshold_sign < 0) {
          df_first_loss[kk] <- df_first_loss[kk] + 1
        } else if (hosp_sign > 0) {
          df_second_win[kk] <- df_second_win[kk] + 1
        } else if (hosp_sign < 0) {
          df_second_loss[kk] <- df_second_loss[kk] + 1
        }

        if (hosp_sign > 0) {
          hf_first_win[kk] <- hf_first_win[kk] + 1
        } else if (hosp_sign < 0) {
          hf_first_loss[kk] <- hf_first_loss[kk] + 1
        } else if (death_threshold_sign > 0) {
          hf_second_win[kk] <- hf_second_win[kk] + 1
        } else if (death_threshold_sign < 0) {
          hf_second_loss[kk] <- hf_second_loss[kk] + 1
        }
      }
    }
  }

  list(
    total_pairs = length(treatment) * length(control),
    t_grid = t_grid,
    df_first_win = df_first_win,
    df_first_loss = df_first_loss,
    df_second_win = df_second_win,
    df_second_loss = df_second_loss,
    hf_first_win = hf_first_win,
    hf_first_loss = hf_first_loss,
    hf_second_win = hf_second_win,
    hf_second_loss = hf_second_loss
  )
}

.wr_two_endpoint_base_counts <- function(ds) {
  ds <- .wr_prepare_recurrent_cache(ds)
  tab <- ds$table.output

  if (exists("cpp_two_endpoint_base_counts", mode = "function")) {
    return(cpp_two_endpoint_base_counts(
      futime = as.numeric(tab$FUTIME),
      cnsr = as.integer(tab$CNSR),
      arm = as.integer(tab$ARM),
      hosp_times = as.numeric(ds$hosp.flat),
      hosp_start = as.integer(ds$hosp.start),
      hosp_len = as.integer(ds$hosp.len)
    ))
  }

  .wr_two_endpoint_base_counts_R(ds)
}

.wr_two_endpoint_threshold_counts <- function(ds, t_grid) {
  ds <- .wr_prepare_recurrent_cache(ds)
  tab <- ds$table.output

  if (exists("cpp_two_endpoint_threshold_counts", mode = "function")) {
    return(cpp_two_endpoint_threshold_counts(
      futime = as.numeric(tab$FUTIME),
      cnsr = as.integer(tab$CNSR),
      arm = as.integer(tab$ARM),
      hosp_times = as.numeric(ds$hosp.flat),
      hosp_start = as.integer(ds$hosp.start),
      hosp_len = as.integer(ds$hosp.len),
      t_grid = as.numeric(t_grid)
    ))
  }

  .wr_two_endpoint_threshold_counts_R(ds, t_grid)
}

.wr_two_endpoint_candidate_row <- function(candidate,
                                           first_win,
                                           first_loss,
                                           second_win,
                                           second_loss,
                                           total_pairs,
                                           eps) {
  p <- as.numeric(candidate$p1[1])
  win_score <- p * first_win + (1 - p) * second_win
  loss_score <- p * first_loss + (1 - p) * second_loss
  win_pairs <- first_win + second_win
  loss_pairs <- first_loss + second_loss
  tie_count <- total_pairs - win_pairs - loss_pairs

  wr <- .wr_safe_wr_eps(
    win_score = win_score,
    loss_score = loss_score,
    total_pairs = total_pairs,
    eps = eps
  )

  wo <- .wr_ratio_safe(
    win_score + 0.5 * tie_count,
    loss_score + 0.5 * tie_count
  )

  out <- candidate
  out$WR_statistic <- wr
  out$WO_statistic <- wo
  out$WR_abslog <- .wr_abslog_safe(wr)
  out$WO_abslog <- .wr_abslog_safe(wo)
  out$weighted_win <- win_score
  out$weighted_loss <- loss_score
  out$tie_count <- tie_count
  out$tie_proportion <- tie_count / total_pairs
  out$n_pairs <- total_pairs
  out$win_rank1 <- first_win
  out$win_rank2 <- second_win
  out$win_rank3 <- NA_real_
  out$loss_rank1 <- first_loss
  out$loss_rank2 <- second_loss
  out$loss_rank3 <- NA_real_
  out
}

evaluate_two_endpoint_recurrent_candidates <- function(ds,
                                                       candidate_grid,
                                                       eps = 1e-8) {
  if (is.null(ds$hosp.times.list)) {
    stop("Recurrent event times are required for pairwise common-follow-up comparison.")
  }

  required <- c(
    "order_key", "candidate_key",
    "p1", "p2",
    "threshold_e1", "threshold_e2"
  )

  missing_cols <- setdiff(required, names(candidate_grid))
  if (length(missing_cols) > 0) {
    stop("Candidate grid is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(abs(candidate_grid$threshold_e2) > 1e-12, na.rm = TRUE)) {
    stop("The specialized two-endpoint recurrent engine applies thresholds only to endpoint 1.")
  }

  base <- .wr_two_endpoint_base_counts(ds)
  positive_thresholds <- sort(unique(as.numeric(candidate_grid$threshold_e1)))
  positive_thresholds <- positive_thresholds[is.finite(positive_thresholds) & positive_thresholds > 0]

  threshold_counts <- if (length(positive_thresholds) > 0) {
    .wr_two_endpoint_threshold_counts(ds, positive_thresholds)
  } else {
    NULL
  }

  rows <- vector("list", nrow(candidate_grid))
  total_pairs <- as.numeric(base$total_pairs)

  for (i in seq_len(nrow(candidate_grid))) {
    candidate <- candidate_grid[i, , drop = FALSE]
    order_key <- as.character(candidate$order_key[1])
    threshold <- as.numeric(candidate$threshold_e1[1])

    if (!is.finite(threshold) || threshold <= 0) {
      if (order_key == "1->2") {
        first_win <- as.numeric(base$D1_win)
        first_loss <- as.numeric(base$D1_loss)
        second_win <- as.numeric(base$H2_win)
        second_loss <- as.numeric(base$H2_loss)
      } else if (order_key == "2->1") {
        first_win <- as.numeric(base$H1_win)
        first_loss <- as.numeric(base$H1_loss)
        second_win <- as.numeric(base$D2_win)
        second_loss <- as.numeric(base$D2_loss)
      } else {
        stop("Specialized two-endpoint recurrent engine supports orders 1->2 and 2->1 only.")
      }
    } else {
      kk <- match(threshold, as.numeric(threshold_counts$t_grid))
      if (is.na(kk)) stop("Threshold was not found in the evaluated threshold grid.")

      if (order_key == "1->2") {
        first_win <- as.numeric(threshold_counts$df_first_win[kk])
        first_loss <- as.numeric(threshold_counts$df_first_loss[kk])
        second_win <- as.numeric(threshold_counts$df_second_win[kk])
        second_loss <- as.numeric(threshold_counts$df_second_loss[kk])
      } else if (order_key == "2->1") {
        first_win <- as.numeric(threshold_counts$hf_first_win[kk])
        first_loss <- as.numeric(threshold_counts$hf_first_loss[kk])
        second_win <- as.numeric(threshold_counts$hf_second_win[kk])
        second_loss <- as.numeric(threshold_counts$hf_second_loss[kk])
      } else {
        stop("Specialized two-endpoint recurrent engine supports orders 1->2 and 2->1 only.")
      }
    }

    rows[[i]] <- .wr_two_endpoint_candidate_row(
      candidate = candidate,
      first_win = first_win,
      first_loss = first_loss,
      second_win = second_win,
      second_loss = second_loss,
      total_pairs = total_pairs,
      eps = eps
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
