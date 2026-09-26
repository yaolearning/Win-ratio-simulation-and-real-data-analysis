.wr_select_best_candidate <- function(counts,
                                      candidate_keys,
                                      measure = c("WR", "WO"),
                                      side = c("one", "two")) {
  measure <- match.arg(measure)
  side <- match.arg(side)

  d <- counts[counts$candidate_key %in% candidate_keys, , drop = FALSE]
  if (nrow(d) == 0) stop("No candidate rows found for selection.")

  stat_col <- if (measure == "WR") "WR_statistic" else "WO_statistic"
  abs_col <- if (measure == "WR") "WR_abslog" else "WO_abslog"

  score <- if (side == "one") d[[stat_col]] else d[[abs_col]]
  score[is.na(score)] <- -Inf

  if (length(score) == 0 || all(score == -Inf)) {
    return(d[FALSE, , drop = FALSE])
  }

  d[which.max(score), , drop = FALSE]
}

.wr_get_stat_value <- function(row,
                               measure = c("WR", "WO"),
                               side = c("one", "two")) {
  measure <- match.arg(measure)
  side <- match.arg(side)

  if (nrow(row) == 0) return(NA_real_)
  if (measure == "WR" && side == "one") return(as.numeric(row$WR_statistic[1]))
  if (measure == "WO" && side == "one") return(as.numeric(row$WO_statistic[1]))
  if (measure == "WR" && side == "two") return(as.numeric(row$WR_abslog[1]))
  as.numeric(row$WO_abslog[1])
}

.wr_flatten_selection <- function(row,
                                  measure,
                                  method_id,
                                  side,
                                  m) {
  if (nrow(row) == 0) {
    out <- data.frame(
      measure = measure,
      method_id = method_id,
      side = side,
      statistic = NA_real_,
      abslog_statistic = NA_real_,
      selected_order = NA_character_,
      tie_count = NA_real_,
      tie_proportion = NA_real_,
      stringsAsFactors = FALSE
    )
    for (j in seq_len(m)) {
      out[[paste0("selected_p", j)]] <- NA_real_
      out[[paste0("selected_t_e", j)]] <- NA_real_
    }
    return(out)
  }

  out <- data.frame(
    measure = measure,
    method_id = method_id,
    side = side,
    statistic = if (measure == "WR") row$WR_statistic[1] else row$WO_statistic[1],
    abslog_statistic = if (measure == "WR") row$WR_abslog[1] else row$WO_abslog[1],
    selected_order = as.character(row$order_key[1]),
    tie_count = as.numeric(row$tie_count[1]),
    tie_proportion = as.numeric(row$tie_proportion[1]),
    stringsAsFactors = FALSE
  )

  for (j in seq_len(m)) {
    out[[paste0("selected_p", j)]] <- as.numeric(row[[paste0("p", j)]][1])
    out[[paste0("selected_t_e", j)]] <- as.numeric(row[[paste0("threshold_e", j)]][1])
  }

  out
}

.wr_measures_to_run <- function(control) {
  c(
    if (isTRUE(control$run_wr)) "WR",
    if (isTRUE(control$run_wo)) "WO"
  )
}

.wr_build_observed_combinations <- function(counts,
                                            methods,
                                            method_control,
                                            m) {
  measures <- .wr_measures_to_run(method_control)
  rows <- list()
  idx <- 1L

  for (measure in measures) {
    for (method in methods) {
      keys <- method$grid$candidate_key

      obs_one <- .wr_select_best_candidate(
        counts = counts,
        candidate_keys = keys,
        measure = measure,
        side = "one"
      )

      obs_two <- .wr_select_best_candidate(
        counts = counts,
        candidate_keys = keys,
        measure = measure,
        side = "two"
      )

      rows[[idx]] <- list(
        combo_id = paste(measure, method$id, sep = "__"),
        measure = measure,
        method_id = method$id,
        notation = method$notation,
        description = method$description,
        n_candidates = nrow(method$grid),
        obs_one_stat = .wr_get_stat_value(obs_one, measure, "one"),
        obs_two_abslog = .wr_get_stat_value(obs_two, measure, "two"),
        obs_one_candidate = if (nrow(obs_one) > 0) obs_one$candidate_key[1] else NA_character_,
        obs_two_candidate = if (nrow(obs_two) > 0) obs_two$candidate_key[1] else NA_character_,
        obs_one = .wr_flatten_selection(obs_one, measure, method$id, "one", m),
        obs_two = .wr_flatten_selection(obs_two, measure, method$id, "two", m)
      )

      idx <- idx + 1L
    }
  }

  rows
}
