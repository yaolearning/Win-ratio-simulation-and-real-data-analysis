.wr_results_wide <- function(observed,
                             B) {
  if (nrow(observed) == 0) return(data.frame())

  method_keys <- unique(observed[, c(
    "measure",
    "method_id",
    "notation",
    "description",
    "n_candidates"
  ), drop = FALSE])

  rows <- vector("list", nrow(method_keys))

  for (i in seq_len(nrow(method_keys))) {
    key <- method_keys[i, , drop = FALSE]

    one <- observed[
      observed$measure == key$measure &
        observed$method_id == key$method_id &
        observed$side == "one",
      ,
      drop = FALSE
    ]

    two <- observed[
      observed$measure == key$measure &
        observed$method_id == key$method_id &
        observed$side == "two",
      ,
      drop = FALSE
    ]

    p_cols <- grep(
      "^selected_p[0-9]+$",
      names(observed),
      value = TRUE
    )

    t_cols <- grep(
      "^selected_t_e[0-9]+$",
      names(observed),
      value = TRUE
    )

    row <- data.frame(
      measure = key$measure,
      method_id = key$method_id,
      notation = key$notation,
      description = key$description,
      n_candidates = key$n_candidates,
      B_perm = B,
      statistic_one = if (nrow(one)) one$selected_value[1] else NA_real_,
      p_one = if (nrow(one)) one$adaptive_p_value[1] else NA_real_,
      fixed_selected_p_one = if (nrow(one)) one$fixed_selected_p_value[1] else NA_real_,
      selected_order_one = if (nrow(one)) one$selected_order[1] else NA_character_,
      tie_count_one = if (nrow(one)) one$tie_count[1] else NA_real_,
      tie_proportion_one = if (nrow(one)) one$tie_proportion[1] else NA_real_,
      mean_perm_tie_count_one = if (nrow(one)) one$mean_perm_tie_count[1] else NA_real_,
      mean_perm_tie_proportion_one = if (nrow(one)) one$mean_perm_tie_proportion[1] else NA_real_,
      mean_fixed_perm_tie_count_one = if (nrow(one)) one$mean_fixed_perm_tie_count[1] else NA_real_,
      mean_fixed_perm_tie_proportion_one = if (nrow(one)) one$mean_fixed_perm_tie_proportion[1] else NA_real_,
      statistic_two = if (nrow(two)) two$selected_value[1] else NA_real_,
      abslog_statistic_two = if (nrow(two)) two$selected_abs_log_value[1] else NA_real_,
      two_sided_direction = if (nrow(two)) two$selected_direction[1] else NA_character_,
      p_two = if (nrow(two)) two$adaptive_p_value[1] else NA_real_,
      fixed_selected_p_two = if (nrow(two)) two$fixed_selected_p_value[1] else NA_real_,
      selected_order_two = if (nrow(two)) two$selected_order[1] else NA_character_,
      tie_count_two = if (nrow(two)) two$tie_count[1] else NA_real_,
      tie_proportion_two = if (nrow(two)) two$tie_proportion[1] else NA_real_,
      mean_perm_tie_count_two = if (nrow(two)) two$mean_perm_tie_count[1] else NA_real_,
      mean_perm_tie_proportion_two = if (nrow(two)) two$mean_perm_tie_proportion[1] else NA_real_,
      mean_fixed_perm_tie_count_two = if (nrow(two)) two$mean_fixed_perm_tie_count[1] else NA_real_,
      mean_fixed_perm_tie_proportion_two = if (nrow(two)) two$mean_fixed_perm_tie_proportion[1] else NA_real_,
      stringsAsFactors = FALSE
    )

    for (nm in p_cols) {
      row[[paste0(nm, "_one")]] <- if (nrow(one)) one[[nm]][1] else NA_real_
      row[[paste0(nm, "_two")]] <- if (nrow(two)) two[[nm]][1] else NA_real_
    }

    for (nm in t_cols) {
      row[[paste0(nm, "_one")]] <- if (nrow(one)) one[[nm]][1] else NA_real_
      row[[paste0(nm, "_two")]] <- if (nrow(two)) two[[nm]][1] else NA_real_
    }

    rows[[i]] <- row
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.wr_endpoint_summary <- function(prepared) {
  if (!is.list(prepared) ||
      !all(c("data", "endpoints") %in% names(prepared))) {
    return(data.frame())
  }

  data <- prepared$data
  rows <- vector("list", length(prepared$endpoints))

  for (j in seq_along(prepared$endpoints)) {
    ep <- prepared$endpoints[[j]]

    if (ep$type == "time") {
      values <- data[[ep$time_col]]
      event <- data[[ep$event_col]]

      rows[[j]] <- data.frame(
        endpoint = j,
        name = ep$name,
        type = ep$type,
        comparison = NA_character_,
        n_events_treatment = sum(
          data$arm == 1 & event == 1,
          na.rm = TRUE
        ),
        n_events_control = sum(
          data$arm == 0 & event == 1,
          na.rm = TRUE
        ),
        mean_treatment = mean(
          values[data$arm == 1],
          na.rm = TRUE
        ),
        mean_control = mean(
          values[data$arm == 0],
          na.rm = TRUE
        ),
        median_treatment = stats::median(
          values[data$arm == 1],
          na.rm = TRUE
        ),
        median_control = stats::median(
          values[data$arm == 0],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    } else {
      value_col <- if (ep$type == "count") {
        ep$count_col
      } else {
        ep$value_col
      }

      values <- data[[value_col]]

      rows[[j]] <- data.frame(
        endpoint = j,
        name = ep$name,
        type = ep$type,
        comparison = if (ep$type == "count") ep$comparison else NA_character_,
        n_events_treatment = if (ep$type %in% c("count", "binary")) {
          sum(
            data$arm == 1 & values > 0,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
        n_events_control = if (ep$type %in% c("count", "binary")) {
          sum(
            data$arm == 0 & values > 0,
            na.rm = TRUE
          )
        } else {
          NA_real_
        },
        mean_treatment = mean(
          values[data$arm == 1],
          na.rm = TRUE
        ),
        mean_control = mean(
          values[data$arm == 0],
          na.rm = TRUE
        ),
        median_treatment = stats::median(
          values[data$arm == 1],
          na.rm = TRUE
        ),
        median_control = stats::median(
          values[data$arm == 0],
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.wr_subject_summary <- function(prepared_general,
                                backend,
                                B,
                                threshold_object,
                                weight_grid) {
  data <- prepared_general$data

  recurrent_ids <- which(vapply(
    prepared_general$endpoints,
    function(ep) {
      identical(ep$type, "count") &&
        identical(ep$comparison, "pairwise_common_followup")
    },
    logical(1)
  ))

  data.frame(
    n_subjects = nrow(data),
    n_treatment = sum(data$arm == 1, na.rm = TRUE),
    n_control = sum(data$arm == 0, na.rm = TRUE),
    n_endpoints = length(prepared_general$endpoints),
    recurrent_common_followup_endpoint_ids = paste(
      recurrent_ids,
      collapse = ","
    ),
    backend = backend,
    B_perm = as.integer(B),
    threshold_time_endpoint_ids = paste(
      threshold_object$active_ids,
      collapse = ","
    ),
    n_threshold_candidate_combinations = nrow(threshold_object$grid),
    n_weight_candidates = nrow(weight_grid),
    stringsAsFactors = FALSE
  )
}

.wr_tie_summary <- function(x) {
  if (!inherits(x, "win_analysis")) stop("x must be a win_analysis object.")
  if (nrow(x$selection) == 0) return(data.frame())

  keep <- c(
    "measure",
    "method_id",
    "notation",
    "side",
    "tie_count",
    "tie_proportion",
    "mean_perm_tie_count",
    "mean_perm_tie_proportion",
    "mean_fixed_perm_tie_count",
    "mean_fixed_perm_tie_proportion"
  )

  x$selection[, keep, drop = FALSE]
}

.wr_settings_table <- function(x) {
  if (!inherits(x, "win_analysis")) stop("x must be a win_analysis object.")

  recurrent_ids <- which(vapply(
    x$endpoint_specs,
    function(ep) {
      identical(ep$type, "count") &&
        identical(ep$comparison, "pairwise_common_followup")
    },
    logical(1)
  ))

  data.frame(
    setting = c(
      "backend",
      "number_of_endpoints",
      "permutation_B",
      "WR_one_sided_statistic",
      "WR_two_sided_statistic",
      "WO_one_sided_statistic",
      "WO_two_sided_statistic",
      "adaptive_permutation_rule",
      "recurrent_common_followup_endpoints",
      "threshold_rule",
      "weight_candidates"
    ),
    value = c(
      x$backend,
      as.character(length(x$endpoint_specs)),
      as.character(x$B),
      "WR; larger values favor treatment",
      "abs(log(WR))",
      "WO; larger values favor treatment",
      "abs(log(WO))",
      "Full order/weight/threshold selection is repeated inside every treatment-label permutation",
      if (length(recurrent_ids) == 0) "none" else paste(recurrent_ids, collapse = ","),
      "Time-endpoint threshold comparison uses strict |Delta time| > t",
      paste(
        apply(
          as.data.frame(x$weight_grid),
          1,
          function(z) .wr_make_weight_string(as.numeric(z))
        ),
        collapse = "; "
      )
    ),
    stringsAsFactors = FALSE
  )
}

.wr_output_readme <- function(x) {
  c(
    "Win ratio / win odds package output",
    "",
    "method_results.csv: one- and two-sided WR/WO results by method.",
    "observed_selection.csv: observed selected candidate for every measure/method/side.",
    "observed_candidates.csv: all unique candidates evaluated in the observed data.",
    "all_candidates.csv: candidate definitions before observed-data evaluation.",
    "permutation_selected.csv: selected adaptive candidate in every permutation.",
    "method_overview.csv: method names, notation, descriptions, and candidate counts.",
    "threshold_info.csv: endpoint-level threshold configuration.",
    "threshold_grid.csv: threshold candidate combinations.",
    "weight_grid.csv: weight candidates.",
    "endpoint_summary.csv: descriptive endpoint summary.",
    "subject_summary.csv: sample size and analysis-engine summary.",
    "tie_summary.csv: observed and permutation-average tie counts/proportions.",
    "logrank_results.csv: time-endpoint log-rank comparator results.",
    "settings.csv: analysis settings recorded in a flat table.",
    "full_results.rds: complete analysis object.",
    "figures/: p-value, tie, and Kaplan-Meier figures.",
    "",
    paste0("Backend: ", x$backend),
    paste0("Permutation B: ", x$B)
  )
}

write_win_results <- function(x,
                              output_dir,
                              prefix = "WIN",
                              save_plots = TRUE) {
  if (!inherits(x, "win_analysis")) {
    stop("x must be a win_analysis object.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  write_csv <- function(obj, suffix) {
    utils::write.csv(
      obj,
      file.path(
        output_dir,
        paste0(prefix, "_", suffix, ".csv")
      ),
      row.names = FALSE
    )
  }

  write_csv(
    x$results,
    "method_results"
  )

  write_csv(
    x$selection,
    "observed_selection"
  )

  write_csv(
    x$observed_candidates,
    "observed_candidates"
  )

  if (!is.null(x$all_candidates)) {
    write_csv(
      x$all_candidates,
      "all_candidates"
    )
  }

  if (nrow(x$permutation_selected) > 0) {
    write_csv(
      x$permutation_selected,
      "permutation_selected"
    )
  }

  write_csv(
    x$method_overview,
    "method_overview"
  )

  write_csv(
    x$threshold_info,
    "threshold_info"
  )

  write_csv(
    x$threshold_grid,
    "threshold_grid"
  )

  write_csv(
    x$weight_grid,
    "weight_grid"
  )

  if (nrow(x$endpoint_summary) > 0) {
    write_csv(
      x$endpoint_summary,
      "endpoint_summary"
    )
  }

  if (!is.null(x$subject_summary) &&
      nrow(x$subject_summary) > 0) {
    write_csv(
      x$subject_summary,
      "subject_summary"
    )
  }

  tie_summary <- .wr_tie_summary(x)

  if (nrow(tie_summary) > 0) {
    write_csv(
      tie_summary,
      "tie_summary"
    )
  }

  if (nrow(x$logrank) > 0) {
    write_csv(
      x$logrank,
      "logrank_results"
    )
  }

  write_csv(
    .wr_settings_table(x),
    "settings"
  )

  writeLines(
    .wr_output_readme(x),
    con = file.path(
      output_dir,
      paste0(prefix, "_OUTPUT_README.txt")
    )
  )

  saveRDS(
    x,
    file.path(
      output_dir,
      paste0(prefix, "_full_results.rds")
    )
  )

  if (isTRUE(save_plots) &&
      exists("save_win_plots", mode = "function")) {
    save_win_plots(
      x = x,
      output_dir = file.path(
        output_dir,
        "figures"
      ),
      prefix = prefix
    )
  }

  manifest <- data.frame(
    file = sort(
      list.files(
        output_dir,
        recursive = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    manifest,
    file.path(
      output_dir,
      paste0(prefix, "_output_manifest.csv")
    ),
    row.names = FALSE
  )

  invisible(output_dir)
}

.wr_simulation_tie_summary <- function(x) {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  d <- x$trial_results
  if (nrow(d) == 0) return(data.frame())

  keys <- unique(d[, c(
    "scenario_index",
    "scenario_id",
    "measure",
    "method_id",
    "notation"
  ), drop = FALSE])

  rows <- vector("list", nrow(keys))

  for (i in seq_len(nrow(keys))) {
    k <- keys[i, , drop = FALSE]

    z <- d[
      d$scenario_id == k$scenario_id &
        d$measure == k$measure &
        d$method_id == k$method_id,
      ,
      drop = FALSE
    ]

    rows[[i]] <- data.frame(
      scenario_index = k$scenario_index,
      scenario_id = k$scenario_id,
      measure = k$measure,
      method_id = k$method_id,
      notation = k$notation,
      mean_observed_tie_proportion_one = mean(
        z$tie_proportion_one,
        na.rm = TRUE
      ),
      mean_observed_tie_proportion_two = mean(
        z$tie_proportion_two,
        na.rm = TRUE
      ),
      mean_perm_tie_proportion_one = mean(
        z$mean_perm_tie_proportion_one,
        na.rm = TRUE
      ),
      mean_perm_tie_proportion_two = mean(
        z$mean_perm_tie_proportion_two,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

write_win_simulation_results <- function(x,
                                         output_dir,
                                         prefix = "GLOBAL",
                                         save_plots = TRUE,
                                         save_scenario_tables = TRUE) {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::write.csv(
    x$scenarios,
    file.path(output_dir, paste0(prefix, "_scenario_config.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    x$trial_results,
    file.path(output_dir, paste0(prefix, "_method_trial_results.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    x$logrank_results,
    file.path(output_dir, paste0(prefix, "_logrank_trial_results.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    x$selection_results,
    file.path(output_dir, paste0(prefix, "_selection_results.csv")),
    row.names = FALSE
  )

  utils::write.csv(
    x$power,
    file.path(output_dir, paste0(prefix, "_power_summary.csv")),
    row.names = FALSE
  )

  tie_summary <- .wr_simulation_tie_summary(x)

  utils::write.csv(
    tie_summary,
    file.path(output_dir, paste0(prefix, "_tie_summary.csv")),
    row.names = FALSE
  )

  settings <- data.frame(
    setting = c(
      "NSIM",
      "B_PERM",
      "MASTER_SEED",
      "ALPHA",
      "number_of_scenarios",
      "WR_two_sided_statistic",
      "WO_two_sided_statistic",
      "permutation_rule",
      "recurrent_count_rule"
    ),
    value = c(
      as.character(x$settings$nsim),
      as.character(x$settings$B),
      as.character(x$settings$seed),
      as.character(x$settings$alpha),
      as.character(nrow(x$scenarios)),
      "abs(log(WR))",
      "abs(log(WO))",
      "Adaptive order/weight/threshold selection is repeated inside every treatment-label permutation",
      "Hospitalization counts are compared within each treatment-control pair's common follow-up time"
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    settings,
    file.path(output_dir, paste0(prefix, "_settings.csv")),
    row.names = FALSE
  )

  saveRDS(
    x,
    file.path(output_dir, paste0(prefix, "_full_simulation_results.rds"))
  )

  if (isTRUE(save_scenario_tables)) {
    for (i in seq_len(nrow(x$scenarios))) {
      scenario_id <- x$scenarios$scenario_id[i]
      scenario_index <- x$scenarios$scenario_index[i]

      scenario_dir <- file.path(
        output_dir,
        paste0(
          sprintf("%02d", scenario_index),
          "_",
          scenario_id
        )
      )

      dir.create(
        scenario_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      utils::write.csv(
        x$trial_results[
          x$trial_results$scenario_id == scenario_id,
          ,
          drop = FALSE
        ],
        file.path(scenario_dir, "method_trial_results.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        x$selection_results[
          x$selection_results$scenario_id == scenario_id,
          ,
          drop = FALSE
        ],
        file.path(scenario_dir, "selection_results.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        x$logrank_results[
          x$logrank_results$scenario_id == scenario_id,
          ,
          drop = FALSE
        ],
        file.path(scenario_dir, "logrank_trial_results.csv"),
        row.names = FALSE
      )

      utils::write.csv(
        x$power[
          x$power$scenario_id == scenario_id,
          ,
          drop = FALSE
        ],
        file.path(scenario_dir, "power_summary.csv"),
        row.names = FALSE
      )
    }
  }

  if (isTRUE(save_plots) &&
      exists("save_win_simulation_plots", mode = "function")) {
    save_win_simulation_plots(
      x = x,
      output_dir = file.path(output_dir, "figures"),
      prefix = prefix
    )
  }

  manifest <- data.frame(
    file = sort(
      list.files(
        output_dir,
        recursive = TRUE
      )
    ),
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    manifest,
    file.path(output_dir, paste0(prefix, "_output_manifest.csv")),
    row.names = FALSE
  )

  invisible(output_dir)
}

print.win_analysis <- function(x, ...) {
  cat("Win ratio / win odds analysis\n")
  cat("Backend:", x$backend, "\n")
  cat("Endpoints:", length(x$endpoint_specs), "\n")
  cat("Permutation B:", x$B, "\n")
  cat("Methods:", nrow(x$method_overview), "\n")
  invisible(x)
}

summary.win_analysis <- function(object, ...) {
  object$results
}

as.data.frame.win_analysis <- function(x, ...) {
  x$results
}
