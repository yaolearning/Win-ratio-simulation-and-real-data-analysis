.wr_simulation_endpoints <- function() {
  list(
    endpoint_time(
      name = "Death",
      time = "FUTIME",
      event = "CNSR",
      unit = "years"
    ),
    endpoint_count(
      name = "Hospitalization",
      count = "NUMHOSP",
      comparison = "pairwise_common_followup",
      recurrent_time = "HOSPTIME",
      recurrent_id = "SUBJID",
      unit = "years"
    )
  )
}

.wr_simulation_trial <- function(scenario_row,
                                 sim_index,
                                 B,
                                 master_seed,
                                 alpha) {
  simulated <- generate_win_scenario_dataset(
    scenario = scenario_row,
    scenario_grid = scenario_row,
    sim_index = sim_index,
    seed = master_seed
  )

  threshold_values <- .wr_simulation_threshold_values(
    simulated$subjects
  )

  input <- win_data(
    data = simulated$subjects,
    recurrent_data = simulated$recurrent_events
  )

  fit <- win_analysis(
    data = input,
    endpoints = .wr_simulation_endpoints(),
    format = "wide",
    id = "SUBJID",
    treatment = "ARM",
    threshold = threshold_control(
      time_endpoints = 1L,
      candidates = threshold_values,
      unit = "years",
      data_probs = numeric(0),
      max_event_times = 800L,
      max_thresholds = max(10L, length(threshold_values)),
      seed = master_seed
    ),
    weight = weight_control(
      mode = "vertices",
      constraint = "ordered"
    ),
    methods = method_control(
      run_wr = TRUE,
      run_wo = TRUE,
      report_all_fixed_orders = TRUE
    ),
    permutation = permutation_control(
      enabled = TRUE,
      B = B,
      seed = simulated$trial_seed + 100000L
    ),
    engine = engine_control(
      backend = "two_endpoint_recurrent"
    ),
    run_logrank = TRUE,
    run_composite_logrank_if_possible = FALSE
  )

  results <- fit$results
  results$scenario_index <- scenario_row$scenario_index[1]
  results$scenario_id <- scenario_row$scenario_id[1]
  results$scenario_type <- scenario_row$scenario_type[1]
  results$description <- scenario_row$description[1]
  results$sim_index <- sim_index
  results$trial_seed <- simulated$trial_seed
  results$reject_one <- results$p_one < alpha
  results$reject_two <- results$p_two < alpha

  logrank <- fit$logrank

  if (nrow(logrank) > 0) {
    logrank$scenario_index <- scenario_row$scenario_index[1]
    logrank$scenario_id <- scenario_row$scenario_id[1]
    logrank$scenario_type <- scenario_row$scenario_type[1]
    logrank$description <- scenario_row$description[1]
    logrank$sim_index <- sim_index
    logrank$trial_seed <- simulated$trial_seed
    logrank$reject_one <- logrank$p_one_sided_benefit < alpha
    logrank$reject_two <- logrank$p_two_sided < alpha
  }

  selection <- fit$selection
  selection$scenario_index <- scenario_row$scenario_index[1]
  selection$scenario_id <- scenario_row$scenario_id[1]
  selection$sim_index <- sim_index

  list(
    results = results,
    logrank = logrank,
    selection = selection
  )
}

.wr_simulation_power_summary <- function(results,
                                         logrank,
                                         alpha) {
  rows <- list()
  idx <- 1L

  if (nrow(results) > 0) {
    keys <- unique(results[, c(
      "scenario_index",
      "scenario_id",
      "scenario_type",
      "description",
      "measure",
      "method_id",
      "notation"
    ), drop = FALSE])

    for (i in seq_len(nrow(keys))) {
      k <- keys[i, , drop = FALSE]

      d <- results[
        results$scenario_id == k$scenario_id &
          results$measure == k$measure &
          results$method_id == k$method_id,
        ,
        drop = FALSE
      ]

      rows[[idx]] <- data.frame(
        scenario_index = k$scenario_index,
        scenario_id = k$scenario_id,
        scenario_type = k$scenario_type,
        description = k$description,
        measure = k$measure,
        method_id = k$method_id,
        notation = k$notation,
        nsim_available = nrow(d),
        rejection_proportion_one_sided = mean(d$p_one < alpha, na.rm = TRUE),
        rejection_proportion_two_sided = mean(d$p_two < alpha, na.rm = TRUE),
        mean_p_one_sided = mean(d$p_one, na.rm = TRUE),
        mean_p_two_sided = mean(d$p_two, na.rm = TRUE),
        mean_fixed_p_one_sided = mean(d$fixed_selected_p_one, na.rm = TRUE),
        mean_fixed_p_two_sided = mean(d$fixed_selected_p_two, na.rm = TRUE),
        mean_tie_proportion_one_sided = mean(d$tie_proportion_one, na.rm = TRUE),
        mean_tie_proportion_two_sided = mean(d$tie_proportion_two, na.rm = TRUE),
        alpha = alpha,
        stringsAsFactors = FALSE
      )

      idx <- idx + 1L
    }
  }

  if (nrow(logrank) > 0) {
    scenario_ids <- unique(logrank$scenario_id)

    for (scenario_id in scenario_ids) {
      d <- logrank[logrank$scenario_id == scenario_id, , drop = FALSE]
      d <- d[d$test == "Log-rank outcome 1", , drop = FALSE]

      if (nrow(d) > 0) {
        rows[[idx]] <- data.frame(
          scenario_index = d$scenario_index[1],
          scenario_id = d$scenario_id[1],
          scenario_type = d$scenario_type[1],
          description = d$description[1],
          measure = "Log-rank",
          method_id = "logrank_outcome1",
          notation = "Log-rank outcome 1",
          nsim_available = nrow(d),
          rejection_proportion_one_sided = mean(d$p_one_sided_benefit < alpha, na.rm = TRUE),
          rejection_proportion_two_sided = mean(d$p_two_sided < alpha, na.rm = TRUE),
          mean_p_one_sided = mean(d$p_one_sided_benefit, na.rm = TRUE),
          mean_p_two_sided = mean(d$p_two_sided, na.rm = TRUE),
          mean_fixed_p_one_sided = NA_real_,
          mean_fixed_p_two_sided = NA_real_,
          mean_tie_proportion_one_sided = NA_real_,
          mean_tie_proportion_two_sided = NA_real_,
          alpha = alpha,
          stringsAsFactors = FALSE
        )

        idx <- idx + 1L
      }
    }
  }

  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Run a repeated MaxWin simulation study.
#' @param scenarios Scenario indices, IDs, a scenario data frame, or `NULL` for all scenarios.
#' @param scenario_grid Scenario definition data frame.
#' @param nsim Number of simulated trials per scenario.
#' @param B Number of treatment-label permutations per simulated trial.
#' @param seed Master random seed.
#' @param alpha Significance level used for empirical rejection proportions.
#' @param output_dir Optional output directory.
#' @param checkpoint_every Frequency of simulation-level checkpoint saves.
#' @param resume Whether to resume a simulation-level checkpoint if available.
#' @param verbose Whether to print progress.
#' @param save_plots Whether to save simulation figures when `output_dir` is supplied.
#' @param save_scenario_tables Whether to save scenario-specific tables.
#' @return A `win_simulation` object.
#' @export
run_win_simulation <- function(scenarios = NULL,
                               scenario_grid = default_win_scenarios(),
                               nsim = 1000L,
                               B = 500L,
                               seed = 2026L,
                               alpha = 0.05,
                               output_dir = NULL,
                               checkpoint_every = 25L,
                               resume = FALSE,
                               verbose = TRUE,
                               save_plots = TRUE,
                               save_scenario_tables = TRUE) {
  nsim <- as.integer(nsim)
  B <- as.integer(B)
  checkpoint_every <- as.integer(checkpoint_every)

  if (nsim < 1L) stop("nsim must be at least 1.")
  if (B < 1L) stop("B must be at least 1.")
  if (checkpoint_every < 1L) stop("checkpoint_every must be at least 1.")

  scenario_set <- .wr_resolve_scenario_set(
    scenarios = scenarios,
    scenario_grid = scenario_grid
  )

  all_results <- data.frame()
  all_logrank <- data.frame()
  all_selection <- data.frame()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
      scenario_set,
      file.path(output_dir, "simulation_scenarios.csv"),
      row.names = FALSE
    )
  }

  for (i in seq_len(nrow(scenario_set))) {
    scenario_row <- scenario_set[i, , drop = FALSE]

    if (isTRUE(verbose)) {
      cat(
        "\nScenario",
        scenario_row$scenario_index,
        ":",
        scenario_row$scenario_id,
        "\n"
      )
    }

    start_sim <- 1L

    if (!is.null(output_dir)) {
      scenario_dir <- file.path(
        output_dir,
        paste0(
          sprintf("%02d", scenario_row$scenario_index),
          "_",
          scenario_row$scenario_id
        )
      )

      dir.create(
        scenario_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      checkpoint_file <- file.path(
        scenario_dir,
        "simulation_checkpoint.rds"
      )

      if (isTRUE(resume) && file.exists(checkpoint_file)) {
        checkpoint <- readRDS(checkpoint_file)
        all_results <- .wr_rbind_fill(all_results, checkpoint$results)
        all_logrank <- .wr_rbind_fill(all_logrank, checkpoint$logrank)
        all_selection <- .wr_rbind_fill(all_selection, checkpoint$selection)
        start_sim <- checkpoint$next_sim
      }
    } else {
      scenario_dir <- NULL
      checkpoint_file <- NULL
    }

    if (start_sim <= nsim) {
      for (sim_index in start_sim:nsim) {
        if (isTRUE(verbose) &&
            (sim_index == start_sim ||
             sim_index == nsim ||
             sim_index %% checkpoint_every == 0)) {
          cat("  simulated trial", sim_index, "of", nsim, "\n")
        }

        trial <- .wr_simulation_trial(
          scenario_row = scenario_row,
          sim_index = sim_index,
          B = B,
          master_seed = seed,
          alpha = alpha
        )

        all_results <- .wr_rbind_fill(
          all_results,
          trial$results
        )

        all_logrank <- .wr_rbind_fill(
          all_logrank,
          trial$logrank
        )

        all_selection <- .wr_rbind_fill(
          all_selection,
          trial$selection
        )

        if (!is.null(checkpoint_file) &&
            (sim_index %% checkpoint_every == 0 ||
             sim_index == nsim)) {
          saveRDS(
            list(
              next_sim = sim_index + 1L,
              results = all_results[
                all_results$scenario_id == scenario_row$scenario_id,
                ,
                drop = FALSE
              ],
              logrank = all_logrank[
                all_logrank$scenario_id == scenario_row$scenario_id,
                ,
                drop = FALSE
              ],
              selection = all_selection[
                all_selection$scenario_id == scenario_row$scenario_id,
                ,
                drop = FALSE
              ]
            ),
            checkpoint_file
          )
        }
      }
    }
  }

  power <- .wr_simulation_power_summary(
    results = all_results,
    logrank = all_logrank,
    alpha = alpha
  )

  out <- list(
    scenarios = scenario_set,
    trial_results = all_results,
    logrank_results = all_logrank,
    selection_results = all_selection,
    power = power,
    settings = list(
      nsim = nsim,
      B = B,
      seed = seed,
      alpha = alpha,
      checkpoint_every = checkpoint_every,
      resume = resume
    )
  )

  class(out) <- c("win_simulation", "list")

  if (!is.null(output_dir)) {
    write_win_simulation_results(
      x = out,
      output_dir = output_dir,
      prefix = "GLOBAL",
      save_plots = save_plots,
      save_scenario_tables = save_scenario_tables
    )
  }

  out
}

#' Print a simulation-study summary.
#' @param x A `win_simulation` object.
#' @param ... Additional arguments.
#' @return `x` invisibly.
#' @export
print.win_simulation <- function(x, ...) {
  cat("Win-ratio simulation study\n")
  cat("Scenarios:", nrow(x$scenarios), "\n")
  cat("NSIM:", x$settings$nsim, "\n")
  cat("Permutation B:", x$settings$B, "\n")
  invisible(x)
}

#' Summarize a simulation study.
#' @param object A `win_simulation` object.
#' @param ... Additional arguments.
#' @return The empirical rejection/power summary.
#' @export
summary.win_simulation <- function(object, ...) {
  object$power
}

#' Convert a simulation study to a data frame.
#' @param x A `win_simulation` object.
#' @param ... Additional arguments.
#' @return The empirical rejection/power summary.
#' @export
as.data.frame.win_simulation <- function(x, ...) {
  x$power
}
