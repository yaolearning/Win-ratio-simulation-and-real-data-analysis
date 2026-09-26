default_win_scenarios <- function() {
  base_mortality <- -log(0.6)

  scenario_id <- c(
    "S01_equivalence_equal_arms",
    "S02_similarity_small_benefit_both",
    "S03_similarity_small_harm_both",
    "S04_superiority_moderate_better_both",
    "S05_superiority_strong_better_both",
    "S06_inferiority_moderate_worse_both",
    "S07_inferiority_strong_worse_both",
    "S08_death_benefit_only",
    "S09_hosp_benefit_only",
    "S10_death_benefit_hosp_harm",
    "S11_death_harm_hosp_benefit",
    "S12_high_random_censoring_superiority",
    "S13_low_censoring_long_followup_superiority",
    "S14_strong_hosp_benefit_only",
    "S15_weak_death_benefit_strong_hosp_benefit"
  )

  hazard_ratio <- c(
    1.00, 0.90, 1.10, 0.75, 0.60,
    1.25, 1.50, 0.60, 1.00, 0.60,
    1.50, 0.60, 0.60, 1.00, 0.90
  )

  hosp_scale_treatment <- c(
    1.00, 0.90, 1.10, 0.75, 0.50,
    1.25, 1.50, 1.00, 0.50, 1.50,
    0.50, 0.50, 0.50, 0.25, 0.25
  )

  followup <- c(rep(1, 12), 3, 1, 1)
  censor_rate <- c(rep(0, 11), 1.50, 0.05, 0, 0)

  scenario_type <- c(
    "null / equivalence",
    "similarity / small benefit",
    "similarity / small harm",
    "superiority / moderate benefit",
    "superiority / strong benefit",
    "inferiority / moderate harm",
    "inferiority / strong harm",
    "death benefit only",
    "hospitalization benefit only",
    "discordant: death benefit, hospitalization harm",
    "discordant: death harm, hospitalization benefit",
    "sensitivity: high random censoring",
    "sensitivity: longer follow-up with light censoring",
    "targeted non-terminal benefit",
    "targeted weak death plus strong non-terminal benefit"
  )

  description <- c(
    "No treatment effect on death or hospitalization.",
    "Small benefit on both death and hospitalization.",
    "Small harm on both death and hospitalization.",
    "Moderate treatment benefit on both endpoints.",
    "Strong treatment benefit on both endpoints.",
    "Moderate treatment harm on both endpoints.",
    "Strong treatment harm on both endpoints.",
    "Treatment improves death only; hospitalization burden is unchanged.",
    "Treatment improves hospitalization only; death is unchanged.",
    "Treatment improves death but worsens hospitalization burden.",
    "Treatment worsens death but improves hospitalization burden.",
    "Same strong benefit pattern as S05 with high independent exponential censoring.",
    "Same strong benefit pattern as S05 with longer follow-up and light censoring.",
    "Treatment has no death effect but strongly reduces hospitalization.",
    "Treatment has weak death benefit and strong hospitalization benefit."
  )

  data.frame(
    scenario_index = seq_along(scenario_id),
    scenario_id = scenario_id,
    scenario_type = scenario_type,
    description = description,
    N0 = rep(50L, length(scenario_id)),
    N1 = rep(50L, length(scenario_id)),
    mort_rate_control = rep(base_mortality, length(scenario_id)),
    hazard_ratio = hazard_ratio,
    hosp_shape = rep(5, length(scenario_id)),
    hosp_scale_control = rep(1, length(scenario_id)),
    hosp_scale_treatment = hosp_scale_treatment,
    followup = followup,
    censor_rate = censor_rate,
    stringsAsFactors = FALSE
  )
}

.wr_resolve_scenario <- function(scenario,
                                 scenario_grid = default_win_scenarios()) {
  if (!is.data.frame(scenario_grid) || nrow(scenario_grid) == 0) {
    stop("scenario_grid must be a non-empty data frame.")
  }

  if (is.data.frame(scenario)) {
    if (nrow(scenario) != 1L) {
      stop("A scenario data frame must contain exactly one row.")
    }
    return(scenario)
  }

  if (is.numeric(scenario) && length(scenario) == 1L) {
    k <- as.integer(scenario)

    if ("scenario_index" %in% names(scenario_grid) &&
        k %in% scenario_grid$scenario_index) {
      return(
        scenario_grid[
          match(k, scenario_grid$scenario_index),
          ,
          drop = FALSE
        ]
      )
    }

    if (k < 1L || k > nrow(scenario_grid)) {
      stop("Scenario index is out of range.")
    }

    return(scenario_grid[k, , drop = FALSE])
  }

  if (is.character(scenario) && length(scenario) == 1L) {
    k <- match(scenario, scenario_grid$scenario_id)
    if (is.na(k)) stop("Unknown scenario_id: ", scenario)
    return(scenario_grid[k, , drop = FALSE])
  }

  stop("scenario must be a scenario index, scenario_id, or one-row data frame.")
}

.wr_resolve_scenario_set <- function(scenarios,
                                     scenario_grid = default_win_scenarios()) {
  if (is.null(scenarios)) return(scenario_grid)
  if (is.data.frame(scenarios)) return(scenarios)

  if (is.numeric(scenarios)) {
    rows <- lapply(
      scenarios,
      .wr_resolve_scenario,
      scenario_grid = scenario_grid
    )
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  if (is.character(scenarios)) {
    rows <- lapply(
      scenarios,
      .wr_resolve_scenario,
      scenario_grid = scenario_grid
    )
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  stop("scenarios must be NULL, scenario indices, scenario IDs, or a data frame.")
}

generate_win_scenario_dataset <- function(scenario = 1L,
                                          scenario_grid = default_win_scenarios(),
                                          sim_index = 1L,
                                          seed = 2026L) {
  row <- .wr_resolve_scenario(
    scenario = scenario,
    scenario_grid = scenario_grid
  )

  sim_index <- as.integer(sim_index)

  if (is.na(sim_index) || sim_index < 1L) {
    stop("sim_index must be at least 1.")
  }

  scenario_index <- if ("scenario_index" %in% names(row)) {
    as.integer(row$scenario_index[1])
  } else {
    1L
  }

  trial_seed <- as.integer(
    seed +
      scenario_index * 1000000L +
      sim_index * 104729L
  )

  out <- generate_win_dataset(
    n_control = row$N0[1],
    n_treatment = row$N1[1],
    mort_rate_control = row$mort_rate_control[1],
    hazard_ratio = row$hazard_ratio[1],
    hosp_shape = row$hosp_shape[1],
    hosp_scale_control = row$hosp_scale_control[1],
    hosp_scale_treatment = row$hosp_scale_treatment[1],
    followup = row$followup[1],
    censor_rate = row$censor_rate[1],
    seed = trial_seed
  )

  out$scenario <- row
  out$sim_index <- sim_index
  out$trial_seed <- trial_seed
  out
}

.wr_simulation_threshold_values <- function(subjects,
                                            clinical_months = c(1, 3, 6, 12, 18, 24),
                                            max_candidates = 10L) {
  max_followup <- max(subjects$FUTIME, na.rm = TRUE)

  clinical <- as.numeric(clinical_months) / 12
  clinical <- clinical[
    is.finite(clinical) &
      clinical > 0 &
      clinical <= max_followup
  ]

  treatment <- subjects$FUTIME[subjects$ARM == 1]
  control <- subjects$FUTIME[subjects$ARM == 0]
  empirical <- numeric(0)

  if (length(treatment) > 0 && length(control) > 0) {
    diffs <- abs(as.vector(outer(treatment, control, "-")))
    diffs <- diffs[
      is.finite(diffs) &
        diffs > 0 &
        diffs <= max_followup
    ]

    if (length(diffs) > 0) {
      empirical <- as.numeric(stats::quantile(
        diffs,
        probs = c(0.25, 0.50, 0.75, 0.90),
        na.rm = TRUE,
        names = FALSE
      ))
    }
  }

  candidates <- sort(unique(round(c(clinical, empirical), 8)))
  candidates <- candidates[
    is.finite(candidates) &
      candidates > 0 &
      candidates <= max_followup
  ]

  if (length(candidates) == 0) {
    candidates <- min(max_followup, 1 / 12)
  }

  if (length(candidates) > max_candidates) {
    idx <- unique(round(seq(
      1,
      length(candidates),
      length.out = max_candidates
    )))
    candidates <- candidates[idx]
  }

  candidates
}
