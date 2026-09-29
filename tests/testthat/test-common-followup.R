test_that("general recurrent count uses pairwise common follow-up", {
  subjects <- data.frame(
    SUBJID = c("T1", "C1"),
    ARM = c(1, 0),
    FUTIME = c(1, 2),
    CNSR = c(0, 0),
    NUMHOSP = c(1, 2),
    SCORE = c(0, 0)
  )

  recurrent <- data.frame(
    SUBJID = c("T1", "C1", "C1"),
    HOSPTIME = c(0.5, 0.75, 1.5)
  )

  dat <- win_data(
    data = subjects,
    recurrent_data = recurrent
  )

  endpoints_common <- list(
    endpoint_time(
      "Death",
      time = "FUTIME",
      event = "CNSR",
      unit = "years"
    ),
    endpoint_count(
      "Hospitalization",
      count = "NUMHOSP",
      comparison = "pairwise_common_followup",
      recurrent_time = "HOSPTIME",
      recurrent_id = "SUBJID",
      followup = "FUTIME",
      unit = "years"
    ),
    endpoint_continuous(
      "Score",
      value = "SCORE"
    )
  )

  fit_common <- win_analysis(
    data = dat,
    endpoints = endpoints_common,
    permutation = permutation_control(enabled = FALSE),
    engine = engine_control(backend = "general"),
    run_logrank = FALSE
  )

  d_common <- fit_common$observed_candidates
  d_common <- d_common[
    d_common$order_key == "1->2->3" &
      abs(d_common$p1 - 1 / 3) < 1e-8 &
      abs(d_common$p2 - 1 / 3) < 1e-8 &
      abs(d_common$p3 - 1 / 3) < 1e-8 &
      abs(d_common$threshold_e1) < 1e-12,
    ,
    drop = FALSE
  ]

  expect_true(nrow(d_common) >= 1)
  expect_equal(d_common$tie_count[1], 1)
  expect_equal(d_common$win_rank2[1], 0)

  endpoints_total <- endpoints_common
  endpoints_total[[2]] <- endpoint_count(
    "Hospitalization",
    count = "NUMHOSP",
    comparison = "observed_total"
  )

  fit_total <- win_analysis(
    data = subjects,
    endpoints = endpoints_total,
    permutation = permutation_control(enabled = FALSE),
    engine = engine_control(backend = "general"),
    run_logrank = FALSE
  )

  d_total <- fit_total$observed_candidates
  d_total <- d_total[
    d_total$order_key == "1->2->3" &
      abs(d_total$p1 - 1 / 3) < 1e-8 &
      abs(d_total$p2 - 1 / 3) < 1e-8 &
      abs(d_total$p3 - 1 / 3) < 1e-8 &
      abs(d_total$threshold_e1) < 1e-12,
    ,
    drop = FALSE
  ]

  expect_true(nrow(d_total) >= 1)
  expect_equal(d_total$tie_count[1], 0)
  expect_equal(d_total$win_rank2[1], 1)
})
