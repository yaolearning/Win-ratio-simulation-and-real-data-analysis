test_that("permutation inference is reproducible", {
  sim <- generate_win_dataset(
    n_control = 4,
    n_treatment = 4,
    seed = 77
  )

  dat <- win_data(
    data = sim$subjects,
    recurrent_data = sim$recurrent_events
  )

  endpoints <- list(
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
    )
  )

  ctl <- permutation_control(
    enabled = TRUE,
    B = 9,
    seed = 101
  )

  a <- win_analysis(
    data = dat,
    endpoints = endpoints,
    permutation = ctl,
    run_logrank = FALSE
  )

  b <- win_analysis(
    data = dat,
    endpoints = endpoints,
    permutation = ctl,
    run_logrank = FALSE
  )

  expect_equal(a$results$p_one, b$results$p_one)
  expect_equal(a$results$p_two, b$results$p_two)
  expect_true(all(a$results$p_one >= 0 & a$results$p_one <= 1, na.rm = TRUE))
  expect_true(all(a$results$p_two >= 0 & a$results$p_two <= 1, na.rm = TRUE))
  expect_true(all(a$results$fixed_selected_p_one >= 0 & a$results$fixed_selected_p_one <= 1, na.rm = TRUE))
  expect_true(all(a$results$fixed_selected_p_two >= 0 & a$results$fixed_selected_p_two <= 1, na.rm = TRUE))
})
