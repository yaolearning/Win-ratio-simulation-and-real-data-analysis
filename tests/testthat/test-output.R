test_that("real-data output writer creates the expected files", {
  sim <- generate_win_dataset(
    n_control = 3,
    n_treatment = 3,
    seed = 11
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

  fit <- win_analysis(
    data = dat,
    endpoints = endpoints,
    permutation = permutation_control(enabled = FALSE),
    run_logrank = FALSE
  )

  outdir <- file.path(tempdir(), paste0("maxwin-test-", Sys.getpid()))
  unlink(outdir, recursive = TRUE)

  write_win_results(
    fit,
    output_dir = outdir,
    prefix = "TEST",
    save_plots = FALSE
  )

  expect_true(file.exists(file.path(outdir, "TEST_method_results.csv")))
  expect_true(file.exists(file.path(outdir, "TEST_observed_candidates.csv")))
  expect_true(file.exists(file.path(outdir, "TEST_tie_summary.csv")))
  expect_true(file.exists(file.path(outdir, "TEST_settings.csv")))
  expect_true(file.exists(file.path(outdir, "TEST_full_results.rds")))
  expect_true(file.exists(file.path(outdir, "TEST_output_manifest.csv")))
})
