test_that("run_win_simulation completes a minimal simulation", {
  sim <- run_win_simulation(
    scenarios = 1,
    nsim = 1,
    B = 2,
    seed = 2026,
    output_dir = NULL,
    verbose = FALSE
  )

  expect_s3_class(sim, "win_simulation")
  expect_equal(nrow(sim$scenarios), 1)
  expect_true(nrow(sim$trial_results) > 0)
  expect_true(nrow(sim$power) > 0)
  expect_equal(sim$settings$nsim, 1L)
  expect_equal(sim$settings$B, 2L)
})
