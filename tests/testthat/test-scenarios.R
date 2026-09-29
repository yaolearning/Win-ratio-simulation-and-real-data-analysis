test_that("default scenarios contain the intended 15 settings", {
  s <- default_win_scenarios()

  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 15)
  expect_equal(s$scenario_index, 1:15)
  expect_true(all(c(
    "scenario_id",
    "hazard_ratio",
    "hosp_scale_treatment",
    "followup",
    "censor_rate"
  ) %in% names(s)))
})

test_that("scenario dataset generation is reproducible", {
  a <- generate_win_scenario_dataset(
    scenario = 8,
    sim_index = 2,
    seed = 99
  )
  b <- generate_win_scenario_dataset(
    scenario = 8,
    sim_index = 2,
    seed = 99
  )

  expect_equal(a$trial_seed, b$trial_seed)
  expect_equal(a$subjects, b$subjects)
  expect_equal(a$scenario$scenario_id, "S08_death_benefit_only")
})
