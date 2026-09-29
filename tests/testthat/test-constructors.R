test_that("endpoint constructors preserve requested settings", {
  e1 <- endpoint_time("Death", time = "FUTIME", event = "CNSR", unit = "years")
  e2 <- endpoint_count(
    "Hospitalization",
    count = "NUMHOSP",
    comparison = "pairwise_common_followup",
    recurrent_time = "HOSPTIME",
    recurrent_id = "SUBJID",
    followup = "FUTIME",
    unit = "years"
  )
  e3 <- endpoint_continuous("Score", value = "SCORE", higher_better = TRUE)

  expect_s3_class(e1, "win_endpoint")
  expect_identical(e1$type, "time")
  expect_identical(e2$comparison, "pairwise_common_followup")
  expect_identical(e2$followup_col, "FUTIME")
  expect_identical(e3$type, "continuous")
})

test_that("control constructors are valid", {
  expect_s3_class(permutation_control(enabled = FALSE), "win_permutation_control")
  expect_s3_class(method_control(), "win_method_control")
  expect_s3_class(engine_control(), "win_engine_control")
  expect_s3_class(threshold_control(), "win_threshold_control")
  expect_s3_class(weight_control(), "win_weight_control")
})
