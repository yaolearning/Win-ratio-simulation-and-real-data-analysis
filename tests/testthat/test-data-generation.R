test_that("generate_win_dataset is reproducible", {
  a <- generate_win_dataset(
    n_control = 5,
    n_treatment = 4,
    seed = 123
  )
  b <- generate_win_dataset(
    n_control = 5,
    n_treatment = 4,
    seed = 123
  )

  expect_s3_class(a, "win_simulated_dataset")
  expect_equal(nrow(a$subjects), 9)
  expect_equal(sum(a$subjects$ARM == 0), 5)
  expect_equal(sum(a$subjects$ARM == 1), 4)
  expect_equal(a$subjects, b$subjects)
  expect_equal(a$recurrent_events, b$recurrent_events)
  expect_true(all(c("SUBJID", "HOSPTIME") %in% names(a$recurrent_events)))
})
