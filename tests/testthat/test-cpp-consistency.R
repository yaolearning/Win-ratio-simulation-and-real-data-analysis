test_that("general R and C++ candidate cores agree", {
  skip_if(!exists("cpp_counts_for_candidate_core", envir = asNamespace("MaxWin"), inherits = FALSE))

  subjects <- data.frame(
    SUBJID = c("T1", "T2", "C1", "C2"),
    ARM = c(1, 1, 0, 0),
    FUTIME = c(1.0, 1.5, 1.2, 1.7),
    CNSR = c(1, 0, 1, 0),
    NUMHOSP = c(1, 1, 2, 0),
    SCORE = c(2, 1, 0, 2)
  )

  recurrent <- data.frame(
    SUBJID = c("T1", "T2", "C1", "C1"),
    HOSPTIME = c(0.4, 0.9, 0.3, 1.1)
  )

  endpoints <- list(
    endpoint_time("Death", "FUTIME", "CNSR", unit = "years"),
    endpoint_count(
      "Hospitalization",
      "NUMHOSP",
      comparison = "pairwise_common_followup",
      recurrent_time = "HOSPTIME",
      recurrent_id = "SUBJID",
      followup = "FUTIME",
      unit = "years"
    ),
    endpoint_continuous("Score", "SCORE")
  )

  prepared <- MaxWin:::prepare_win_data(
    data = win_data(subjects, recurrent_data = recurrent),
    endpoints = endpoints,
    format = "wide",
    id = "SUBJID",
    treatment = "ARM"
  )

  pairs <- MaxWin:::.wr_get_pair_indices(prepared$data$arm)
  cache <- MaxWin:::.wr_prepare_endpoint_cache(
    prepared$data,
    prepared$endpoints,
    prepared$recurrent_cache
  )

  order_vec <- 1:3
  weights <- rep(1 / 3, 3)
  thresholds <- c(0, 0, 0)

  r <- MaxWin:::.wr_counts_for_candidate_R(
    data = prepared$data,
    arm = prepared$data$arm,
    endpoint_specs = prepared$endpoints,
    recurrent_cache = prepared$recurrent_cache,
    order_vec = order_vec,
    weights = weights,
    thresholds_by_endpoint = thresholds,
    pairs = pairs
  )

  cpp <- MaxWin:::cpp_counts_for_candidate_core(
    idx_treatment = as.integer(pairs$treatment),
    idx_control = as.integer(pairs$control),
    order_vec = as.integer(order_vec),
    weights = as.numeric(weights),
    thresholds_by_endpoint = as.numeric(thresholds),
    type_code = as.integer(cache$type_code),
    time_mat = cache$time_mat,
    event_mat = cache$event_mat,
    value_mat = cache$value_mat,
    followup_mat = cache$followup_mat,
    recurrent_times = cache$recurrent_times,
    recurrent_start = cache$recurrent_start,
    recurrent_len = cache$recurrent_len
  )

  expect_equal(r$WR_statistic, cpp$WR_statistic)
  expect_equal(r$WO_statistic, cpp$WO_statistic)
  expect_equal(r$tie_count, cpp$tie_count)
  expect_equal(r$wins_by_rank, as.numeric(cpp$wins_by_rank))
  expect_equal(r$losses_by_rank, as.numeric(cpp$losses_by_rank))
})
