.wr_resolve_backend <- function(data,
                                endpoints,
                                format,
                                engine) {
  if (engine$backend != "auto") return(engine$backend)

  if (format == "wide" &&
      length(endpoints) == 2L &&
      identical(endpoints[[1]]$type, "time") &&
      identical(endpoints[[2]]$type, "count") &&
      identical(endpoints[[2]]$comparison, "pairwise_common_followup")) {
    return("two_endpoint_recurrent")
  }

  "general"
}

.wr_prepare_specialized_two_endpoint <- function(data,
                                                 endpoints,
                                                 id,
                                                 treatment) {
  ep1 <- endpoints[[1]]
  ep2 <- endpoints[[2]]
  recurrent_data <- .wr_recurrent_data(data)

  recurrent_id <- if (!is.null(ep2$recurrent_id_col)) {
    ep2$recurrent_id_col
  } else {
    id
  }

  recurrent_event_time <- if (!is.null(ep2$recurrent_time_col) &&
                              !is.null(recurrent_data) &&
                              is.data.frame(recurrent_data) &&
                              ep2$recurrent_time_col %in% names(recurrent_data)) {
    ep2$recurrent_time_col
  } else {
    "HOSPTIME"
  }

  subject_recurrent_time <- if (!is.null(ep2$recurrent_time_col)) {
    ep2$recurrent_time_col
  } else {
    "HOSP_TIMES"
  }

  prepare_two_endpoint_recurrent_data(
    data = data,
    id = id,
    treatment = treatment,
    followup = ep1$time_col,
    event = ep1$event_col,
    count = ep2$count_col,
    recurrent_time = subject_recurrent_time,
    time_unit = ep1$unit,
    recurrent_time_unit = ep2$unit,
    recurrent_time_type = "absolute",
    recurrent_id = recurrent_id,
    recurrent_event_time = recurrent_event_time,
    legacy_count_approximation = FALSE
  )
}

preview_win_plan <- function(data,
                             endpoints,
                             format = c("wide", "event_long"),
                             id = "SUBJID",
                             treatment = "ARM",
                             event_time = NULL,
                             event_status = NULL,
                             event_time_unit = "years",
                             threshold = threshold_control(),
                             weight = weight_control(),
                             methods = method_control(),
                             engine = engine_control()) {
  format <- match.arg(format)

  prepared_general <- prepare_win_data(
    data = data,
    endpoints = endpoints,
    format = format,
    id = id,
    treatment = treatment,
    event_time = event_time,
    event_status = event_status,
    event_time_unit = event_time_unit
  )

  plan <- build_win_method_plan(
    endpoint_specs = prepared_general$endpoints,
    threshold = threshold,
    weight = weight,
    methods = methods,
    data = prepared_general$data
  )

  backend <- .wr_resolve_backend(
    data = data,
    endpoints = endpoints,
    format = format,
    engine = engine
  )

  list(
    backend = backend,
    endpoint_specs = prepared_general$endpoints,
    method_overview = plan$method_overview,
    threshold_info = plan$threshold$info,
    threshold_grid = plan$threshold$grid,
    weight_grid = plan$weights,
    n_unique_candidates = nrow(plan$all_candidates)
  )
}

win_analysis <- function(data,
                         endpoints,
                         format = c("wide", "event_long"),
                         id = "SUBJID",
                         treatment = "ARM",
                         event_time = NULL,
                         event_status = NULL,
                         event_time_unit = "years",
                         threshold = threshold_control(),
                         weight = weight_control(),
                         methods = method_control(),
                         permutation = permutation_control(),
                         engine = engine_control(),
                         run_logrank = TRUE,
                         run_composite_logrank_if_possible = TRUE,
                         output_dir = NULL,
                         output_prefix = "WIN") {
  format <- match.arg(format)

  prepared_general <- prepare_win_data(
    data = data,
    endpoints = endpoints,
    format = format,
    id = id,
    treatment = treatment,
    event_time = event_time,
    event_status = event_status,
    event_time_unit = event_time_unit
  )

  plan <- build_win_method_plan(
    endpoint_specs = prepared_general$endpoints,
    threshold = threshold,
    weight = weight,
    methods = methods,
    data = prepared_general$data
  )

  backend <- .wr_resolve_backend(
    data = data,
    endpoints = endpoints,
    format = format,
    engine = engine
  )

  engine_resolved <- engine
  engine_resolved$backend <- backend

  prepared_engine <- if (backend == "two_endpoint_recurrent") {
    .wr_prepare_specialized_two_endpoint(
      data = data,
      endpoints = endpoints,
      id = id,
      treatment = treatment
    )
  } else {
    prepared_general
  }

  permutation_out <- run_win_permutation(
    prepared = prepared_engine,
    candidate_grid = plan$all_candidates,
    methods = plan$methods,
    method_control = methods,
    permutation = permutation,
    engine = engine_resolved
  )

  results <- .wr_results_wide(
    observed = permutation_out$observed,
    B = permutation_out$B
  )

  logrank <- run_win_logrank(
    prepared = prepared_general,
    run = run_logrank,
    composite_if_possible = run_composite_logrank_if_possible
  )

  out <- list(
    results = results,
    selection = permutation_out$observed,
    permutation_selected = permutation_out$perm_selected,
    observed_candidates = permutation_out$observed_counts,
    all_candidates = plan$all_candidates,
    method_overview = plan$method_overview,
    threshold_info = plan$threshold$info,
    threshold_grid = plan$threshold$grid,
    weight_grid = plan$weights,
    endpoint_summary = .wr_endpoint_summary(prepared_general),
    logrank = logrank,
    prepared_data = prepared_general$data,
    prepared_recurrent_cache = prepared_general$recurrent_cache,
    endpoint_specs = prepared_general$endpoints,
    backend = backend,
    B = permutation_out$B,
    subject_summary = .wr_subject_summary(
      prepared_general = prepared_general,
      backend = backend,
      B = permutation_out$B,
      threshold_object = plan$threshold,
      weight_grid = plan$weights
    ),
    controls = list(
      threshold = threshold,
      weight = weight,
      methods = methods,
      permutation = permutation,
      engine = engine_resolved
    )
  )

  class(out) <- c("win_analysis", "list")

  if (!is.null(output_dir)) {
    write_win_results(
      x = out,
      output_dir = output_dir,
      prefix = output_prefix
    )
  }

  out
}
