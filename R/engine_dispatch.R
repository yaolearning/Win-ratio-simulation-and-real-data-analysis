evaluate_win_candidates <- function(prepared,
                                    candidate_grid,
                                    engine = engine_control()) {
  backend <- engine$backend

  if (backend == "auto") {
    if (is.list(prepared) &&
        all(c("table.output", "hosp.times.list") %in% names(prepared)) &&
        !is.null(prepared$hosp.times.list)) {
      backend <- "two_endpoint_recurrent"
    } else {
      backend <- "general"
    }
  }

  if (backend == "two_endpoint_recurrent") {
    return(evaluate_two_endpoint_recurrent_candidates(
      ds = prepared,
      candidate_grid = candidate_grid,
      eps = engine$eps
    ))
  }

  if (!is.list(prepared) ||
      !all(c("data", "endpoints") %in% names(prepared))) {
    stop("The general engine requires prepared$data and prepared$endpoints.")
  }

  evaluate_general_candidates(
    data = prepared$data,
    arm = prepared$data$arm,
    endpoint_specs = prepared$endpoints,
    candidate_grid = candidate_grid,
    recurrent_cache = prepared$recurrent_cache
  )
}
