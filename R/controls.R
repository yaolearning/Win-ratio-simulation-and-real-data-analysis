permutation_control <- function(enabled = TRUE,
                                B = 500L,
                                seed = 2026L) {
  B <- as.integer(B)

  if (isTRUE(enabled) && (is.na(B) || B < 1L)) {
    stop("B must be at least 1 when permutation is enabled.")
  }

  structure(
    list(
      enabled = isTRUE(enabled),
      B = if (isTRUE(enabled)) B else 0L,
      seed = as.integer(seed)
    ),
    class = c("win_permutation_control", "list")
  )
}

method_control <- function(run_wr = TRUE,
                           run_wo = TRUE,
                           report_all_fixed_orders = TRUE) {
  if (!isTRUE(run_wr) && !isTRUE(run_wo)) {
    stop("At least one of run_wr or run_wo must be TRUE.")
  }

  structure(
    list(
      run_wr = isTRUE(run_wr),
      run_wo = isTRUE(run_wo),
      report_all_fixed_orders = isTRUE(report_all_fixed_orders)
    ),
    class = c("win_method_control", "list")
  )
}

engine_control <- function(backend = c("auto", "general", "two_endpoint_recurrent"),
                           eps = 1e-8) {
  backend <- match.arg(backend)
  eps <- as.numeric(eps)[1]

  if (!is.finite(eps) || eps <= 0) {
    stop("eps must be positive.")
  }

  structure(
    list(
      backend = backend,
      eps = eps
    ),
    class = c("win_engine_control", "list")
  )
}
