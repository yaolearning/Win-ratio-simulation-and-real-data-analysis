.wr_zero_threshold_grid <- function(m) {
  m <- as.integer(m)
  out <- as.data.frame(
    setNames(
      replicate(m, 0, simplify = FALSE),
      paste0("threshold_e", seq_len(m))
    ),
    stringsAsFactors = FALSE
  )
  out$threshold_key <- paste(rep("0.00000000", m), collapse = "|")
  out
}

.wr_candidate_key <- function(df, m) {
  p_cols <- paste0("p", seq_len(m))
  t_cols <- paste0("threshold_e", seq_len(m))

  p_part <- apply(
    df[, p_cols, drop = FALSE],
    1,
    function(z) paste(sprintf("%.8f", as.numeric(z)), collapse = "|")
  )

  t_part <- apply(
    df[, t_cols, drop = FALSE],
    1,
    function(z) paste(sprintf("%.8f", as.numeric(z)), collapse = "|")
  )

  paste(df$order_key, p_part, t_part, sep = "__")
}

.wr_make_candidate_grid <- function(orders,
                                    weights,
                                    thresholds,
                                    m) {
  m <- as.integer(m)
  if (!is.data.frame(weights)) weights <- as.data.frame(weights)
  if (!is.data.frame(thresholds)) thresholds <- as.data.frame(thresholds)

  p_cols <- paste0("p", seq_len(m))
  t_cols <- paste0("threshold_e", seq_len(m))

  if (!all(p_cols %in% names(weights))) {
    stop("Weight grid is missing columns: ", paste(setdiff(p_cols, names(weights)), collapse = ", "))
  }

  if (!all(t_cols %in% names(thresholds))) {
    stop("Threshold grid is missing columns: ", paste(setdiff(t_cols, names(thresholds)), collapse = ", "))
  }

  rows <- vector("list", length(orders) * nrow(weights) * nrow(thresholds))
  idx <- 1L

  for (ord in orders) {
    ord_key <- .wr_format_order(ord)

    for (i in seq_len(nrow(weights))) {
      w <- as.numeric(weights[i, p_cols, drop = TRUE])

      for (j in seq_len(nrow(thresholds))) {
        row <- data.frame(
          order_key = ord_key,
          stringsAsFactors = FALSE
        )

        for (k in seq_len(m)) row[[paste0("p", k)]] <- w[k]
        for (k in seq_len(m)) {
          row[[paste0("threshold_e", k)]] <- as.numeric(thresholds[[paste0("threshold_e", k)]][j])
        }

        rows[[idx]] <- row
        idx <- idx + 1L
      }
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$candidate_key <- .wr_candidate_key(out, m)
  out[!duplicated(out$candidate_key), , drop = FALSE]
}

.wr_build_methods <- function(m,
                              weight_grid,
                              threshold_grid,
                              control = method_control()) {
  m <- as.integer(m)
  natural_order <- .wr_natural_order(m)
  all_orders <- .wr_all_orders(m)
  equal_weights <- .wr_equal_weights(m)

  equal_df <- as.data.frame(as.list(equal_weights))
  names(equal_df) <- paste0("p", seq_len(m))

  zero_thresholds <- .wr_zero_threshold_grid(m)
  methods <- list()

  natural_grid <- .wr_make_candidate_grid(
    orders = list(natural_order),
    weights = equal_df,
    thresholds = zero_thresholds,
    m = m
  )

  methods[[length(methods) + 1L]] <- list(
    id = "original",
    notation = if (m == 2L) "WR/WO" else "WR/WO^{(3)}",
    description = "Fixed natural order, equal weights, no threshold",
    grid = natural_grid
  )

  if (isTRUE(control$report_all_fixed_orders)) {
    for (ord in all_orders) {
      if (identical(as.integer(ord), as.integer(natural_order))) next
      ord_key <- .wr_format_order(ord)

      methods[[length(methods) + 1L]] <- list(
        id = paste0("fixed_order_", gsub("->", "_", ord_key, fixed = TRUE)),
        notation = paste0("Fixed order ", ord_key),
        description = paste0("Fixed order ", ord_key, ", equal weights, no threshold"),
        grid = .wr_make_candidate_grid(
          orders = list(ord),
          weights = equal_df,
          thresholds = zero_thresholds,
          m = m
        )
      )
    }
  }

  methods[[length(methods) + 1L]] <- list(
    id = "M_ord",
    notation = if (m == 2L) "M_ord" else "M_ord^{(3)}",
    description = "Maximize over endpoint order; equal weights; no threshold",
    grid = .wr_make_candidate_grid(
      orders = all_orders,
      weights = equal_df,
      thresholds = zero_thresholds,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_wt",
    notation = if (m == 2L) "M_wt^0.5" else "M_wt^{(3)}",
    description = "Fixed natural order; maximize over ordered weights; no threshold",
    grid = .wr_make_candidate_grid(
      orders = list(natural_order),
      weights = weight_grid,
      thresholds = zero_thresholds,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt",
    notation = if (m == 2L) "M_ord,wt^0.5" else "M_ord,wt^{(3)}",
    description = "Maximize over endpoint order and ordered weights; no threshold",
    grid = .wr_make_candidate_grid(
      orders = all_orders,
      weights = weight_grid,
      thresholds = zero_thresholds,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_thresh",
    notation = if (m == 2L) "M_thresh" else "M_thresh^{(3)}",
    description = "Fixed natural order; equal weights; maximize threshold(s) on time-to-event endpoint(s)",
    grid = .wr_make_candidate_grid(
      orders = list(natural_order),
      weights = equal_df,
      thresholds = threshold_grid,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_wt_thresh",
    notation = if (m == 2L) "M_wt,thresh^0.5" else "M_wt,thresh^{(3)}",
    description = "Fixed natural order; maximize ordered weights and threshold(s) on time-to-event endpoint(s)",
    grid = .wr_make_candidate_grid(
      orders = list(natural_order),
      weights = weight_grid,
      thresholds = threshold_grid,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_thresh",
    notation = if (m == 2L) "M_ord,thresh" else "M_ord,thresh^{(3)}",
    description = "Maximize over endpoint order and threshold(s) on time-to-event endpoint(s); equal weights",
    grid = .wr_make_candidate_grid(
      orders = all_orders,
      weights = equal_df,
      thresholds = threshold_grid,
      m = m
    )
  )

  methods[[length(methods) + 1L]] <- list(
    id = "M_ord_wt_thresh",
    notation = if (m == 2L) "M_ord,wt,thresh^0.5" else "M_ord,wt,thresh^{(3)}",
    description = "Full adaptive: maximize over order, ordered weights, and time-endpoint threshold(s)",
    grid = .wr_make_candidate_grid(
      orders = all_orders,
      weights = weight_grid,
      thresholds = threshold_grid,
      m = m
    )
  )

  methods
}

.wr_method_overview <- function(methods) {
  rows <- lapply(methods, function(x) {
    data.frame(
      method_id = x$id,
      notation = x$notation,
      description = x$description,
      n_candidates = nrow(x$grid),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.wr_all_unique_candidates <- function(methods) {
  out <- .wr_rbind_fill(lapply(methods, function(x) x$grid))
  if (nrow(out) == 0) return(out)
  out[!duplicated(out$candidate_key), , drop = FALSE]
}

build_win_method_plan <- function(endpoint_specs,
                                  threshold = threshold_control(),
                                  weight = weight_control(),
                                  methods = method_control(),
                                  data = NULL) {
  m <- length(endpoint_specs)
  if (!(m %in% c(2L, 3L))) stop("Only 2 or 3 endpoints are supported.")

  if (is.null(data)) {
    stop("Prepared subject-level data are required to build data-driven threshold candidates.")
  }

  threshold_object <- .wr_make_threshold_grid(
    data = data,
    endpoint_specs = endpoint_specs,
    control = threshold
  )

  weight_grid <- .wr_make_weight_grid(
    m = m,
    control = weight
  )

  method_list <- .wr_build_methods(
    m = m,
    weight_grid = weight_grid,
    threshold_grid = threshold_object$grid,
    control = methods
  )

  list(
    methods = method_list,
    method_overview = .wr_method_overview(method_list),
    all_candidates = .wr_all_unique_candidates(method_list),
    threshold = threshold_object,
    weights = weight_grid
  )
}
