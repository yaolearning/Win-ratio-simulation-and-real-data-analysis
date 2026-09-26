plot_win_pvalues <- function(x,
                             measure = c("WR", "WO"),
                             side = c("one", "two"),
                             alpha = 0.05,
                             main = NULL) {
  if (!inherits(x, "win_analysis")) {
    stop("x must be a win_analysis object.")
  }

  measure <- match.arg(measure)
  side <- match.arg(side)

  d <- x$results[
    x$results$measure == measure,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    stop("No results found for measure: ", measure)
  }

  p_col <- if (side == "one") "p_one" else "p_two"
  values <- as.numeric(d[[p_col]])
  labels <- paste(measure, d$notation)

  if (all(is.na(values))) {
    graphics::plot.new()
    graphics::title(
      main = if (is.null(main)) {
        paste0(measure, ": permutation not run")
      } else {
        main
      }
    )
    graphics::text(
      0.5,
      0.5,
      "No permutation p-values available"
    )
    return(invisible(d))
  }

  ord <- order(
    values,
    decreasing = TRUE,
    na.last = TRUE
  )

  values <- values[ord]
  labels <- labels[ord]

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  graphics::par(
    mar = c(5.2, 13, 3.2, 1.2)
  )

  bp <- graphics::barplot(
    values,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    xlim = c(0, 1),
    xlab = "Permutation p-value",
    main = if (is.null(main)) {
      paste0(
        measure,
        ": ",
        side,
        "-sided permutation p-values"
      )
    } else {
      main
    },
    cex.names = 0.8
  )

  graphics::abline(
    v = alpha,
    lty = 2
  )

  graphics::text(
    x = pmin(values + 0.035, 0.97),
    y = bp,
    labels = ifelse(
      is.na(values),
      "",
      sprintf("%.3f", values)
    ),
    cex = 0.75
  )

  invisible(
    d[ord, , drop = FALSE]
  )
}

plot_win_km <- function(x,
                        endpoint = 1L,
                        main = NULL) {
  if (!inherits(x, "win_analysis")) {
    stop("x must be a win_analysis object.")
  }

  endpoint <- as.integer(endpoint)

  if (endpoint < 1L ||
      endpoint > length(x$endpoint_specs)) {
    stop("endpoint is out of range.")
  }

  ep <- x$endpoint_specs[[endpoint]]

  if (!identical(ep$type, "time")) {
    stop("Selected endpoint is not time-to-event.")
  }

  if (is.null(x$prepared_data)) {
    stop("Prepared subject-level data are not stored in this analysis object.")
  }

  data <- x$prepared_data

  surv <- survival::Surv(
    data[[ep$time_col]],
    data[[ep$event_col]] == 1
  )

  fit <- survival::survfit(
    surv ~ arm,
    data = data
  )

  graphics::plot(
    fit,
    col = c("blue", "red"),
    lwd = 2,
    xlab = "Time (years)",
    ylab = "Event-free probability",
    main = if (is.null(main)) {
      paste0(
        "Kaplan-Meier curve: ",
        ep$name
      )
    } else {
      main
    }
  )

  graphics::legend(
    "bottomleft",
    legend = c(
      "Control",
      "Treatment"
    ),
    col = c(
      "blue",
      "red"
    ),
    lwd = 2,
    bty = "n"
  )

  invisible(fit)
}

plot_win_ties <- function(x,
                          measure = c("WR", "WO"),
                          side = c("one", "two"),
                          main = NULL) {
  if (!inherits(x, "win_analysis")) {
    stop("x must be a win_analysis object.")
  }

  measure <- match.arg(measure)
  side <- match.arg(side)

  d <- x$results[
    x$results$measure == measure,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    stop("No results found for measure: ", measure)
  }

  tie_col <- if (side == "one") {
    "tie_proportion_one"
  } else {
    "tie_proportion_two"
  }

  perm_col <- if (side == "one") {
    "mean_perm_tie_proportion_one"
  } else {
    "mean_perm_tie_proportion_two"
  }

  mat <- rbind(
    observed = as.numeric(d[[tie_col]]),
    permutation_mean = as.numeric(d[[perm_col]])
  )

  colnames(mat) <- d$notation

  ymax <- suppressWarnings(
    max(mat, na.rm = TRUE)
  )

  if (!is.finite(ymax)) ymax <- 1
  ymax <- max(1, ymax)

  old_par <- graphics::par(no.readonly = TRUE)

  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  graphics::par(
    mar = c(9, 5, 3.2, 1.2)
  )

  graphics::barplot(
    mat,
    beside = TRUE,
    las = 2,
    ylim = c(0, ymax),
    ylab = "Tie proportion",
    main = if (is.null(main)) {
      paste0(
        measure,
        ": ",
        side,
        "-sided tie summary"
      )
    } else {
      main
    },
    legend.text = TRUE,
    args.legend = list(
      x = "topright",
      bty = "n"
    ),
    cex.names = 0.8
  )

  invisible(d)
}

plot.win_analysis <- function(x,
                              type = c("pvalue", "km", "ties"),
                              measure = c("WR", "WO"),
                              side = c("one", "two"),
                              endpoint = 1L,
                              alpha = 0.05,
                              ...) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  side <- match.arg(side)

  if (type == "pvalue") {
    return(
      plot_win_pvalues(
        x = x,
        measure = measure,
        side = side,
        alpha = alpha,
        ...
      )
    )
  }

  if (type == "km") {
    return(
      plot_win_km(
        x = x,
        endpoint = endpoint,
        ...
      )
    )
  }

  plot_win_ties(
    x = x,
    measure = measure,
    side = side,
    ...
  )
}

.wr_open_png <- function(file,
                         width = 2200,
                         height = 1400,
                         res = 180) {
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  grDevices::png(
    file,
    width = width,
    height = height,
    res = res
  )
}

.wr_open_pdf <- function(file,
                         width = 10,
                         height = 7) {
  dir.create(
    dirname(file),
    recursive = TRUE,
    showWarnings = FALSE
  )

  grDevices::pdf(
    file,
    width = width,
    height = height
  )
}

.wr_save_dual_plot <- function(png_file,
                               pdf_file,
                               plot_fun,
                               png_width = 2200,
                               png_height = 1400,
                               pdf_width = 10,
                               pdf_height = 7) {
  .wr_open_png(
    png_file,
    width = png_width,
    height = png_height
  )

  try(
    plot_fun(),
    silent = TRUE
  )

  grDevices::dev.off()

  .wr_open_pdf(
    pdf_file,
    width = pdf_width,
    height = pdf_height
  )

  try(
    plot_fun(),
    silent = TRUE
  )

  grDevices::dev.off()

  invisible(NULL)
}

save_win_plots <- function(x,
                           output_dir,
                           prefix = "WIN",
                           alpha = 0.05) {
  if (!inherits(x, "win_analysis")) {
    stop("x must be a win_analysis object.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (measure in unique(x$results$measure)) {
    if (!(measure %in% c("WR", "WO"))) next

    for (side in c("one", "two")) {
      .wr_save_dual_plot(
        png_file = file.path(
          output_dir,
          paste0(
            prefix,
            "_",
            measure,
            "_",
            side,
            "_pvalues.png"
          )
        ),
        pdf_file = file.path(
          output_dir,
          paste0(
            prefix,
            "_",
            measure,
            "_",
            side,
            "_pvalues.pdf"
          )
        ),
        plot_fun = function() {
          plot_win_pvalues(
            x = x,
            measure = measure,
            side = side,
            alpha = alpha
          )
        }
      )

      .wr_save_dual_plot(
        png_file = file.path(
          output_dir,
          paste0(
            prefix,
            "_",
            measure,
            "_",
            side,
            "_ties.png"
          )
        ),
        pdf_file = file.path(
          output_dir,
          paste0(
            prefix,
            "_",
            measure,
            "_",
            side,
            "_ties.pdf"
          )
        ),
        plot_fun = function() {
          plot_win_ties(
            x = x,
            measure = measure,
            side = side
          )
        }
      )
    }
  }

  time_ids <- which(vapply(
    x$endpoint_specs,
    function(ep) identical(ep$type, "time"),
    logical(1)
  ))

  for (endpoint in time_ids) {
    .wr_save_dual_plot(
      png_file = file.path(
        output_dir,
        paste0(
          prefix,
          "_KM_outcome",
          endpoint,
          ".png"
        )
      ),
      pdf_file = file.path(
        output_dir,
        paste0(
          prefix,
          "_KM_outcome",
          endpoint,
          ".pdf"
        )
      ),
      plot_fun = function() {
        plot_win_km(
          x = x,
          endpoint = endpoint
        )
      }
    )
  }

  invisible(output_dir)
}

.wr_simulation_power_column <- function(side) {
  if (side == "one") {
    "rejection_proportion_one_sided"
  } else {
    "rejection_proportion_two_sided"
  }
}

plot_win_simulation_power <- function(x,
                                      measure = c("WR", "WO", "Log-rank"),
                                      side = c("one", "two"),
                                      main = NULL) {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  measure <- match.arg(measure)
  side <- match.arg(side)

  d <- x$power[
    x$power$measure == measure,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    stop("No simulation results found for measure: ", measure)
  }

  power_col <- .wr_simulation_power_column(side)
  scenarios <- unique(d$scenario_id)
  methods <- unique(d$notation)

  mat <- matrix(
    NA_real_,
    nrow = length(methods),
    ncol = length(scenarios),
    dimnames = list(
      methods,
      scenarios
    )
  )

  for (i in seq_len(nrow(d))) {
    mat[
      match(d$notation[i], methods),
      match(d$scenario_id[i], scenarios)
    ] <- d[[power_col]][i]
  }

  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat),
    zlim = c(0, 1),
    xaxt = "n",
    yaxt = "n",
    xlab = "Scenario",
    ylab = "Method",
    main = if (is.null(main)) {
      paste0(
        measure,
        ": ",
        side,
        "-sided rejection proportion"
      )
    } else {
      main
    }
  )

  graphics::axis(
    1,
    at = seq_len(ncol(mat)),
    labels = colnames(mat),
    las = 2,
    cex.axis = 0.65
  )

  graphics::axis(
    2,
    at = seq_len(nrow(mat)),
    labels = rownames(mat),
    las = 2,
    cex.axis = 0.7
  )

  invisible(mat)
}

plot_win_simulation_gain <- function(x,
                                     measure = c("WR", "WO"),
                                     side = c("one", "two"),
                                     baseline_method = "original",
                                     main = NULL) {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  measure <- match.arg(measure)
  side <- match.arg(side)

  d <- x$power[
    x$power$measure == measure,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    stop("No simulation results found for measure: ", measure)
  }

  power_col <- .wr_simulation_power_column(side)

  baseline <- d[
    d$method_id == baseline_method,
    c(
      "scenario_id",
      power_col
    ),
    drop = FALSE
  ]

  names(baseline)[2] <- "baseline_power"

  d <- merge(
    d,
    baseline,
    by = "scenario_id",
    all.x = TRUE,
    sort = FALSE
  )

  d$gain <- d[[power_col]] - d$baseline_power

  d <- d[
    d$method_id != baseline_method,
    ,
    drop = FALSE
  ]

  scenarios <- unique(d$scenario_id)
  methods <- unique(d$notation)

  mat <- matrix(
    NA_real_,
    nrow = length(methods),
    ncol = length(scenarios),
    dimnames = list(
      methods,
      scenarios
    )
  )

  for (i in seq_len(nrow(d))) {
    mat[
      match(d$notation[i], methods),
      match(d$scenario_id[i], scenarios)
    ] <- d$gain[i]
  }

  limit <- suppressWarnings(
    max(abs(mat), na.rm = TRUE)
  )

  if (!is.finite(limit) || limit == 0) {
    limit <- 1
  }

  graphics::image(
    x = seq_len(ncol(mat)),
    y = seq_len(nrow(mat)),
    z = t(mat),
    zlim = c(-limit, limit),
    xaxt = "n",
    yaxt = "n",
    xlab = "Scenario",
    ylab = "Method",
    main = if (is.null(main)) {
      paste0(
        measure,
        ": ",
        side,
        "-sided gain vs ",
        baseline_method
      )
    } else {
      main
    }
  )

  graphics::axis(
    1,
    at = seq_len(ncol(mat)),
    labels = colnames(mat),
    las = 2,
    cex.axis = 0.65
  )

  graphics::axis(
    2,
    at = seq_len(nrow(mat)),
    labels = rownames(mat),
    las = 2,
    cex.axis = 0.7
  )

  invisible(mat)
}

plot_win_simulation_scenario_gain <- function(x,
                                              scenario,
                                              measure = c("WR", "WO"),
                                              side = c("one", "two"),
                                              baseline_method = "original",
                                              main = NULL) {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  measure <- match.arg(measure)
  side <- match.arg(side)

  scenario_row <- .wr_resolve_scenario(
    scenario = scenario,
    scenario_grid = x$scenarios
  )

  scenario_id <- scenario_row$scenario_id[1]
  power_col <- .wr_simulation_power_column(side)

  d <- x$power[
    x$power$scenario_id == scenario_id &
      x$power$measure == measure,
    ,
    drop = FALSE
  ]

  if (nrow(d) == 0) {
    stop("No simulation results found for the requested scenario and measure.")
  }

  baseline <- d[
    d$method_id == baseline_method,
    power_col,
    drop = TRUE
  ]

  if (length(baseline) != 1L || !is.finite(baseline)) {
    stop("Baseline method is unavailable for this scenario.")
  }

  d <- d[
    d$method_id != baseline_method,
    ,
    drop = FALSE
  ]

  gain <- d[[power_col]] - baseline

  ord <- order(
    gain,
    decreasing = FALSE,
    na.last = TRUE
  )

  gain <- gain[ord]
  labels <- d$notation[ord]

  old_par <- graphics::par(no.readonly = TRUE)

  on.exit(
    graphics::par(old_par),
    add = TRUE
  )

  graphics::par(
    mar = c(5.2, 13, 3.2, 1.2)
  )

  graphics::barplot(
    gain,
    names.arg = labels,
    horiz = TRUE,
    las = 1,
    xlab = paste0(
      "Rejection-proportion gain vs ",
      baseline_method
    ),
    main = if (is.null(main)) {
      paste0(
        scenario_id,
        ": ",
        measure,
        " ",
        side,
        "-sided gain"
      )
    } else {
      main
    },
    cex.names = 0.8
  )

  graphics::abline(
    v = 0,
    lty = 2
  )

  invisible(d[ord, , drop = FALSE])
}

plot.win_simulation <- function(x,
                                type = c("power", "gain", "scenario_gain"),
                                measure = c("WR", "WO", "Log-rank"),
                                side = c("one", "two"),
                                baseline_method = "original",
                                scenario = 1L,
                                ...) {
  type <- match.arg(type)
  measure <- match.arg(measure)
  side <- match.arg(side)

  if (type == "power") {
    return(
      plot_win_simulation_power(
        x = x,
        measure = measure,
        side = side,
        ...
      )
    )
  }

  if (measure == "Log-rank") {
    stop("Gain plots require WR or WO.")
  }

  if (type == "scenario_gain") {
    return(
      plot_win_simulation_scenario_gain(
        x = x,
        scenario = scenario,
        measure = measure,
        side = side,
        baseline_method = baseline_method,
        ...
      )
    )
  }

  plot_win_simulation_gain(
    x = x,
    measure = measure,
    side = side,
    baseline_method = baseline_method,
    ...
  )
}

save_win_simulation_plots <- function(x,
                                      output_dir,
                                      prefix = "SIM") {
  if (!inherits(x, "win_simulation")) {
    stop("x must be a win_simulation object.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (measure in c("WR", "WO", "Log-rank")) {
    if (!any(x$power$measure == measure)) next

    for (side in c("one", "two")) {
      file <- file.path(
        output_dir,
        paste0(
          prefix,
          "_",
          gsub(
            "[^A-Za-z0-9]+",
            "_",
            measure
          ),
          "_",
          side,
          "_power.png"
        )
      )

      .wr_open_png(
        file,
        width = 3000,
        height = 1800
      )

      try(
        plot_win_simulation_power(
          x = x,
          measure = measure,
          side = side
        ),
        silent = TRUE
      )

      grDevices::dev.off()
    }
  }

  for (measure in c("WR", "WO")) {
    for (side in c("one", "two")) {
      file <- file.path(
        output_dir,
        paste0(
          prefix,
          "_",
          measure,
          "_",
          side,
          "_gain_heatmap.png"
        )
      )

      .wr_open_png(
        file,
        width = 3000,
        height = 1800
      )

      try(
        plot_win_simulation_gain(
          x = x,
          measure = measure,
          side = side
        ),
        silent = TRUE
      )

      grDevices::dev.off()
    }
  }

  scenario_gain_dir <- file.path(
    output_dir,
    "scenario_gain"
  )

  dir.create(
    scenario_gain_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in seq_len(nrow(x$scenarios))) {
    scenario_id <- x$scenarios$scenario_id[i]

    for (measure in c("WR", "WO")) {
      for (side in c("one", "two")) {
        file <- file.path(
          scenario_gain_dir,
          paste0(
            sprintf("%02d", x$scenarios$scenario_index[i]),
            "_",
            scenario_id,
            "_",
            measure,
            "_",
            side,
            "_gain.png"
          )
        )

        .wr_open_png(
          file,
          width = 2400,
          height = 1500
        )

        try(
          plot_win_simulation_scenario_gain(
            x = x,
            scenario = scenario_id,
            measure = measure,
            side = side
          ),
          silent = TRUE
        )

        grDevices::dev.off()
      }
    }
  }

  invisible(output_dir)
}
