.wr_read_data_file <- function(path, object_name = NULL) {
  if (is.null(path) || length(path) == 0 || !nzchar(as.character(path)[1])) {
    stop("A file path is required.")
  }
  path <- as.character(path)[1]
  if (!file.exists(path)) stop("File not found: ", path)

  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("csv", "txt")) {
    return(utils::read.csv(path, stringsAsFactors = FALSE))
  }

  if (ext == "rds") {
    return(readRDS(path))
  }

  if (ext %in% c("rdata", "rda")) {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)

    if (!is.null(object_name)) {
      object_name <- as.character(object_name)[1]
      if (!exists(object_name, envir = env, inherits = FALSE)) {
        stop("Object not found in RData file: ", object_name)
      }
      return(get(object_name, envir = env, inherits = FALSE))
    }

    if (length(loaded) == 1L) {
      return(get(loaded[1], envir = env, inherits = FALSE))
    }

    data_objects <- loaded[vapply(
      loaded,
      function(nm) {
        obj <- get(nm, envir = env, inherits = FALSE)
        is.data.frame(obj) || is.list(obj)
      },
      logical(1)
    )]

    if (length(data_objects) != 1L) {
      stop("RData file contains multiple candidate objects. Supply object_name.")
    }

    return(get(data_objects[1], envir = env, inherits = FALSE))
  }

  stop("Unsupported file extension: ", ext)
}

.wr_load_package_data <- function(package, dataset) {
  package <- trimws(as.character(package)[1])
  dataset <- trimws(as.character(dataset)[1])

  if (!nzchar(package) || !nzchar(dataset)) {
    stop("Both package and dataset are required.")
  }

  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Package is not installed: ", package)
  }

  env <- new.env(parent = emptyenv())
  try(utils::data(list = dataset, package = package, envir = env), silent = TRUE)

  if (!exists(dataset, envir = env, inherits = FALSE)) {
    stop("Dataset not found in package: ", package, "::", dataset)
  }

  get(dataset, envir = env, inherits = FALSE)
}

.wr_resolve_one_source <- function(data = NULL,
                                   file = NULL,
                                   package = NULL,
                                   dataset = NULL,
                                   object_name = NULL,
                                   label = "data") {
  has_data <- !is.null(data)
  has_file <- !is.null(file) && length(file) > 0 && nzchar(as.character(file)[1])
  has_package <- !is.null(package) && length(package) > 0 && nzchar(as.character(package)[1])

  n_sources <- sum(c(has_data, has_file, has_package))

  if (n_sources == 0L) return(NULL)
  if (n_sources > 1L) stop("Supply only one source for ", label, ".")

  if (has_data) return(data)
  if (has_file) return(.wr_read_data_file(file, object_name = object_name))

  .wr_load_package_data(package, dataset)
}

win_data <- function(data = NULL,
                     file = NULL,
                     package = NULL,
                     dataset = NULL,
                     object_name = NULL,
                     recurrent_data = NULL,
                     recurrent_file = NULL,
                     recurrent_package = NULL,
                     recurrent_dataset = NULL,
                     recurrent_object_name = NULL) {
  primary <- .wr_resolve_one_source(
    data = data,
    file = file,
    package = package,
    dataset = dataset,
    object_name = object_name,
    label = "primary data"
  )

  if (is.null(primary)) stop("A primary data source is required.")

  recurrent <- .wr_resolve_one_source(
    data = recurrent_data,
    file = recurrent_file,
    package = recurrent_package,
    dataset = recurrent_dataset,
    object_name = recurrent_object_name,
    label = "recurrent-event data"
  )

  structure(
    list(
      data = primary,
      recurrent = recurrent,
      source = list(
        primary = if (!is.null(data)) "object" else if (!is.null(file)) "file" else "package",
        recurrent = if (!is.null(recurrent_data)) {
          "object"
        } else if (!is.null(recurrent_file)) {
          "file"
        } else if (!is.null(recurrent_package)) {
          "package"
        } else {
          "none"
        }
      )
    ),
    class = c("win_data", "list")
  )
}

.wr_primary_data <- function(x) {
  if (inherits(x, "win_data")) return(x$data)
  x
}

.wr_recurrent_data <- function(x) {
  if (inherits(x, "win_data")) return(x$recurrent)
  NULL
}

prepare_win_data <- function(data,
                             endpoints,
                             format = c("wide", "event_long"),
                             id,
                             treatment,
                             event_time = NULL,
                             event_status = NULL,
                             event_time_unit = "years") {
  format <- match.arg(format)
  raw <- .wr_primary_data(data)

  if (!is.data.frame(raw)) {
    stop("The primary input must be a data frame for the general endpoint engine.")
  }

  recurrent <- .wr_recurrent_data(data)

  if (format == "wide") {
    return(.wr_prepare_wide_data(
      raw = raw,
      id_col = id,
      arm_col = treatment,
      endpoints = endpoints,
      recurrent_data = recurrent
    ))
  }

  if (is.null(event_time) || is.null(event_status)) {
    stop("event_time and event_status are required for event_long data.")
  }

  .wr_prepare_event_long_data(
    raw = raw,
    id_col = id,
    arm_col = treatment,
    time_col = event_time,
    status_col = event_status,
    endpoints = endpoints,
    time_unit = event_time_unit
  )
}

.wr_normalize_arm <- function(x,
                              treatment_values = c("1", "treatment", "treated", "trt", "active", "intervention", "digoxin"),
                              control_values = c("0", "control", "ctrl", "placebo", "standard", "usual care")) {
  if (is.numeric(x) || is.integer(x)) {
    z <- suppressWarnings(as.integer(x))
    if (all(is.na(z) | z %in% c(0L, 1L))) return(z)
  }

  sx <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))
  out[sx %in% tolower(treatment_values)] <- 1L
  out[sx %in% tolower(control_values)] <- 0L

  if (any(is.na(out) & !is.na(sx) & nzchar(sx))) {
    stop("Treatment column could not be converted to 0/1 labels.")
  }

  out
}

.wr_parse_numeric_times <- function(x) {
  if (is.null(x) || length(x) == 0) return(numeric(0))
  if (is.list(x)) x <- x[[1]]
  if (length(x) > 1 && is.numeric(x)) {
    z <- as.numeric(x)
    return(z[is.finite(z)])
  }
  if (all(is.na(x))) return(numeric(0))
  if (is.numeric(x)) {
    z <- as.numeric(x)
    return(z[is.finite(z)])
  }

  s <- as.character(x[1])
  if (is.na(s) || !nzchar(trimws(s))) return(numeric(0))
  parts <- unlist(strsplit(s, split = "[;,|[:space:]]+"))
  parts <- parts[nzchar(parts)]
  z <- suppressWarnings(as.numeric(parts))
  z[is.finite(z)]
}

.wr_absolute_times_to_gaps <- function(times, followup) {
  times <- sort(as.numeric(times[is.finite(times)]))
  times <- times[times > 0 & times <= followup]
  if (length(times) == 0) return(numeric(0))
  diff(c(0, times))
}

.wr_gap_times_to_gaps <- function(gaps, followup) {
  gaps <- as.numeric(gaps[is.finite(gaps) & gaps > 0])
  if (length(gaps) == 0) return(numeric(0))
  abs_times <- cumsum(gaps)
  keep <- abs_times <= followup
  abs_keep <- abs_times[keep]
  if (length(abs_keep) == 0) return(numeric(0))
  diff(c(0, abs_keep))
}

.wr_evenly_spaced_gaps_from_count <- function(n, followup) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0 || is.na(followup) || followup <= 0) return(numeric(0))
  abs_times <- seq(
    followup / (n + 1),
    followup * n / (n + 1),
    length.out = n
  )
  diff(c(0, abs_times))
}

.wr_make_recurrent_gap_list <- function(subject_table,
                                        recurrent_data = NULL,
                                        id_col = "SUBJID",
                                        recurrent_id_col = "SUBJID",
                                        recurrent_time_col = "HOSPTIME",
                                        recurrent_time_unit = "years",
                                        recurrent_time_type = c("absolute", "gap"),
                                        subject_times_col = NULL,
                                        legacy_count_approximation = FALSE,
                                        count_col = "NUMHOSP") {
  recurrent_time_type <- match.arg(recurrent_time_type)
  n <- nrow(subject_table)

  if (!is.null(recurrent_data)) {
    if (!is.data.frame(recurrent_data)) recurrent_data <- as.data.frame(recurrent_data)
    if (!(recurrent_id_col %in% names(recurrent_data))) {
      stop("Recurrent subject ID column not found: ", recurrent_id_col)
    }
    if (!(recurrent_time_col %in% names(recurrent_data))) {
      stop("Recurrent event-time column not found: ", recurrent_time_col)
    }

    subject_ids <- as.character(subject_table[[id_col]])
    recurrent_ids <- as.character(recurrent_data[[recurrent_id_col]])
    recurrent_times <- .wr_time_to_years(recurrent_data[[recurrent_time_col]], recurrent_time_unit)

    out <- vector("list", n)
    for (i in seq_len(n)) {
      z <- recurrent_times[recurrent_ids == subject_ids[i]]
      if (recurrent_time_type == "absolute") {
        out[[i]] <- .wr_absolute_times_to_gaps(z, subject_table$FUTIME[i])
      } else {
        out[[i]] <- .wr_gap_times_to_gaps(z, subject_table$FUTIME[i])
      }
    }
    return(out)
  }

  if (!is.null(subject_times_col) && subject_times_col %in% names(subject_table)) {
    out <- vector("list", n)
    for (i in seq_len(n)) {
      z <- .wr_parse_numeric_times(subject_table[[subject_times_col]][i])
      z <- .wr_time_to_years(z, recurrent_time_unit)
      if (recurrent_time_type == "absolute") {
        out[[i]] <- .wr_absolute_times_to_gaps(z, subject_table$FUTIME[i])
      } else {
        out[[i]] <- .wr_gap_times_to_gaps(z, subject_table$FUTIME[i])
      }
    }
    return(out)
  }

  if (isTRUE(legacy_count_approximation) && count_col %in% names(subject_table)) {
    out <- vector("list", n)
    for (i in seq_len(n)) {
      out[[i]] <- .wr_evenly_spaced_gaps_from_count(
        subject_table[[count_col]][i],
        subject_table$FUTIME[i]
      )
    }
    return(out)
  }

  NULL
}

prepare_two_endpoint_recurrent_data <- function(data,
                                                id = "SUBJID",
                                                treatment = "ARM",
                                                followup = "FUTIME",
                                                event = "CNSR",
                                                count = "NUMHOSP",
                                                recurrent_time = "HOSP_TIMES",
                                                time_unit = "years",
                                                recurrent_time_unit = "years",
                                                recurrent_time_type = c("absolute", "gap"),
                                                recurrent_id = "SUBJID",
                                                recurrent_event_time = "HOSPTIME",
                                                legacy_count_approximation = FALSE) {
  recurrent_time_type <- match.arg(recurrent_time_type)
  raw <- .wr_primary_data(data)
  recurrent <- .wr_recurrent_data(data)

  if (is.list(raw) && all(c("table.output", "hosp.times.list") %in% names(raw))) {
    ds <- raw
    ds$table.output$ARM <- .wr_normalize_arm(ds$table.output$ARM)
    ds$table.output$FUTIME <- as.numeric(ds$table.output$FUTIME)
    ds$table.output$CNSR <- as.integer(ds$table.output$CNSR)
    if (!("SUBJID" %in% names(ds$table.output))) {
      ds$table.output$SUBJID <- seq_len(nrow(ds$table.output))
    }
    ds$table.output$NUMHOSP <- as.integer(lengths(ds$hosp.times.list))
    ds$has_recurrent_times <- TRUE
    return(ds)
  }

  if (!is.data.frame(raw)) raw <- as.data.frame(raw)

  required <- c(treatment, followup, event)
  missing_cols <- setdiff(required, names(raw))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  ids <- if (!is.null(id) && id %in% names(raw)) raw[[id]] else seq_len(nrow(raw))
  arm <- .wr_normalize_arm(raw[[treatment]])
  futime <- .wr_time_to_years(raw[[followup]], time_unit)
  cnsr <- as.integer(.wr_safe_num(raw[[event]]))

  if (any(!is.na(cnsr) & !(cnsr %in% c(0L, 1L)))) {
    stop("Event indicator must be coded 1 = event and 0 = censored.")
  }

  num_hosp <- if (!is.null(count) && count %in% names(raw)) {
    as.integer(.wr_safe_num(raw[[count]]))
  } else {
    rep(NA_integer_, nrow(raw))
  }

  table_output <- data.frame(
    SUBJID = as.character(ids),
    ARM = arm,
    FUTIME = as.numeric(futime),
    CNSR = cnsr,
    NUMHOSP = num_hosp,
    stringsAsFactors = FALSE
  )

  if (any(is.na(table_output$ARM))) stop("Treatment assignment contains missing values.")
  if (any(is.na(table_output$FUTIME)) || any(table_output$FUTIME < 0)) {
    stop("Follow-up time must be non-missing and non-negative.")
  }

  subject_times_col <- if (!is.null(recurrent_time) && recurrent_time %in% names(raw)) recurrent_time else NULL

  subject_with_times <- table_output
  if (!is.null(subject_times_col)) {
    subject_with_times[[subject_times_col]] <- raw[[subject_times_col]]
  }

  hosp_times_list <- .wr_make_recurrent_gap_list(
    subject_table = subject_with_times,
    recurrent_data = recurrent,
    id_col = "SUBJID",
    recurrent_id_col = recurrent_id,
    recurrent_time_col = recurrent_event_time,
    recurrent_time_unit = recurrent_time_unit,
    recurrent_time_type = recurrent_time_type,
    subject_times_col = subject_times_col,
    legacy_count_approximation = legacy_count_approximation,
    count_col = "NUMHOSP"
  )

  has_recurrent_times <- !is.null(hosp_times_list)

  if (has_recurrent_times) {
    table_output$NUMHOSP <- as.integer(lengths(hosp_times_list))
  }

  list(
    table.output = table_output,
    hosp.times.list = hosp_times_list,
    has_recurrent_times = has_recurrent_times,
    count_only_available = any(!is.na(table_output$NUMHOSP))
  )
}

.wr_prepare_recurrent_cache <- function(ds) {
  if (is.null(ds$hosp.times.list)) {
    stop("Recurrent event times are required for pairwise common-follow-up comparison.")
  }

  n <- nrow(ds$table.output)
  hosp_abs_times_list <- vector("list", n)

  for (i in seq_len(n)) {
    x <- ds$hosp.times.list[[i]]
    if (length(x) == 0 || all(is.na(x))) {
      hosp_abs_times_list[[i]] <- numeric(0)
    } else {
      x <- as.numeric(x[!is.na(x)])
      hosp_abs_times_list[[i]] <- cumsum(x)
    }
  }

  lens <- as.integer(lengths(hosp_abs_times_list))
  starts <- if (n == 0) integer(0) else as.integer(cumsum(c(0L, lens[-n])))
  flat <- as.numeric(unlist(hosp_abs_times_list, use.names = FALSE))
  if (length(flat) == 0) flat <- numeric(0)

  ds$hosp.abs.times.list <- hosp_abs_times_list
  ds$hosp.flat <- flat
  ds$hosp.start <- starts
  ds$hosp.len <- lens
  ds
}
