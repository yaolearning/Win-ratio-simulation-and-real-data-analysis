install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_if_missing("survival")
install_if_missing("Rcpp")

library(survival)
library(Rcpp)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

OUTDIR <- "WR outputs"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Main real-data controls.
if (!exists("B_REAL")) B_REAL <- 2000
if (!exists("BATCH_SIZE_REAL")) BATCH_SIZE_REAL <- 50
if (!exists("MASTER_SEED")) MASTER_SEED <- 2026
if (!exists("ALPHA")) ALPHA <- 0.05
if (!exists("RESUME_IF_EXISTS")) RESUME_IF_EXISTS <- TRUE
if (!exists("VERBOSE")) VERBOSE <- TRUE
if (!exists("AUTO_RUN")) AUTO_RUN <- TRUE

# Weighted p grids.
# p is the weight assigned to the first endpoint in the selected hierarchy.
# For death-first order: p weights death/survival; 1-p weights hospitalization.
# For hospitalization-first order: p weights hospitalization; 1-p weights death.
# No p = 0 value is used. Values below 0.5 are exploratory only;
# the primary weighted analysis remains p in [0.5, 1].
P_GRID_EXPLORATORY <- seq(0.01, 0.49, by = 0.01)
P_GRID_LOW <- P_GRID_EXPLORATORY
P_GRID_PRIMARY <- seq(0.50, 1.00, by = 0.01)
P_GRID_ALL <- sort(unique(c(P_GRID_EXPLORATORY, P_GRID_PRIMARY)))

# Clinically interpretable backup candidates for a time-difference threshold t.
# These candidates are combined with first-endpoint data-informed candidates.
# They are sensitivity candidates, not universal margins; update them for the
# disease area and statistical analysis plan if stronger disease-specific
# candidates are available.
CLINICAL_T_MONTHS <- c(1, 3, 6, 12, 18, 24)
MAX_T_GRID_SIZE <- 10
EPS_WR <- 1e-8

#data import blocks

# Subject-level file format expected by default:
#   SUBJID    subject ID
#   ARM       0 = control, 1 = treatment; text labels are also accepted
#   FUTIME    follow-up time
#   CNSR      death/event indicator for the first endpoint, 1 = death, 0 = censored
#   NUMHOSP   optional hospitalization count
#   HOSP_TIMES optional hospitalization times, e.g. "0.2;0.6;1.1"
#
# Recurrent hospitalization file is optional. If supplied, it should contain:
#   SUBJID    subject ID
#   HOSPTIME  hospitalization time
#
# Times are converted to years using time_unit / hosp_time_unit.
# hosp_time_type = "absolute" means HOSP_TIMES/HOSPTIME are times from baseline.
# hosp_time_type = "gap" means they are gap times between hospitalizations.

if (!exists("REAL_DATASETS")) {
  REAL_DATASETS <- data.frame(
    dataset_id = c("realdata1"),
    dataset_label = c("Real data 1"),

    # source_type = "file" reads subject_file / recurrent_file.
    # source_type = "package" reads data_object from package_name using data().
    source_type = c("file"),

    # File input fields. Leave blank when source_type = "package".
    subject_file = c("realdata_subjects.csv"),
    recurrent_file = c(""),

    # R-package input fields. Example:
    # source_type = "package", package_name = "asympTest", data_object = "DIGdata"
    package_name = c(""),
    data_object = c(""),

    # Optional recurrent hospitalization data from an R package.
    # Leave recurrent_source_type = "file" if recurrent_file is a CSV/RDS path.
    # Use recurrent_source_type = "package" if recurrent_data_object is stored in an R package.
    recurrent_source_type = c("file"),
    recurrent_package_name = c(""),
    recurrent_data_object = c(""),

    time_unit = c("years"),
    hosp_time_unit = c("years"),
    hosp_time_type = c("absolute"),
    stringsAsFactors = FALSE
  )
}

# Column maps can be edited once if all imported datasets share the same names.
# If a column is unavailable, set it to "". The first four columns are required.
if (!exists("COLUMN_MAP")) {
  COLUMN_MAP <- list(
    id = "SUBJID",
    arm = "ARM",
    followup = "FUTIME",
    death_event = "CNSR",
    num_hosp = "NUMHOSP",
    hosp_times = "HOSP_TIMES",
    survival_time = "SURVTIME",
    censor_time = "CNSRTIME",
    freq_hosp = "FREQHOSP"
  )
}

if (!exists("RECURRENT_COLUMN_MAP")) {
  RECURRENT_COLUMN_MAP <- list(
    id = "SUBJID",
    hosp_time = "HOSPTIME"
  )
}

# Text values interpreted as treatment/control if ARM is not already coded 0/1.
ARM_TREATMENT_VALUES <- c("1", "treatment", "treated", "trt", "active", "intervention", "digoxin")
ARM_CONTROL_VALUES <- c("0", "control", "ctrl", "placebo", "standard", "usual care")

#data import helpers
sanitize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nchar(x) == 0, "dataset", x)
}

row_value <- function(row, name, default = "") {
  if (!(name %in% names(row))) return(default)
  value <- row[[name]][1]
  if (length(value) == 0 || is.na(value)) return(default)
  as.character(value)
}

time_to_years <- function(x, unit = "years") {
  unit <- tolower(as.character(unit))
  x <- as.numeric(x)
  if (unit %in% c("year", "years", "yr", "yrs")) return(x)
  if (unit %in% c("month", "months", "mo", "mos")) return(x / 12)
  if (unit %in% c("day", "days", "d")) return(x / 365.25)
  stop("Unsupported time unit: ", unit)
}

read_table_or_rds <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("csv", "txt")) return(utils::read.csv(path, stringsAsFactors = FALSE))
  stop("Unsupported file type for ", path, ". Use .csv, .txt, or .rds.")
}

load_package_data_object <- function(package.name, object.name, dataset.id = "realdata") {
  package.name <- trimws(as.character(package.name))
  object.name <- trimws(as.character(object.name))
  if (package.name == "" || object.name == "") {
    stop("Package input for ", dataset.id, " requires both package_name and data_object.")
  }
  if (!requireNamespace(package.name, quietly = TRUE)) {
    stop("Package not installed for ", dataset.id, ": ", package.name)
  }
  env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    utils::data(list = object.name, package = package.name, envir = env),
    error = function(e) e
  )
  if (inherits(loaded, "error") || !exists(object.name, envir = env, inherits = FALSE)) {
    stop("Could not load data object '", object.name, "' from package '", package.name,
         "' for ", dataset.id, ".")
  }
  get(object.name, envir = env, inherits = FALSE)
}

col_or_null <- function(df, colname) {
  if (is.null(colname) || is.na(colname) || colname == "") return(NULL)
  if (!(colname %in% names(df))) return(NULL)
  df[[colname]]
}

require_col <- function(df, colname, role) {
  if (is.null(colname) || is.na(colname) || colname == "" || !(colname %in% names(df))) {
    stop("Missing required ", role, " column: ", colname)
  }
  df[[colname]]
}

parse_numeric_times <- function(x) {
  if (is.null(x) || length(x) == 0) return(numeric(0))
  if (is.list(x)) x <- x[[1]]
  if (length(x) > 1 && is.numeric(x)) return(as.numeric(x[is.finite(x)]))
  if (all(is.na(x))) return(numeric(0))
  if (is.numeric(x)) {
    x <- as.numeric(x)
    return(x[is.finite(x)])
  }
  s <- as.character(x[1])
  if (is.na(s) || trimws(s) == "") return(numeric(0))
  parts <- unlist(strsplit(s, split = "[;,|[:space:]]+"))
  parts <- parts[parts != ""]
  z <- suppressWarnings(as.numeric(parts))
  z <- z[is.finite(z)]
  as.numeric(z)
}

normalize_arm <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    z <- as.integer(x)
    if (all(z %in% c(0L, 1L))) return(z)
  }
  sx <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(sx))
  out[sx %in% ARM_TREATMENT_VALUES] <- 1L
  out[sx %in% ARM_CONTROL_VALUES] <- 0L
  if (any(is.na(out))) {
    lev <- unique(sx[is.na(out)])
    lev <- lev[!is.na(lev) & lev != ""]
    if (length(lev) == 2) {
      sorted <- sort(lev)
      out[sx == sorted[1]] <- 0L
      out[sx == sorted[2]] <- 1L
      warning("ARM values were not recognized as 0/1 or common labels. Using alphabetical coding: ",
              sorted[1], " = control, ", sorted[2], " = treatment.")
    }
  }
  if (any(is.na(out))) stop("ARM column could not be converted to 0/1 treatment labels.")
  out
}

absolute_times_to_gaps <- function(times, followup) {
  times <- sort(as.numeric(times[is.finite(times)]))
  times <- times[times > 0 & times <= followup]
  if (length(times) == 0) return(numeric(0))
  diff(c(0, times))
}

gap_times_to_gaps <- function(gaps, followup) {
  gaps <- as.numeric(gaps[is.finite(gaps) & gaps > 0])
  if (length(gaps) == 0) return(numeric(0))
  abs.times <- cumsum(gaps)
  keep <- abs.times <= followup
  abs.keep <- abs.times[keep]
  if (length(abs.keep) == 0) return(numeric(0))
  diff(c(0, abs.keep))
}

evenly_spaced_gaps_from_count <- function(nh, followup) {
  nh <- as.integer(nh)
  if (is.na(nh) || nh <= 0 || is.na(followup) || followup <= 0) return(numeric(0))
  abs.times <- seq(followup / (nh + 1), followup * nh / (nh + 1), length.out = nh)
  diff(c(0, abs.times))
}

make_hosp_list_from_recurrent <- function(tab, recurrent.df, id.col, time.col,
                                          time.unit = "years", time.type = "absolute") {
  if (is.null(recurrent.df) || nrow(recurrent.df) == 0) return(NULL)
  ids <- as.character(tab$SUBJID)
  rid <- as.character(require_col(recurrent.df, id.col, "recurrent subject ID"))
  rtime <- time_to_years(require_col(recurrent.df, time.col, "recurrent hospitalization time"), time.unit)
  hosp.list <- vector("list", nrow(tab))
  time.type <- tolower(time.type)
  for (i in seq_len(nrow(tab))) {
    z <- rtime[rid == ids[i]]
    if (time.type == "absolute") hosp.list[[i]] <- absolute_times_to_gaps(z, tab$FUTIME[i])
    else if (time.type == "gap") hosp.list[[i]] <- gap_times_to_gaps(z, tab$FUTIME[i])
    else stop("hosp_time_type must be 'absolute' or 'gap'.")
  }
  hosp.list
}

standardize.real.subjects <- function(subject.df,
                                      recurrent.df = NULL,
                                      column.map = COLUMN_MAP,
                                      recurrent.column.map = RECURRENT_COLUMN_MAP,
                                      time.unit = "years",
                                      hosp.time.unit = "years",
                                      hosp.time.type = "absolute",
                                      dataset.id = "realdata") {
  if (!is.data.frame(subject.df)) subject.df <- as.data.frame(subject.df)

  # Subject ID is useful but not required. If the package data has no ID column,
  # set COLUMN_MAP$id <- "" and row numbers will be used as SUBJID.
  if (is.null(column.map$id) || is.na(column.map$id) || column.map$id == "" ||
      !(column.map$id %in% names(subject.df))) {
    id <- seq_len(nrow(subject.df))
    warning("No subject ID column supplied/found for ", dataset.id,
            "; using row number as SUBJID.")
  } else {
    id <- require_col(subject.df, column.map$id, "subject ID")
  }

  arm <- normalize_arm(require_col(subject.df, column.map$arm, "arm"))
  futime <- time_to_years(require_col(subject.df, column.map$followup, "follow-up"), time.unit)
  cnsr <- as.integer(require_col(subject.df, column.map$death_event, "death/event indicator"))
  if (!all(cnsr %in% c(0L, 1L))) stop("Death/event indicator must be coded 1 = event/death, 0 = censored.")

  n <- length(id)
  num.hosp.input <- col_or_null(subject.df, column.map$num_hosp)
  hosp.times.input <- col_or_null(subject.df, column.map$hosp_times)
  surv.time.input <- col_or_null(subject.df, column.map$survival_time)
  cens.time.input <- col_or_null(subject.df, column.map$censor_time)
  freq.hosp.input <- col_or_null(subject.df, column.map$freq_hosp)

  table.output <- data.frame(
    SUBJID = as.character(id),
    ARM = as.integer(arm),
    FUTIME = as.numeric(futime),
    CNSR = as.integer(cnsr),
    SURVTIME = if (!is.null(surv.time.input)) time_to_years(surv.time.input, time.unit) else as.numeric(futime),
    CNSRTIME = if (!is.null(cens.time.input)) time_to_years(cens.time.input, time.unit) else as.numeric(futime),
    FREQHOSP = if (!is.null(freq.hosp.input)) as.numeric(freq.hosp.input) else NA_real_,
    NUMHOSP = if (!is.null(num.hosp.input)) as.integer(num.hosp.input) else 0L,
    stringsAsFactors = FALSE
  )

  if (any(is.na(table.output$FUTIME)) || any(table.output$FUTIME < 0)) {
    stop("FUTIME must be non-missing and non-negative after time-unit conversion.")
  }

  hosp.times.list <- NULL
  if (!is.null(recurrent.df)) {
    hosp.times.list <- make_hosp_list_from_recurrent(
      tab = table.output,
      recurrent.df = recurrent.df,
      id.col = recurrent.column.map$id,
      time.col = recurrent.column.map$hosp_time,
      time.unit = hosp.time.unit,
      time.type = hosp.time.type
    )
  }

  if (is.null(hosp.times.list) && !is.null(hosp.times.input)) {
    hosp.times.list <- vector("list", n)
    for (i in seq_len(n)) {
      z <- parse_numeric_times(hosp.times.input[i])
      z <- time_to_years(z, hosp.time.unit)
      if (tolower(hosp.time.type) == "absolute") hosp.times.list[[i]] <- absolute_times_to_gaps(z, table.output$FUTIME[i])
      else if (tolower(hosp.time.type) == "gap") hosp.times.list[[i]] <- gap_times_to_gaps(z, table.output$FUTIME[i])
      else stop("hosp_time_type must be 'absolute' or 'gap'.")
    }
  }

  if (is.null(hosp.times.list)) {
    if (!is.null(num.hosp.input)) {
      warning("No recurrent hospitalization times were supplied for ", dataset.id,
              ". Hospitalization times are evenly spaced from NUMHOSP as an approximation. ",
              "Provide recurrent_file or HOSP_TIMES for exact pair-specific truncation.")
      hosp.times.list <- vector("list", n)
      for (i in seq_len(n)) {
        hosp.times.list[[i]] <- evenly_spaced_gaps_from_count(table.output$NUMHOSP[i], table.output$FUTIME[i])
      }
    } else {
      hosp.times.list <- replicate(n, numeric(0), simplify = FALSE)
    }
  }

  table.output$NUMHOSP <- as.integer(lengths(hosp.times.list))

  list(
    table.output = table.output,
    hosp.times.list = hosp.times.list,
    import.note = data.frame(
      dataset_id = dataset.id,
      n_subjects = n,
      n_control = sum(table.output$ARM == 0),
      n_treatment = sum(table.output$ARM == 1),
      total_hospitalizations = sum(table.output$NUMHOSP, na.rm = TRUE),
      time_unit = time.unit,
      hosp_time_unit = hosp.time.unit,
      hosp_time_type = hosp.time.type,
      stringsAsFactors = FALSE
    )
  )
}

load.real.dataset <- function(registry.row,
                              column.map = COLUMN_MAP,
                              recurrent.column.map = RECURRENT_COLUMN_MAP) {
  dataset.id <- sanitize_id(row_value(registry.row, "dataset_id", "realdata"))
  source.type <- tolower(row_value(registry.row, "source_type", "file"))
  subject.file <- row_value(registry.row, "subject_file", "")
  recurrent.file <- row_value(registry.row, "recurrent_file", "")
  package.name <- row_value(registry.row, "package_name", "")
  data.object <- row_value(registry.row, "data_object", "")
  recurrent.source.type <- tolower(row_value(registry.row, "recurrent_source_type", "file"))
  recurrent.package.name <- row_value(registry.row, "recurrent_package_name", package.name)
  recurrent.data.object <- row_value(registry.row, "recurrent_data_object", "")
  time.unit <- row_value(registry.row, "time_unit", "years")
  hosp.time.unit <- row_value(registry.row, "hosp_time_unit", time.unit)
  hosp.time.type <- row_value(registry.row, "hosp_time_type", "absolute")

  if (source.type %in% c("package", "r_package", "data")) {
    subject.obj <- load_package_data_object(package.name, data.object, dataset.id = dataset.id)
  } else {
    if (subject.file == "" || !file.exists(subject.file)) {
      stop("Subject-level file not found for ", dataset.id, ": ", subject.file)
    }
    subject.obj <- read_table_or_rds(subject.file)
  }

  if (is.list(subject.obj) && all(c("table.output", "hosp.times.list") %in% names(subject.obj))) {
    ds <- subject.obj
    ds$table.output$ARM <- normalize_arm(ds$table.output$ARM)
    ds$table.output$FUTIME <- as.numeric(ds$table.output$FUTIME)
    ds$table.output$CNSR <- as.integer(ds$table.output$CNSR)
    if (!("SUBJID" %in% names(ds$table.output))) ds$table.output$SUBJID <- seq_len(nrow(ds$table.output))
    if (!("SURVTIME" %in% names(ds$table.output))) ds$table.output$SURVTIME <- ds$table.output$FUTIME
    if (!("CNSRTIME" %in% names(ds$table.output))) ds$table.output$CNSRTIME <- ds$table.output$FUTIME
    if (!("FREQHOSP" %in% names(ds$table.output))) ds$table.output$FREQHOSP <- NA_real_
    ds$table.output$NUMHOSP <- as.integer(lengths(ds$hosp.times.list))
    ds$import.note <- data.frame(
      dataset_id = dataset.id,
      n_subjects = nrow(ds$table.output),
      n_control = sum(ds$table.output$ARM == 0),
      n_treatment = sum(ds$table.output$ARM == 1),
      total_hospitalizations = sum(ds$table.output$NUMHOSP, na.rm = TRUE),
      time_unit = "already standardized to years",
      hosp_time_unit = "already standardized to years",
      hosp_time_type = "gap list supplied",
      stringsAsFactors = FALSE
    )
    return(ds)
  }

  recurrent.df <- NULL
  if (recurrent.source.type %in% c("package", "r_package", "data") && recurrent.data.object != "") {
    if (recurrent.package.name == "") recurrent.package.name <- package.name
    recurrent.df <- load_package_data_object(
      recurrent.package.name,
      recurrent.data.object,
      dataset.id = paste0(dataset.id, " recurrent data")
    )
  } else if (recurrent.file != "" && file.exists(recurrent.file)) {
    recurrent.df <- read_table_or_rds(recurrent.file)
  }

  standardize.real.subjects(
    subject.df = subject.obj,
    recurrent.df = recurrent.df,
    column.map = column.map,
    recurrent.column.map = recurrent.column.map,
    time.unit = time.unit,
    hosp.time.unit = hosp.time.unit,
    hosp.time.type = hosp.time.type,
    dataset.id = dataset.id
  )
}

write.realdata.templates <- function(outdir = OUTDIR) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  subject.template <- data.frame(
    SUBJID = c("001", "002"),
    ARM = c(0, 1),
    FUTIME = c(1.0, 1.0),
    CNSR = c(0, 1),
    NUMHOSP = c(1, 2),
    HOSP_TIMES = c("0.4", "0.2;0.7"),
    stringsAsFactors = FALSE
  )
  recurrent.template <- data.frame(
    SUBJID = c("001", "002", "002"),
    HOSPTIME = c(0.4, 0.2, 0.7),
    stringsAsFactors = FALSE
  )
  write.csv(subject.template, file.path(outdir, "TEMPLATE_subject_level_file.csv"), row.names = FALSE)
  write.csv(recurrent.template, file.path(outdir, "TEMPLATE_recurrent_hospitalization_file.csv"), row.names = FALSE)
}

#Preparation
prepare.ds.fast <- function(ds) {
  # Convert hospitalization gap times to absolute hospitalization times.
  # This avoids repeated cumsum() inside pairwise loops.
  n <- nrow(ds$table.output)
  hosp.abs.times.list <- vector("list", n)

  for (i in seq_len(n)) {
    x <- ds$hosp.times.list[[i]]
    if (length(x) == 0 || all(is.na(x))) {
      hosp.abs.times.list[[i]] <- numeric(0)
    } else {
      x <- as.numeric(x[!is.na(x)])
      hosp.abs.times.list[[i]] <- cumsum(x)
    }
  }

  lens <- as.integer(lengths(hosp.abs.times.list))
  starts <- if (n == 0) integer(0) else as.integer(cumsum(c(0L, lens[-n])))
  flat <- as.numeric(unlist(hosp.abs.times.list, use.names = FALSE))
  if (length(flat) == 0) flat <- numeric(0)

  ds$hosp.abs.times.list <- hosp.abs.times.list
  ds$hosp.flat <- flat
  ds$hosp.start <- starts  # zero-based starts for C++
  ds$hosp.len <- lens
  ds
}

safe_mean <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}

safe_quantile <- function(x, p) {
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

safe_wr <- function(win.score, loss.score, total.pairs, eps = EPS_WR) {
  ((win.score / total.pairs) + eps) / ((loss.score / total.pairs) + eps)
}

extract_nearest <- function(x, value) {
  which.min(abs(x - value))
}


#fast Rcpp for U matrix

Rcpp::sourceCpp(code = '
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
using namespace Rcpp;

int count_hosp_until_cpp(const NumericVector& hosp_times,
                         const IntegerVector& hosp_start,
                         const IntegerVector& hosp_len,
                         int idx,
                         double t) {
  int len = hosp_len[idx];
  if (len <= 0) return 0;
  int start = hosp_start[idx];
  int lo = 0;
  int hi = len;
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    if (hosp_times[start + mid] <= t) {
      lo = mid + 1;
    } else {
      hi = mid;
    }
  }
  return lo;
}

double safe_wr_cpp(double win_score, double loss_score, double total_pairs, double eps) {
  return ((win_score / total_pairs) + eps) / ((loss_score / total_pairs) + eps);
}

// [[Rcpp::export]]
List fast_wr_core_revised_cpp(NumericVector futime,
                              IntegerVector cnsr,
                              IntegerVector arm,
                              NumericVector hosp_times,
                              IntegerVector hosp_start,
                              IntegerVector hosp_len,
                              NumericVector p_grid,
                              NumericVector t_grid,
                              double eps = 1e-8) {

  int n = futime.size();
  std::vector<int> trt;
  std::vector<int> ctrl;
  trt.reserve(n);
  ctrl.reserve(n);

  for (int i = 0; i < n; ++i) {
    if (arm[i] == 1) trt.push_back(i);
    if (arm[i] == 0) ctrl.push_back(i);
  }

  int n_trt = trt.size();
  int n_ctrl = ctrl.size();
  double total_pairs = static_cast<double>(n_trt) * static_cast<double>(n_ctrl);
  if (total_pairs <= 0) {
    stop("Need at least one treatment patient and one control patient.");
  }

  int K = t_grid.size();
  int P = p_grid.size();

  // Death-first hierarchy counts.
  double D1_win = 0.0;
  double D1_loss = 0.0;
  double H2_win = 0.0;
  double H2_loss = 0.0;

  // Hospitalization-first hierarchy counts.
  double H1_win = 0.0;
  double H1_loss = 0.0;
  double D2_win = 0.0;
  double D2_loss = 0.0;

  // Endpoint signs and true unresolved ties.
  double death_win_pairs = 0.0;
  double death_loss_pairs = 0.0;
  double hosp_win_pairs = 0.0;
  double hosp_loss_pairs = 0.0;
  double true_tie_pairs = 0.0;

  std::vector<double> win_diff_t(K + 1, 0.0);
  std::vector<double> loss_diff_t(K + 1, 0.0);

  auto add_range = [&](std::vector<double>& diff, int start, int end, double value) {
    if (start > end) return;
    diff[start] += value;
    if (end + 1 < static_cast<int>(diff.size())) diff[end + 1] -= value;
  };

  for (int ii = 0; ii < n_trt; ++ii) {
    int pt1 = trt[ii];
    double fu1 = futime[pt1];
    int c1_original = cnsr[pt1];

    for (int jj = 0; jj < n_ctrl; ++jj) {
      int pt2 = ctrl[jj];
      double fu2 = futime[pt2];
      int c2_original = cnsr[pt2];
      double min_fu = fu1 < fu2 ? fu1 : fu2;

      int hosp1 = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, pt1, min_fu);
      int hosp2 = count_hosp_until_cpp(hosp_times, hosp_start, hosp_len, pt2, min_fu);

      int hosp_sign = 0;
      if (hosp2 > hosp1) hosp_sign = 1;       // treatment has fewer hospitalizations
      else if (hosp2 < hosp1) hosp_sign = -1; // treatment has more hospitalizations

      int c1 = c1_original;
      int c2 = c2_original;
      if (fu1 < fu2) c2 = 0;
      if (fu2 < fu1) c1 = 0;

      int death_sign = 0;
      if (c1 == 0 && c2 == 1) death_sign = 1;       // treatment survives longer
      else if (c1 == 1 && c2 == 0) death_sign = -1; // treatment dies earlier

      if (death_sign > 0) death_win_pairs += 1.0;
      if (death_sign < 0) death_loss_pairs += 1.0;
      if (hosp_sign > 0) hosp_win_pairs += 1.0;
      if (hosp_sign < 0) hosp_loss_pairs += 1.0;
      if (death_sign == 0 && hosp_sign == 0) true_tie_pairs += 1.0;

      // Death-first hierarchy.
      if (death_sign > 0) {
        D1_win += 1.0;
      } else if (death_sign < 0) {
        D1_loss += 1.0;
      } else {
        if (hosp_sign > 0) H2_win += 1.0;
        else if (hosp_sign < 0) H2_loss += 1.0;
      }

      // Hospitalization-first hierarchy.
      if (hosp_sign > 0) {
        H1_win += 1.0;
      } else if (hosp_sign < 0) {
        H1_loss += 1.0;
      } else {
        if (death_sign > 0) D2_win += 1.0;
        else if (death_sign < 0) D2_loss += 1.0;
      }

      // Threshold pathway: endpoint 1 is survival/death time difference.
      // If |survival difference| > t, use survival difference;
      // otherwise use hospitalization count.
      double diff_surv = fu1 - fu2;
      double abs_diff = diff_surv >= 0.0 ? diff_surv : -diff_surv;

      int k = 0;
      while (k < K && t_grid[k] < abs_diff) {
        k++;
      }

      if (k > 0) {
        if (diff_surv > 0.0) {
          add_range(win_diff_t, 0, k - 1, 1.0);
        } else if (diff_surv < 0.0) {
          add_range(loss_diff_t, 0, k - 1, 1.0);
        }
      }

      if (k < K) {
        if (hosp_sign > 0) {
          add_range(win_diff_t, k, K - 1, 1.0);
        } else if (hosp_sign < 0) {
          add_range(loss_diff_t, k, K - 1, 1.0);
        }
      }
    }
  }

  NumericVector WR_death_first(P), win_death_first(P), loss_death_first(P);
  NumericVector WR_hosp_first(P), win_hosp_first(P), loss_hosp_first(P);

  for (int a = 0; a < P; ++a) {
    double p = p_grid[a];

    win_death_first[a] = p * D1_win + (1.0 - p) * H2_win;
    loss_death_first[a] = p * D1_loss + (1.0 - p) * H2_loss;
    WR_death_first[a] = safe_wr_cpp(win_death_first[a], loss_death_first[a], total_pairs, eps);

    win_hosp_first[a] = p * H1_win + (1.0 - p) * D2_win;
    loss_hosp_first[a] = p * H1_loss + (1.0 - p) * D2_loss;
    WR_hosp_first[a] = safe_wr_cpp(win_hosp_first[a], loss_hosp_first[a], total_pairs, eps);
  }

  double ordinary_win_score = 0.5 * D1_win + 0.5 * H2_win;
  double ordinary_loss_score = 0.5 * D1_loss + 0.5 * H2_loss;
  double ordinaryWR = safe_wr_cpp(ordinary_win_score, ordinary_loss_score, total_pairs, eps);

  NumericVector t_wins(K), t_losses(K), WRt(K), pr_win_t(K), pr_loss_t(K), pr_tie_t(K);
  double cur_win = 0.0;
  double cur_loss = 0.0;
  for (int kk = 0; kk < K; ++kk) {
    cur_win += win_diff_t[kk];
    cur_loss += loss_diff_t[kk];
    t_wins[kk] = cur_win;
    t_losses[kk] = cur_loss;
    pr_win_t[kk] = cur_win / total_pairs;
    pr_loss_t[kk] = cur_loss / total_pairs;
    pr_tie_t[kk] = 1.0 - pr_win_t[kk] - pr_loss_t[kk];
    WRt[kk] = safe_wr_cpp(cur_win, cur_loss, total_pairs, eps);
  }

  return List::create(
    Named("ordinaryWR") = ordinaryWR,
    Named("ordinary_win_score") = ordinary_win_score,
    Named("ordinary_loss_score") = ordinary_loss_score,

    Named("total_pairs") = total_pairs,
    Named("D1_win") = D1_win,
    Named("D1_loss") = D1_loss,
    Named("H2_win") = H2_win,
    Named("H2_loss") = H2_loss,
    Named("H1_win") = H1_win,
    Named("H1_loss") = H1_loss,
    Named("D2_win") = D2_win,
    Named("D2_loss") = D2_loss,
    Named("death_win_pairs") = death_win_pairs,
    Named("death_loss_pairs") = death_loss_pairs,
    Named("hosp_win_pairs") = hosp_win_pairs,
    Named("hosp_loss_pairs") = hosp_loss_pairs,
    Named("true_tie_pairs") = true_tie_pairs,

    Named("p_grid") = p_grid,
    Named("WR_death_first") = WR_death_first,
    Named("win_death_first") = win_death_first,
    Named("loss_death_first") = loss_death_first,
    Named("WR_hosp_first") = WR_hosp_first,
    Named("win_hosp_first") = win_hosp_first,
    Named("loss_hosp_first") = loss_hosp_first,

    Named("t_grid") = t_grid,
    Named("WRt") = WRt,
    Named("t_wins") = t_wins,
    Named("t_losses") = t_losses,
    Named("pr_win_t") = pr_win_t,
    Named("pr_loss_t") = pr_loss_t,
    Named("pr_tie_t") = pr_tie_t
  );
}
')

max_from_curve <- function(curve, method, lower, upper) {
  ix <- which(curve$p >= lower & curve$p <= upper)
  if (length(ix) == 0) {
    return(data.frame(
      method = method,
      order = NA_character_,
      p = NA_real_,
      WR = NA_real_,
      win.score = NA_real_,
      loss.score = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  k <- ix[which.max(curve$WR[ix])]
  data.frame(
    method = method,
    order = as.character(curve$order[k]),
    p = curve$p[k],
    WR = curve$WR[k],
    win.score = curve$win.score[k],
    loss.score = curve$loss.score[k],
    stringsAsFactors = FALSE
  )
}

fast.wr.engine.revised <- function(ds,
                                   p.grid = P_GRID_ALL,
                                   t.grid = CLINICAL_T_MONTHS / 12,
                                   eps = EPS_WR) {
  if (is.null(ds$hosp.flat)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output

  core <- fast_wr_core_revised_cpp(
    futime = as.numeric(tab$FUTIME),
    cnsr = as.integer(tab$CNSR),
    arm = as.integer(tab$ARM),
    hosp_times = as.numeric(ds$hosp.flat),
    hosp_start = as.integer(ds$hosp.start),
    hosp_len = as.integer(ds$hosp.len),
    p_grid = as.numeric(p.grid),
    t_grid = as.numeric(t.grid),
    eps = eps
  )

  total.pairs <- as.numeric(core$total_pairs)
  D1.win <- as.numeric(core$D1_win)
  D1.loss <- as.numeric(core$D1_loss)
  H2.win <- as.numeric(core$H2_win)
  H2.loss <- as.numeric(core$H2_loss)
  H1.win <- as.numeric(core$H1_win)
  H1.loss <- as.numeric(core$H1_loss)
  D2.win <- as.numeric(core$D2_win)
  D2.loss <- as.numeric(core$D2_loss)
  true.tie.pairs <- as.numeric(core$true_tie_pairs)

  weighted.death.first <- data.frame(
    order = "death_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_death_first),
    win.score = as.numeric(core$win_death_first),
    loss.score = as.numeric(core$loss_death_first),
    stringsAsFactors = FALSE
  )

  weighted.death.first$win.pairs <- ifelse(
    weighted.death.first$p >= 1,
    D1.win,
    D1.win + H2.win
  )
  weighted.death.first$loss.pairs <- ifelse(
    weighted.death.first$p >= 1,
    D1.loss,
    D1.loss + H2.loss
  )
  weighted.death.first$tie.count <- total.pairs -
    weighted.death.first$win.pairs - weighted.death.first$loss.pairs
  weighted.death.first$tie.pr <- weighted.death.first$tie.count / total.pairs

  weighted.hosp.first <- data.frame(
    order = "hospitalization_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_hosp_first),
    win.score = as.numeric(core$win_hosp_first),
    loss.score = as.numeric(core$loss_hosp_first),
    stringsAsFactors = FALSE
  )

  weighted.hosp.first$win.pairs <- ifelse(
    weighted.hosp.first$p >= 1,
    H1.win,
    H1.win + D2.win
  )
  weighted.hosp.first$loss.pairs <- ifelse(
    weighted.hosp.first$p >= 1,
    H1.loss,
    H1.loss + D2.loss
  )
  weighted.hosp.first$tie.count <- total.pairs -
    weighted.hosp.first$win.pairs - weighted.hosp.first$loss.pairs
  weighted.hosp.first$tie.pr <- weighted.hosp.first$tie.count / total.pairs

  weighted.all.orders <- rbind(weighted.death.first, weighted.hosp.first)

  threshold <- data.frame(
    t = as.numeric(core$t_grid),
    t.months = as.numeric(core$t_grid) * 12,
    WR = as.numeric(core$WRt),
    wins = as.numeric(core$t_wins),
    losses = as.numeric(core$t_losses),
    pr.win = as.numeric(core$pr_win_t),
    pr.loss = as.numeric(core$pr_loss_t),
    pr.tie = as.numeric(core$pr_tie_t),
    stringsAsFactors = FALSE
  )
  threshold$tie.count <- total.pairs - threshold$wins - threshold$losses
  best.t.idx <- which.max(threshold$WR)

  p.low.values <- p.grid[p.grid < 0.5]
  p.primary.values <- p.grid[p.grid >= 0.5]
  p.low.lower <- if (length(p.low.values) > 0) min(p.low.values) else NA_real_
  p.low.upper <- if (length(p.low.values) > 0) max(p.low.values) else NA_real_
  p.primary.lower <- if (length(p.primary.values) > 0) min(p.primary.values) else NA_real_
  p.primary.upper <- if (length(p.primary.values) > 0) max(p.primary.values) else NA_real_

  list(
    ordinaryWR = as.numeric(core$ordinaryWR),
    ordinary.win.score = as.numeric(core$ordinary_win_score),
    ordinary.loss.score = as.numeric(core$ordinary_loss_score),

    weighted.death.first = weighted.death.first,
    weighted.hosp.first = weighted.hosp.first,
    weighted.all.orders = weighted.all.orders,

    max.weighted.primary = max_from_curve(weighted.death.first, "maxWRp_primary_death_first", p.primary.lower, p.primary.upper),
    max.weighted.low = max_from_curve(weighted.death.first, "maxWRp_low_death_first", p.low.lower, p.low.upper),
    max.weighted.full = max_from_curve(weighted.death.first, "maxWRp_full_death_first", min(p.grid), max(p.grid)),

    max.order.primary = max_from_curve(weighted.all.orders, "maxOrderWR_primary", p.primary.lower, p.primary.upper),
    max.order.full = max_from_curve(weighted.all.orders, "maxOrderWR_full", min(p.grid), max(p.grid)),

    threshold = threshold,
    max.threshold = threshold[best.t.idx, , drop = FALSE],

    counts = list(
      total.pairs = as.numeric(core$total_pairs),
      D1.win = as.numeric(core$D1_win),
      D1.loss = as.numeric(core$D1_loss),
      H2.win = as.numeric(core$H2_win),
      H2.loss = as.numeric(core$H2_loss),
      H1.win = as.numeric(core$H1_win),
      H1.loss = as.numeric(core$H1_loss),
      D2.win = as.numeric(core$D2_win),
      D2.loss = as.numeric(core$D2_loss),
      death.win.pairs = as.numeric(core$death_win_pairs),
      death.loss.pairs = as.numeric(core$death_loss_pairs),
      hosp.win.pairs = as.numeric(core$hosp_win_pairs),
      hosp.loss.pairs = as.numeric(core$hosp_loss_pairs),
      true.tie.pairs = as.numeric(core$true_tie_pairs),
      true.tie.pr = as.numeric(core$true_tie_pairs) / as.numeric(core$total_pairs)
    )
  )
}

wr_at_order_p <- function(engine.out, order, p) {
  if (is.na(order) || is.na(p)) return(NA_real_)
  curve <- if (order == "death_first") engine.out$weighted.death.first else engine.out$weighted.hosp.first
  idx <- extract_nearest(curve$p, p)
  curve$WR[idx]
}

threshold_at_t <- function(engine.out, t) {
  if (is.na(t)) return(NA_real_)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold$WR[idx]
}

threshold_tie_at_t <- function(engine.out, t) {
  if (is.na(t)) return(NA_real_)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold$pr.tie[idx]
}


#Composite endpoint and threshold-grid helpers
first_hosp_time_vec <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  sapply(ds$hosp.abs.times.list, function(x) if (length(x) > 0) x[1] else Inf)
}

compute.composite.endpoint.data <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  first.hosp <- first_hosp_time_vec(ds)

  death.time <- ifelse(tab$CNSR == 1, tab$FUTIME, Inf)
  comp.time <- pmin(death.time, first.hosp)
  comp.event <- is.finite(comp.time) & comp.time <= tab$FUTIME

  comp.type <- rep("censored", nrow(tab))
  comp.type[comp.event & death.time <= first.hosp] <- "death"
  comp.type[comp.event & first.hosp < death.time] <- "hospitalization"

  data.frame(
    SUBJID = tab$SUBJID,
    ARM = tab$ARM,
    time = ifelse(comp.event, comp.time, tab$FUTIME),
    event = as.integer(comp.event),
    event.type = comp.type,
    first.hosp.time = first.hosp,
    death.time = death.time,
    stringsAsFactors = FALSE
  )
}

composite.statistics <- function(ds) {
  if (is.null(ds$hosp.abs.times.list)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output
  comp <- compute.composite.endpoint.data(ds)

  dat <- data.frame(
    ARM = tab$ARM,
    FUTIME = tab$FUTIME,
    CNSR = tab$CNSR,
    NUMHOSP = tab$NUMHOSP,
    FREQHOSP = tab$FREQHOSP,
    comp.time = comp$time,
    comp.event = comp$event,
    comp.type = comp$event.type,
    stringsAsFactors = FALSE
  )

  out <- do.call(rbind, lapply(split(dat, dat$ARM), function(d) {
    data.frame(
      ARM = unique(d$ARM),
      arm.label = ifelse(unique(d$ARM) == 0, "control", "treatment"),
      n = nrow(d),
      death.events = sum(d$CNSR == 1, na.rm = TRUE),
      death.event.rate = mean(d$CNSR == 1, na.rm = TRUE),
      subjects.with.hosp = sum(d$NUMHOSP > 0, na.rm = TRUE),
      hosp.subject.rate = mean(d$NUMHOSP > 0, na.rm = TRUE),
      total.hosp.events = sum(d$NUMHOSP, na.rm = TRUE),
      mean.num.hosp = mean(d$NUMHOSP, na.rm = TRUE),
      median.num.hosp = median(d$NUMHOSP, na.rm = TRUE),
      mean.followup = mean(d$FUTIME, na.rm = TRUE),
      composite.events = sum(d$comp.event == 1, na.rm = TRUE),
      composite.event.rate = mean(d$comp.event == 1, na.rm = TRUE),
      composite.death.events = sum(d$comp.type == "death", na.rm = TRUE),
      composite.hosp.events = sum(d$comp.type == "hospitalization", na.rm = TRUE),
      mean.composite.time = mean(d$comp.time, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

run.logrank.tests <- function(ds) {
  tab <- ds$table.output

  lr.death <- survival::survdiff(
    survival::Surv(FUTIME, CNSR) ~ ARM,
    data = tab
  )
  p.death <- pchisq(lr.death$chisq, df = length(lr.death$n) - 1, lower.tail = FALSE)

  comp <- compute.composite.endpoint.data(ds)
  lr.comp <- survival::survdiff(
    survival::Surv(time, event) ~ ARM,
    data = comp
  )
  p.comp <- pchisq(lr.comp$chisq, df = length(lr.comp$n) - 1, lower.tail = FALSE)

  data.frame(
    method = c("Log-rank death endpoint", "Log-rank composite endpoint"),
    statistic = c(lr.death$chisq, lr.comp$chisq),
    p.value = c(p.death, p.comp),
    note = c(
      "Death / survival endpoint only",
      "Time to earliest of death or first hospitalization"
    ),
    stringsAsFactors = FALSE
  )
}

choose.threshold.grid.primary <- function(ds,
                                          clinical.months = CLINICAL_T_MONTHS,
                                          max.K = MAX_T_GRID_SIZE) {
  # Only first-endpoint information is used here.
  # No hospitalization time, recurrent-event count, or composite endpoint time is used.
  tab <- ds$table.output
  max.months <- max(tab$FUTIME, na.rm = TRUE) * 12

  primary.months <- tab$FUTIME * 12
  death.months <- tab$FUTIME[tab$CNSR == 1] * 12

  pair.abs.diff <- numeric(0)
  if (length(primary.months) >= 2) {
    if (length(primary.months) <= 300) {
      pair.abs.diff <- as.vector(stats::dist(primary.months, method = "manhattan"))
    } else {
      set.seed(999)
      id1 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      id2 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      pair.abs.diff <- abs(primary.months[id1] - primary.months[id2])
    }
  }

  diagnostics <- data.frame(
    quantity = c(
      "first endpoint follow-up mean",
      "first endpoint follow-up median",
      "observed death time mean among deaths",
      "observed death time median among deaths",
      "absolute pairwise first-endpoint difference mean",
      "absolute pairwise first-endpoint difference median",
      "absolute pairwise first-endpoint difference q25",
      "absolute pairwise first-endpoint difference q75"
    ),
    months = c(
      safe_mean(primary.months),
      safe_quantile(primary.months, 0.50),
      safe_mean(death.months),
      safe_quantile(death.months, 0.50),
      safe_mean(pair.abs.diff),
      safe_quantile(pair.abs.diff, 0.50),
      safe_quantile(pair.abs.diff, 0.25),
      safe_quantile(pair.abs.diff, 0.75)
    ),
    source = "first endpoint only",
    stringsAsFactors = FALSE
  )

  clinical <- clinical.months[clinical.months > 0 & clinical.months <= max.months]

  priority <- c(
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference mean"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference median"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q25"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q75"],
    0.25 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    0.50 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    diagnostics$months[diagnostics$quantity == "observed death time median among deaths"]
  )

  priority <- round(priority, digits = 1)
  priority <- priority[is.finite(priority) & !is.na(priority)]
  priority <- priority[priority > 0 & priority <= max.months]

  selected <- unique(sort(round(clinical, digits = 1)))
  for (x in priority) {
    if (length(selected) >= max.K) break
    selected <- unique(sort(c(selected, x)))
  }

  if (length(selected) == 0) {
    selected <- round(seq(max.months / 4, max.months, length.out = min(4, max.K)), 1)
    selected <- selected[selected > 0]
  }

  if (length(selected) > max.K) {
    keep.idx <- unique(round(seq(1, length(selected), length.out = max.K)))
    selected <- selected[keep.idx]
  }

  selected <- unique(sort(selected))

  threshold.table <- data.frame(
    t.months = selected,
    t.years = selected / 12,
    source = ifelse(selected %in% round(clinical, 1), "clinical candidate", "first-endpoint data-informed"),
    stringsAsFactors = FALSE
  )

  list(
    t.grid = selected / 12,
    threshold.table = threshold.table,
    diagnostics = diagnostics
  )
}

#Permutation tests: max-statistic and fixed-selected p-values
perm.test.revised <- function(ds,
                              B = B_PERM,
                              seed = MASTER_SEED,
                              p.grid = P_GRID_ALL,
                              t.grid = CLINICAL_T_MONTHS / 12,
                              verbose = FALSE) {
  ds <- prepare.ds.fast(ds)
  set.seed(seed)

  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR",
    "maxWRp_primary",
    "maxWRp_low",
    "maxWRp_full",
    "maxOrderWR_primary",
    "maxOrderWR_full",
    "maxWRt"
  )

  fixed.names <- c(
    "ordinaryWR",
    "fixedWRp_primary",
    "fixedWRp_low",
    "fixedWRp_full",
    "fixedOrder_primary",
    "fixedOrder_full",
    "fixedWRt"
  )

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )

  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names))
  colnames(perm.max) <- max.names
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names))
  colnames(perm.fixed) <- fixed.names

  # Store pointwise permutation curves so the script can report
  # average p-values and average tie counts at each candidate p and t.
  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.p.primary = NA_real_,
    selected.p.low = NA_real_,
    selected.p.full = NA_real_,
    selected.order.primary = NA_character_,
    selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_,
    selected.order.p.full = NA_real_,
    selected.t = NA_real_,
    stringsAsFactors = FALSE
  )

  perm.seeds <- seed + seq_len(B) * 1009L

  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b == B || b %% 50 == 0)) {
      cat("  Permutation", b, "of", B, "\n")
    }

    set.seed(perm.seeds[b])
    ds.b <- ds
    ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)

    out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)

    perm.pointwise.death.first[b, ] <- out.b$weighted.death.first$WR
    perm.pointwise.hosp.first[b, ] <- out.b$weighted.hosp.first$WR
    perm.pointwise.threshold[b, ] <- out.b$threshold$WR

    perm.max[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.max[b, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
    perm.max[b, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
    perm.max[b, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
    perm.max[b, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
    perm.max[b, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
    perm.max[b, "maxWRt"] <- out.b$max.threshold$WR[1]

    perm.fixed[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.fixed[b, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
    perm.fixed[b, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
    perm.fixed[b, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
    perm.fixed[b, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
    perm.fixed[b, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
    perm.fixed[b, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

    selected.perm$selected.p.primary[b] <- out.b$max.weighted.primary$p[1]
    selected.perm$selected.p.low[b] <- out.b$max.weighted.low$p[1]
    selected.perm$selected.p.full[b] <- out.b$max.weighted.full$p[1]
    selected.perm$selected.order.primary[b] <- out.b$max.order.primary$order[1]
    selected.perm$selected.order.p.primary[b] <- out.b$max.order.primary$p[1]
    selected.perm$selected.order.full[b] <- out.b$max.order.full$order[1]
    selected.perm$selected.order.p.full[b] <- out.b$max.order.full$p[1]
    selected.perm$selected.t[b] <- out.b$max.threshold$t[1]
  }

  calc_perm_p <- function(perm.vec, obs.value) {
    (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  }

  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))

  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k])
  })
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k])
  })
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) {
    calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k])
  })

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(
      death_first = perm.pointwise.death.first,
      hospitalization_first = perm.pointwise.hosp.first,
      threshold = perm.pointwise.threshold
    ),
    p.value.pointwise = list(
      death_first = p.value.pointwise.death.first,
      hospitalization_first = p.value.pointwise.hosp.first,
      threshold = p.value.pointwise.threshold
    ),
    selected.fixed = list(
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    p.grid = p.grid,
    t.grid = t.grid
  )
}


#Output tables and plotting 
fmt_p <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x < 0.001, "<0.001", sprintf("%.4f", x)))
}

save_png <- function(outdir, filename, width = 1800, height = 1100, res = 160) {
  png(file.path(outdir, filename), width = width, height = height, res = res)
}

make.pvalue.table <- function(test.out, logrank.results = NULL) {
  max.tab <- data.frame(
    method = names(test.out$T.obs.max),
    statistic.type = "max-statistic permutation p-value",
    observed.statistic = as.numeric(test.out$T.obs.max),
    p.value = as.numeric(test.out$p.value.max),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  fixed.tab <- data.frame(
    method = names(test.out$T.obs.fixed),
    statistic.type = "fixed selected-parameter permutation p-value",
    observed.statistic = as.numeric(test.out$T.obs.fixed),
    p.value = as.numeric(test.out$p.value.fixed),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  out <- rbind(max.tab, fixed.tab)

  if (!is.null(logrank.results)) {
    lr.tab <- data.frame(
      method = logrank.results$method,
      statistic.type = "log-rank chi-square p-value",
      observed.statistic = logrank.results$statistic,
      p.value = logrank.results$p.value,
      B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.tab)
  }

  out
}

make.max.vs.fixed.table <- function(test.out) {
  data.frame(
    method = c(
      "Traditional WR",
      "Weighted p: death first, p in [0.5, 1]",
      "Weighted p: death first, p < 0.5",
      "Weighted p: death first, full p grid",
      "Maximum order WR: p in [0.5, 1]",
      "Maximum order WR: full p grid",
      "Threshold WR(t)"
    ),
    selected.parameter = c(
      "p = 0.5, death-first",
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), ", death-first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), ", death-first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), ", death-first"),
      paste0("order = ", test.out$selected.fixed$order.primary,
             ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", test.out$selected.fixed$order.full,
             ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    max.statistic = as.numeric(test.out$T.obs.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    fixed.parameter.statistic = as.numeric(test.out$T.obs.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    stringsAsFactors = FALSE
  )
}


get_selected_weighted_row <- function(engine.out, order, p) {
  if (is.na(order) || is.na(p)) return(NULL)
  curve <- if (order == "death_first") engine.out$weighted.death.first else engine.out$weighted.hosp.first
  idx <- extract_nearest(curve$p, p)
  curve[idx, , drop = FALSE]
}

get_selected_threshold_row <- function(engine.out, t) {
  if (is.na(t)) return(NULL)
  idx <- extract_nearest(engine.out$threshold$t, t)
  engine.out$threshold[idx, , drop = FALSE]
}

make.selected.max.table <- function(test.out) {
  obs <- test.out$observed

  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)

  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])

  data.frame(
    method = c(
      "Traditional WR",
      "Weighted p: death first, p in [0.5, 1]",
      "Weighted p: death first, p < 0.5",
      "Weighted p: death first, full p grid",
      "Maximum order WR: p in [0.5, 1]",
      "Maximum order WR: full p grid",
      "Threshold WR(t)"
    ),
    selected.order = c(
      "death_first",
      "death_first",
      "death_first",
      "death_first",
      test.out$selected.fixed$order.primary,
      test.out$selected.fixed$order.full,
      NA_character_
    ),
    selected.p = c(
      0.5,
      test.out$selected.fixed$p.primary,
      test.out$selected.fixed$p.low,
      test.out$selected.fixed$p.full,
      test.out$selected.fixed$order.p.primary,
      test.out$selected.fixed$order.p.full,
      NA_real_
    ),
    selected.t.months = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      test.out$selected.fixed$t.months
    ),
    observed.max.statistic = as.numeric(test.out$T.obs.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[c(
      "ordinaryWR",
      "maxWRp_primary",
      "maxWRp_low",
      "maxWRp_full",
      "maxOrderWR_primary",
      "maxOrderWR_full",
      "maxWRt"
    )]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR",
      "fixedWRp_primary",
      "fixedWRp_low",
      "fixedWRp_full",
      "fixedOrder_primary",
      "fixedOrder_full",
      "fixedWRt"
    )]),
    selected.tie.count = c(
      get_tie_count(w.traditional),
      get_tie_count(w.primary),
      get_tie_count(w.low),
      get_tie_count(w.full),
      get_tie_count(o.primary),
      get_tie_count(o.full),
      get_tie_count(t.row)
    ),
    selected.tie.pr = c(
      get_tie_pr(w.traditional),
      get_tie_pr(w.primary),
      get_tie_pr(w.low),
      get_tie_pr(w.full),
      get_tie_pr(o.primary),
      get_tie_pr(o.full),
      get_tie_pr(t.row)
    ),
    stringsAsFactors = FALSE
  )
}

make.pointwise.table <- function(test.out) {
  obs <- test.out$observed

  death <- obs$weighted.death.first
  death$family <- "weighted_p"
  death$parameter.name <- "p"
  death$parameter <- death$p
  death$parameter.months <- NA_real_
  death$parameter.label <- paste0("p = ", sprintf("%.2f", death$p))
  death$observed.WR <- death$WR
  death$pointwise.p.value <- as.numeric(test.out$p.value.pointwise$death_first)

  hosp <- obs$weighted.hosp.first
  hosp$family <- "weighted_p"
  hosp$parameter.name <- "p"
  hosp$parameter <- hosp$p
  hosp$parameter.months <- NA_real_
  hosp$parameter.label <- paste0("p = ", sprintf("%.2f", hosp$p))
  hosp$observed.WR <- hosp$WR
  hosp$pointwise.p.value <- as.numeric(test.out$p.value.pointwise$hospitalization_first)

  thr <- obs$threshold
  thr.out <- data.frame(
    order = NA_character_,
    p = NA_real_,
    WR = thr$WR,
    win.score = thr$wins,
    loss.score = thr$losses,
    win.pairs = thr$wins,
    loss.pairs = thr$losses,
    tie.count = thr$tie.count,
    tie.pr = thr$pr.tie,
    family = "threshold_t",
    parameter.name = "t",
    parameter = thr$t,
    parameter.months = thr$t.months,
    parameter.label = paste0("t = ", sprintf("%.1f", thr$t.months), " months"),
    observed.WR = thr$WR,
    pointwise.p.value = as.numeric(test.out$p.value.pointwise$threshold),
    stringsAsFactors = FALSE
  )

  keep.cols <- c(
    "family", "order", "parameter.name", "parameter", "parameter.months",
    "parameter.label", "observed.WR", "pointwise.p.value", "win.score",
    "loss.score", "win.pairs", "loss.pairs", "tie.count", "tie.pr"
  )

  out <- rbind(
    death[, keep.cols, drop = FALSE],
    hosp[, keep.cols, drop = FALSE],
    thr.out[, keep.cols, drop = FALSE]
  )
  out$total.pairs <- test.out$observed$counts$total.pairs
  out
}

plot.weighted.death.first <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_weighted_death_first_WRp.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  tab <- engine.out$weighted.death.first
  plot(tab$p, tab$WR, type = "l", lwd = 2,
       xlab = "p: death/survival weight",
       ylab = "WR(p)",
       main = "Death-first weighted WR(p)")
  abline(v = 0.5, lty = 3, lwd = 2)
  points(engine.out$max.weighted.primary$p, engine.out$max.weighted.primary$WR, pch = 19)
  points(engine.out$max.weighted.low$p, engine.out$max.weighted.low$WR, pch = 17)
  legend("topright",
         legend = c(
           paste0("primary max p = ", sprintf("%.2f", engine.out$max.weighted.primary$p),
                  ", WR = ", sprintf("%.3f", engine.out$max.weighted.primary$WR)),
           paste0("low-p max p = ", sprintf("%.2f", engine.out$max.weighted.low$p),
                  ", WR = ", sprintf("%.3f", engine.out$max.weighted.low$WR))
         ),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.order.curves <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_maximum_order_WRp.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  dtab <- engine.out$weighted.death.first
  htab <- engine.out$weighted.hosp.first
  yr <- range(c(dtab$WR, htab$WR), finite = TRUE)
  plot(dtab$p, dtab$WR, type = "l", lwd = 2, ylim = yr,
       xlab = "p: first-endpoint weight",
       ylab = "WR(p)",
       main = "Effect of hierarchy order on maximum WR")
  lines(htab$p, htab$WR, lwd = 2, lty = 2)
  points(engine.out$max.order.primary$p, engine.out$max.order.primary$WR, pch = 19)
  abline(v = 0.5, lty = 3, lwd = 2)
  legend("topright",
         legend = c(
           "death first",
           "hospitalization first",
           paste0("selected order = ", engine.out$max.order.primary$order),
           paste0("selected p = ", sprintf("%.2f", engine.out$max.order.primary$p)),
           paste0("max WR = ", sprintf("%.3f", engine.out$max.order.primary$WR))
         ),
         lty = c(1, 2, NA, NA, NA),
         pch = c(NA, NA, NA, NA, NA),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.threshold.curve <- function(engine.out, outdir, prefix, p.value = NA_real_) {
  save_png(outdir, paste0(prefix, "_threshold_WRt.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  tab <- engine.out$threshold
  plot(tab$t.months, tab$WR, type = "b", pch = 19, lwd = 2,
       xlab = "threshold t (months)",
       ylab = "WR(t)",
       main = "Threshold WR(t)")
  abline(v = engine.out$max.threshold$t.months, lty = 2, lwd = 2)
  legend("topright",
         legend = c(
           paste0("selected t = ", sprintf("%.1f", engine.out$max.threshold$t.months), " months"),
           paste0("max WR = ", sprintf("%.3f", engine.out$max.threshold$WR)),
           paste0("max-stat p = ", fmt_p(p.value))
         ),
         bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}

plot.max.vs.fixed.pvalues <- function(max.vs.fixed.table, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_max_vs_fixed_pvalues.png"), width = 2100, height = 1150)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(7.5, 5.2, 3.2, 1.2))
  mat <- rbind(
    "Max-statistic" = max.vs.fixed.table$max.statistic.p.value,
    "Fixed-selected" = max.vs.fixed.table$fixed.parameter.p.value
  )
  bp <- barplot(mat, beside = TRUE,
                names.arg = c("Trad WR", "p>=0.5", "p<0.5", "full p", "order p>=0.5", "order full", "threshold"),
                ylim = c(0, 1),
                ylab = "Permutation p-value",
                main = "Max-statistic vs fixed-selected p-values",
                legend.text = TRUE,
                args.legend = list(x = "topright", bty = "n"),
                las = 2, cex.names = 0.75)
  abline(h = ALPHA, lty = 2, lwd = 2)
  text(x = bp, y = pmin(mat + 0.035, 0.98), labels = sprintf("%.3f", mat), cex = 0.65)
  par(old.par)
  dev.off()
}

plot.km.with.threshold <- function(ds, engine.out, outdir, prefix) {
  tab <- ds$table.output
  comp <- compute.composite.endpoint.data(ds)
  selected.t <- engine.out$max.threshold$t[1]

  save_png(outdir, paste0(prefix, "_KM_death_with_selected_t.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  fit.death <- survival::survfit(survival::Surv(FUTIME, CNSR) ~ ARM, data = tab)
  plot(fit.death, lwd = 2, mark.time = TRUE,
       xlab = "Time (years)", ylab = "Survival probability",
       main = "Kaplan-Meier curve: death endpoint")
  abline(v = selected.t, lty = 2, lwd = 2)
  legend("bottomleft",
         legend = c("ARM 0: control", "ARM 1: treatment",
                    paste0("selected t = ", sprintf("%.1f", selected.t * 12), " months")),
         lty = c(1, 1, 2), lwd = c(2, 2, 2), bty = "n", cex = 0.85)
  par(old.par)
  dev.off()

  save_png(outdir, paste0(prefix, "_KM_composite_with_selected_t.png"))
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.4, 1.2))
  fit.comp <- survival::survfit(survival::Surv(time, event) ~ ARM, data = comp)
  plot(fit.comp, lwd = 2, mark.time = TRUE,
       xlab = "Time (years)", ylab = "Composite-event-free probability",
       main = "Kaplan-Meier curve: composite endpoint")
  abline(v = selected.t, lty = 2, lwd = 2)
  legend("bottomleft",
         legend = c("ARM 0: control", "ARM 1: treatment",
                    paste0("selected t = ", sprintf("%.1f", selected.t * 12), " months")),
         lty = c(1, 1, 2), lwd = c(2, 2, 2), bty = "n", cex = 0.85)
  par(old.par)
  dev.off()
}



# final reporting labels, scenario comparisons,
scenario.display.label <- function(x) {
  map <- c(
    S01_equivalence_equal_arms = "S01 Equal",
    S02_similarity_small_benefit_both = "S02 Similar +",
    S03_similarity_small_harm_both = "S03 Similar -",
    S04_superiority_moderate_better_both = "S04 Sup mod",
    S05_superiority_strong_better_both = "S05 Sup strong",
    S06_inferiority_moderate_worse_both = "S06 Inf mod",
    S07_inferiority_strong_worse_both = "S07 Inf strong",
    S08_death_benefit_only = "S08 Death +",
    S09_hosp_benefit_only = "S09 Hosp +",
    S10_death_benefit_hosp_harm = "S10 Death + / Hosp -",
    S11_death_harm_hosp_benefit = "S11 Death - / Hosp +",
    S12_high_random_censoring_superiority = "S12 High censor",
    S13_low_censoring_long_followup_superiority = "S13 Long FU"
  )
  out <- as.character(x)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

method.long.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional WR",
    maxWRp_primary = "Weighted max WR, p >= 0.5",
    maxWRp_low = "Weighted max WR, p < 0.5",
    maxWRp_full = "Weighted max WR, full p grid",
    maxOrderWR_primary = "Maximum-order WR, p >= 0.5",
    maxOrderWR_full = "Maximum-order WR, full p grid",
    maxWRt = "Threshold max WR(t)",
    logrank_death = "Log-rank death",
    logrank_composite = "Log-rank composite",
    fixedWRp_primary = "Weighted fixed-selected WR, p >= 0.5",
    fixedWRp_low = "Weighted fixed-selected WR, p < 0.5",
    fixedWRp_full = "Weighted fixed-selected WR, full p grid",
    fixedOrderWR_primary = "Maximum-order fixed-selected WR, p >= 0.5",
    fixedOrderWR_full = "Maximum-order fixed-selected WR, full p grid",
    fixedWRt = "Threshold fixed-selected WR(t)"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

method.display.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional\nWR",
    maxWRp_primary = "Weighted\nmax",
    maxWRp_low = "Low-p\nmax",
    maxWRp_full = "Full-p\nmax",
    maxOrderWR_primary = "Order\nmax",
    maxOrderWR_full = "Order full\nmax",
    maxWRt = "Threshold\nmax",
    logrank_death = "LR\ndeath",
    logrank_composite = "LR\ncomposite",
    fixedWRp_primary = "Weighted\nfixed",
    fixedWRp_low = "Low-p\nfixed",
    fixedWRp_full = "Full-p\nfixed",
    fixedOrderWR_primary = "Order\nfixed",
    fixedOrderWR_full = "Order full\nfixed",
    fixedWRt = "Threshold\nfixed"
  )
  out <- as.character(method)
  hit <- out %in% names(map)
  out[hit] <- unname(map[out[hit]])
  out[!hit] <- gsub("_", " ", out[!hit])
  out
}

clean.order.label <- function(x) {
  out <- as.character(x)
  out <- gsub("death_first", "death first", out)
  out <- gsub("hospitalization_first", "hospitalization first", out)
  out <- gsub("_", " ", out)
  out
}

safe.mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe.median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe.col <- function(df, col) {
  if (is.na(col) || !(col %in% names(df))) return(rep(NA_real_, nrow(df)))
  as.numeric(df[[col]])
}

mode.string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

get.power.value <- function(power.table, method) {
  if (is.na(method) || is.null(power.table) || nrow(power.table) == 0) return(NA_real_)
  z <- power.table$rejection.proportion[power.table$method == method]
  if (length(z) == 0) return(NA_real_)
  as.numeric(z[1])
}


#result tables and figures
make.pathway.results.table <- function(test.out, logrank.results = NULL) {
  obs <- test.out$observed

  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)

  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])

  method.id <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )

  out <- data.frame(
    method = method.id,
    method_label = method.long.label(method.id),
    statistic_type = "Win ratio",
    statistic_name = c(
      "WR(p = 0.50), death first",
      "max WR(p), death first, p in [0.50, 1.00]",
      "max WR(p), death first, p in [0.01, 0.49]",
      "max WR(p), death first, p in [0.01, 1.00]",
      "max WR(p), selected order, p in [0.50, 1.00]",
      "max WR(p), selected order, p in [0.01, 1.00]",
      "max WR(t)"
    ),
    observed_statistic = as.numeric(test.out$T.obs.max[method.id]),
    selected_order = c(
      "death_first", "death_first", "death_first", "death_first",
      test.out$selected.fixed$order.primary,
      test.out$selected.fixed$order.full,
      NA_character_
    ),
    selected_p = c(
      0.5,
      test.out$selected.fixed$p.primary,
      test.out$selected.fixed$p.low,
      test.out$selected.fixed$p.full,
      test.out$selected.fixed$order.p.primary,
      test.out$selected.fixed$order.p.full,
      NA_real_
    ),
    selected_t_months = c(
      NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
      test.out$selected.fixed$t.months
    ),
    selected_value = c(
      "p = 0.50; death first",
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), "; death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary),
             "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full),
             "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    permutation_p_value = as.numeric(test.out$p.value.max[method.id]),
    selected_parameter_p_value = as.numeric(test.out$p.value.fixed[c(
      "ordinaryWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
      "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
    )]),
    tie_count = c(
      get_tie_count(w.traditional), get_tie_count(w.primary), get_tie_count(w.low),
      get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full),
      get_tie_count(t.row)
    ),
    tie_proportion = c(
      get_tie_pr(w.traditional), get_tie_pr(w.primary), get_tie_pr(w.low),
      get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full),
      get_tie_pr(t.row)
    ),
    B = test.out$B,
    stringsAsFactors = FALSE
  )

  if (!is.null(logrank.results)) {
    lr.death <- logrank.results[logrank.results$method == "Log-rank death endpoint", , drop = FALSE]
    lr.comp <- logrank.results[logrank.results$method == "Log-rank composite endpoint", , drop = FALSE]
    lr.rows <- data.frame(
      method = c("logrank_death", "logrank_composite"),
      method_label = method.long.label(c("logrank_death", "logrank_composite")),
      statistic_type = "Log-rank chi-square",
      statistic_name = c("Death endpoint", "Composite endpoint"),
      observed_statistic = c(lr.death$statistic[1], lr.comp$statistic[1]),
      selected_order = NA_character_,
      selected_p = NA_real_,
      selected_t_months = NA_real_,
      selected_value = NA_character_,
      permutation_p_value = c(lr.death$p.value[1], lr.comp$p.value[1]),
      selected_parameter_p_value = NA_real_,
      tie_count = NA_real_,
      tie_proportion = NA_real_,
      B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.rows)
  }

  out
}

plot.pvalue.bar <- function(tab, outdir, filename, title, p.col = "permutation_p_value",
                            method.order = NULL, show.alpha = TRUE) {
  if (is.null(tab) || nrow(tab) == 0 || !(p.col %in% names(tab))) return(invisible(NULL))
  d <- tab
  if (!is.null(method.order)) d <- d[d$method %in% method.order, , drop = FALSE]
  d <- d[!is.na(d[[p.col]]), , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  if (!is.null(method.order)) {
    d$method <- factor(d$method, levels = method.order)
    d <- d[order(d$method), , drop = FALSE]
  }
  d <- d[order(d[[p.col]], decreasing = TRUE), , drop = FALSE]

  save_png(outdir, filename, width = 2100, height = 1150)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 13.0, 3.2, 1.2))
  bp <- barplot(d[[p.col]], horiz = TRUE,
                names.arg = method.long.label(as.character(d$method)),
                xlim = c(0, 1), las = 1,
                xlab = "P-value",
                main = title,
                cex.names = 0.86)
  if (show.alpha) abline(v = ALPHA, lty = 2, lwd = 2)
  text(x = pmin(d[[p.col]] + 0.035, 0.96), y = bp,
       labels = fmt_p(d[[p.col]]), cex = 0.82)
  if (show.alpha) legend("bottomright", legend = paste0("alpha = ", ALPHA), bty = "n")
  par(old.par)
  dev.off()
}

plot.example.pvalue.figures <- function(pathway.tab, outdir, prefix) {
  all.methods <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt",
    "logrank_death", "logrank_composite"
  )
  wr.methods <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  weighted.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full")
  logrank.methods <- c("logrank_death", "logrank_composite")
  threshold.methods <- c("maxWRt")

  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_all_methods_with_logrank.png"),
                  "All method p-values", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_WR_pathways_only_no_logrank.png"),
                  "WR pathway p-values", "permutation_p_value", wr.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_traditional_and_weighted_only.png"),
                  "Traditional and weighted WR p-values", "permutation_p_value", weighted.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_logrank_pvalues_separate.png"),
                  "Log-rank p-values", "permutation_p_value", logrank.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_threshold_pvalue_separate.png"),
                  "Threshold pathway p-value", "permutation_p_value", threshold.methods)
}

plot.permutation.null.revised <- function(test.out, method, outdir, prefix) {
  if (!(method %in% colnames(test.out$T.perm.max))) return(invisible(NULL))
  x <- as.numeric(test.out$T.perm.max[, method])
  obs.value <- as.numeric(test.out$T.obs.max[method])
  p.value <- as.numeric(test.out$p.value.max[method])

  save_png(outdir, paste0(prefix, "_perm_null_", method, ".png"), width = 1600, height = 1000)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(5.2, 5.2, 3.2, 1.2))
  hist(x, breaks = 35,
       main = paste0("Permutation null: ", method.long.label(method)),
       xlab = "Permutation statistic",
       ylab = "Frequency")
  abline(v = obs.value, lty = 2, lwd = 2)
  legend("topright",
         legend = c(
           paste0("observed = ", sprintf("%.3f", obs.value)),
           paste0("p = ", fmt_p(p.value)),
           paste0("B = ", test.out$B)
         ), bty = "n")
  par(old.par)
  dev.off()
}

plot.example.permutation.nulls <- function(test.out, outdir, prefix) {
  for (m in c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")) {
    plot.permutation.null.revised(test.out, m, outdir, prefix)
  }
}

# Override one-example output to include the DIG-style plots plus the new pathway table.
save.one.dataset.outputs <- function(ds, test.out, threshold.info, logrank.results,
                                     comp.stats, outdir, prefix = "example") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  obs <- test.out$observed

  write.csv(obs$weighted.death.first, file.path(outdir, paste0(prefix, "_weighted_death_first_curve.csv")), row.names = FALSE)
  write.csv(obs$weighted.hosp.first, file.path(outdir, paste0(prefix, "_weighted_hospitalization_first_curve.csv")), row.names = FALSE)
  write.csv(obs$threshold, file.path(outdir, paste0(prefix, "_threshold_curve.csv")), row.names = FALSE)
  write.csv(threshold.info$diagnostics, file.path(outdir, paste0(prefix, "_threshold_diagnostics_first_endpoint_only.csv")), row.names = FALSE)
  write.csv(threshold.info$threshold.table, file.path(outdir, paste0(prefix, "_threshold_grid.csv")), row.names = FALSE)
  write.csv(logrank.results, file.path(outdir, paste0(prefix, "_logrank_results.csv")), row.names = FALSE)
  write.csv(comp.stats, file.path(outdir, paste0(prefix, "_composite_statistics.csv")), row.names = FALSE)
  write.csv(make.pvalue.table(test.out, logrank.results), file.path(outdir, paste0(prefix, "_pvalue_table.csv")), row.names = FALSE)
  write.csv(make.selected.max.table(test.out), file.path(outdir, paste0(prefix, "_selected_max_parameter_table.csv")), row.names = FALSE)
  write.csv(make.pointwise.table(test.out), file.path(outdir, paste0(prefix, "_pointwise_pvalue_and_tie_table.csv")), row.names = FALSE)

  pathway.tab <- make.pathway.results.table(test.out, logrank.results)
  write.csv(pathway.tab, file.path(outdir, paste0(prefix, "_pathway_results_table.csv")), row.names = FALSE)

  max.vs.fixed.table <- make.max.vs.fixed.table(test.out)
  write.csv(max.vs.fixed.table, file.path(outdir, paste0(prefix, "_max_vs_fixed_parameter_pvalues.csv")), row.names = FALSE)

  plot.weighted.death.first(obs, outdir, prefix)
  plot.order.curves(obs, outdir, prefix)
  plot.threshold.curve(obs, outdir, prefix, p.value = test.out$p.value.max["maxWRt"])
  plot.max.vs.fixed.pvalues(max.vs.fixed.table, outdir, prefix)
  plot.km.with.threshold(ds, obs, outdir, prefix)
  plot.example.pvalue.figures(pathway.tab, outdir, prefix)
  plot.example.permutation.nulls(test.out, outdir, prefix)
}



# Base-R row bind that tolerates tables with different columns.
rbind_fill_base <- function(a, b) {
  if (is.null(a) || nrow(a) == 0) return(b)
  if (is.null(b) || nrow(b) == 0) return(a)
  cols <- union(names(a), names(b))
  for (cc in setdiff(cols, names(a))) a[[cc]] <- NA
  for (cc in setdiff(cols, names(b))) b[[cc]] <- NA
  rbind(a[, cols, drop = FALSE], b[, cols, drop = FALSE])
}


#Threshold source table and permutation wrapper
threshold.candidate.source.table <- function(clinical.months = CLINICAL_T_MONTHS) {
  data.frame(
    candidate_months = paste(clinical.months, collapse = ", "),
    candidate_type = "clinical sensitivity candidates",
    how_used = paste(
      "These time-difference candidates are pre-specified sensitivity candidates.",
      "They are combined with first-endpoint data-informed candidates.",
      "Final clinical margins should be justified for the disease area and SAP."
    ),
    source_key = c(
      "Pocock 2012 win-ratio hierarchy; FDA multiple-endpoints guidance; ICH E9(R1) estimand framework; RMST literature for time-scale interpretation"
    ),
    stringsAsFactors = FALSE
  )
}

# Override the threshold-grid helper to add source fields while keeping the
# first-endpoint-only rule. No hospitalization or composite endpoint data are
# used to choose data-informed t candidates.
choose.threshold.grid.primary <- function(ds,
                                          clinical.months = CLINICAL_T_MONTHS,
                                          max.K = MAX_T_GRID_SIZE) {
  tab <- ds$table.output
  max.months <- max(tab$FUTIME, na.rm = TRUE) * 12

  primary.months <- tab$FUTIME * 12
  death.months <- tab$FUTIME[tab$CNSR == 1] * 12

  pair.abs.diff <- numeric(0)
  if (length(primary.months) >= 2) {
    if (length(primary.months) <= 300) {
      pair.abs.diff <- as.vector(stats::dist(primary.months, method = "manhattan"))
    } else {
      set.seed(999)
      id1 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      id2 <- sample(seq_along(primary.months), 5000, replace = TRUE)
      pair.abs.diff <- abs(primary.months[id1] - primary.months[id2])
    }
  }

  diagnostics <- data.frame(
    quantity = c(
      "first endpoint follow-up mean",
      "first endpoint follow-up median",
      "observed death time mean among deaths",
      "observed death time median among deaths",
      "absolute pairwise first-endpoint difference mean",
      "absolute pairwise first-endpoint difference median",
      "absolute pairwise first-endpoint difference q25",
      "absolute pairwise first-endpoint difference q75"
    ),
    months = c(
      safe_mean(primary.months),
      safe_quantile(primary.months, 0.50),
      safe_mean(death.months),
      safe_quantile(death.months, 0.50),
      safe_mean(pair.abs.diff),
      safe_quantile(pair.abs.diff, 0.50),
      safe_quantile(pair.abs.diff, 0.25),
      safe_quantile(pair.abs.diff, 0.75)
    ),
    source = "first endpoint only",
    stringsAsFactors = FALSE
  )

  clinical <- clinical.months[clinical.months > 0 & clinical.months <= max.months]

  priority <- c(
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference mean"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference median"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q25"],
    diagnostics$months[diagnostics$quantity == "absolute pairwise first-endpoint difference q75"],
    0.25 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    0.50 * diagnostics$months[diagnostics$quantity == "first endpoint follow-up mean"],
    diagnostics$months[diagnostics$quantity == "observed death time median among deaths"]
  )

  priority <- round(priority, digits = 1)
  priority <- priority[is.finite(priority) & !is.na(priority)]
  priority <- priority[priority > 0 & priority <= max.months]

  selected <- unique(sort(round(clinical, digits = 1)))
  for (x in priority) {
    if (length(selected) >= max.K) break
    selected <- unique(sort(c(selected, x)))
  }

  if (length(selected) == 0) {
    selected <- round(seq(max.months / 4, max.months, length.out = min(4, max.K)), 1)
    selected <- selected[selected > 0]
  }

  if (length(selected) > max.K) {
    keep.idx <- unique(round(seq(1, length(selected), length.out = max.K)))
    selected <- selected[keep.idx]
  }

  selected <- unique(sort(selected))
  clinical.round <- round(clinical, 1)

  threshold.table <- data.frame(
    t.months = selected,
    t.years = selected / 12,
    source = ifelse(selected %in% clinical.round, "clinical candidate", "first-endpoint data-informed"),
    source_detail = ifelse(
      selected %in% clinical.round,
      "Pre-specified clinical sensitivity candidate; justify or replace using disease-specific clinical input.",
      "Selected from first-endpoint follow-up/death-time/pairwise-difference summaries only."
    ),
    reference_key = ifelse(
      selected %in% clinical.round,
      "Pocock2012/FDA-MultipleEndpoints/ICH-E9R1/RMST-time-scale",
      "First-endpoint data only; no recurrent hospitalization or composite endpoint information used"
    ),
    stringsAsFactors = FALSE
  )

  list(
    t.grid = selected / 12,
    threshold.table = threshold.table,
    diagnostics = diagnostics,
    source.table = threshold.candidate.source.table(clinical.months)
  )
}

perm.test.revised.batch <- function(ds,
                                    B = B_REAL,
                                    batch.size = BATCH_SIZE_REAL,
                                    seed = MASTER_SEED,
                                    p.grid = P_GRID_ALL,
                                    t.grid = CLINICAL_T_MONTHS / 12,
                                    outdir = OUTDIR,
                                    cache.prefix = "realdata",
                                    resume = RESUME_IF_EXISTS,
                                    verbose = VERBOSE) {
  ds <- prepare.ds.fast(ds)
  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.names <- c(
    "ordinaryWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )

  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names), dimnames = list(NULL, max.names))
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names), dimnames = list(NULL, fixed.names))
  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.p.primary = NA_real_,
    selected.p.low = NA_real_,
    selected.p.full = NA_real_,
    selected.order.primary = NA_character_,
    selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_,
    selected.order.p.full = NA_real_,
    selected.t = NA_real_,
    stringsAsFactors = FALSE
  )

  batch.dir <- file.path(outdir, paste0(cache.prefix, "_permutation_batches"))
  if (!dir.exists(batch.dir)) dir.create(batch.dir, recursive = TRUE)
  n.batch <- ceiling(B / batch.size)
  perm.seeds <- seed + seq_len(B) * 1009L

  for (bb in seq_len(n.batch)) {
    start.b <- (bb - 1L) * batch.size + 1L
    end.b <- min(bb * batch.size, B)
    idx <- start.b:end.b
    batch.file <- file.path(batch.dir, sprintf("perm_batch_%04d_%04d_to_%04d.rds", bb, start.b, end.b))

    if (resume && file.exists(batch.file)) {
      if (verbose) cat("Loading existing batch", bb, "of", n.batch, "\n")
      batch.out <- readRDS(batch.file)
    } else {
      if (verbose) cat("Running batch", bb, "of", n.batch, ": permutations", start.b, "to", end.b, "\n")

      n.idx <- length(idx)
      batch.max <- matrix(NA_real_, nrow = n.idx, ncol = length(max.names), dimnames = list(NULL, max.names))
      batch.fixed <- matrix(NA_real_, nrow = n.idx, ncol = length(fixed.names), dimnames = list(NULL, fixed.names))
      batch.point.death <- matrix(NA_real_, nrow = n.idx, ncol = length(p.grid))
      batch.point.hosp <- matrix(NA_real_, nrow = n.idx, ncol = length(p.grid))
      batch.point.thr <- matrix(NA_real_, nrow = n.idx, ncol = length(t.grid))
      batch.selected <- selected.perm[idx, , drop = FALSE]

      for (r in seq_along(idx)) {
        b <- idx[r]
        set.seed(perm.seeds[b])
        ds.b <- ds
        ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)
        out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)

        batch.point.death[r, ] <- out.b$weighted.death.first$WR
        batch.point.hosp[r, ] <- out.b$weighted.hosp.first$WR
        batch.point.thr[r, ] <- out.b$threshold$WR

        batch.max[r, "ordinaryWR"] <- out.b$ordinaryWR
        batch.max[r, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
        batch.max[r, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
        batch.max[r, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
        batch.max[r, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
        batch.max[r, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
        batch.max[r, "maxWRt"] <- out.b$max.threshold$WR[1]

        batch.fixed[r, "ordinaryWR"] <- out.b$ordinaryWR
        batch.fixed[r, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
        batch.fixed[r, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
        batch.fixed[r, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
        batch.fixed[r, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
        batch.fixed[r, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
        batch.fixed[r, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

        batch.selected$selected.p.primary[r] <- out.b$max.weighted.primary$p[1]
        batch.selected$selected.p.low[r] <- out.b$max.weighted.low$p[1]
        batch.selected$selected.p.full[r] <- out.b$max.weighted.full$p[1]
        batch.selected$selected.order.primary[r] <- out.b$max.order.primary$order[1]
        batch.selected$selected.order.p.primary[r] <- out.b$max.order.primary$p[1]
        batch.selected$selected.order.full[r] <- out.b$max.order.full$order[1]
        batch.selected$selected.order.p.full[r] <- out.b$max.order.full$p[1]
        batch.selected$selected.t[r] <- out.b$max.threshold$t[1]
      }

      batch.out <- list(
        idx = idx,
        perm.max = batch.max,
        perm.fixed = batch.fixed,
        point.death = batch.point.death,
        point.hosp = batch.point.hosp,
        point.threshold = batch.point.thr,
        selected = batch.selected
      )
      saveRDS(batch.out, batch.file)
    }

    perm.max[idx, ] <- batch.out$perm.max
    perm.fixed[idx, ] <- batch.out$perm.fixed
    perm.pointwise.death.first[idx, ] <- batch.out$point.death
    perm.pointwise.hosp.first[idx, ] <- batch.out$point.hosp
    perm.pointwise.threshold[idx, ] <- batch.out$point.threshold
    selected.perm[idx, ] <- batch.out$selected
  }

  calc_perm_p <- function(perm.vec, obs.value) {
    (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  }

  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))

  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k])
  })
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) {
    calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k])
  })
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) {
    calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k])
  })

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(
      death_first = perm.pointwise.death.first,
      hospitalization_first = perm.pointwise.hosp.first,
      threshold = perm.pointwise.threshold
    ),
    p.value.pointwise = list(
      death_first = p.value.pointwise.death.first,
      hospitalization_first = p.value.pointwise.hosp.first,
      threshold = p.value.pointwise.threshold
    ),
    selected.fixed = list(
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    batch.size = batch.size,
    p.grid = p.grid,
    t.grid = t.grid,
    batch.dir = batch.dir
  )
}

#data summary tables and figures
dataset.display.label <- function(id, label = NULL) {
  if (!is.null(label) && length(label) > 0 && !is.na(label) && label != "") return(as.character(label))
  gsub("_", " ", as.character(id))
}

get_comp_arm_value <- function(comp.stats, arm, varname) {
  z <- comp.stats[comp.stats$ARM == arm, varname]
  if (length(z) == 0) return(NA_real_)
  as.numeric(z[1])
}

make.realdata.summary.row <- function(dataset.id, dataset.label, test.out, logrank.results,
                                      comp.stats, threshold.info, import.note) {
  obs <- test.out$observed
  death.first.primary <- obs$max.weighted.primary
  hosp.first.primary <- max_from_curve(obs$weighted.hosp.first, "hosp_first_primary", 0.50, 1.00)
  order.diff <- hosp.first.primary$WR[1] - death.first.primary$WR[1]

  w.traditional.row <- get_selected_weighted_row(obs, "death_first", 0.5)
  w.primary.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.primary$p[1])
  w.low.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.low$p[1])
  w.full.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.full$p[1])
  o.primary.row <- get_selected_weighted_row(obs, obs$max.order.primary$order[1], obs$max.order.primary$p[1])
  o.full.row <- get_selected_weighted_row(obs, obs$max.order.full$order[1], obs$max.order.full$p[1])

  data.frame(
    dataset_id = dataset.id,
    dataset_label = dataset.label,
    n_subjects = import.note$n_subjects[1],
    n_control = import.note$n_control[1],
    n_treatment = import.note$n_treatment[1],

    ordinaryWR = test.out$T.obs.max["ordinaryWR"],
    traditional.tie.count = w.traditional.row$tie.count[1],
    traditional.tie.pr = w.traditional.row$tie.pr[1],
    max.pvalue.ordinaryWR = test.out$p.value.max["ordinaryWR"],
    fixed.pvalue.ordinaryWR = test.out$p.value.fixed["ordinaryWR"],

    maxWRp_primary = test.out$T.obs.max["maxWRp_primary"],
    selected.p.primary = obs$max.weighted.primary$p[1],
    weighted.primary.tie.count = w.primary.row$tie.count[1],
    weighted.primary.tie.pr = w.primary.row$tie.pr[1],
    max.pvalue.maxWRp_primary = test.out$p.value.max["maxWRp_primary"],
    fixed.pvalue.maxWRp_primary = test.out$p.value.fixed["fixedWRp_primary"],

    maxWRp_low = test.out$T.obs.max["maxWRp_low"],
    selected.p.low = obs$max.weighted.low$p[1],
    weighted.low.tie.count = w.low.row$tie.count[1],
    weighted.low.tie.pr = w.low.row$tie.pr[1],
    max.pvalue.maxWRp_low = test.out$p.value.max["maxWRp_low"],
    fixed.pvalue.maxWRp_low = test.out$p.value.fixed["fixedWRp_low"],

    maxWRp_full = test.out$T.obs.max["maxWRp_full"],
    selected.p.full = obs$max.weighted.full$p[1],
    weighted.full.tie.count = w.full.row$tie.count[1],
    weighted.full.tie.pr = w.full.row$tie.pr[1],
    selected.full.is.low.p = as.integer(obs$max.weighted.full$p[1] < 0.5),
    max.pvalue.maxWRp_full = test.out$p.value.max["maxWRp_full"],
    fixed.pvalue.maxWRp_full = test.out$p.value.fixed["fixedWRp_full"],

    death.first.primary.WR = death.first.primary$WR[1],
    hosp.first.primary.WR = hosp.first.primary$WR[1],
    order.diff.hosp.minus.death = order.diff,
    order.changed.primary = as.integer(obs$max.order.primary$order[1] == "hospitalization_first"),
    maxOrderWR_primary = test.out$T.obs.max["maxOrderWR_primary"],
    selected.order.primary = obs$max.order.primary$order[1],
    selected.order.p.primary = obs$max.order.primary$p[1],
    order.primary.tie.count = o.primary.row$tie.count[1],
    order.primary.tie.pr = o.primary.row$tie.pr[1],
    max.pvalue.maxOrderWR_primary = test.out$p.value.max["maxOrderWR_primary"],
    fixed.pvalue.maxOrderWR_primary = test.out$p.value.fixed["fixedOrder_primary"],

    maxOrderWR_full = test.out$T.obs.max["maxOrderWR_full"],
    selected.order.full = obs$max.order.full$order[1],
    selected.order.p.full = obs$max.order.full$p[1],
    order.full.tie.count = o.full.row$tie.count[1],
    order.full.tie.pr = o.full.row$tie.pr[1],
    selected.order.full.is.low.p = as.integer(obs$max.order.full$p[1] < 0.5),
    max.pvalue.maxOrderWR_full = test.out$p.value.max["maxOrderWR_full"],
    fixed.pvalue.maxOrderWR_full = test.out$p.value.fixed["fixedOrder_full"],

    maxWRt = test.out$T.obs.max["maxWRt"],
    selected.t.years = obs$max.threshold$t[1],
    selected.t.months = obs$max.threshold$t.months[1],
    threshold.tie.count = obs$max.threshold$tie.count[1],
    threshold.tie.pr = obs$max.threshold$pr.tie[1],
    max.pvalue.maxWRt = test.out$p.value.max["maxWRt"],
    fixed.pvalue.maxWRt = test.out$p.value.fixed["fixedWRt"],

    true.hierarchical.tie.count = obs$counts$true.tie.pairs,
    true.hierarchical.tie.pr = obs$counts$true.tie.pr,
    total.pairs = obs$counts$total.pairs,

    logrank.death.statistic = logrank.results$statistic[logrank.results$method == "Log-rank death endpoint"],
    logrank.death.p = logrank.results$p.value[logrank.results$method == "Log-rank death endpoint"],
    logrank.composite.statistic = logrank.results$statistic[logrank.results$method == "Log-rank composite endpoint"],
    logrank.composite.p = logrank.results$p.value[logrank.results$method == "Log-rank composite endpoint"],

    death.event.rate.control = get_comp_arm_value(comp.stats, 0, "death.event.rate"),
    death.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "death.event.rate"),
    composite.event.rate.control = get_comp_arm_value(comp.stats, 0, "composite.event.rate"),
    composite.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "composite.event.rate"),
    mean.num.hosp.control = get_comp_arm_value(comp.stats, 0, "mean.num.hosp"),
    mean.num.hosp.treatment = get_comp_arm_value(comp.stats, 1, "mean.num.hosp"),

    threshold.grid.months = paste(sprintf("%.1f", threshold.info$threshold.table$t.months), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

make.realdata.significance.table <- function(pathway.tab, alpha = ALPHA) {
  out <- pathway.tab
  out$significant_by_permutation_p_value <- ifelse(is.na(out$permutation_p_value), NA_integer_, as.integer(out$permutation_p_value < alpha))
  out$significant_by_selected_parameter_p_value <- ifelse(is.na(out$selected_parameter_p_value), NA_integer_, as.integer(out$selected_parameter_p_value < alpha))
  out$alpha <- alpha
  out
}

plot.realdata.pathway.comparison <- function(pathway.tab, outdir, prefix, title.label) {
  all.methods <- c(
    "ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt",
    "logrank_death", "logrank_composite"
  )
  wr.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")

  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_permutation_pvalue_comparison.png"),
                  paste0(title.label, ": permutation/log-rank p-value comparison"),
                  "permutation_p_value", all.methods)

  fixed.tab <- pathway.tab[pathway.tab$method %in% wr.methods, , drop = FALSE]
  fixed.tab$method <- factor(fixed.tab$method, levels = wr.methods)
  fixed.tab <- fixed.tab[order(fixed.tab$method), , drop = FALSE]
  mat <- rbind(
    "Max-statistic" = fixed.tab$permutation_p_value,
    "Fixed-selected" = fixed.tab$selected_parameter_p_value
  )
  colnames(mat) <- method.display.label(as.character(fixed.tab$method))
  save_png(outdir, paste0(prefix, "_max_vs_fixed_selected_pvalues.png"), width = 2300, height = 1200)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  bp <- barplot(mat, beside = TRUE, ylim = c(0, 1), las = 1,
                ylab = "Permutation p-value",
                main = paste0(title.label, ": max-statistic vs fixed-selected p-values"),
                legend.text = TRUE,
                args.legend = list(x = "topright", bty = "n"),
                cex.names = 0.82)
  abline(h = ALPHA, lty = 2, lwd = 2)
  text(bp, pmin(mat + 0.04, 0.98), labels = sprintf("%.3f", mat), cex = 0.60)
  par(old.par)
  dev.off()
}

plot.global.realdata.pvalues <- function(pathway.table, outdir) {
  if (is.null(pathway.table) || nrow(pathway.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite")
  tab <- pathway.table[pathway.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  datasets <- unique(tab$dataset_id)
  dataset.labels <- sapply(datasets, function(x) tab$dataset_label[match(x, tab$dataset_id)])
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(datasets),
                dimnames = list(method.display.label(key.methods), dataset.labels))
  for (i in seq_len(nrow(tab))) {
    rr <- method.display.label(tab$method[i])
    cc <- tab$dataset_label[i]
    mat[rr, cc] <- tab$permutation_p_value[i]
  }
  save_png(outdir, "GLOBAL_permutation_pvalue_comparison_by_dataset.png", width = 3000, height = 1450)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2,
          ylab = "Permutation / log-rank p-value",
          main = "Permutation p-value comparison across real datasets",
          legend.text = TRUE,
          args.legend = list(x = "topright", bty = "n", cex = 0.72),
          cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2)
  par(old.par)
  dev.off()
}

plot.global.realdata.fixed.pvalues <- function(pathway.table, outdir) {
  if (is.null(pathway.table) || nrow(pathway.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")
  tab <- pathway.table[pathway.table$method %in% key.methods, , drop = FALSE]
  if (nrow(tab) == 0) return(invisible(NULL))
  datasets <- unique(tab$dataset_id)
  dataset.labels <- sapply(datasets, function(x) tab$dataset_label[match(x, tab$dataset_id)])
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(datasets),
                dimnames = list(method.display.label(key.methods), dataset.labels))
  for (i in seq_len(nrow(tab))) {
    rr <- method.display.label(tab$method[i])
    cc <- tab$dataset_label[i]
    mat[rr, cc] <- tab$selected_parameter_p_value[i]
  }
  save_png(outdir, "GLOBAL_fixed_selected_pvalue_comparison_by_dataset.png", width = 2600, height = 1350)
  old.par <- par(no.readonly = TRUE)
  par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2,
          ylab = "Fixed-selected permutation p-value",
          main = "Fixed-selected p-value comparison across real datasets",
          legend.text = TRUE,
          args.legend = list(x = "topright", bty = "n", cex = 0.75),
          cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2)
  par(old.par)
  dev.off()
}

summarise.average.across.datasets <- function(raw.table) {
  if (is.null(raw.table) || nrow(raw.table) == 0) return(data.frame())
  numeric.cols <- names(raw.table)[sapply(raw.table, is.numeric)]
  means <- sapply(numeric.cols, function(cc) safe.mean(raw.table[[cc]]))
  data.frame(
    n_datasets = nrow(raw.table),
    statistic = names(means),
    mean_across_datasets = as.numeric(means),
    stringsAsFactors = FALSE
  )
}

run.one.real.dataset <- function(registry.row,
                                 B = B_REAL,
                                 batch.size = BATCH_SIZE_REAL,
                                 outdir = OUTDIR,
                                 seed = MASTER_SEED,
                                 verbose = VERBOSE) {
  dataset.id <- sanitize_id(row_value(registry.row, "dataset_id", "realdata"))
  dataset.label <- dataset.display.label(dataset.id, row_value(registry.row, "dataset_label", dataset.id))
  dataset.outdir <- file.path(outdir, dataset.id)
  if (!dir.exists(dataset.outdir)) dir.create(dataset.outdir, recursive = TRUE)

  cat("\n===== Analyzing real dataset:", dataset.label, "=====\n")

  ds <- load.real.dataset(registry.row)
  ds <- prepare.ds.fast(ds)
  write.csv(ds$import.note, file.path(dataset.outdir, paste0(dataset.id, "_import_summary.csv")), row.names = FALSE)

  threshold.info <- choose.threshold.grid.primary(ds, clinical.months = CLINICAL_T_MONTHS, max.K = MAX_T_GRID_SIZE)
  write.csv(threshold.info$diagnostics, file.path(dataset.outdir, paste0(dataset.id, "_threshold_diagnostics_first_endpoint_only.csv")), row.names = FALSE)
  write.csv(threshold.info$threshold.table, file.path(dataset.outdir, paste0(dataset.id, "_threshold_grid_selected.csv")), row.names = FALSE)
  write.csv(threshold.info$source.table, file.path(dataset.outdir, paste0(dataset.id, "_threshold_candidate_sources.csv")), row.names = FALSE)

  test.out <- perm.test.revised.batch(
    ds = ds,
    B = B,
    batch.size = batch.size,
    seed = seed,
    p.grid = P_GRID_ALL,
    t.grid = threshold.info$t.grid,
    outdir = dataset.outdir,
    cache.prefix = dataset.id,
    resume = RESUME_IF_EXISTS,
    verbose = verbose
  )

  logrank.results <- run.logrank.tests(ds)
  comp.stats <- composite.statistics(ds)

  save.one.dataset.outputs(
    ds = ds,
    test.out = test.out,
    threshold.info = threshold.info,
    logrank.results = logrank.results,
    comp.stats = comp.stats,
    outdir = dataset.outdir,
    prefix = dataset.id
  )

  pathway.tab <- make.pathway.results.table(test.out, logrank.results)
  pathway.tab$dataset_id <- dataset.id
  pathway.tab$dataset_label <- dataset.label
  pathway.tab <- pathway.tab[, c("dataset_id", "dataset_label", setdiff(names(pathway.tab), c("dataset_id", "dataset_label"))), drop = FALSE]
  write.csv(pathway.tab, file.path(dataset.outdir, paste0(dataset.id, "_pathway_method_comparison.csv")), row.names = FALSE)

  significance.tab <- make.realdata.significance.table(pathway.tab, alpha = ALPHA)
  write.csv(significance.tab, file.path(dataset.outdir, paste0(dataset.id, "_significance_indicators.csv")), row.names = FALSE)

  summary.row <- make.realdata.summary.row(
    dataset.id = dataset.id,
    dataset.label = dataset.label,
    test.out = test.out,
    logrank.results = logrank.results,
    comp.stats = comp.stats,
    threshold.info = threshold.info,
    import.note = ds$import.note
  )
  write.csv(summary.row, file.path(dataset.outdir, paste0(dataset.id, "_observed_summary_statistics.csv")), row.names = FALSE)

  # Keep the same naming convention as the simulation files, but in real data this
  # is a single-trial observed table, not an average over simulated trials.
  write.csv(summary.row, file.path(dataset.outdir, paste0(dataset.id, "_average_statistics_single_real_trial.csv")), row.names = FALSE)

  pointwise.tab <- make.pointwise.table(test.out)
  pointwise.tab$dataset_id <- dataset.id
  pointwise.tab$dataset_label <- dataset.label
  pointwise.tab <- pointwise.tab[, c("dataset_id", "dataset_label", setdiff(names(pointwise.tab), c("dataset_id", "dataset_label"))), drop = FALSE]
  write.csv(pointwise.tab, file.path(dataset.outdir, paste0(dataset.id, "_pointwise_results.csv")), row.names = FALSE)

  plot.realdata.pathway.comparison(pathway.tab, dataset.outdir, dataset.id, dataset.label)

  saveRDS(
    list(
      ds = ds,
      threshold.info = threshold.info,
      test.out = test.out,
      logrank.results = logrank.results,
      composite.statistics = comp.stats,
      pathway.table = pathway.tab,
      significance.table = significance.tab,
      summary.row = summary.row,
      pointwise.table = pointwise.tab
    ),
    file.path(dataset.outdir, paste0(dataset.id, "_full_results_bundle.rds"))
  )

  list(
    dataset_id = dataset.id,
    dataset_label = dataset.label,
    raw = summary.row,
    pathway = pathway.tab,
    significance = significance.tab,
    pointwise = pointwise.tab,
    threshold.info = threshold.info
  )
}

run.all.real.datasets <- function(dataset.registry = REAL_DATASETS,
                                  B = B_REAL,
                                  batch.size = BATCH_SIZE_REAL,
                                  outdir = OUTDIR,
                                  master.seed = MASTER_SEED) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  write.realdata.templates(outdir)

  dataset.registry$dataset_id <- sanitize_id(dataset.registry$dataset_id)
  if (!("dataset_label" %in% names(dataset.registry))) dataset.registry$dataset_label <- dataset.registry$dataset_id
  write.csv(dataset.registry, file.path(outdir, "GLOBAL_realdata_registry.csv"), row.names = FALSE)
  write.csv(threshold.candidate.source.table(CLINICAL_T_MONTHS), file.path(outdir, "GLOBAL_threshold_candidate_sources.csv"), row.names = FALSE)

  settings.table <- data.frame(
    setting = c(
      "B_REAL", "BATCH_SIZE_REAL", "MASTER_SEED", "ALPHA", "number_of_datasets",
      "p_grid_low_exploratory", "p_grid_primary", "p_grid_full", "clinical_t_months"
    ),
    value = c(
      as.character(B),
      as.character(batch.size),
      as.character(master.seed),
      as.character(ALPHA),
      as.character(nrow(dataset.registry)),
      paste0(sprintf("%.2f", min(P_GRID_EXPLORATORY)), " to ", sprintf("%.2f", max(P_GRID_EXPLORATORY))),
      paste0(sprintf("%.2f", min(P_GRID_PRIMARY)), " to ", sprintf("%.2f", max(P_GRID_PRIMARY))),
      paste0(sprintf("%.2f", min(P_GRID_ALL)), " to ", sprintf("%.2f", max(P_GRID_ALL)), "; no p = 0"),
      paste(CLINICAL_T_MONTHS, collapse = ", ")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(settings.table, file.path(outdir, "GLOBAL_settings.csv"), row.names = FALSE)

  cat("\n===== Real-data analysis settings =====\n")
  cat("Output folder:", outdir, "\n")
  cat("B_REAL =", B, "; batch size =", batch.size, "; datasets =", nrow(dataset.registry), "\n")
  cat("Primary weighted p grid: [", min(P_GRID_PRIMARY), ", ", max(P_GRID_PRIMARY), "]\n", sep = "")
  cat("Exploratory low-p grid: [", min(P_GRID_EXPLORATORY), ", ", max(P_GRID_EXPLORATORY), "]\n", sep = "")

  all.raw <- data.frame()
  all.pathway <- data.frame()
  all.significance <- data.frame()
  all.pointwise <- data.frame()

  for (i in seq_len(nrow(dataset.registry))) {
    out.i <- run.one.real.dataset(
      registry.row = dataset.registry[i, , drop = FALSE],
      B = B,
      batch.size = batch.size,
      outdir = outdir,
      seed = master.seed + i * 10000L,
      verbose = VERBOSE
    )
    all.raw <- rbind_fill_base(all.raw, out.i$raw)
    all.pathway <- rbind_fill_base(all.pathway, out.i$pathway)
    all.significance <- rbind_fill_base(all.significance, out.i$significance)
    all.pointwise <- rbind_fill_base(all.pointwise, out.i$pointwise)

    write.csv(all.raw, file.path(outdir, "GLOBAL_all_datasets_observed_summary_statistics.csv"), row.names = FALSE)
    write.csv(all.pathway, file.path(outdir, "GLOBAL_all_datasets_pathway_method_comparison.csv"), row.names = FALSE)
    write.csv(all.significance, file.path(outdir, "GLOBAL_all_datasets_significance_indicators.csv"), row.names = FALSE)
    write.csv(all.pointwise, file.path(outdir, "GLOBAL_all_datasets_pointwise_results.csv"), row.names = FALSE)
  }

  avg.table <- summarise.average.across.datasets(all.raw)
  write.csv(avg.table, file.path(outdir, "GLOBAL_average_across_datasets_table_only.csv"), row.names = FALSE)

  plot.global.realdata.pvalues(all.pathway, outdir)
  plot.global.realdata.fixed.pvalues(all.pathway, outdir)

  manifest <- data.frame(file = list.files(outdir, recursive = TRUE), stringsAsFactors = FALSE)
  write.csv(manifest, file.path(outdir, "GLOBAL_output_manifest.csv"), row.names = FALSE)

  saveRDS(
    list(
      registry = dataset.registry,
      raw = all.raw,
      pathway = all.pathway,
      significance = all.significance,
      pointwise = all.pointwise,
      average.across.datasets = avg.table,
      settings = list(
        B_REAL = B,
        BATCH_SIZE_REAL = batch.size,
        MASTER_SEED = master.seed,
        ALPHA = ALPHA,
        P_GRID_EXPLORATORY = P_GRID_EXPLORATORY,
        P_GRID_PRIMARY = P_GRID_PRIMARY,
        P_GRID_ALL = P_GRID_ALL,
        CLINICAL_T_MONTHS = CLINICAL_T_MONTHS
      )
    ),
    file.path(outdir, "GLOBAL_full_results_bundle.rds")
  )

  cat("\n===== Real-data analysis complete =====\n")
  cat("Output folder:", outdir, "\n")
  cat("Main tables:\n")
  cat("  GLOBAL_all_datasets_pathway_method_comparison.csv\n")
  cat("  GLOBAL_all_datasets_observed_summary_statistics.csv\n")
  cat("  GLOBAL_all_datasets_significance_indicators.csv\n")
  cat("  GLOBAL_average_across_datasets_table_only.csv\n")
  cat("  GLOBAL_all_datasets_pointwise_results.csv\n")
  cat("Main global figures:\n")
  cat("  GLOBAL_permutation_pvalue_comparison_by_dataset.png\n")
  cat("  GLOBAL_fixed_selected_pvalue_comparison_by_dataset.png\n")
  cat("Each dataset folder contains DIG-style curves, KM curves, permutation null plots, and pathway p-value comparisons.\n")

  invisible(list(
    raw = all.raw,
    pathway = all.pathway,
    significance = all.significance,
    pointwise = all.pointwise,
    average.across.datasets = avg.table
  ))
}


#Traditional order-sensitivity
#Fixed p = 0.50, compare death-first vs hospitalization-first.
#This block intentionally overrides selected functions above while
#keeping the original DIG-style workflow and output organization.
fast.wr.engine.revised <- function(ds,
                                   p.grid = P_GRID_ALL,
                                   t.grid = CLINICAL_T_MONTHS / 12,
                                   eps = EPS_WR) {
  if (is.null(ds$hosp.flat)) ds <- prepare.ds.fast(ds)
  tab <- ds$table.output

  core <- fast_wr_core_revised_cpp(
    futime = as.numeric(tab$FUTIME),
    cnsr = as.integer(tab$CNSR),
    arm = as.integer(tab$ARM),
    hosp_times = as.numeric(ds$hosp.flat),
    hosp_start = as.integer(ds$hosp.start),
    hosp_len = as.integer(ds$hosp.len),
    p_grid = as.numeric(p.grid),
    t_grid = as.numeric(t.grid),
    eps = eps
  )

  total.pairs <- as.numeric(core$total_pairs)
  D1.win <- as.numeric(core$D1_win)
  D1.loss <- as.numeric(core$D1_loss)
  H2.win <- as.numeric(core$H2_win)
  H2.loss <- as.numeric(core$H2_loss)
  H1.win <- as.numeric(core$H1_win)
  H1.loss <- as.numeric(core$H1_loss)
  D2.win <- as.numeric(core$D2_win)
  D2.loss <- as.numeric(core$D2_loss)

  weighted.death.first <- data.frame(
    order = "death_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_death_first),
    win.score = as.numeric(core$win_death_first),
    loss.score = as.numeric(core$loss_death_first),
    stringsAsFactors = FALSE
  )
  weighted.death.first$win.pairs <- ifelse(weighted.death.first$p >= 1, D1.win, D1.win + H2.win)
  weighted.death.first$loss.pairs <- ifelse(weighted.death.first$p >= 1, D1.loss, D1.loss + H2.loss)
  weighted.death.first$tie.count <- total.pairs - weighted.death.first$win.pairs - weighted.death.first$loss.pairs
  weighted.death.first$tie.pr <- weighted.death.first$tie.count / total.pairs

  weighted.hosp.first <- data.frame(
    order = "hospitalization_first",
    p = as.numeric(core$p_grid),
    WR = as.numeric(core$WR_hosp_first),
    win.score = as.numeric(core$win_hosp_first),
    loss.score = as.numeric(core$loss_hosp_first),
    stringsAsFactors = FALSE
  )
  weighted.hosp.first$win.pairs <- ifelse(weighted.hosp.first$p >= 1, H1.win, H1.win + D2.win)
  weighted.hosp.first$loss.pairs <- ifelse(weighted.hosp.first$p >= 1, H1.loss, H1.loss + D2.loss)
  weighted.hosp.first$tie.count <- total.pairs - weighted.hosp.first$win.pairs - weighted.hosp.first$loss.pairs
  weighted.hosp.first$tie.pr <- weighted.hosp.first$tie.count / total.pairs

  weighted.all.orders <- rbind(weighted.death.first, weighted.hosp.first)

  threshold <- data.frame(
    t = as.numeric(core$t_grid),
    t.months = as.numeric(core$t_grid) * 12,
    WR = as.numeric(core$WRt),
    wins = as.numeric(core$t_wins),
    losses = as.numeric(core$t_losses),
    pr.win = as.numeric(core$pr_win_t),
    pr.loss = as.numeric(core$pr_loss_t),
    pr.tie = as.numeric(core$pr_tie_t),
    stringsAsFactors = FALSE
  )
  threshold$tie.count <- total.pairs - threshold$wins - threshold$losses
  best.t.idx <- which.max(threshold$WR)

  p.low.values <- p.grid[p.grid < 0.5]
  p.primary.values <- p.grid[p.grid >= 0.5]
  p.low.lower <- if (length(p.low.values) > 0) min(p.low.values) else NA_real_
  p.low.upper <- if (length(p.low.values) > 0) max(p.low.values) else NA_real_
  p.primary.lower <- if (length(p.primary.values) > 0) min(p.primary.values) else NA_real_
  p.primary.upper <- if (length(p.primary.values) > 0) max(p.primary.values) else NA_real_

  traditional.death.first <- max_from_curve(weighted.death.first, "ordinaryWR_death_first", 0.50, 0.50)
  traditional.hosp.first <- max_from_curve(weighted.hosp.first, "ordinaryWR_hosp_first", 0.50, 0.50)
  max.traditional.order <- max_from_curve(weighted.all.orders, "traditionalOrderWR", 0.50, 0.50)

  list(
    ordinaryWR = as.numeric(core$ordinaryWR),
    ordinary.win.score = as.numeric(core$ordinary_win_score),
    ordinary.loss.score = as.numeric(core$ordinary_loss_score),
    traditional.death.first = traditional.death.first,
    traditional.hosp.first = traditional.hosp.first,
    max.traditional.order = max.traditional.order,
    weighted.death.first = weighted.death.first,
    weighted.hosp.first = weighted.hosp.first,
    weighted.all.orders = weighted.all.orders,
    max.weighted.primary = max_from_curve(weighted.death.first, "maxWRp_primary_death_first", p.primary.lower, p.primary.upper),
    max.weighted.low = max_from_curve(weighted.death.first, "maxWRp_low_death_first", p.low.lower, p.low.upper),
    max.weighted.full = max_from_curve(weighted.death.first, "maxWRp_full_death_first", min(p.grid), max(p.grid)),
    max.order.primary = max_from_curve(weighted.all.orders, "maxOrderWR_primary", p.primary.lower, p.primary.upper),
    max.order.full = max_from_curve(weighted.all.orders, "maxOrderWR_full", min(p.grid), max(p.grid)),
    threshold = threshold,
    max.threshold = threshold[best.t.idx, , drop = FALSE],
    counts = list(
      total.pairs = as.numeric(core$total_pairs),
      D1.win = as.numeric(core$D1_win),
      D1.loss = as.numeric(core$D1_loss),
      H2.win = as.numeric(core$H2_win),
      H2.loss = as.numeric(core$H2_loss),
      H1.win = as.numeric(core$H1_win),
      H1.loss = as.numeric(core$H1_loss),
      D2.win = as.numeric(core$D2_win),
      D2.loss = as.numeric(core$D2_loss),
      death.win.pairs = as.numeric(core$death_win_pairs),
      death.loss.pairs = as.numeric(core$death_loss_pairs),
      hosp.win.pairs = as.numeric(core$hosp_win_pairs),
      hosp.loss.pairs = as.numeric(core$hosp_loss_pairs),
      true.tie.pairs = as.numeric(core$true_tie_pairs),
      true.tie.pr = as.numeric(core$true_tie_pairs) / as.numeric(core$total_pairs)
    )
  )
}

traditional_order_stat <- function(engine.out) {
  death.wr <- wr_at_order_p(engine.out, "death_first", 0.50)
  hosp.wr <- wr_at_order_p(engine.out, "hospitalization_first", 0.50)
  selected.order <- if (!is.na(hosp.wr) && (is.na(death.wr) || hosp.wr > death.wr)) {
    "hospitalization_first"
  } else {
    "death_first"
  }
  data.frame(
    WR = max(c(death.wr, hosp.wr), na.rm = TRUE),
    order = selected.order,
    p = 0.50,
    death.first.WR = death.wr,
    hosp.first.WR = hosp.wr,
    stringsAsFactors = FALSE
  )
}

perm.test.revised <- function(ds,
                              B = B_REAL,
                              seed = MASTER_SEED,
                              p.grid = P_GRID_ALL,
                              t.grid = CLINICAL_T_MONTHS / 12,
                              verbose = FALSE) {
  ds <- prepare.ds.fast(ds)
  set.seed(seed)
  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR",
    "traditionalWR_hosp_first",
    "traditionalOrderWR",
    "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.names <- c(
    "ordinaryWR",
    "fixedTraditionalHospFirst",
    "fixedTraditionalOrderWR",
    "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )

  fixed.traditional.order <- obs$max.traditional.order$order[1]
  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    traditionalWR_hosp_first = wr_at_order_p(obs, "hospitalization_first", 0.50),
    traditionalOrderWR = obs$max.traditional.order$WR[1],
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )

  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedTraditionalHospFirst = wr_at_order_p(obs, "hospitalization_first", 0.50),
    fixedTraditionalOrderWR = wr_at_order_p(obs, fixed.traditional.order, 0.50),
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names)); colnames(perm.max) <- max.names
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names)); colnames(perm.fixed) <- fixed.names

  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.traditional.order = NA_character_,
    selected.p.primary = NA_real_, selected.p.low = NA_real_, selected.p.full = NA_real_,
    selected.order.primary = NA_character_, selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_, selected.order.p.full = NA_real_,
    selected.t = NA_real_, stringsAsFactors = FALSE
  )

  perm.seeds <- seed + seq_len(B) * 1009L
  for (b in seq_len(B)) {
    if (verbose && (b == 1 || b == B || b %% 50 == 0)) cat("  Permutation", b, "of", B, "\n")
    set.seed(perm.seeds[b])
    ds.b <- ds
    ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)
    out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)

    perm.pointwise.death.first[b, ] <- out.b$weighted.death.first$WR
    perm.pointwise.hosp.first[b, ] <- out.b$weighted.hosp.first$WR
    perm.pointwise.threshold[b, ] <- out.b$threshold$WR

    perm.max[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.max[b, "traditionalWR_hosp_first"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
    perm.max[b, "traditionalOrderWR"] <- out.b$max.traditional.order$WR[1]
    perm.max[b, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
    perm.max[b, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
    perm.max[b, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
    perm.max[b, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
    perm.max[b, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
    perm.max[b, "maxWRt"] <- out.b$max.threshold$WR[1]

    perm.fixed[b, "ordinaryWR"] <- out.b$ordinaryWR
    perm.fixed[b, "fixedTraditionalHospFirst"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
    perm.fixed[b, "fixedTraditionalOrderWR"] <- wr_at_order_p(out.b, fixed.traditional.order, 0.50)
    perm.fixed[b, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
    perm.fixed[b, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
    perm.fixed[b, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
    perm.fixed[b, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
    perm.fixed[b, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
    perm.fixed[b, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

    selected.perm$selected.traditional.order[b] <- out.b$max.traditional.order$order[1]
    selected.perm$selected.p.primary[b] <- out.b$max.weighted.primary$p[1]
    selected.perm$selected.p.low[b] <- out.b$max.weighted.low$p[1]
    selected.perm$selected.p.full[b] <- out.b$max.weighted.full$p[1]
    selected.perm$selected.order.primary[b] <- out.b$max.order.primary$order[1]
    selected.perm$selected.order.p.primary[b] <- out.b$max.order.primary$p[1]
    selected.perm$selected.order.full[b] <- out.b$max.order.full$order[1]
    selected.perm$selected.order.p.full[b] <- out.b$max.order.full$p[1]
    selected.perm$selected.t[b] <- out.b$max.threshold$t[1]
  }

  calc_perm_p <- function(perm.vec, obs.value) (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))
  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k]))
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k]))
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k]))

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(death_first = perm.pointwise.death.first, hospitalization_first = perm.pointwise.hosp.first, threshold = perm.pointwise.threshold),
    p.value.pointwise = list(death_first = p.value.pointwise.death.first, hospitalization_first = p.value.pointwise.hosp.first, threshold = p.value.pointwise.threshold),
    selected.fixed = list(
      traditional.order = fixed.traditional.order,
      traditional.p = 0.50,
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    p.grid = p.grid,
    t.grid = t.grid
  )
}

method.long.label <- function(method) {
  map <- c(
    ordinaryWR = "Traditional WR, death first",
    traditionalWR_hosp_first = "Traditional WR, hospitalization first",
    traditionalOrderWR = "Traditional order-selected WR",
    maxWRp_primary = "Weighted max WR, p >= 0.5",
    maxWRp_low = "Weighted max WR, p < 0.5",
    maxWRp_full = "Weighted max WR, full p grid",
    maxOrderWR_primary = "Maximum-order weighted WR, p >= 0.5",
    maxOrderWR_full = "Maximum-order weighted WR, full p grid",
    maxWRt = "Threshold max WR(t)",
    logrank_death = "Log-rank death",
    logrank_composite = "Log-rank composite",
    fixedTraditionalHospFirst = "Traditional fixed WR, hospitalization first",
    fixedTraditionalOrderWR = "Traditional fixed-selected order WR",
    fixedWRp_primary = "Weighted fixed-selected WR, p >= 0.5",
    fixedWRp_low = "Weighted fixed-selected WR, p < 0.5",
    fixedWRp_full = "Weighted fixed-selected WR, full p grid",
    fixedOrderWR_primary = "Maximum-order fixed-selected WR, p >= 0.5",
    fixedOrderWR_full = "Maximum-order fixed-selected WR, full p grid",
    fixedWRt = "Threshold fixed-selected WR(t)"
  )
  out <- as.character(method); hit <- out %in% names(map); out[hit] <- unname(map[out[hit]]); out[!hit] <- gsub("_", " ", out[!hit]); out
}

method.display.label <- function(method) {
  map <- c(
    ordinaryWR = "Trad death\nWR",
    traditionalWR_hosp_first = "Trad hosp\nWR",
    traditionalOrderWR = "Trad order\nmax",
    maxWRp_primary = "Weighted\nmax",
    maxWRp_low = "Low-p\nmax",
    maxWRp_full = "Full-p\nmax",
    maxOrderWR_primary = "Weighted order\nmax",
    maxOrderWR_full = "Order full\nmax",
    maxWRt = "Threshold\nmax",
    logrank_death = "LR\ndeath",
    logrank_composite = "LR\ncomposite",
    fixedTraditionalHospFirst = "Trad hosp\nfixed",
    fixedTraditionalOrderWR = "Trad order\nfixed",
    fixedWRp_primary = "Weighted\nfixed",
    fixedWRp_low = "Low-p\nfixed",
    fixedWRp_full = "Full-p\nfixed",
    fixedOrderWR_primary = "Weighted order\nfixed",
    fixedOrderWR_full = "Order full\nfixed",
    fixedWRt = "Threshold\nfixed"
  )
  out <- as.character(method); hit <- out %in% names(map); out[hit] <- unname(map[out[hit]]); out[!hit] <- gsub("_", " ", out[!hit]); out
}

make.max.vs.fixed.table <- function(test.out) {
  method.ids <- c(
    "ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR",
    "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.ids <- c(
    "ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR",
    "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )
  data.frame(
    method.id = method.ids,
    method = method.long.label(method.ids),
    selected.parameter = c(
      "p = 0.50, death first",
      "p = 0.50, hospitalization first",
      paste0("selected order = ", clean.order.label(test.out$selected.fixed$traditional.order), "; p = 0.50"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), ", death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), ", death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), ", death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary), ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full), ", p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    max.statistic = as.numeric(test.out$T.obs.max[method.ids]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[method.ids]),
    fixed.parameter.statistic = as.numeric(test.out$T.obs.fixed[fixed.ids]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[fixed.ids]),
    stringsAsFactors = FALSE
  )
}

make.selected.max.table <- function(test.out) {
  obs <- test.out$observed
  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional <- get_selected_weighted_row(obs, test.out$selected.fixed$traditional.order, 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)
  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])
  method.ids <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  fixed.ids <- c("ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full", "fixedOrder_primary", "fixedOrder_full", "fixedWRt")
  data.frame(
    method = method.long.label(method.ids),
    method.id = method.ids,
    selected.order = c("death_first", "hospitalization_first", test.out$selected.fixed$traditional.order, "death_first", "death_first", "death_first", test.out$selected.fixed$order.primary, test.out$selected.fixed$order.full, NA_character_),
    selected.p = c(0.5, 0.5, 0.5, test.out$selected.fixed$p.primary, test.out$selected.fixed$p.low, test.out$selected.fixed$p.full, test.out$selected.fixed$order.p.primary, test.out$selected.fixed$order.p.full, NA_real_),
    selected.t.months = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, test.out$selected.fixed$t.months),
    observed.max.statistic = as.numeric(test.out$T.obs.max[method.ids]),
    max.statistic.p.value = as.numeric(test.out$p.value.max[method.ids]),
    fixed.parameter.p.value = as.numeric(test.out$p.value.fixed[fixed.ids]),
    selected.tie.count = c(get_tie_count(w.traditional), get_tie_count(h.traditional), get_tie_count(to.traditional), get_tie_count(w.primary), get_tie_count(w.low), get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full), get_tie_count(t.row)),
    selected.tie.pr = c(get_tie_pr(w.traditional), get_tie_pr(h.traditional), get_tie_pr(to.traditional), get_tie_pr(w.primary), get_tie_pr(w.low), get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full), get_tie_pr(t.row)),
    stringsAsFactors = FALSE
  )
}

make.pathway.results.table <- function(test.out, logrank.results = NULL) {
  obs <- test.out$observed
  w.traditional <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional <- get_selected_weighted_row(obs, test.out$selected.fixed$traditional.order, 0.5)
  w.primary <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.primary)
  w.low <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.low)
  w.full <- get_selected_weighted_row(obs, "death_first", test.out$selected.fixed$p.full)
  o.primary <- get_selected_weighted_row(obs, test.out$selected.fixed$order.primary, test.out$selected.fixed$order.p.primary)
  o.full <- get_selected_weighted_row(obs, test.out$selected.fixed$order.full, test.out$selected.fixed$order.p.full)
  t.row <- get_selected_threshold_row(obs, test.out$selected.fixed$t)
  get_tie_count <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.count[1])
  get_tie_pr <- function(x) if (is.null(x)) NA_real_ else as.numeric(x$tie.pr[1])
  method.id <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  fixed.id <- c("ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR", "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full", "fixedOrder_primary", "fixedOrder_full", "fixedWRt")
  out <- data.frame(
    method = method.id,
    method_label = method.long.label(method.id),
    statistic_type = "Win ratio",
    statistic_name = c(
      "WR(p = 0.50), death first",
      "WR(p = 0.50), hospitalization first",
      "max WR(order), p = 0.50",
      "max WR(p), death first, p in [0.50, 1.00]",
      "max WR(p), death first, p in [0.01, 0.49]",
      "max WR(p), death first, p in [0.01, 1.00]",
      "max WR(order, p), p in [0.50, 1.00]",
      "max WR(order, p), p in [0.01, 1.00]",
      "max WR(t)"
    ),
    observed_statistic = as.numeric(test.out$T.obs.max[method.id]),
    selected_order = c("death_first", "hospitalization_first", test.out$selected.fixed$traditional.order, "death_first", "death_first", "death_first", test.out$selected.fixed$order.primary, test.out$selected.fixed$order.full, NA_character_),
    selected_p = c(0.5, 0.5, 0.5, test.out$selected.fixed$p.primary, test.out$selected.fixed$p.low, test.out$selected.fixed$p.full, test.out$selected.fixed$order.p.primary, test.out$selected.fixed$order.p.full, NA_real_),
    selected_t_months = c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, test.out$selected.fixed$t.months),
    selected_value = c(
      "p = 0.50; death first",
      "p = 0.50; hospitalization first",
      paste0("selected order = ", clean.order.label(test.out$selected.fixed$traditional.order), "; p = 0.50"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.primary), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.low), "; death first"),
      paste0("p = ", sprintf("%.2f", test.out$selected.fixed$p.full), "; death first"),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.primary), "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.primary)),
      paste0("order = ", clean.order.label(test.out$selected.fixed$order.full), "; p = ", sprintf("%.2f", test.out$selected.fixed$order.p.full)),
      paste0("t = ", sprintf("%.1f", test.out$selected.fixed$t.months), " months")
    ),
    permutation_p_value = as.numeric(test.out$p.value.max[method.id]),
    selected_parameter_p_value = as.numeric(test.out$p.value.fixed[fixed.id]),
    tie_count = c(get_tie_count(w.traditional), get_tie_count(h.traditional), get_tie_count(to.traditional), get_tie_count(w.primary), get_tie_count(w.low), get_tie_count(w.full), get_tie_count(o.primary), get_tie_count(o.full), get_tie_count(t.row)),
    tie_proportion = c(get_tie_pr(w.traditional), get_tie_pr(h.traditional), get_tie_pr(to.traditional), get_tie_pr(w.primary), get_tie_pr(w.low), get_tie_pr(w.full), get_tie_pr(o.primary), get_tie_pr(o.full), get_tie_pr(t.row)),
    B = test.out$B,
    stringsAsFactors = FALSE
  )
  if (!is.null(logrank.results)) {
    lr.death <- logrank.results[logrank.results$method == "Log-rank death endpoint", , drop = FALSE]
    lr.comp <- logrank.results[logrank.results$method == "Log-rank composite endpoint", , drop = FALSE]
    lr.rows <- data.frame(
      method = c("logrank_death", "logrank_composite"),
      method_label = method.long.label(c("logrank_death", "logrank_composite")),
      statistic_type = "Log-rank chi-square",
      statistic_name = c("Death endpoint", "Composite endpoint"),
      observed_statistic = c(lr.death$statistic[1], lr.comp$statistic[1]),
      selected_order = NA_character_, selected_p = NA_real_, selected_t_months = NA_real_, selected_value = NA_character_,
      permutation_p_value = c(lr.death$p.value[1], lr.comp$p.value[1]),
      selected_parameter_p_value = NA_real_, tie_count = NA_real_, tie_proportion = NA_real_, B = NA_integer_,
      stringsAsFactors = FALSE
    )
    out <- rbind(out, lr.rows)
  }
  out
}

plot.order.curves <- function(engine.out, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_maximum_order_WRp.png"), width = 1900, height = 1150)
  old.par <- par(no.readonly = TRUE); par(mar = c(5.2, 5.2, 3.4, 1.2))
  dtab <- engine.out$weighted.death.first; htab <- engine.out$weighted.hosp.first
  yr <- range(c(dtab$WR, htab$WR), finite = TRUE)
  plot(dtab$p, dtab$WR, type = "l", lwd = 2, ylim = yr,
       xlab = "p: first-endpoint weight", ylab = "WR(p)", main = "Hierarchy order comparison")
  lines(htab$p, htab$WR, lwd = 2, lty = 2)
  points(0.50, engine.out$traditional.death.first$WR[1], pch = 1, cex = 1.5)
  points(0.50, engine.out$traditional.hosp.first$WR[1], pch = 2, cex = 1.5)
  points(engine.out$max.traditional.order$p[1], engine.out$max.traditional.order$WR[1], pch = 17, cex = 1.3)
  points(engine.out$max.order.primary$p, engine.out$max.order.primary$WR, pch = 19)
  abline(v = 0.5, lty = 3, lwd = 2)
  legend("topright",
         legend = c("death first", "hospitalization first",
                    paste0("traditional selected = ", clean.order.label(engine.out$max.traditional.order$order[1])),
                    paste0("traditional order WR = ", sprintf("%.3f", engine.out$max.traditional.order$WR[1])),
                    paste0("weighted selected = ", clean.order.label(engine.out$max.order.primary$order[1])),
                    paste0("weighted p = ", sprintf("%.2f", engine.out$max.order.primary$p[1]))),
         lty = c(1, 2, NA, NA, NA, NA), pch = c(NA, NA, 17, NA, 19, NA), bty = "n", cex = 0.82)
  par(old.par); dev.off()
}

plot.max.vs.fixed.pvalues <- function(max.vs.fixed.table, outdir, prefix) {
  save_png(outdir, paste0(prefix, "_max_vs_fixed_pvalues.png"), width = 2500, height = 1250)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  mat <- rbind("Max-statistic" = max.vs.fixed.table$max.statistic.p.value, "Fixed-selected" = max.vs.fixed.table$fixed.parameter.p.value)
  labels <- if ("method.id" %in% names(max.vs.fixed.table)) method.display.label(max.vs.fixed.table$method.id) else max.vs.fixed.table$method
  bp <- barplot(mat, beside = TRUE, names.arg = labels, ylim = c(0, 1), ylab = "Permutation p-value", main = "Max-statistic vs fixed-selected p-values", legend.text = TRUE, args.legend = list(x = "topright", bty = "n"), las = 2, cex.names = 0.74)
  abline(h = ALPHA, lty = 2, lwd = 2)
  text(x = bp, y = pmin(mat + 0.035, 0.98), labels = sprintf("%.3f", mat), cex = 0.58)
  par(old.par); dev.off()
}

plot.example.pvalue.figures <- function(pathway.tab, outdir, prefix) {
  all.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite")
  wr.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  weighted.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full")
  logrank.methods <- c("logrank_death", "logrank_composite")
  threshold.methods <- c("maxWRt")
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_all_methods_with_logrank.png"), "All method p-values", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_WR_pathways_only_no_logrank.png"), "WR pathway p-values", "permutation_p_value", wr.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_traditional_and_weighted_only.png"), "Traditional/order and weighted WR p-values", "permutation_p_value", weighted.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_logrank_pvalues_separate.png"), "Log-rank p-values", "permutation_p_value", logrank.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_threshold_pvalue_separate.png"), "Threshold pathway p-value", "permutation_p_value", threshold.methods)
}

plot.example.permutation.nulls <- function(test.out, outdir, prefix) {
  for (m in c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")) {
    plot.permutation.null.revised(test.out, m, outdir, prefix)
  }
}


perm.test.revised.batch <- function(ds,
                                    B = B_REAL,
                                    batch.size = BATCH_SIZE_REAL,
                                    seed = MASTER_SEED,
                                    p.grid = P_GRID_ALL,
                                    t.grid = CLINICAL_T_MONTHS / 12,
                                    outdir = OUTDIR,
                                    cache.prefix = "realdata",
                                    resume = RESUME_IF_EXISTS,
                                    verbose = VERBOSE) {
  ds <- prepare.ds.fast(ds)
  obs <- fast.wr.engine.revised(ds, p.grid = p.grid, t.grid = t.grid)

  max.names <- c(
    "ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR",
    "maxWRp_primary", "maxWRp_low", "maxWRp_full",
    "maxOrderWR_primary", "maxOrderWR_full", "maxWRt"
  )
  fixed.names <- c(
    "ordinaryWR", "fixedTraditionalHospFirst", "fixedTraditionalOrderWR",
    "fixedWRp_primary", "fixedWRp_low", "fixedWRp_full",
    "fixedOrder_primary", "fixedOrder_full", "fixedWRt"
  )

  fixed.traditional.order <- obs$max.traditional.order$order[1]
  fixed.p.primary <- obs$max.weighted.primary$p[1]
  fixed.p.low <- obs$max.weighted.low$p[1]
  fixed.p.full <- obs$max.weighted.full$p[1]
  fixed.order.primary <- obs$max.order.primary$order[1]
  fixed.order.p.primary <- obs$max.order.primary$p[1]
  fixed.order.full <- obs$max.order.full$order[1]
  fixed.order.p.full <- obs$max.order.full$p[1]
  fixed.t <- obs$max.threshold$t[1]

  T.obs.max <- c(
    ordinaryWR = obs$ordinaryWR,
    traditionalWR_hosp_first = wr_at_order_p(obs, "hospitalization_first", 0.50),
    traditionalOrderWR = obs$max.traditional.order$WR[1],
    maxWRp_primary = obs$max.weighted.primary$WR[1],
    maxWRp_low = obs$max.weighted.low$WR[1],
    maxWRp_full = obs$max.weighted.full$WR[1],
    maxOrderWR_primary = obs$max.order.primary$WR[1],
    maxOrderWR_full = obs$max.order.full$WR[1],
    maxWRt = obs$max.threshold$WR[1]
  )
  T.obs.fixed <- c(
    ordinaryWR = obs$ordinaryWR,
    fixedTraditionalHospFirst = wr_at_order_p(obs, "hospitalization_first", 0.50),
    fixedTraditionalOrderWR = wr_at_order_p(obs, fixed.traditional.order, 0.50),
    fixedWRp_primary = wr_at_order_p(obs, "death_first", fixed.p.primary),
    fixedWRp_low = wr_at_order_p(obs, "death_first", fixed.p.low),
    fixedWRp_full = wr_at_order_p(obs, "death_first", fixed.p.full),
    fixedOrder_primary = wr_at_order_p(obs, fixed.order.primary, fixed.order.p.primary),
    fixedOrder_full = wr_at_order_p(obs, fixed.order.full, fixed.order.p.full),
    fixedWRt = threshold_at_t(obs, fixed.t)
  )

  perm.max <- matrix(NA_real_, nrow = B, ncol = length(max.names)); colnames(perm.max) <- max.names
  perm.fixed <- matrix(NA_real_, nrow = B, ncol = length(fixed.names)); colnames(perm.fixed) <- fixed.names
  perm.pointwise.death.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.hosp.first <- matrix(NA_real_, nrow = B, ncol = length(p.grid))
  perm.pointwise.threshold <- matrix(NA_real_, nrow = B, ncol = length(t.grid))
  colnames(perm.pointwise.death.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.hosp.first) <- paste0("p_", sprintf("%.2f", p.grid))
  colnames(perm.pointwise.threshold) <- paste0("t_months_", sprintf("%.1f", t.grid * 12))

  selected.perm <- data.frame(
    b = seq_len(B),
    selected.traditional.order = NA_character_,
    selected.p.primary = NA_real_,
    selected.p.low = NA_real_,
    selected.p.full = NA_real_,
    selected.order.primary = NA_character_,
    selected.order.p.primary = NA_real_,
    selected.order.full = NA_character_,
    selected.order.p.full = NA_real_,
    selected.t = NA_real_,
    stringsAsFactors = FALSE
  )

  perm.seeds <- seed + seq_len(B) * 1009L
  batch.dir <- file.path(outdir, paste0(cache.prefix, "_tradorder_permutation_batches"))
  if (!dir.exists(batch.dir)) dir.create(batch.dir, recursive = TRUE)
  n.batch <- ceiling(B / batch.size)

  for (bb in seq_len(n.batch)) {
    start.b <- (bb - 1L) * batch.size + 1L
    end.b <- min(bb * batch.size, B)
    idx <- start.b:end.b
    batch.file <- file.path(batch.dir, sprintf("perm_batch_%04d_%04d_to_%04d.rds", bb, start.b, end.b))

    if (resume && file.exists(batch.file)) {
      if (verbose) cat("Loading existing batch", bb, "of", n.batch, "\n")
      batch.out <- readRDS(batch.file)
    } else {
      if (verbose) cat("Running batch", bb, "of", n.batch, ": permutations", start.b, "to", end.b, "\n")
      batch.max <- matrix(NA_real_, nrow = length(idx), ncol = length(max.names)); colnames(batch.max) <- max.names
      batch.fixed <- matrix(NA_real_, nrow = length(idx), ncol = length(fixed.names)); colnames(batch.fixed) <- fixed.names
      batch.death <- matrix(NA_real_, nrow = length(idx), ncol = length(p.grid))
      batch.hosp <- matrix(NA_real_, nrow = length(idx), ncol = length(p.grid))
      batch.thr <- matrix(NA_real_, nrow = length(idx), ncol = length(t.grid))
      batch.selected <- selected.perm[idx, , drop = FALSE]

      for (r in seq_along(idx)) {
        b <- idx[r]
        set.seed(perm.seeds[b])
        ds.b <- ds
        ds.b$table.output$ARM <- sample(ds$table.output$ARM, replace = FALSE)
        out.b <- fast.wr.engine.revised(ds.b, p.grid = p.grid, t.grid = t.grid)

        batch.death[r, ] <- out.b$weighted.death.first$WR
        batch.hosp[r, ] <- out.b$weighted.hosp.first$WR
        batch.thr[r, ] <- out.b$threshold$WR
        batch.max[r, "ordinaryWR"] <- out.b$ordinaryWR
        batch.max[r, "traditionalWR_hosp_first"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
        batch.max[r, "traditionalOrderWR"] <- out.b$max.traditional.order$WR[1]
        batch.max[r, "maxWRp_primary"] <- out.b$max.weighted.primary$WR[1]
        batch.max[r, "maxWRp_low"] <- out.b$max.weighted.low$WR[1]
        batch.max[r, "maxWRp_full"] <- out.b$max.weighted.full$WR[1]
        batch.max[r, "maxOrderWR_primary"] <- out.b$max.order.primary$WR[1]
        batch.max[r, "maxOrderWR_full"] <- out.b$max.order.full$WR[1]
        batch.max[r, "maxWRt"] <- out.b$max.threshold$WR[1]

        batch.fixed[r, "ordinaryWR"] <- out.b$ordinaryWR
        batch.fixed[r, "fixedTraditionalHospFirst"] <- wr_at_order_p(out.b, "hospitalization_first", 0.50)
        batch.fixed[r, "fixedTraditionalOrderWR"] <- wr_at_order_p(out.b, fixed.traditional.order, 0.50)
        batch.fixed[r, "fixedWRp_primary"] <- wr_at_order_p(out.b, "death_first", fixed.p.primary)
        batch.fixed[r, "fixedWRp_low"] <- wr_at_order_p(out.b, "death_first", fixed.p.low)
        batch.fixed[r, "fixedWRp_full"] <- wr_at_order_p(out.b, "death_first", fixed.p.full)
        batch.fixed[r, "fixedOrder_primary"] <- wr_at_order_p(out.b, fixed.order.primary, fixed.order.p.primary)
        batch.fixed[r, "fixedOrder_full"] <- wr_at_order_p(out.b, fixed.order.full, fixed.order.p.full)
        batch.fixed[r, "fixedWRt"] <- threshold_at_t(out.b, fixed.t)

        batch.selected$selected.traditional.order[r] <- out.b$max.traditional.order$order[1]
        batch.selected$selected.p.primary[r] <- out.b$max.weighted.primary$p[1]
        batch.selected$selected.p.low[r] <- out.b$max.weighted.low$p[1]
        batch.selected$selected.p.full[r] <- out.b$max.weighted.full$p[1]
        batch.selected$selected.order.primary[r] <- out.b$max.order.primary$order[1]
        batch.selected$selected.order.p.primary[r] <- out.b$max.order.primary$p[1]
        batch.selected$selected.order.full[r] <- out.b$max.order.full$order[1]
        batch.selected$selected.order.p.full[r] <- out.b$max.order.full$p[1]
        batch.selected$selected.t[r] <- out.b$max.threshold$t[1]
      }
      batch.out <- list(max = batch.max, fixed = batch.fixed, death = batch.death, hosp = batch.hosp, threshold = batch.thr, selected = batch.selected)
      saveRDS(batch.out, batch.file)
    }

    perm.max[idx, colnames(batch.out$max)] <- batch.out$max
    perm.fixed[idx, colnames(batch.out$fixed)] <- batch.out$fixed
    perm.pointwise.death.first[idx, ] <- batch.out$death
    perm.pointwise.hosp.first[idx, ] <- batch.out$hosp
    perm.pointwise.threshold[idx, ] <- batch.out$threshold
    selected.perm[idx, names(batch.out$selected)] <- batch.out$selected
  }

  calc_perm_p <- function(perm.vec, obs.value) {
    (1 + sum(perm.vec >= obs.value, na.rm = TRUE)) / (sum(!is.na(perm.vec)) + 1)
  }
  p.value.max <- sapply(names(T.obs.max), function(nm) calc_perm_p(perm.max[, nm], T.obs.max[nm]))
  p.value.fixed <- sapply(names(T.obs.fixed), function(nm) calc_perm_p(perm.fixed[, nm], T.obs.fixed[nm]))
  p.value.pointwise.death.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.death.first[, k], obs$weighted.death.first$WR[k]))
  p.value.pointwise.hosp.first <- sapply(seq_along(p.grid), function(k) calc_perm_p(perm.pointwise.hosp.first[, k], obs$weighted.hosp.first$WR[k]))
  p.value.pointwise.threshold <- sapply(seq_along(t.grid), function(k) calc_perm_p(perm.pointwise.threshold[, k], obs$threshold$WR[k]))

  list(
    observed = obs,
    T.obs.max = T.obs.max,
    T.perm.max = perm.max,
    p.value.max = p.value.max,
    T.obs.fixed = T.obs.fixed,
    T.perm.fixed = perm.fixed,
    p.value.fixed = p.value.fixed,
    T.perm.pointwise = list(death_first = perm.pointwise.death.first, hospitalization_first = perm.pointwise.hosp.first, threshold = perm.pointwise.threshold),
    p.value.pointwise = list(death_first = p.value.pointwise.death.first, hospitalization_first = p.value.pointwise.hosp.first, threshold = p.value.pointwise.threshold),
    selected.fixed = list(
      traditional.order = fixed.traditional.order,
      traditional.p = 0.50,
      p.primary = fixed.p.primary,
      p.low = fixed.p.low,
      p.full = fixed.p.full,
      order.primary = fixed.order.primary,
      order.p.primary = fixed.order.p.primary,
      order.full = fixed.order.full,
      order.p.full = fixed.order.p.full,
      t = fixed.t,
      t.months = fixed.t * 12
    ),
    selected.perm = selected.perm,
    B = B,
    p.grid = p.grid,
    t.grid = t.grid,
    batch.dir = batch.dir
  )
}

make.realdata.summary.row <- function(dataset.id, dataset.label, test.out, logrank.results,
                                      comp.stats, threshold.info, import.note) {
  obs <- test.out$observed
  death.first.primary <- obs$max.weighted.primary
  hosp.first.primary <- max_from_curve(obs$weighted.hosp.first, "hosp_first_primary", 0.50, 1.00)
  order.diff <- hosp.first.primary$WR[1] - death.first.primary$WR[1]

  w.traditional.row <- get_selected_weighted_row(obs, "death_first", 0.5)
  h.traditional.row <- get_selected_weighted_row(obs, "hospitalization_first", 0.5)
  to.traditional.row <- get_selected_weighted_row(obs, obs$max.traditional.order$order[1], 0.5)
  w.primary.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.primary$p[1])
  w.low.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.low$p[1])
  w.full.row <- get_selected_weighted_row(obs, "death_first", obs$max.weighted.full$p[1])
  o.primary.row <- get_selected_weighted_row(obs, obs$max.order.primary$order[1], obs$max.order.primary$p[1])
  o.full.row <- get_selected_weighted_row(obs, obs$max.order.full$order[1], obs$max.order.full$p[1])

  data.frame(
    dataset_id = dataset.id,
    dataset_label = dataset.label,
    n_subjects = import.note$n_subjects[1],
    n_control = import.note$n_control[1],
    n_treatment = import.note$n_treatment[1],

    ordinaryWR = test.out$T.obs.max["ordinaryWR"],
    traditional.tie.count = w.traditional.row$tie.count[1],
    traditional.tie.pr = w.traditional.row$tie.pr[1],
    max.pvalue.ordinaryWR = test.out$p.value.max["ordinaryWR"],
    fixed.pvalue.ordinaryWR = test.out$p.value.fixed["ordinaryWR"],

    traditionalWR_hosp_first = test.out$T.obs.max["traditionalWR_hosp_first"],
    traditional.hosp.first.tie.count = h.traditional.row$tie.count[1],
    traditional.hosp.first.tie.pr = h.traditional.row$tie.pr[1],
    max.pvalue.traditionalWR_hosp_first = test.out$p.value.max["traditionalWR_hosp_first"],
    fixed.pvalue.traditionalWR_hosp_first = test.out$p.value.fixed["fixedTraditionalHospFirst"],

    traditional.death.first.WR = w.traditional.row$WR[1],
    traditional.hosp.first.WR = h.traditional.row$WR[1],
    traditional.order.diff.hosp.minus.death = h.traditional.row$WR[1] - w.traditional.row$WR[1],
    traditional.order.changed = as.integer(obs$max.traditional.order$order[1] == "hospitalization_first"),
    traditionalOrderWR = test.out$T.obs.max["traditionalOrderWR"],
    selected.traditional.order = obs$max.traditional.order$order[1],
    traditional.order.tie.count = to.traditional.row$tie.count[1],
    traditional.order.tie.pr = to.traditional.row$tie.pr[1],
    max.pvalue.traditionalOrderWR = test.out$p.value.max["traditionalOrderWR"],
    fixed.pvalue.traditionalOrderWR = test.out$p.value.fixed["fixedTraditionalOrderWR"],

    maxWRp_primary = test.out$T.obs.max["maxWRp_primary"],
    selected.p.primary = obs$max.weighted.primary$p[1],
    weighted.primary.tie.count = w.primary.row$tie.count[1],
    weighted.primary.tie.pr = w.primary.row$tie.pr[1],
    max.pvalue.maxWRp_primary = test.out$p.value.max["maxWRp_primary"],
    fixed.pvalue.maxWRp_primary = test.out$p.value.fixed["fixedWRp_primary"],

    maxWRp_low = test.out$T.obs.max["maxWRp_low"],
    selected.p.low = obs$max.weighted.low$p[1],
    weighted.low.tie.count = w.low.row$tie.count[1],
    weighted.low.tie.pr = w.low.row$tie.pr[1],
    max.pvalue.maxWRp_low = test.out$p.value.max["maxWRp_low"],
    fixed.pvalue.maxWRp_low = test.out$p.value.fixed["fixedWRp_low"],

    maxWRp_full = test.out$T.obs.max["maxWRp_full"],
    selected.p.full = obs$max.weighted.full$p[1],
    weighted.full.tie.count = w.full.row$tie.count[1],
    weighted.full.tie.pr = w.full.row$tie.pr[1],
    selected.full.is.low.p = as.integer(obs$max.weighted.full$p[1] < 0.5),
    max.pvalue.maxWRp_full = test.out$p.value.max["maxWRp_full"],
    fixed.pvalue.maxWRp_full = test.out$p.value.fixed["fixedWRp_full"],

    death.first.primary.WR = death.first.primary$WR[1],
    hosp.first.primary.WR = hosp.first.primary$WR[1],
    order.diff.hosp.minus.death = order.diff,
    order.changed.primary = as.integer(obs$max.order.primary$order[1] == "hospitalization_first"),
    maxOrderWR_primary = test.out$T.obs.max["maxOrderWR_primary"],
    selected.order.primary = obs$max.order.primary$order[1],
    selected.order.p.primary = obs$max.order.primary$p[1],
    order.primary.tie.count = o.primary.row$tie.count[1],
    order.primary.tie.pr = o.primary.row$tie.pr[1],
    max.pvalue.maxOrderWR_primary = test.out$p.value.max["maxOrderWR_primary"],
    fixed.pvalue.maxOrderWR_primary = test.out$p.value.fixed["fixedOrder_primary"],

    maxOrderWR_full = test.out$T.obs.max["maxOrderWR_full"],
    selected.order.full = obs$max.order.full$order[1],
    selected.order.p.full = obs$max.order.full$p[1],
    order.full.tie.count = o.full.row$tie.count[1],
    order.full.tie.pr = o.full.row$tie.pr[1],
    selected.order.full.is.low.p = as.integer(obs$max.order.full$p[1] < 0.5),
    max.pvalue.maxOrderWR_full = test.out$p.value.max["maxOrderWR_full"],
    fixed.pvalue.maxOrderWR_full = test.out$p.value.fixed["fixedOrder_full"],

    maxWRt = test.out$T.obs.max["maxWRt"],
    selected.t.years = obs$max.threshold$t[1],
    selected.t.months = obs$max.threshold$t.months[1],
    threshold.tie.count = obs$max.threshold$tie.count[1],
    threshold.tie.pr = obs$max.threshold$pr.tie[1],
    max.pvalue.maxWRt = test.out$p.value.max["maxWRt"],
    fixed.pvalue.maxWRt = test.out$p.value.fixed["fixedWRt"],

    true.hierarchical.tie.count = obs$counts$true.tie.pairs,
    true.hierarchical.tie.pr = obs$counts$true.tie.pr,
    total.pairs = obs$counts$total.pairs,

    logrank.death.statistic = logrank.results$statistic[logrank.results$method == "Log-rank death endpoint"],
    logrank.death.p = logrank.results$p.value[logrank.results$method == "Log-rank death endpoint"],
    logrank.composite.statistic = logrank.results$statistic[logrank.results$method == "Log-rank composite endpoint"],
    logrank.composite.p = logrank.results$p.value[logrank.results$method == "Log-rank composite endpoint"],

    death.event.rate.control = get_comp_arm_value(comp.stats, 0, "death.event.rate"),
    death.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "death.event.rate"),
    composite.event.rate.control = get_comp_arm_value(comp.stats, 0, "composite.event.rate"),
    composite.event.rate.treatment = get_comp_arm_value(comp.stats, 1, "composite.event.rate"),
    mean.num.hosp.control = get_comp_arm_value(comp.stats, 0, "mean.num.hosp"),
    mean.num.hosp.treatment = get_comp_arm_value(comp.stats, 1, "mean.num.hosp"),

    threshold.grid.months = paste(sprintf("%.1f", threshold.info$threshold.table$t.months), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

plot.realdata.pathway.comparison <- function(pathway.tab, outdir, prefix, title.label = prefix) {
  all.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt", "logrank_death", "logrank_composite")
  wr.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full", "maxOrderWR_primary", "maxOrderWR_full", "maxWRt")
  weighted.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxWRp_low", "maxWRp_full")
  logrank.methods <- c("logrank_death", "logrank_composite")
  threshold.methods <- c("maxWRt")
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_all_methods_with_logrank.png"), "All method p-values", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_WR_pathways_only_no_logrank.png"), "WR pathway p-values", "permutation_p_value", wr.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_pvalue_traditional_and_weighted_only.png"), "Traditional/order and weighted WR p-values", "permutation_p_value", weighted.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_logrank_pvalues_separate.png"), "Log-rank p-values", "permutation_p_value", logrank.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_threshold_pvalue_separate.png"), "Threshold pathway p-value", "permutation_p_value", threshold.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_permutation_pvalue_comparison.png"), "Permutation p-value comparison", "permutation_p_value", all.methods)
  plot.pvalue.bar(pathway.tab, outdir, paste0(prefix, "_fixed_selected_pvalue_comparison.png"), "Fixed-selected p-value comparison", "selected_parameter_p_value", wr.methods)
}

plot.global.realdata.pvalues <- function(pathway.table, outdir) {
  if (nrow(pathway.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt", "logrank_death", "logrank_composite")
  tab <- pathway.table[pathway.table$method %in% key.methods, , drop = FALSE]
  datasets <- unique(tab$dataset_label)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(datasets), dimnames = list(method.display.label(key.methods), datasets))
  for (i in seq_len(nrow(tab))) mat[method.display.label(tab$method[i]), tab$dataset_label[i]] <- tab$permutation_p_value[i]
  save_png(outdir, "GLOBAL_permutation_pvalue_comparison_by_dataset.png", width = 3000, height = 1450)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2, ylab = "Permutation p-value", main = "Permutation p-values by dataset", legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.68), cex.names = 0.82)
  abline(h = ALPHA, lty = 2, lwd = 2); par(old.par); dev.off()
}

plot.global.realdata.fixed.pvalues <- function(pathway.table, outdir) {
  if (nrow(pathway.table) == 0) return(invisible(NULL))
  key.methods <- c("ordinaryWR", "traditionalWR_hosp_first", "traditionalOrderWR", "maxWRp_primary", "maxOrderWR_primary", "maxWRt")
  tab <- pathway.table[pathway.table$method %in% key.methods, , drop = FALSE]
  datasets <- unique(tab$dataset_label)
  mat <- matrix(NA_real_, nrow = length(key.methods), ncol = length(datasets), dimnames = list(method.display.label(key.methods), datasets))
  for (i in seq_len(nrow(tab))) mat[method.display.label(tab$method[i]), tab$dataset_label[i]] <- tab$selected_parameter_p_value[i]
  save_png(outdir, "GLOBAL_fixed_selected_pvalue_comparison_by_dataset.png", width = 2800, height = 1400)
  old.par <- par(no.readonly = TRUE); par(mar = c(8.8, 5.2, 3.2, 1.2))
  barplot(mat, beside = TRUE, ylim = c(0, 1), las = 2, ylab = "Fixed-selected p-value", main = "Fixed-selected p-values by dataset", legend.text = TRUE, args.legend = list(x = "topright", bty = "n", cex = 0.72), cex.names = 0.84)
  abline(h = ALPHA, lty = 2, lwd = 2); par(old.par); dev.off()
}


#run block
if (AUTO_RUN) {
  source.types <- if ("source_type" %in% names(REAL_DATASETS)) {
    tolower(as.character(REAL_DATASETS$source_type))
  } else {
    rep("file", nrow(REAL_DATASETS))
  }
  file.rows <- which(!(source.types %in% c("package", "r_package", "data")))
  missing.files <- character(0)
  if (length(file.rows) > 0 && "subject_file" %in% names(REAL_DATASETS)) {
    candidate.files <- as.character(REAL_DATASETS$subject_file[file.rows])
    missing.files <- candidate.files[candidate.files == "" | !file.exists(candidate.files)]
  }

  if (length(missing.files) > 0) {
    write.realdata.templates(OUTDIR)
    cat("\nReal-data script loaded, but subject_file path(s) were not found.\n")
    cat("Edit REAL_DATASETS at the top of the script, set source_type = 'package', or set REAL_DATASETS before sourcing.\n")
    cat("Template files were written to:", OUTDIR, "\n")
    cat("Missing subject_file values:\n")
    print(missing.files)
    cat("\nAfter editing, run:\n")
    cat("out <- run.all.real.datasets()\n")
  } else {
    out <- run.all.real.datasets(
      dataset.registry = REAL_DATASETS,
      B = B_REAL,
      batch.size = BATCH_SIZE_REAL,
      outdir = OUTDIR,
      master.seed = MASTER_SEED
    )
  }
}
