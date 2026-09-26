.wr_safe_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

.wr_normalize_time_unit <- function(unit) {
  unit <- tolower(trimws(as.character(unit)[1]))
  if (unit %in% c("year", "years", "yr", "yrs")) return("years")
  if (unit %in% c("month", "months", "mo", "mos")) return("months")
  if (unit %in% c("day", "days", "d")) return("days")
  stop("Unsupported time unit: ", unit)
}

.wr_time_to_years <- function(x, unit = "years") {
  unit <- .wr_normalize_time_unit(unit)
  x <- .wr_safe_num(x)
  if (unit == "years") return(x)
  if (unit == "months") return(x / 12)
  x / 365.25
}

.wr_years_to_unit <- function(x, unit = "years") {
  unit <- .wr_normalize_time_unit(unit)
  x <- as.numeric(x)
  if (unit == "years") return(x)
  if (unit == "months") return(x * 12)
  x * 365.25
}

.wr_first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[1]
}

.wr_ratio_safe <- function(num, den) {
  if (is.na(num) || is.na(den)) return(NA_real_)
  if (den == 0 && num == 0) return(1)
  if (den == 0 && num > 0) return(Inf)
  if (num == 0 && den > 0) return(0)
  num / den
}

.wr_log_ratio <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  out[is.infinite(x) & x > 0] <- Inf
  ok <- is.finite(x) & x > 0
  out[ok] <- log(x[ok])
  out
}

.wr_abslog_safe <- function(x) {
  abs(.wr_log_ratio(x))
}

.wr_right_tail_perm_p <- function(perm_stat, obs_stat) {
  perm_stat <- as.numeric(perm_stat)
  obs_stat <- as.numeric(obs_stat)[1]
  ok <- is.finite(perm_stat) | is.infinite(perm_stat)
  if (is.na(obs_stat) || sum(ok) == 0) return(NA_real_)
  (1 + sum(perm_stat[ok] >= obs_stat)) / (sum(ok) + 1)
}

.wr_safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

.wr_safe_sd <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(NA_real_)
  stats::sd(x)
}

.wr_safe_median <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

.wr_mode_string <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

.wr_format_order <- function(order_vec) {
  paste(as.integer(order_vec), collapse = "->")
}

.wr_parse_order <- function(order_key) {
  as.integer(strsplit(as.character(order_key), "->", fixed = TRUE)[[1]])
}

.wr_all_permutations <- function(x) {
  x <- as.integer(x)
  if (length(x) == 1) return(list(x))
  out <- list()
  idx <- 1L
  for (i in seq_along(x)) {
    rest <- x[-i]
    sub <- .wr_all_permutations(rest)
    for (s in sub) {
      out[[idx]] <- c(x[i], s)
      idx <- idx + 1L
    }
  }
  out
}

.wr_rbind_fill <- function(...) {
  xs <- list(...)
  flatten <- function(z) {
    out <- list()
    for (item in z) {
      if (is.null(item)) {
        next
      } else if (is.data.frame(item)) {
        out <- c(out, list(item))
      } else if (is.list(item)) {
        out <- c(out, flatten(item))
      }
    }
    out
  }
  xs <- flatten(xs)
  xs <- xs[vapply(xs, function(x) nrow(x) > 0, logical(1))]
  if (length(xs) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}

.wr_clean_file_name <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x)
}

.wr_make_weight_string <- function(w) {
  paste0("(", paste(sprintf("%.3f", as.numeric(w)), collapse = ","), ")")
}

.wr_make_threshold_string <- function(thr_vec) {
  paste0("(", paste(sprintf("%.4f", as.numeric(thr_vec)), collapse = ","), ")")
}
