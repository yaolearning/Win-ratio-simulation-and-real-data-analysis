weight_control <- function(mode = c("vertices", "grid"),
                           step = 0.05,
                           constraint = c("ordered", "simplex")) {
  mode <- match.arg(mode)
  constraint <- match.arg(constraint)
  step <- as.numeric(step)[1]
  if (!is.finite(step) || step <= 0 || step > 1) stop("step must be in (0, 1].")

  structure(
    list(
      mode = mode,
      step = step,
      constraint = constraint
    ),
    class = c("win_weight_control", "list")
  )
}

.wr_make_weight_grid <- function(m,
                                 control = weight_control()) {
  m <- as.integer(m)[1]
  if (!(m %in% c(2L, 3L))) stop("Only 2 or 3 endpoints are supported.")

  mode <- control$mode
  step <- control$step
  constraint <- control$constraint

  if (m == 2L) {
    if (mode == "vertices") {
      if (constraint == "ordered") {
        w <- rbind(c(1, 0), c(0.5, 0.5))
      } else {
        w <- rbind(c(1, 0), c(0, 1), c(0.5, 0.5))
      }
    } else {
      p1 <- seq(0, 1, by = step)
      w <- cbind(p1, 1 - p1)
      if (constraint == "ordered") {
        w <- w[w[, 1] >= w[, 2] - 1e-9, , drop = FALSE]
      }
    }
  } else {
    if (mode == "vertices") {
      if (constraint == "ordered") {
        w <- rbind(
          c(1, 0, 0),
          c(0.5, 0.5, 0),
          c(1 / 3, 1 / 3, 1 / 3)
        )
      } else {
        w <- rbind(
          c(1, 0, 0),
          c(0, 1, 0),
          c(0, 0, 1),
          c(1 / 3, 1 / 3, 1 / 3)
        )
      }
    } else {
      vals <- seq(0, 1, by = step)
      rows <- list()
      idx <- 1L

      for (p1 in vals) {
        for (p2 in vals) {
          p3 <- 1 - p1 - p2
          if (p3 < -1e-9) next
          if (p3 < 0) p3 <- 0

          ww <- c(p1, p2, p3)
          if (abs(sum(ww) - 1) > 1e-8) next

          if (constraint == "ordered") {
            if (!(ww[1] >= ww[2] - 1e-9 && ww[2] >= ww[3] - 1e-9)) next
          }

          rows[[idx]] <- ww
          idx <- idx + 1L
        }
      }

      if (length(rows) == 0) stop("No valid weight combinations were generated.")
      w <- do.call(rbind, rows)
    }
  }

  w <- unique(round(w, 8))
  colnames(w) <- paste0("p", seq_len(ncol(w)))
  as.data.frame(w)
}

.wr_equal_weights <- function(m) {
  rep(1 / as.integer(m), as.integer(m))
}

.wr_natural_order <- function(m) {
  seq_len(as.integer(m))
}

.wr_all_orders <- function(m) {
  .wr_all_permutations(seq_len(as.integer(m)))
}
