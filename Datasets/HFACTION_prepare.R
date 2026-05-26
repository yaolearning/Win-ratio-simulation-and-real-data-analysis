# =====================================================================
# HF-ACTION dataset preparation for adaptive win ratio analysis
# =====================================================================
# Purpose
#   Load hfaction_cpx9 from the WR package and reshape it into the full
#   data structure expected by the adaptive win ratio pipeline. Save the
#   COMPLETE ds.hf object (patient-level table + per-patient
#   hospitalization timestamps) as an .RData file so the analysis script
#   can just load() it without redoing any reshape work.
#
# Outputs (to ./hfaction_prepared/)
#   hfaction_reshaped.RData  -- complete ds.hf list:
#                                 ds.hf$table.output     data frame
#                                 ds.hf$hosp.times.list  per-patient hosp times
#                               (use this in the analysis script)
#
#   hfaction_reshaped.csv    -- patient-level table only (for inspection
#                               in Excel; cannot hold hosp times because
#                               they're a variable-length list per row)
#
# Source
#   hfaction_cpx9 from the WR package on CRAN.
#   Original trial: HF-ACTION (NCT00047437, O'Connor et al. JAMA 2009).
#   Subgroup: 426 non-ischemic patients with baseline CPX <= 9 minutes.
# =====================================================================


# ---------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

install_if_missing("WR")
library(WR)


# ---------------------------------------------------------------------
# 1. Load hfaction_cpx9 (long format)
# ---------------------------------------------------------------------

data(hfaction_cpx9, package = "WR")
hf.long <- hfaction_cpx9

cat("\n===== HF-ACTION cpx9 dataset loaded =====\n")
cat("Dimensions (long format):", paste(dim(hf.long), collapse = " x "), "\n")
cat("Variables:", paste(names(hf.long), collapse = ", "), "\n\n")

cat("Status codes: 0 = censored, 1 = death, 2 = hospitalization\n")
cat("Status frequencies (event rows):\n")
print(table(hf.long$status))
cat("\nUnique patients:", length(unique(hf.long$patid)), "\n")


# ---------------------------------------------------------------------
# 2. Reshape long -> ds list (table.output + hosp.times.list)
# ---------------------------------------------------------------------
# For each patient:
#   - max(time) is the follow-up endpoint (terminal row is status==0
#     for censored or status==1 for death; both are terminal).
#   - DEATH (CNSR) = 1 if any row has status == 1.
#   - Hospitalization times: ALL rows with status == 2 (exact times,
#     sorted ascending) -> converted to years -> stored as GAPS
#     (deltas), so prepare.ds.fast() in the analysis script can
#     cumsum() them back to absolute times.
# ---------------------------------------------------------------------

reshape.hfaction.to.ds <- function(hf.long, max.followup.years = NULL) {
  hf.by.patient <- split(hf.long, hf.long$patid)
  n <- length(hf.by.patient)
  
  arm <- integer(n)
  fu.years <- numeric(n)
  death <- integer(n)
  numhosp <- integer(n)
  age60 <- integer(n)
  patid.char <- character(n)
  hosp.times.list <- vector("list", n)
  
  for (i in seq_len(n)) {
    rows <- hf.by.patient[[i]]
    patid.char[i] <- as.character(rows$patid[1])
    arm[i] <- as.integer(rows$trt_ab[1])
    age60[i] <- as.integer(rows$age60[1])
    
    # Follow-up endpoint: months -> years.
    fu.years[i] <- max(rows$time) / 12
    
    # Death indicator.
    death[i] <- as.integer(any(rows$status == 1))
    
    # Exact hospitalization times in years, sorted ascending.
    hosp.y <- sort(rows$time[rows$status == 2]) / 12
    
    # Optional follow-up cap.
    if (!is.null(max.followup.years)) {
      if (fu.years[i] > max.followup.years) {
        fu.years[i] <- max.followup.years
        death[i] <- 0
      }
      hosp.y <- hosp.y[hosp.y <= max.followup.years]
    }
    
    numhosp[i] <- length(hosp.y)
    
    # Store as GAPS (deltas) so prepare.ds.fast() can cumsum() to abs times.
    if (length(hosp.y) == 0) {
      hosp.times.list[[i]] <- numeric(0)
    } else {
      hosp.times.list[[i]] <- diff(c(0, hosp.y))
    }
  }
  
  table.output <- data.frame(
    SUBJID   = seq_len(n),
    ARM      = arm,                   # 0 = usual care, 1 = exercise training
    FUTIME   = fu.years,
    CNSR     = death,                 # 1 = died, 0 = censored/alive
    SURVTIME = fu.years,
    CNSRTIME = fu.years,
    FREQHOSP = NA_real_,
    NUMHOSP  = numhosp,
    AGE60    = age60,
    PATID    = patid.char,
    stringsAsFactors = FALSE
  )
  
  list(table.output = table.output, hosp.times.list = hosp.times.list)
}

ds.hf <- reshape.hfaction.to.ds(hf.long)


# ---------------------------------------------------------------------
# 3. Sanity checks / summary
# ---------------------------------------------------------------------

tab <- ds.hf$table.output

cat("\n===== Reshaped HF-ACTION dataset =====\n")
cat("Number of patients:", nrow(tab), "\n\n")

cat("Treatment counts, ARM: 0 = usual care, 1 = exercise training\n")
print(table(tab$ARM))

cat("\nDeath counts by treatment arm (CNSR = 1 means death):\n")
print(table(tab$ARM, tab$CNSR))

cat("\nHospitalization count (NUMHOSP) summary:\n")
print(summary(tab$NUMHOSP))

cat("\nFollow-up time in years (FUTIME) summary:\n")
print(summary(tab$FUTIME))

cat("\nPer-arm summary:\n")
arm.summary <- data.frame(
  ARM         = c(0, 1),
  n           = c(sum(tab$ARM == 0), sum(tab$ARM == 1)),
  deaths      = c(sum(tab$CNSR[tab$ARM == 0]), sum(tab$CNSR[tab$ARM == 1])),
  death.pct   = round(c(100 * mean(tab$CNSR[tab$ARM == 0]),
                        100 * mean(tab$CNSR[tab$ARM == 1])), 1),
  total.hosps = c(sum(tab$NUMHOSP[tab$ARM == 0]), sum(tab$NUMHOSP[tab$ARM == 1])),
  mean.hosps  = round(c(mean(tab$NUMHOSP[tab$ARM == 0]),
                        mean(tab$NUMHOSP[tab$ARM == 1])), 2),
  mean.fu.yrs = round(c(mean(tab$FUTIME[tab$ARM == 0]),
                        mean(tab$FUTIME[tab$ARM == 1])), 2)
)
print(arm.summary)

cat("\nFirst 5 patients - hospitalization gaps (years):\n")
for (i in seq_len(min(5, nrow(tab)))) {
  cat(sprintf("  Pt %d (PATID=%s, ARM=%d): NUMHOSP=%d, gaps=[%s]\n",
              i, tab$PATID[i], tab$ARM[i], tab$NUMHOSP[i],
              paste(sprintf("%.3f", ds.hf$hosp.times.list[[i]]), collapse = ", ")))
}


# ---------------------------------------------------------------------
# 4. Save outputs
# ---------------------------------------------------------------------

OUTDIR <- "hfaction_prepared"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Complete ds.hf object (table.output + hosp.times.list) for the
# analysis pipeline. The analysis script loads this with one line:
#   load("hfaction_prepared/hfaction_reshaped.RData")  # provides ds.hf
save(ds.hf, file = file.path(OUTDIR, "hfaction_reshaped.RData"))

# Patient-level table only, for inspection. CSV cannot hold the
# variable-length hosp.times.list per row.
write.csv(tab,
          file      = file.path(OUTDIR, "hfaction_reshaped.csv"),
          row.names = FALSE)

cat("\n===== Saved =====\n")
cat("  ", file.path(OUTDIR, "hfaction_reshaped.RData"),
    "  <-- complete ds.hf (use this for analysis)\n")
cat("  ", file.path(OUTDIR, "hfaction_reshaped.csv"),
    "    <-- patient-level table only (for inspection)\n\n")

cat("To use in analysis script:\n")
cat("  load(\"hfaction_prepared/hfaction_reshaped.RData\")\n")
cat("  # 'ds.hf' is now in your environment with table.output + hosp.times.list\n")