# =====================================================================
# PBC data preparation for adaptive win ratio analysis (order-swap demo)
# =====================================================================
# Purpose
#   Load the Mayo Clinic Primary Biliary Cirrhosis (PBC) trial dataset
#   from the survival package and reshape it so the adaptive WR pipeline
#   can analyze it with the ENDPOINT ORDER PATHWAY (Section 1.3 of the
#   methodology).
#
#   Why PBC for the order-swap demonstration?
#     Both events recorded in PBC are clinically severe:
#       - status == 2  ->  death from liver disease
#       - status == 1  ->  liver transplant
#     A liver transplant is not "obviously less bad" than death:
#       major surgery, lifelong immunosuppression, operative mortality.
#     Reasonable hepatologists could argue either endpoint should be
#     prioritized first. This is exactly the clinical context where
#     T_ord, T^primary_{ord,p}, and T^full_{ord,p} become meaningful.
#
# Mapping onto the existing pipeline schema
#   CNSR     ->  death indicator (status == 2)
#   NUMHOSP  ->  transplant indicator (status == 1), 0 or 1 per patient
#   hosp.times.list[[i]]  ->  transplant time in years if status == 1,
#                             else numeric(0)
#
#   This lets the same downstream code analyze D -> NT (death first)
#   and NT -> D (transplant first) by simply swapping which column is
#   treated as the primary endpoint.
#
# Outputs (to ./pbc_prepared/)
#   pbc_reshaped.RData  -- complete ds.pbc list:
#                            ds.pbc$table.output     data frame
#                            ds.pbc$hosp.times.list  per-patient transplant times
#   pbc_reshaped.csv    -- patient-level table only, for Excel inspection
#
# Source
#   survival::pbc, Mayo Clinic D-penicillamine PBC trial.
#   Therneau & Grambsch (2000), Modeling Survival Data.
#   https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbc.html
# =====================================================================


# ---------------------------------------------------------------------
# 0. Setup
# ---------------------------------------------------------------------

# survival is part of base R, so no install_if_missing wrapper needed.
library(survival)


# ---------------------------------------------------------------------
# 1. Load PBC dataset
# ---------------------------------------------------------------------

data(pbc, package = "survival")
pbc.raw <- pbc

cat("\n===== PBC dataset loaded =====\n")
cat("Total dimensions:", paste(dim(pbc.raw), collapse = " x "), "\n")
cat("Variables:", paste(names(pbc.raw), collapse = ", "), "\n\n")

cat("First 5 rows (selected columns):\n")
print(pbc.raw[1:5, c("id", "time", "status", "trt", "age", "sex", "stage")])

cat("\nOriginal status codes (0 = censored, 1 = transplant, 2 = death):\n")
print(table(pbc.raw$status, useNA = "ifany"))

cat("\nTreatment codes (1 = D-penicillamine, 2 = placebo, NA = not in RCT):\n")
print(table(pbc.raw$trt, useNA = "ifany"))


# ---------------------------------------------------------------------
# 2. Reshape PBC -> ds list (RCT subset only)
# ---------------------------------------------------------------------
# Steps:
#   - Restrict to the 312 patients who were in the randomized trial
#     (those with non-NA trt).
#   - Convert time from days to years (divide by 365.25).
#   - Map status to two indicator columns (death, transplant).
#   - Build the hospitalization times list with transplant time
#     when applicable, empty otherwise.
# ---------------------------------------------------------------------

reshape.pbc.to.ds <- function(pbc.raw, rct.only = TRUE) {

  # Restrict to RCT participants (those with assigned treatment).
  if (rct.only) {
    pbc.raw <- pbc.raw[!is.na(pbc.raw$trt), ]
  }

  n <- nrow(pbc.raw)

  # Treatment indicator (our pipeline convention):
  #   1 = D-penicillamine (active treatment)
  #   0 = placebo (control)
  arm <- as.integer(pbc.raw$trt == 1)

  # Death indicator (status == 2 means death from liver disease).
  death <- as.integer(pbc.raw$status == 2)

  # Transplant indicator (status == 1 means liver transplant).
  transplant <- as.integer(pbc.raw$status == 1)

  # Time to event/censoring, in years (PBC times are in days).
  fu.years <- pbc.raw$time / 365.25

  # Transplant times list. Each patient has either 0 transplants
  # (numeric(0)) or 1 transplant at fu.years[i]. Stored as a single-
  # element vector (gap from 0 = absolute time, so cumsum recovers it).
  hosp.times.list <- vector("list", n)
  for (i in seq_len(n)) {
    if (transplant[i] == 1L) {
      hosp.times.list[[i]] <- fu.years[i]
    } else {
      hosp.times.list[[i]] <- numeric(0)
    }
  }

  # Patient-level wide-format table.
  table.output <- data.frame(
    SUBJID      = seq_len(n),
    ARM         = arm,                  # 0 = placebo, 1 = D-penicillamine
    FUTIME      = fu.years,
    CNSR        = death,                # 1 = died from liver disease, 0 = no
    SURVTIME    = fu.years,
    CNSRTIME    = fu.years,
    FREQHOSP    = NA_real_,
    NUMHOSP     = transplant,           # 1 = transplant, 0 = no
    PBC_STATUS  = pbc.raw$status,       # original 0/1/2 for traceability
    PBC_AGE     = pbc.raw$age,          # age in years at registration
    PBC_SEX     = as.character(pbc.raw$sex),
    PBC_STAGE   = pbc.raw$stage,        # histologic stage 1-4
    PBC_ORIG_ID = pbc.raw$id,
    stringsAsFactors = FALSE
  )

  list(table.output = table.output, hosp.times.list = hosp.times.list)
}

ds.pbc <- reshape.pbc.to.ds(pbc.raw, rct.only = TRUE)


# ---------------------------------------------------------------------
# 3. Sanity checks / summary
# ---------------------------------------------------------------------

tab <- ds.pbc$table.output

cat("\n===== Reshaped PBC dataset (RCT subset) =====\n")
cat("Number of patients (RCT only):", nrow(tab), "\n\n")

cat("Treatment counts, ARM: 0 = placebo, 1 = D-penicillamine\n")
print(table(tab$ARM))

cat("\nDeath counts by treatment arm (CNSR = 1 means died from liver disease):\n")
print(table(tab$ARM, tab$CNSR))

cat("\nTransplant counts by treatment arm (NUMHOSP = 1 means received transplant):\n")
print(table(tab$ARM, tab$NUMHOSP))

cat("\nJoint distribution of original PBC status by treatment arm:\n")
cat("  status: 0 = censored, 1 = transplant, 2 = death\n")
print(table(tab$ARM, tab$PBC_STATUS))

cat("\nFollow-up time in years (FUTIME) summary:\n")
print(summary(tab$FUTIME))

cat("\nPer-arm summary:\n")
arm.summary <- data.frame(
  ARM             = c(0, 1),
  n               = c(sum(tab$ARM == 0), sum(tab$ARM == 1)),
  deaths          = c(sum(tab$CNSR[tab$ARM == 0]),
                      sum(tab$CNSR[tab$ARM == 1])),
  transplants     = c(sum(tab$NUMHOSP[tab$ARM == 0]),
                      sum(tab$NUMHOSP[tab$ARM == 1])),
  censored        = c(sum(tab$PBC_STATUS[tab$ARM == 0] == 0),
                      sum(tab$PBC_STATUS[tab$ARM == 1] == 0)),
  mean_fu_years   = round(c(mean(tab$FUTIME[tab$ARM == 0]),
                            mean(tab$FUTIME[tab$ARM == 1])), 2),
  median_fu_years = round(c(median(tab$FUTIME[tab$ARM == 0]),
                            median(tab$FUTIME[tab$ARM == 1])), 2)
)
print(arm.summary)

cat("\nFirst 5 patients - transplant times (years):\n")
for (i in seq_len(min(5, nrow(tab)))) {
  cat(sprintf("  Pt %d (ARM=%d, status=%d): transplant_time=[%s]\n",
              i, tab$ARM[i], tab$PBC_STATUS[i],
              paste(sprintf("%.3f", ds.pbc$hosp.times.list[[i]]),
                    collapse = ", ")))
}


# ---------------------------------------------------------------------
# 4. How the order-swap pathway uses this data (for reference)
# ---------------------------------------------------------------------

cat("\n===== How to use this for the order-swap analysis =====\n")
cat("\n")
cat("Two endpoints are available, both clinically severe:\n")
cat("  CNSR    = 1  ->  death from liver disease\n")
cat("  NUMHOSP = 1  ->  liver transplant\n")
cat("\n")
cat("D -> NT (death-first, traditional priority):\n")
cat("  Each pair: compare death first; if undecided, compare transplant.\n")
cat("\n")
cat("NT -> D (transplant-first, swapped priority):\n")
cat("  Each pair: compare transplant first; if undecided, compare death.\n")
cat("\n")
cat("T_ord  = max{ WR_o(0.5) : o in O }  where O = {D->NT, NT->D}.\n")
cat("T_ord_p_primary = max{ WR_o(p) : o in O, p in [0.5, 1] }.\n")
cat("\n")


# ---------------------------------------------------------------------
# 5. Save outputs
# ---------------------------------------------------------------------

OUTDIR <- "pbc_prepared"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Complete ds.pbc object for the analysis pipeline.
# The analysis script loads this with:
#   load("pbc_prepared/pbc_reshaped.RData")
save(ds.pbc, file = file.path(OUTDIR, "pbc_reshaped.RData"))

# Patient-level table only, for inspection.
write.csv(tab,
          file      = file.path(OUTDIR, "pbc_reshaped.csv"),
          row.names = FALSE)

cat("===== Saved =====\n")
cat("  ", file.path(OUTDIR, "pbc_reshaped.RData"),
    "  <-- complete ds.pbc (use this for analysis)\n", sep = "")
cat("  ", file.path(OUTDIR, "pbc_reshaped.csv"),
    "    <-- patient-level table only (for inspection)\n\n", sep = "")


# ---------------------------------------------------------------------
# 6. Instructions for teammates
# ---------------------------------------------------------------------

cat("===== How teammates can use the output =====\n")
cat("\n")
cat("Two lines of R (NO survival package needed if you have the RData):\n")
cat("\n")
cat('  load("pbc_prepared/pbc_reshaped.RData")\n')
cat('  # ds.pbc is now available with:\n')
cat('  #   ds.pbc$table.output    (data frame, one row per patient)\n')
cat('  #   ds.pbc$hosp.times.list (per-patient transplant times)\n')
cat("\n")
cat("For Excel inspection: open pbc_reshaped.csv directly.\n\n")
