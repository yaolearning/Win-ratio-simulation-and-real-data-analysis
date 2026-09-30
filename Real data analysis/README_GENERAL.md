# DIG 3-Endpoint WR/WO Analysis

This folder contains the current generalized GitHub-ready implementation of the DIG three-endpoint analysis.

## Current analysis settings

- Three three-endpoint sets, all including death
- WR and WO
- Optimization over endpoint ordering, weighting, thresholds, and joint combinations
- Death-only thresholds: 0, 7, 14, 21 days
- Adaptive statistic is re-maximized inside every permutation
- Default B = 10,000 permutations
- Default OpenMP threads = 8
- Two-sided results are exported separately for manuscript use

## Files

- `RUN_DIG_3SETS_THRESHOLD_WR_WO_OPENMP8_B10000_GENERAL.R`
- `dig3_3sets_threshold_openmp_engine_GENERAL.cpp`

For easiest use, rename the C++ file to:

`dig3_3sets_threshold_openmp_engine.cpp`

or update `CPP_FILE` in the R script.

## Endpoint sets

1. Death + first all-cause hospitalization + recurrent hospitalization
2. Death + first all-cause hospitalization + SVA hospitalization
3. Death + recurrent hospitalization + SVA hospitalization

## Notes

The script loads `DIGdata` from the `asympTest` R package. Output paths are relative to the current working directory so the code is portable across Windows, macOS, and Linux.

Manuscript-facing three-endpoint results should use the two-sided outputs.
