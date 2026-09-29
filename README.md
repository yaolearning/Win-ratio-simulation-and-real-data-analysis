# MaxWin

**MaxWin** is an R package for Win Ratio (WR) and Win Odds (WO) analyses of hierarchical composite endpoints. It supports standard and adaptive methods based on endpoint ordering, endpoint weighting, and time-endpoint thresholds, together with treatment-label permutation inference, recurrent-event comparisons, real-data analysis, and simulation studies.

## Installation

MaxWin can be installed directly from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "yaolearning/Win-ratio-simulation-and-real-data-analysis"
)

library(MaxWin)
```

Check the installed version with:

```r
packageVersion("MaxWin")
```

For the current release:

```text
0.1.0
```

Because MaxWin contains C++ code through Rcpp, Windows users installing from source may need a version of Rtools compatible with their R installation.

## Main features

MaxWin currently provides:

Win Ratio and Win Odds analysis.
Two- and three-endpoint hierarchical composite outcomes.
Time-to-event, count, recurrent-count, binary, and continuous endpoints.
Recurrent-event count comparisons using pairwise common follow-up.
Fixed and adaptive endpoint ordering.
Endpoint weighting.
Time-endpoint threshold methods.
Combined order, weight, and threshold selection.
One-sided inference using the raw WR or WO statistic.
Two-sided inference using `abs(log(WR))` or `abs(log(WO))`.
Treatment-label permutation testing.
Full adaptive re-selection inside each permutation.
Fixed-selected permutation p-values.
Tie counts and permutation-average tie summaries.
Log-rank comparator analyses for time-to-event endpoints.
Real-data tables, CSV/RDS export, Kaplan-Meier plots, p-value plots, and tie plots.
Simulation data generation.
Scenario-based data generation.
Fifteen predefined simulation scenarios.
Repeated simulation studies with empirical rejection/power summaries and gain plots.

## Quick start: generate a simulated dataset

```r
library(MaxWin)

dat <- generate_win_dataset(
  n_control = 50,
  n_treatment = 50,
  seed = 2026
)

dat
```

The returned object contains:

```r
dat$subjects
dat$recurrent_events
dat$analysis_data
dat$settings
```

## Define hierarchical endpoints

For a two-endpoint death and recurrent-hospitalization analysis:

```r
endpoints <- list(
  endpoint_time(
    name = "Death",
    time = "FUTIME",
    event = "CNSR",
    unit = "years"
  ),
  endpoint_count(
    name = "Hospitalization",
    count = "NUMHOSP",
    comparison = "pairwise_common_followup",
    recurrent_time = "HOSPTIME",
    recurrent_id = "SUBJID",
    followup = "FUTIME",
    unit = "years"
  )
)
```

For recurrent counts, the comparison can be restricted to the pair-specific common follow-up:

\[
F_{ij} = \min(F_i, F_j).
\]

Only recurrent events occurring by \(F_{ij}\) are counted for that treatment-control pair.

## Real-data analysis

If subject-level and recurrent-event data are stored separately:

```r
input <- win_data(
  data = subject_data,
  recurrent_data = recurrent_event_data
)
```

Run the analysis:

```r
fit <- win_analysis(
  data = input,
  endpoints = endpoints,
  id = "SUBJID",
  treatment = "ARM",
  permutation = permutation_control(
    enabled = TRUE,
    B = 500,
    seed = 2026
  )
)
```

View the main results:

```r
summary(fit)
```

or:

```r
as.data.frame(fit)
```

## Adaptive WR/WO methods

MaxWin can evaluate methods based on:

the original endpoint order,
fixed alternative orders,
adaptive order selection,
adaptive endpoint weighting,
adaptive time thresholds,
order + weight,
weight + threshold,
order + threshold,
full order + weight + threshold selection.

Threshold comparisons for time-to-event endpoints use the strict rule:

\[
|\Delta T| > t.
\]

For adaptive methods, the complete selection procedure is repeated inside each treatment-label permutation.

## One-sided and two-sided inference

For one-sided benefit-oriented inference, larger WR or WO values favor treatment.

For two-sided inference, MaxWin uses distance from the null value 1:

\[
|\log(WR)|
\]

and

\[
|\log(WO)|.
\]

## Plotting

Examples:

```r
plot(
  fit,
  type = "pvalue",
  measure = "WR",
  side = "two"
)
```

```r
plot(
  fit,
  type = "km",
  endpoint = 1
)
```

```r
plot(
  fit,
  type = "ties",
  measure = "WR",
  side = "one"
)
```

## Export real-data results

```r
write_win_results(
  fit,
  output_dir = "MaxWin_results",
  prefix = "analysis"
)
```

The output can include:

method results,
observed selections,
all candidate definitions,
observed candidate statistics,
permutation selections,
threshold and weight grids,
endpoint summaries,
tie summaries,
log-rank results,
settings,
figures,
an RDS object containing the complete analysis,
an output manifest.

## Default simulation scenarios

The package includes 15 predefined scenarios:

```r
scenarios <- default_win_scenarios()

nrow(scenarios)
head(scenarios)
```

These include null, benefit, harm, discordant endpoint effects, censoring sensitivity, follow-up sensitivity, and scenarios emphasizing non-terminal outcomes.

## Generate one dataset from a scenario

```r
scenario_dat <- generate_win_scenario_dataset(
  scenario = 8,
  sim_index = 1,
  seed = 2026
)
```

A scenario can be selected by index, scenario ID, or a one-row custom scenario data frame.

## Run a simulation study

A small example:

```r
sim <- run_win_simulation(
  scenarios = 1:3,
  nsim = 20,
  B = 50,
  seed = 2026,
  output_dir = "simulation_output"
)
```

For a full simulation study, increase `nsim` and `B` as needed. Large values can be computationally intensive.

View the empirical rejection/power summary:

```r
summary(sim)
```

Simulation plots include:

```r
plot(
  sim,
  type = "power",
  measure = "WR",
  side = "two"
)
```

```r
plot(
  sim,
  type = "gain",
  measure = "WR",
  side = "two"
)
```

```r
plot(
  sim,
  type = "scenario_gain",
  scenario = 5,
  measure = "WR",
  side = "two"
)
```

## Help

After installation:

```r
?win_analysis
?endpoint_time
?endpoint_count
?generate_win_dataset
?generate_win_scenario_dataset
?run_win_simulation
```

## Development status

MaxWin version 0.1.0 has been locally validated with:

```text
R CMD check: 0 errors, 0 warnings, 0 notes
```

The package includes unit tests for data generation, scenarios, recurrent-event common-follow-up comparisons, permutation reproducibility, output generation, and core package functionality.

## Repository

Source code and development history:

```text
https://github.com/yaolearning/Win-ratio-simulation-and-real-data-analysis
```
## Contact

For questions, feedback, or collaboration inquiries:

Contact: yy5933@nyu.edu
Alternative contact: yaoyifan64@gmail.com