# Thesis code — Comparing UVB and ADR audit sampling methods in EU subsidy audits

Andre Blokland (s4113675), MSc Statistics and Data Science, Leiden University.

## How to reproduce the results

Set the R working directory to this folder, then run the scripts in this order:

1. `create_dbase_and_templates.R` — reads `BA23-25.xlsx` and writes the `.rds`
   data objects used by the simulation.
2. `run_simulation_parallel.R` — runs the Monte Carlo simulation and writes
   `pc_results_<n>sim_<timestamp>.rds`.
3. `script_section5_results_analysis_clean.R` — produces the tables and
   figures for Chapter 5.
4. `analyse_hypotheses.R` — evaluates hypotheses H1–H3.

Each script has a header explaining what it does and what it depends on.

## R packages

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "patchwork",
  "readxl", "flextable", "officer"
))
```
