# Overview

This repository provides an analytical pipeline for calibrated projections of SARS-CoV-2 transmission and COVID19 burden for specified regions.

# Analysis Engineering

This analysis is designed to be run from the command line, specifically in an high performance computing setup.

The core model is [covidm](https://github.com/nicholasdavies/covidm). This repository does not manage installation of that library and it's dependencies.

This repository does define all other dependencies (*i.e.*, R packages) and inputs (*e.g.* case and death data), and provides targets for their download and cleaning.

The analyses defined herein are written in an embarassingly parallel fashion. Inputs and outputs are segregated by file paths; within each source (inputs) and sink (outputs) subfolder, there are identical files. Generic inputs and outputs are defined at the source and sink folder level, while locale-specific ones are defined in subfolders.

The `Makefile` defines the relationships between inputs and outputs, in terms of various scripts in the `analysis` folder. Assorted top level variables, *e.g.* path to the source and sink folders, can be defined in a `local.makefile`; these variables are distinguishable by the `?=` assignment in the `Makefile`. The `example.makefile` shows the structure for a `local.makefile`, and `make local.makefile` will create a copy for you to edit.

In general, each target represents a step in the analysis, typically accomplished by executing a single script with dependencies and the target as arguments, *i.e.*:

```
%/output.last: somescript.R %/someinput %/moreinput
  Rscript $^ $@
```

A `somescript.R` generally begins:

```
suppressPackageWarnings({
#' require assorted packages
})

.debug <- "PATTERN" #' an example % match from Makefile

#'
.args <- if (interactive()) sprintf(c(
  "%s/someinput", "%s/moreinput",
  "output.last"
), .debug) else commandArgs(trailingOnly = FALSE)

#' assign .args into sensible variables
```

## Naming

Scripts follow these naming conventions:
 - `get_...`: gets data from a remote source; requires open internet access
 - `gen_...`: synthesizes local data into new local data; only re-organization, no analysis
 - `est_...`: imputes some result from local data; maybe either deterministic or stochastic
 - `sim_...`: simulates (projects) some process using assumptions in local data
 - `dig_...`: digests simulation products into summaries
 - `fig_...`: visualizes some element of the analysis

# Setup

`make setup` will attempt to install dependencies and gather inputs.

# Data Sources

## Reported Cases and Deaths

We use the European Centre for Disease Prevention and Control (ECDC) [time series for cases and deaths](https://opendata.ecdc.europa.eu/covid19/casedistribution/csv). Aside from superficial organizational changes, we perform one cleaning step to address negative case counts. We assume negative reports represent a correction to earlier data. We distribute the negative cases proportionally into all previous daily case counts.

# Analysis Steps

 - determine when pre-/post- [and optionally, modification] intervention periods occur
 - calculate static Rt for these periods
 - use the static Rt as an approximate R
  * this implies pre-R == R0 with no interventions
  * post-R == R0 with interventions in place
  * optional modification-R == Reff (i.e. with both susceptible depletion and some intervention)
 - using R values, fit assorted intervention scenarios
  * for the modification-R, this means fitting the pre-/post-R first, then projecting (to deplete susceptibles)
 - with intervention fits, run projections