# Overview

This repository provides an analytical pipeline for calibrated projections of SARS-CoV-2 transmission and COVID19 burden for specified regions.

# Analysis Engineering

This analysis is designed to be run from the command line, specifically in an high performance computing setup, but with allowances for development on a personal machine.

The core model is [covidm](https://github.com/nicholasdavies/covidm). There is a convenience target for cloning that repository, but this does not manage installation of covidm dependencies. The `.install` target will attempt to compile covidm; subsequent uses of covidm will *not* recompile that library.

This repository does define other dependencies (*i.e.*, R packages) and inputs (*e.g.* case and death data), and provides targets for their download and cleaning.

This analyses is written in "map-reduce" paradigm, with the parallel steps defined by different ISO country codes. Inputs and outputs are segregated by file paths; within each source (inputs) and sink (outputs) subfolder, there are identical files. Generic inputs and outputs are defined at the source and sink folder level, while locale-specific ones are defined in subfolders or files named by ISO code.

The `Makefile` defines the relationships between inputs and outputs, in terms of various scripts. Assorted top level variables, *e.g.* path to the source and sink folders, can be defined in a `local.makefile`; these variables are distinguishable by the `?=` assignment in the `Makefile`. The `example.makefile` shows the structure for a `local.makefile`, and `make local.makefile` will create a copy for you to edit.

In general, each target represents a step in the analysis, typically accomplished by executing a single script with dependencies and the target as arguments, *e.g.*:

```
%/output.last: somescript.R %/someinput %/moreinput
  Rscript $^ $@
```

A `somescript.R` generally follows this pattern:

```
suppressPackageWarnings({
#' require assorted packages
})

.debug <- c("DATADIR", "PATTERN") #' an example % match from Makefile

#' this structure will use `.debug` to create the relevant arguments
#' when working interactively with the script
.args <- if (interactive()) sprintf(c(
  "%s/input/population/%s.csv", "%s/output/intermediate_result/%s.rds",
  "%s/output/next_result_step/%s.rds"
), .debug[1], .debug[2]) else commandArgs(trailingOnly = FALSE)

#' assign `.args` into sensible variables
#' note that typically the target is the final item in `.args`
tarfile <- tail(.args, 1)

#' PERFORM VARIOUS CLEANING / ANALYSIS / PLOTTING / ETC

saveRDS(result, tarfile)
#' or save, qsave, ggsave, etc as appropriate to the script target
```

## Naming

Scripts follow these naming conventions:

 - `get_...`: gets data from a remote source; requires open internet access; typically puts products in SOURCES (*i.e.* inputs) directory
 - `gen_...`: synthesizes local data into new local data; only re-organization, no analysis; typically puts products in SOURCES (*i.e.* inputs) directory, though to SINKS (*i.e.* outputs) if it's about reorganizing structure of outputs
 - `est_...`: imputes some result from local data; maybe either deterministic or stochastic; products to SINKS
 - `sim_...`: simulates (projects) some process using assumptions in local data; products in SINKS
 - `fig_...`: visualizes some element of the analysis; products in SINKS

# Setup

`make .install` will attempt to install dependencies. `make allinputs` will attempt to get raw inputs.

# Data Sources

## Reported Cases and Deaths

We use the JHU Center for Systems Science and Engineering (CSSE) time series for cases and deaths [(repository link)](https://github.com/CSSEGISandData/COVID-19). Aside from superficial organizational changes, we perform one cleaning step to address negative case counts (if necessary; the process emits a warning message when this occurs). We assume negative reports represent a correction to earlier data, thus we distribute negative incidence into previous daily counts. See `get_epi_data.R` for specifics.

## Sequence Data

We obtain sequence metadata from the (Network for Genomics Surveillance in South Africa)[https://www.krisp.org.za/ngs-sa/index.php], via [NextStrain](https://nextstrain.org/groups/ngs-sa/COVID19-ZA-2021.01.18).

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