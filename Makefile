
default: setup generate

#################### SETUP ###################################################

export
# export all variable definitions except those tagged for "unexport"
# provides all these top-level variables to sub-make invocations

-include local.makefile
# if present, include local.makefile to provide alternative definitions of
# elements assigned via `?=` (e.g. default variables, paths)

PROJRT := $(shell pwd)
# root filesystem location for inputs & outputs
DATART ?= ${PROJRT}/analysis
# example local.makefile overrides this to point to a Dropbox folder

# directories for the various tasks in the analysis
# grouped by kind, not order, though order is generally:
# inputs -> generation -> estimation -> figures
SOURCE := ${DATART}/ins
GEND   := ${DATART}/gen
ESTD   := ${DATART}/est
FIGS   := ${DATART}/fig

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${FIGS} ${GEND} ${ESTD}

africaisos.txt: gen_isos.R ${SOURCE}/epi_data.rds
	${R}

ISOS = NGA GHA PAK
# ISOS ?= $(shell cat africaisos.txt) PAK

# provides non-analysis support
# n.b. this is reloaded in sub-make files
include support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# performs compile step for covidm repo
# writes .install if it succeeds otherwise errors uninformatively
.install: get_install.R rpack.txt | ${COVIDM}
	${Rpipe}


SM = ${MAKE} -wC $(notdir $|)

${SOURCE}/%: | ${SOURCE}
	${SM} $@

setup: | ${SOURCE}
	${SM}

${GEND}/%: | ${GEND}
	${SM} $@

generate: | ${GEND}
	${SM}

${ESTD}/%: | ${ESTD}
	${SM} $@

estimate: | ${ESTD}
	${SM}

${FIGS}/%: | ${FIGS}
	${SM} $@

figures: | ${FIGS}
	${SM}

