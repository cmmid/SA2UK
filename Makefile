
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
ACCT ?= TESTTESTTEST
# for HPC submission script generation


# directories for the various tasks in the analysis
# grouped by kind, not order, though order is generally:
# inputs -> generation -> estimation -> figures
SOURCE := ${DATART}/ins
GEND   := ${DATART}/gen
ESTD   := ${DATART}/est
FIGS   := ${DATART}/fig

RAWDATA := ${SOURCE}/epi_data.rds
EPIDATA := ${SOURCE}/adj_data.rds

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${FIGS} ${GEND} ${ESTD}

# override with local.makefile
NSAMPLES ?= 4000

africaisos.txt: gen_isos.R ${SOURCE}/epi_data.rds
	${R}

ISOS = NGA GHA PAK
# ISOS ?= $(shell cat africaisos.txt) PAK

# provides non-analysis support
# n.b. this is reloaded in sub-make files
include support.makefile

${SOURCE}/%: | ${SOURCE}
	${MAKE} -wC $(notdir $|) $@

setup: | ${SOURCE}
	${MAKE} -wC $(notdir $|)

${GEND}/%: | ${GEND}
	${MAKE} -wC $(notdir $|) $@

generate: | ${GEND}
	${MAKE} -wC $(notdir $|)

${ESTD}/%: | ${ESTD}
	${MAKE} -wC $(notdir $|) $@

estimate r0 intros params: | ${ESTD}
	${MAKE} -wC $(notdir $|) $@

${FIGS}/%: | ${FIGS}
	${SM} $@

figures: | ${FIGS}
	${SM}

hpcclean:
	rm *.out *.err
