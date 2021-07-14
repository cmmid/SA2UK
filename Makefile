
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
SIMD   := ${DATART}/sim
FIGS   := ${DATART}/fig

RAWDATA := ${SOURCE}/epi_data.rds
EPIDATA := ${SOURCE}/adj_data.rds

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${GEND} ${ESTD} ${SIMD} ${FIGS}

# override with local.makefile
NSAMPLES ?= 4000

africaisos.txt: gen_isos.R ${SOURCE}/epi_data.rds
	${R}

ISOS = NGA GHA PAK
# ISOS ?= $(shell cat africaisos.txt) PAK

# provides non-analysis support
# n.b. this is reloaded in sub-make files
include support.makefile

# n.b. while the following all have the same recipe structure
# can't define a variable like, e.g., ${R}, because make explicitly
# looks for ${MAKE} when following special sub-make capabilities
# e.g. continuing with a `-n` option

${SOURCE}/%: | ${SOURCE}
	${MAKE} -wC $(notdir $|) $@

setup: | ${SOURCE}
	${MAKE} -wC $(notdir $|)

install: | ${SOURCE}
	${MAKE} -wC $(notdir $|) ${SOURCE}/.$@

${GEND}/%: | ${GEND}
	${MAKE} -wC $(notdir $|) $@

generate: | ${GEND}
	${MAKE} -wC $(notdir $|)

${ESTD}/%: | ${ESTD}
	${MAKE} -wC $(notdir $|) $@

estimate r0 intros params clean_params: | ${ESTD}
	${MAKE} -wC $(notdir $|) $@

${SIMD}/%: | ${SIMD}
	${MAKE} -wC $(notdir $|) $@

${FIGS}/%: | ${FIGS}
	${MAKE} -wC $(notdir $|) $@

figures: | ${FIGS}
	${MAKE} -wC $(notdir $|) $@

hpcclean:
	rm *.out *.err
