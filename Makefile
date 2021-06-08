
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

SOURCE := ${DATART}/inputs
SINK   := ${DATART}/outputs

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${SINK} ${SINK}/intervention_timing

africaisos.txt: gen_isos.R ${SOURCE}/epi_data.rds
	${R}

ISOS = NGA GHA ETH PAK
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

# get + subset the JHU data
# was ECDC data, but now that's only weekly
RAWDATA := ${SOURCE}/epi_data.rds
EPIDATA := ${SINK}/adj_data.rds

${RAWDATA}: get_epi_data.R | ${SOURCE}
	${R}

${EPIDATA}: est_imputed_data.R ${RAWDATA} | ${SINK}
	${R}

datasetup: ${RAWDATA} ${EPIDATA}

allfigs epireview adjreview timing: | allgen
	${MAKE} -wC figs $@

allgen ${SINK}/intervention_timing/%.rds:
	${MAKE} -wC gen $@

#' assumes the replacement trend of beta in ZAF
${SINK}/phylo.rds: est_phylo_share.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv
	${R}

