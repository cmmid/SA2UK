
insdefault: ${SOURCE}/.install ${RAWDATA} ${EPIDATA} $(addprefix ${SOURCE}/,\
	phylo.rds sampling.json google_mobility.csv ox_si_schools.csv ox_si_schools_flag.csv lifetables.csv\
	sampling.json scenario_reference.rds\
)

# provides non-analysis support
# n.b. this is reloaded in sub-make files
include ../support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# performs compile step for covidm repo
# writes .install if it succeeds otherwise errors uninformatively
${SOURCE}/.install: install.R ${PROJRT}/ins/rpack.csv | ${COVIDM}
	${Rpipe}

${SOURCE}/shared_funs.rda: shared_funs.R
	${R}

${RAWDATA}: epi_data.R | ${SOURCE}
	${R}

${EPIDATA}: impute_data.R ${RAWDATA} | ${SOURCE}
	${R}

# assumes the replacement trend of beta in ZAF
# requires nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv be manually added
${SOURCE}/phylo.rds: phylo_share.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv
	${R}

GOOGURL := https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv
OXURLS := $(patsubst %,https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/c1_%.csv,school_closing flag)

${SOURCE}/google_mobility.csv:
	${WGET} ${GOOGURL}

${SOURCE}/ox_si_schools.csv:
	${WGET} $(firstword ${OXURLS})

${SOURCE}/ox_si_schools_flag.csv:
	${WGET} $(lastword ${OXURLS})

MORTURL := https://population.un.org/wpp/Download/Files/1_Indicators (Standard)/CSV_FILES/WPP2019_Life_Table_Medium.csv
LIFEURL := https://population.un.org/wpp/Download/Files/1_Indicators (Standard)/CSV_FILES/WPP2019_Period_Indicators_Medium.csv

${SOURCE}/lifetables.csv:
	${WGET} "${MORTURL}"

${SOURCE}/fertility.csv:
	${WGET} "${LIFEURL}"

${SOURCE}/sampling.json: config.R
	${RSCRIPT} $^ ${NSAMPLES} $@

#' this sets up the intervention scenarios to be considered

#' these provide tables of vaccine efficacy, by product, dose[, and variant?]
${SOURCE}/vaccines.rds: vaccines.R
	${R}

#' this is a function for generating the campaign scenarios in covidm-ready format,
#' from covidm parameter object + campaign number
${SOURCE}/campaign_function.rds: campaign_function.R
	${R}

${SOURCE}/scenario_reference.rds: generic_scenarios.R $(patsubst %, ${SOURCE}/%.rds,vaccines campaign_function)
	${R}