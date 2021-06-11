
insdefault: ${SOURCE}/.install ${RAWDATA} ${EPIDATA} $(addprefix ${SOURCE}/,\
	phylo.rds sampling.json google_mobility.csv ox_si_schools.csv ox_si_schools_flag.csv lifetables.csv\
	sampling.json\
)

# provides non-analysis support
# n.b. this is reloaded in sub-make files
include ../support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# performs compile step for covidm repo
# writes .install if it succeeds otherwise errors uninformatively
${SOURCE}/.install: get_install.R ${PROJRT}/ins/rpack.txt | ${COVIDM}
	${Rpipe}


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

WPPURLS := $(patsubst %,"https://population.un.org/wpp/Download/Files/1_Indicators (Standard)/CSV_FILES/WPP2019_%_Medium.csv",Life_Table Period_Indicators)
MORTURL := $(firstword ${WPPURLS})
LIFEURL := $(lastword ${WPPURLS})

${SOURCE}/lifetables.csv:
	${WGET} ${MORTURL}

${SOURCE}/fertility.csv:
	${WGET} ${LIFEURL}

${SOURCE}/sampling.json: config.R
	${RSCRIPT} $^ ${NSAMPLES} $@
