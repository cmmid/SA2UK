
# if present, include local.makefile
-include local.makefile

# default assumes subfolders of this repository;
# example local.makefile overrides these to Dropbox folder
DATART ?= .

SOURCE := ${DATART}/inputs
SINK := ${DATART}/outputs

# for running EpiNow2; should override in local.makefile
NCORES ?= 4
NSAMPS ?= 8e3

# TODO define parallel path, url lists; feed to something in support.makefile
# TODO that support.makefile should also do existence checks, and pull if repo
#   already exists

# default assumption: covidm is a sibling to this repository
COVIDM ?= ../covidm
COVIDMGIT := https://github.com/nicholasdavies/covidm

# TODO correct these - seems to clone in this folder instead?
${COVIDM}:
	cd $(dir $@); git clone ${COVIDMGIT} $(notdir $@)

GITLIBS := ${COVIDM}

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${SINK} $(addprefix ${SINK}/,intervention_timing r0 fits introductions scenarios mod_scenarios projections figs relaxation) $(addprefix ${SOURCE}/,pops yuqs) ${MIRDIR}

# provides non-analysis support
include support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# writes .install if it succeeds otherwise errors uninformatively
.install: get_install.R rpack.txt
	${R}

# get + subset the JHU data
# was ECDC data, but now that's only weekly
${SOURCE}/epi_data.rds: get_epi_data.R | ${SOURCE}
	${R}

${SINK}/phylo.rds: est_phylo_share.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv
	${R}

# only works for ZAF, placeholder that just sets by hand values identified in covidLMIC
${SINK}/intervention_timing/%.rds: gen_r0_est_timing.R | ${SINK}/intervention_timing
	${Rstar}

#${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds ${SINK}/introductions/%.rds | ${SINK}/intervention_timing
${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/%.rds | ${SINK}/intervention_timing
	${Rstar}

${SOURCE}/populations.rds: gen_populations.R | ${SOURCE}
	${R}

WBURL := http://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv

${SOURCE}/urbanization.rds: gen_urbanization.R | ${SOURCE}
	curl ${WBURL} --output tmp.urban.zip
	unzip tmp.urban.zip API_SP.URB.TOTL.IN.ZS_DS2_*.csv
	mv API_SP.URB.TOTL.IN.ZS_DS2_*.csv tmp.urban.csv
	Rscript $^ tmp.urban.csv $@

${SOURCE}/pops/%.rds: gen_covidm_pop.R | ${COVIDM} ${SOURCE}/pops
	${RSCRIPT} $^ $* ${COVIDM} $@

NGM.rda: NGM.R
	${R}

${SOURCE}/yuqs/%.rds: gen_reference_qs.R ${SOURCE}/covidm_fit_yu.qs ${SOURCE}/pops/%.rds | ${SOURCE}/yuqs NGM.rda
	${R}

#' TODO strip relaxation calculation
${SINK}/r0/%.rds: est_r0.R ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/yuqs/%.rds | ${SINK}/r0
	${RSCRIPT} $^ ${NCORES} ${NSAMPS} $* $@

${SINK}/relaxation/%.rds: est_relaxation_rt.R ${SOURCE}/epi_data.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs | ${SINK}/relaxation NGM.rda
	Rscript $^ 2020-06-10 2020-10-15 $* $@

relax: ${SINK}/relaxation/ZAF.rds

int_r0: ${SINK}/r0/ZAF.rds

${SINK}/introductions/%.rds: est_introductions.R ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${SOURCE}/epi_data.rds ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

${SINK}/fits/%.rds: est_fits.R ${SINK}/r0/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs | ${SINK}/fits NGM.rda
	${Rstar}

${SINK}/scenarios/%.rds: gen_scenarios.R ${SINK}/fits/%.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds | ${SINK}/scenarios
	${Rstar}

${SINK}/projections/%.qs: sim_scenarios.R ${SINK}/scenarios/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs ${SINK}/r0/%.rds ${SINK}/introductions/%.rds ${SOURCE}/urbanization.rds | ${SINK}/projections
	${RSCRIPT} $^ $* ${COVIDM} $@

${SINK}/mod_scenarios/%.rds: est_mod_r0.R ${SINK}/scenarios/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs ${SINK}/r0/%.rds ${SINK}/introductions/%.rds ${SOURCE}/urbanization.rds | ${SINK}/mod_scenarios
	${RSCRIPT} $^ $* ${COVIDM} $@

default: ${SINK}/projections/ZAF.qs ${SINK}/mod_scenarios/ZAF.rds

${SINK}/projections/%.png: fig_projection.R ${SINK}/projections/%.qs ${SOURCE}/pops/%.rds ${SINK}/introductions/%.rds ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/urbanization.rds
	${Rstar}

${SINK}/figs/phylo.rds: fig_phylo_share_ts.R ${SINK}/phylo.rds | ${SINK}/figs
	${R}

${SINK}/figs/cfr.rds: fig_crude_cfr.R ${SOURCE}/epi_data.rds | ${SINK}/figs
	${R}

${SINK}/figs/relaxed_rt.rds: fig_relaxed_rt.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv | ${SINK}/figs
	${R}

figpieces: $(patsubst %,${SINK}/figs/%.rds,phylo cfr)