
# if present, include local.makefile
# a local.makefile should provide alternative definitions
# of default variables, paths, etc
-include local.makefile

# default assumes subfolders of this repository;
# example local.makefile overrides this to point to Dropbox folder
DATART ?= .

SOURCE := ${DATART}/inputs
SINK   := ${DATART}/outputs

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
MKDIRS := ${SOURCE} ${SINK} $(addprefix ${SINK}/,intervention_timing r0 introductions sample params projections figs variant) $(addprefix ${SOURCE}/,pops yuqs) ${MIRDIR}

# provides non-analysis support
include support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# writes .install if it succeeds otherwise errors uninformatively
.install: get_install.R rpack.txt ${COVIDM}
	${R}

# get + subset the JHU data
# was ECDC data, but now that's only weekly
${SOURCE}/epi_data.rds: get_epi_data.R | ${SOURCE}
	${R}

# get + subset the JHU data
# was ECDC data, but now that's only weekly
${SOURCE}/prov_data.rds: get_prov_data.R | ${SOURCE}
	${R}

${SINK}/cfrs.rds: est_scale_prov.R ${SOURCE}/prov_data.rds | ${SINK}
	${R}

epi: ${SOURCE}/epi_data.rds ${SOURCE}/prov_data.rds ${SINK}/cfrs.rds

resetepi:
	rm ${SOURCE}/epi_data.rds
	make ${SOURCE}/epi_data.rds

${SINK}/phylo.rds: est_phylo_share.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv
	${R}

# only works for ZAF, placeholder that just sets by hand values identified in covidLMIC
${SINK}/intervention_timing/%.rds: gen_r0_est_timing.R | ${SINK}/intervention_timing
	${Rstar}

#${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds ${SINK}/introductions/%.rds | ${SINK}/intervention_timing
${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/%.rds ${SINK}/phylo.rds | ${SINK}/intervention_timing
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

int_r0: ${SINK}/r0/ZAF.rds

${SINK}/introductions/%.rds: est_introductions.R ${SINK}/yuqs/%.rds ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${SOURCE}/epi_data.rds ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

STARTID ?= 0001

${SINK}/sample/%.rds: gen_sample.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds | ${SINK}/sample
	${R}

tarsample: ${SINK}/sample/ZAF.rds

${SINK}/params/%.rds: est_parameters.R ${SOURCE}/pops/%.rds ${SOURCE}/urbanization.rds ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${SINK}/sample/%.rds | ${SINK}/params NGM.rda ${COVIDM}
	Rscript $^ $* ${STARTID} ${COVIDM} $(subst $*,$*_${STARTID},$@)

.SECONDEXPANSION:
${SINK}/params/%_consolidated.rds: gen_consolidate.R $$(wildcard ${SINK}/params/%_*.rds)
	Rscript $< $(@D) $* $@

testpars: ${SINK}/params/ZAF.rds
conspars: ${SINK}/params/ZAF_consolidated.rds

${SINK}/scenarios/%.rds: gen_scenarios.R ${SINK}/fits/%.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds | ${SINK}/scenarios
	${Rstar}

${SINK}/projections/%.rds: sim_relax.R ${SINK}/params/%_consolidated.rds ${SOURCE}/pops/%.rds ${SINK}/introductions/%.rds ${SOURCE}/urbanization.rds ${SINK}/intervention_timing/%.rds | ${SINK}/projections
	${RSCRIPT} $^ $* ${COVIDM} $@

${SINK}/variant/%.rds: est_variant.R ${SINK}/params/%_consolidated.rds ${SOURCE}/pops/%.rds ${SOURCE}/urbanization.rds ${SINK}/projections/%.rds ${SINK}/sample/%.rds | ${SINK}/variant
	${Rstar}

int_scen: ${SINK}/scenarios/ZAF.rds
int_proj: ${SINK}/projections/ZAF.rds
int_var: ${SINK}/variant/ZAF.rds

${SINK}/figs/timeseries.rds: fig_relax_proj.R ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/ZAF.rds ${SINK}/phylo.rds ${SINK}/projections/ZAF.rds | ${SINK}/figs
	${Rstar}

${SINK}/figs/phylo.rds: fig_phylo_share_ts.R ${SINK}/phylo.rds | ${SINK}/figs
	${R}

${SINK}/figs/ccfr.rds: fig_scaled_cfr.R ${SINK}/cfrs.rds | ${SINK}/figs
	${R}

${SINK}/figs/AR.rds: fig_relax_AR.R ${SOURCE}/urbanization.rds ${SOURCE}/pops/ZAF.rds ${SINK}/projections/ZAF.rds | ${SINK}/figs
	Rscript $^ ZAF $@

${SINK}/figs/relaxed_rt.rds: fig_relaxed_rt.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv | ${SINK}/figs
	${R}

figpieces: $(patsubst %,${SINK}/figs/%.rds,phylo cfr timeseries AR)