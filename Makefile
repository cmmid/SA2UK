
#################### SETUP ###################################################

-include local.makefile
# if present, include local.makefile to provide alternative definitions of
# elements assigned via `?=` (e.g. default variables, paths)

# root filesystem location for inputs & outputs
DATART ?= .
# example local.makefile overrides this to point to a Dropbox folder

SOURCE := ${DATART}/inputs
SINK   := ${DATART}/outputs

NCORES ?= 4
NSAMPS ?= 4e3
# sets # of samples and non-embarassingly parallel (cores) populations
# primarily for running EpiNow2

# default assumption: covidm is a sibling to this repository
COVIDMGIT := https://github.com/nicholasdavies/covidm
COVIDM ?= ../covidm

# TODO correct these - seems to clone in this folder instead?
${COVIDM}:
	cd $(dir $@); git clone ${COVIDMGIT} $(notdir $@)

GITLIBS := ${COVIDM}

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${SINK} \
	$(addprefix ${SINK}/,intervention_timing r0 introductions sample params projections figs variant) \
	$(addprefix ${SOURCE}/,pops yuqs figs/epi) ${MIRDIR}

africaisos.txt: gen_isos.R ${SOURCE}/epi_data.rds
	${R}

ISOS = PAK
# ISOS ?= $(shell cat africaisos.txt) PAK

# provides non-analysis support
include support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# performs compile step for covidm repo
# writes .install if it succeeds otherwise errors uninformatively
.install: get_install.R rpack.txt ${COVIDM}
	${R}

# get + subset the JHU data
# was ECDC data, but now that's only weekly
RAWDATA := ${SOURCE}/epi_data.rds

rawdata: ${RAWDATA}

${RAWDATA}: get_epi_data.R | ${SOURCE}
	${R}

# overview of all raw data
${SOURCE}/figs/epi/%.png: fig_epi_overview.R ${RAWDATA} | ${SOURCE}/figs/epi
	${Rstar}

epireview: $(patsubst %,${SOURCE}/figs/epi/%.png,${ISOS})

EPIDATA := ${SINK}/adj_data.rds

${EPIDATA}: est_imputed_data.R ${RAWDATA}
	${R}

adjdata: ${EPIDATA}

${SOURCE}/figs/adjusted/%.png: fig_epi_adjusted.R ${EPIDATA} | ${SOURCE}/figs/adjusted
	${Rstar}

adjreview: $(patsubst %,${SOURCE}/figs/adjusted/%.png,${ISOS})

# get + subset the JHU data
# was ECDC data, but now that's only weekly
${SOURCE}/prov_data.rds: get_prov_data.R | ${SOURCE}
	${R}

${SINK}/cfrs.rds: est_scale_prov.R ${SOURCE}/prov_data.rds | ${SINK}
	${R}

${SOURCE}/urbanization.rds: gen_urbanization.R | ${SOURCE}
	${R}

${SOURCE}/matrices.rds: get_matrices.R | ${SOURCE}
	${R}

${SOURCE}/mortality.rds: gen_mortality.R | ${SOURCE}
	${R}

${SOURCE}/fertility.rds: gen_fertility.R | ${SOURCE}
	${R}

MOBURLS := https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv $(patsubst %,https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/c1_%.csv,school_closing flag)

${SOURCE}/mobility.rds: get_mobility.R | ${SOURCE}
	${RSCRIPT} $^ ${MOBURLS} $@ 

mob: ${SOURCE}/mobility.rds

${SOURCE}/pops/%.rds: gen_covidm_pop.R $(patsubst %,${SOURCE}/%.rds,mortality fertility urbanization matrices) | ${COVIDM} ${SOURCE}/pops
	${RSCRIPT} $^ $* ${COVIDM} $@

allpops: $(patsubst %,${SOURCE}/pops/%.rds,${ISOS})

${SOURCE}/populations.rds: gen_populations.R | ${SOURCE}
	${R}

${SOURCE}/ox_timings.rds: gen_ox_si_timing.R | ${SOURCE}
	${R}

allinputs: ${SOURCE}/epi_data.rds ${SOURCE}/urbanization.rds ${SOURCE}/populations.rds

resetepi:
	rm ${SOURCE}/epi_data.rds
	make ${SOURCE}/epi_data.rds

${SINK}/phylo.rds: est_phylo_share.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv
	${R}

# only works for ZAF, placeholder that just sets by hand values identified in covidLMIC
${SINK}/intervention_timing/%.rds: gen_r0_est_timing.R interventions.csv | ${SINK}/intervention_timing
	${Rstar}

#${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds ${SINK}/introductions/%.rds | ${SINK}/intervention_timing
${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SINK}/phylo.rds | ${SINK}/intervention_timing
	${Rstar}

timing: $(patsubst %,${SINK}/intervention_timing/%.png,${ISOS})

NGM.rda: NGM.R
	${R}

${SOURCE}/yuqs/%.rds: gen_reference_qs.R ${SOURCE}/covidm_fit_yu.qs ${SOURCE}/pops/%.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/mobility.rds | ${SOURCE}/yuqs ${COVIDM}
	${RSCRIPT} $^ $* ${COVIDM} $@

.PRECIOUS: ${SINK}/intervention_timing/%.rds ${SOURCE}/yuqs/%.rds ${SOURCE}/pops/%.rds ${SINK}/introductions/%.rds

#' TODO strip relaxation calculation
${SINK}/r0/%.rds: est_r0.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SOURCE}/yuqs/%.rds ${SOURCE}/pops/%.rds | ${SINK}/r0
	${RSCRIPT} $^ ${NCORES} ${NSAMPS} $* $@

int_r0: $(patsubst %,${SINK}/r0/%.rds,PAK)

${SINK}/introductions/%.rds: est_introductions.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${EPIDATA} ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

intros: $(patsubst %,${SINK}/introductions/%.rds,PAK)

STARTID ?= 0001

${SINK}/sample/%.rds: gen_sample.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds | ${SINK}/sample
	${R}

samples: $(patsubst %,${SINK}/sample/%.rds,PAK)

${SINK}/params/%.rds: est_parameters.R ${SOURCE}/pops/%.rds ${SINK}/r0/%.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${SINK}/sample/%.rds | ${SINK}/params ${COVIDM}
	Rscript $^ $* ${STARTID} ${COVIDM} $(subst $*,$*_${STARTID},$@)

pars:
	make $(patsubst %,${SINK}/params/%.rds,PAK) STARTID=0051
	make $(patsubst %,${SINK}/params/%.rds,PAK) STARTID=0056

.SECONDEXPANSION:
${SINK}/params/%_consolidated.rds: gen_consolidate.R $$(wildcard ${SINK}/params/%_*.rds)
	Rscript $< $(@D) $* $@

testpars: ${SINK}/params/ZAF.rds
conspars: ${SINK}/params/ZAF_consolidated.rds

${SINK}/scenarios.rds: gen_scenarios.R ${SINK}/fits/%.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds | ${SINK}/scenarios
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

########################## CONVENIENCE TARGETS ###########################################

popprep: $(patsubst %,${SOURCE}/%.rds,mortality fertility urbanization matrices)