
# if present, include local.makefile
-include local.makefile

# default assumes subfolders of this repository;
# example local.makefile overrides these to Dropbox folder
# n.b.: these paths cannot be identical, but this is *not*
# enforced by the Makefile
SOURCE ?= ./inputs
SINK ?= ./outputs

# for running EpiNow2; should override in local.makefile
NCORES ?= 4
NSAMPS ?= 8e3

# TODO define parallel path, url lists; feed to something in support.makefile
# TODO that support.makefile should also do existence checks, and pull if repo
#   already exists

# default assumption: covidm is a sibling to this repository
COVIDM ?= ../covidm
COVIDMGIT := https://github.com/nicholasdavies/covidm

# default assumption: Oxford Policy Tracker is a sibling to this repository
OXDATA ?= ../covid-policy-tracker
OXGIT := git@github.com:OxCGRT/covid-policy-tracker.git

# TODO correct these - seems to clone in this folder instead?
${COVIDM}:
	cd $(dir $@); git clone ${COVIDMGIT} $(notdir $@)

${OXDATA}:
	cd $(dir $@); git clone --depth 1 ${OXGIT}

GITLIBS := ${COVIDM} ${OXDATA}

MIRDIR := ~/Dropbox/Covid_LMIC/All_Africa_paper/projections

# support.makefile will provide a directory target for all of these
MKDIRS := ${SOURCE} ${SINK} $(addprefix ${SINK}/,intervention_timing r0 fits introductions scenarios mod_scenarios projections) ${SOURCE}/isos ${SOURCE}/pops ${MIRDIR}

# provides non-analysis support
include support.makefile

# very basic R package installation; THIS MAY NOT "JUST WORK"
# checks availability of packages (via require)
# attempts to install any not available
# writes .install if it succeeds otherwise errors uninformatively
.install: get_install.R rpack.txt
	${R}

# get + clean the ECDC data
${SOURCE}/ecdc_data.rds: get_ecdc_data.R | ${SOURCE}
	${R}

# re-analyze the Oxford CPT data
${SOURCE}/ox_si_timing.csv: est_ox_si_timing.R | ${OXDATA}
	${Rpipe}

${SINK}/interventions.rds: est_interventions.R ${SOURCE}/ox_si_timing.csv ${SOURCE}/isos/all.iso timing_adjust.csv | ${SINK}
	${R}

${SINK}/intervention_timing/%.rds: gen_r0_est_timing.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds | ${SINK}/intervention_timing
	${Rstar}

#${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds ${SINK}/introductions/%.rds | ${SINK}/intervention_timing
${SINK}/intervention_timing/%.png: fig_assess_interventions.R ${SINK}/interventions.rds ${SOURCE}/ecdc_data.rds | ${SINK}/intervention_timing
	${Rstar}

${SOURCE}/populations.rds: gen_populations.R | ${SOURCE}
	${R}

WBURL := http://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv

${SOURCE}/urbanization.rds: gen_urbanization.R | ${SOURCE}
	curl ${WBURL} --output tmp.urban.zip
	unzip tmp.urban.zip API_SP.URB.TOTL.IN.ZS_DS2_*.csv
	mv API_SP.URB.TOTL.IN.ZS_DS2_*.csv tmp.urban.csv
	Rscript $^ tmp.urban.csv $@

${SOURCE}/pops/%.rds: gen_covidm_pop.R ${SOURCE}/contact_substitutes.csv | ${COVIDM} ${SOURCE}/pops
	${RSCRIPT} $^ $* ${COVIDM} $@

NGM.rda: NGM.R
	${R}

${SINK}/r0/%.rds: est_r0.R ${SOURCE}/ecdc_data.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs | ${SINK}/r0 NGM.rda
	${RSCRIPT} $^ ${NCORES} ${NSAMPS} $* $@

${SINK}/introductions/%.rds: est_introductions.R ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${SOURCE}/ecdc_data.rds ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

${SINK}/fits/%.rds: est_fits.R ${SINK}/r0/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs | ${SINK}/r0 NGM.rda
	${Rstar}

${SINK}/scenarios/%.rds: gen_scenarios.R ${SINK}/fits/%.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds | ${SINK}/scenarios
	${Rstar}

${SINK}/mod_scenarios/%.rds: est_mod_r0.R ${SINK}/scenarios/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs ${SINK}/r0/%.rds ${SINK}/introductions/%.rds ${SOURCE}/urbanization.rds | ${SINK}/mod_scenarios
	${RSCRIPT} $^ $* ${COVIDM} $@

${SINK}/projections/%.qs: sim_scenarios.R ${SINK}/mod_scenarios/%.rds ${SOURCE}/pops/%.rds ${SOURCE}/covidm_fit_yu.qs ${SINK}/r0/%.rds ${SINK}/introductions/%.rds ${SOURCE}/urbanization.rds | ${SINK}/projections
	${RSCRIPT} $^ $* ${COVIDM} $@

${SINK}/projections/%.png: fig_projection.R ${SINK}/projections/%.qs ${SOURCE}/pops/%.rds ${SINK}/introductions/%.rds ${SOURCE}/ecdc_data.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/urbanization.rds
	${Rstar}

ALLISOQS := $(shell ls ${SINK}/projections/*.qs)
ALLISOS := $(shell more ${SOURCE}/isos/africa.iso) PAK
#ALLISOS := $(subst .qs,,$(shell cd ${SINK}/projections; ls *.qs))
#$(info ${ALLISOS})
#$(shell cat ${SOURCE}/isos/africa.iso)

testiso: ${SINK}/projections/CAF.png

setup: ${COVIDM} .install $(patsubst %,${SOURCE}/%.rds,ecdc_data populations urbanization) NGM.rda

estimate: $($(patsubst %,${SINK}/%.rds,interventions))

${SOURCE}/isos/all.iso: ${SOURCE}/ox_si_timing.csv | ${SOURCE}/isos
	tail -n+2 $^ | cut -d, -f1 > $@

${SOURCE}/isos/%.iso: gen_iso_files.R ${SOURCE}/isos/all.iso
	${Rstar}

EMPTY :=
SPACE := ${EMPTY} ${EMPTY}
OTHISOS := PAK AFG

${SOURCE}/isos/other.iso:
	echo "$(subst ${SPACE},\n,${OTHISOS})" > $@

isos: $(patsubst %,${SOURCE}/isos/%.iso,all africa other)

evaluate: $(patsubst %,${SINK}/intervention_timing/%.png,${ALLISOS})

pops: $(patsubst %,${SOURCE}/pops/%.rds,${ALLISOS})

r0: $(patsubst %,${SINK}/r0/%.rds,${ALLISOS})

fits: $(patsubst %,${SINK}/fits/%.rds,${ALLISOS})

intros: $(patsubst %,${SINK}/introductions/%.rds,${ALLISOS})

scenarios: $(patsubst %,${SINK}/scenarios/%.rds,${ALLISOS})

mod_scenarios: $(patsubst %,${SINK}/mod_scenarios/%.rds,${ALLISOS})

projections: $(patsubst %,${SINK}/projections/%.qs,${ALLISOS})

int_figs: $(patsubst %,${SINK}/intervention_timing/%.png,${ALLISOS})

clean_int_figs:
	rm ${SINK}/intervention_timing/*.png

figs: $(patsubst %,${SINK}/projections/%.png,${ALLISOS})

.PRECIOUS: ${SINK}/scenarios/%.rds ${SINK}/mod_scenarios/%.rds ${SINK}/introductions/%.rds ${SINK}/fits/%.rds ${SINK}/r0/%.rds ${SOURCE}/pops/%.rds ${SINK}/intervention_timing/%.rds

mirror: | ${MIRDIR}
	cd ${SINK}/projections; rsync -t *.png $|/ 