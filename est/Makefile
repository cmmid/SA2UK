
allest: allr0 allintros

MKDIRS := $(addprefix ${SINK}/,r0 introductions sample)

include ../support.makefile

#' TODO strip relaxation calculation
${SINK}/r0/%.rds: r0.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SOURCE}/yuqs/%.rds ${SOURCE}/pops/%.rds | ${SINK}/r0
	${RSCRIPT} $^ ${NCORES} ${NSAMPS} $* $@

${SINK}/introductions/%.rds: introductions.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${EPIDATA} ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

${SINK}/sample/%.rds: sample.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds | ${SINK}/sample
	${R}

allr0: $(call alliso,${SINK}/r0,rds)
allintros: $(call alliso,${SINK}/r0,rds)
allsamples: $(call alliso,${SINK}/sample,rds)

${SINK}/params/%.rds: est_parameters.R ${SOURCE}/pops/%.rds ${SINK}/r0/%.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${SINK}/sample/%.rds | ${SINK}/params ${COVIDM}
	Rscript $^ $* ${STARTID} ${COVIDM} $(subst $*,$*_${STARTID},$@)

.SECONDEXPANSION:
${SINK}/params/%_consolidated.rds: gen_consolidate.R $$(wildcard ${SINK}/params/%_*.rds)
	Rscript $< $(@D) $* $@

testpars: ${SINK}/params/ZAF.rds
conspars: ${SINK}/params/PAK_consolidated.rds

${SINK}/params/%.png: fig_params_overview.R ${SINK}/params/%_consolidated.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${EPIDATA} | ${SINK}/params ${COVIDM}
	Rscript $^ $* ${COVIDM} $@

${SINK}/scenarios.rds: gen_scenarios.R ${SINK}/fits/%.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds | ${SINK}/scenarios
	${Rstar}

${SINK}/projections/%.rds: sim_relax.R ${SINK}/params/%_consolidated.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${EPIDATA} | ${SINK}/projections ${COVIDM}
	${RSCRIPT} $^ $* ${COVIDM} $@

${SINK}/variant/%.rds: est_variant.R ${SINK}/params/%_consolidated.rds ${SOURCE}/pops/%.rds ${SOURCE}/urbanization.rds ${SINK}/projections/%.rds ${SINK}/sample/%.rds | ${SINK}/variant
	${Rstar}

${SINK}/scenario/%.rds: sim_scenarios.R ${SINK}/pops/%.rds ${SINK}/params/%_consolidated.rds ${SINK}/projections/%.rds ${SINK}/sample/%.rds | ${SINK}/variant
	${RSCRIPT} $^ $* ${COVIDM} $@

int_scen: ${SINK}/scenarios/ZAF.rds
int_proj: ${SINK}/projections/ZAF.rds
int_var: ${SINK}/variant/ZAF.rds

${SINK}/figs/timeseries.rds: fig_relax_proj.R ${SOURCE}/epi_data.rds ${SINK}/intervention_timing/ZAF.rds ${SINK}/phylo.rds ${SINK}/projections/ZAF.rds | ${SINK}/figs
	${Rstar}

${SINK}/figs/phylo.rds: fig_phylo_share_ts.R ${SINK}/phylo.rds | ${SINK}/figs
	${R}

#${SINK}/figs/ccfr.rds: fig_scaled_cfr.R ${SINK}/cfrs.rds | ${SINK}/figs
#	${R}

${SINK}/figs/AR.rds: fig_relax_AR.R ${SOURCE}/urbanization.rds ${SOURCE}/pops/ZAF.rds ${SINK}/projections/ZAF.rds | ${SINK}/figs
	Rscript $^ ZAF $@

${SINK}/figs/relaxed_rt.rds: fig_relaxed_rt.R ${SOURCE}/nextstrain_groups_ngs-sa_COVID19-ZA-2020.12.17_metadata.tsv | ${SINK}/figs
	${R}

figpieces: $(patsubst %,${SINK}/figs/%.rds,phylo cfr timeseries AR)

########################## CONVENIENCE TARGETS ###########################################

popprep: $(patsubst %,${SOURCE}/%.rds,mortality fertility urbanization matrices)