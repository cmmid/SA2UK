
defaultgen: allpops allyuqs allints

MKDIRS := $(addprefix ${GEND}/,pops yuqs intervention_timing)

include ../support.makefile

${GEND}/urbanization.rds: urbanization.R | ${GEND}
	${R}

${GEND}/matrices.rds: matrices.R | ${GEND}
	${R}

${GEND}/mortality.rds: mortality.R ${SOURCE}/lifetables.csv | ${GEND}
	${R}

${GEND}/fertility.rds: fertility.R ${SOURCE}/fertility.csv | ${GEND}
	${R}

MOBDATA := $(patsubst %,${SOURCE}/%.csv,google_mobility ox_si_schools ox_si_schools_flag)

${GEND}/mobility.rds: mobility.R ${MOBDATA} | ${GEND}
	${RSCRIPT} $^ ${MOBURLS} $@

REFDATA := $(patsubst %,${GEND}/%.rds,mortality fertility urbanization matrices mobility)

${GEND}/intervention_timing/%.rds: r0_est_timing.R ${PROJRT}/interventions.csv | ${GEND}/intervention_timing
	${Rstar}

${GEND}/pops/%.rds: covidm_pop.R ${REFDATA} | ${COVIDM} ${GEND}/pops
	${RSCRIPT} $^ $* $(firstword $|) $@

${GEND}/yuqs/%.rds: reference_qs.R ${SOURCE}/covidm_fit_yu.qs ${GEND}/pops/%.rds ${GEND}/intervention_timing/%.rds ${GEND}/mobility.rds | ${COVIDM} ${GEND}/yuqs
	${RSCRIPT} $^ $* $(firstword $|) $@

allints: $(call alliso,${GEND}/intervention_timing,rds)
allpops: $(call alliso,${GEND}/pops,rds)
allyuqs: $(call alliso,${GEND}/yuqs,rds)

${SOURCE}/% ${ESTD}/% ${FIGS}/%:
	${MAKE} -wC .. $@