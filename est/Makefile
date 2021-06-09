
allest: allr0 allintros allsamples allparams

MKDIRS := $(addprefix ${SINK}/,r0 introductions sample params)

include ../support.makefile

#' TODO strip relaxation calculation
${SINK}/r0/%.rds: r0.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SOURCE}/yuqs/%.rds ${SOURCE}/pops/%.rds | ${SINK}/r0
	${RSCRIPT} $^ ${NCORES} ${NSAMPS} $* $@

${SINK}/introductions/%.rds: introductions.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds ${SOURCE}/populations.rds ${SOURCE}/pops/%.rds ${EPIDATA} ene-ifr.csv ${SINK}/intervention_timing/%.rds | ${SINK}/introductions
	${Rstar}

${SINK}/sample/%.rds: sample.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds | ${SINK}/sample
	${R}

${SINK}/params/%.rds: est_parameters.R ${SOURCE}/pops/%.rds ${SINK}/r0/%.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${SINK}/sample/%.rds | ${COVIDM} ${SINK}/params
	Rscript $^ $* $(firstword $|) $@

allr0: $(call alliso,${SINK}/r0,rds)
allintros: $(call alliso,${SINK}/r0,rds)
allsamples: $(call alliso,${SINK}/sample,rds)
allparams: $(call alliso,${SINK}/params,rds)