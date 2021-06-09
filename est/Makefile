
allest: allr0 allintros allsamples allparams

MKDIRS := $(addprefix ${ESTD}/,r0 introductions sample params)

include ../support.makefile

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/yuqs/%.rds ${GEND}/pops/%.rds

#' TODO strip relaxation calculation
${ESTD}/r0/%.rds: r0.R ${EPIDATA} ${SOURCE}/sampling.json ${GENITEMS} | ${ESTD}/r0
	${RSCRIPT} $^ $* $@

${ESTD}/introductions/%.rds: introductions.R ${GENITEMS} ${ESTD}/r0/%.rds ${EPIDATA} ${PROJRT}/ene-ifr.csv | ${SINK}/introductions
	${Rstar}

${SINK}/sample/%.rds: sample.R ${SOURCE}/yuqs/%.rds ${SINK}/r0/%.rds | ${SINK}/sample
	${R}

${SINK}/params/%.rds: parameters.R ${SOURCE}/pops/%.rds ${SINK}/r0/%.rds ${SOURCE}/mobility.rds ${SINK}/intervention_timing/%.rds ${SINK}/introductions/%.rds ${SINK}/sample/%.rds | ${COVIDM} ${SINK}/params
	Rscript $^ $* $(firstword $|) $@

allr0: $(call alliso,${ESTD}/r0,rds)
allintros: $(call alliso,${ESTD}/r0,rds)
allsamples: $(call alliso,${ESTD}/sample,rds)
allparams: $(call alliso,${ESTD}/params,rds)