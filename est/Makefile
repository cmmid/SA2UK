
TARISO ?= PAK

defest: ${ESTD}/sample/${TARISO}.rds

MKDIRS := $(addprefix ${ESTD}/,r0 introductions sample params)

include ../support.makefile

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds

#' TODO strip relaxation calculation
${ESTD}/r0/%.rds: r0.R ${EPIDATA} ${SOURCE}/sampling.json ${GENITEMS} | ${ESTD}/r0
	${RSCRIPT} $^ $* $@

${ESTD}/introductions/%.rds: introductions.R ${GENITEMS} ${ESTD}/r0/%.rds ${EPIDATA} ${PROJRT}/ene-ifr.csv | ${ESTD}/introductions
	${Rstar}

${ESTD}/sample/%.rds: sample.R ${GEND}/yuqs/%.rds ${ESTD}/r0/%.rds | ${ESTD}/sample
	${R}

.PRECIOUS: ${ESTD}/r0/%.rds ${ESTD}/introductions/%.rds ${ESTD}/sample/%.rds

${ESTD}/params/%.rds: parameters.R ${SOURCE}/pops/%.rds ${ESTD}/r0/%.rds ${ESTD}/intervention_timing/%.rds ${ESTD}/introductions/%.rds ${ESTD}/sample/%.rds ${GEND}/mobility.rds | ${COVIDM} ${ESTD}/params
	Rscript $^ $* $(firstword $|) $@

allr0: $(call alliso,${ESTD}/r0,rds)
allintros: $(call alliso,${ESTD}/r0,rds)
allsamples: $(call alliso,${ESTD}/sample,rds)
allparams: $(call alliso,${ESTD}/params,rds)