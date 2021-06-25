
TARISO ?= PAK

estimate: params
r0: ${ESTD}/r0/${TARISO}.rds
intros: ${ESTD}/introductions/${TARISO}.rds
params: ${ESTD}/params/${TARISO}.rds
variant: ${ESTD}/variant/${TARISO}.rds

MKDIRS := $(addprefix ${ESTD}/,r0 introductions sample params variant)

include ../support.makefile

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds

#' TODO strip relaxation calculation
${ESTD}/r0/%.rds: r0.R ${EPIDATA} ${SOURCE}/sampling.json ${GENITEMS} | ${ESTD}/r0
	${RSCRIPT} $^ $* $@

${ESTD}/sample/%.rds: sample.R ${GEND}/yuqs/%.rds ${ESTD}/r0/%.rds | ${ESTD}/sample
	${R}

${ESTD}/introductions/%.rds: introductions.R ${GENITEMS} ${ESTD}/sample/%.rds ${EPIDATA} ${PROJRT}/ene-ifr.csv | ${COVIDM} ${ESTD}/introductions
	${RSCRIPT} $^ $* $(firstword $|) $@

.PRECIOUS: ${ESTD}/r0/%.rds ${ESTD}/introductions/%.rds ${ESTD}/sample/%.rds

${ESTD}/params/%.rds: parameters.R ${GENITEMS} ${GEND}/mobility.rds ${ESTD}/r0/%.rds ${ESTD}/sample/%.rds ${ESTD}/introductions/%.rds ${EPIDATA} | ${COVIDM} ${ESTD}/params
	${Rscript} $^ $* $(firstword $|) $@

${ESTD}/variant/%.rds: variant.R ${GENITEMS} ${ESTD}/params/%.rds ${SIMD}/history/%.rds ${ESTD}/sample/%.rds | ${ESTD}/variant
	${Rstar}

clean_params:
	rm ${ESTD}/params/*.rds