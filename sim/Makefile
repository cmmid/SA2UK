
TARISO ?= PAK

MKDIRS := $(addprefix ${SIMD}/,history)

include ../support.makefile

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds
ESTITEMS := ${ESTD}/introductions/%.rds ${ESTD}/params/%s.rds

#' TODO strip relaxation calculation
${SIMD}/history/%.rds: history.R ${GENITEMS} ${GEND}/mobility.rds ${ESTITEMS} | ${SIMD}/history
	${RSCRIPT} $^ $* $@

${GEND}/%:
	${MAKE} -wC .. $@

${ESTD}/%:
	${MAKE} -wC .. $@