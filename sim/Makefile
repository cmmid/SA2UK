
TARISO ?= PAK

MKDIRS := $(addprefix ${SIMD}/,history)

include ../support.makefile

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds

#' TODO strip relaxation calculation
${SIMD}/history/%.rds: param_review.R ${GENITEMS} | ${SIMD}/history
	${RSCRIPT} $^ $* $@

${GEND}/%:
	${MAKE} -wC .. $@

${ESTD}/%:
	${MAKE} -wC .. $@