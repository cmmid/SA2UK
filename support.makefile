
# from https://stackoverflow.com/questions/10858261/abort-makefile-if-variable-not-set
check_defined = $(strip $(foreach 1,$1,$(call __check_defined,$1,$(strip $(value 2)))))
__check_defined = $(if $(value $1),, $(error Undefined $1$(if $2, ($2))))

${MKDIRS}:
	mkdir -p $@

alldirs: ${MKDIRS}

# will copy example.makefile for a locally edited copy
local.makefile: | example.makefile
	cp $| $@

# ASSORTED CONVENIENCE DEFINITIONS FOR R-RELATED RULES

# get path to Rscript; complain if undefined
RSCRIPT := $(shell which Rscript)
$(call check_defined, RSCRIPT)

# the following define common rule arrangements
# assumes first dependency is something.R, called with other dependencies ($^) &
# the target ($@) as arguments
R = ${RSCRIPT} $^ $@

# for pattern matching rules, to pass the match (rather than extracting it from $@)
Rstar = ${RSCRIPT} $^ $* $@

# for with pipe dependencies
Rpipe = ${RSCRIPT} $^ $| $@

# star + pipe
Rsp = ${RSCRIPT} $^ $* $| $@

FROMDIR := override/path/on/commandline

xfer: | ${FROMDIR}
	rsync -rv ${FROMDIR}/ ${SINK}
