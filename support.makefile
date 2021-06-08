
${MKDIRS}:
	mkdir -p $@

alldirs: ${MKDIRS}

alllibs: ${GITLIBS}

# will copy example.makefile for a locally edited copy
egcopy: | example.makefile
	cp $| local.makefile

# ASSORTED CONVENIENCE DEFINITIONS FOR R-RELATED RULES

# get path to Rscript; complain if undefined
RSCRIPT := $(shell which Rscript)

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

FROMDIR ?= override/path/on/commandline

xfer: | ${FROMDIR}
	rsync -rv $|/ ${SINK}

# intended to be used as $(call alliso,PATH,MIME)
alliso = $(patsubst %,$(1)/%.$(2),${ISOS})