
${MKDIRS}:
	mkdir -p $@

alldirs: ${MKDIRS}

alllibs: ${GITLIBS}

# will copy example.makefile for a locally edited copy
egcopy: | example.makefile
	cp $| local.makefile

# DATA FETCHING DEFINITION; USE AS ${WGET} ${THEURL}
WGET = wget -c -O $@

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

# default assumption: covidm is the version shipped with this repository
# but could be defined to live elsewhere in local.makefile
COVIDM ?= ${PROJRT}/covidm
COVIDMGIT := git@github.com:nicholasdavies/covidm.git

# TODO correct these - seems to clone in this folder instead?
${COVIDM}:
	cd $(dir $@); git clone --single-branch --branch ngmupdate ${COVIDMGIT} $(notdir $@)

rsync: | ${DATART}
	rsync -a ${HPCDIR} ${DATART}

rsyncvn: | ${DATART}
	rsync -avn ${HPCDIR} ${DATART}