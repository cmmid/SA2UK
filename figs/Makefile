# organization file for making figures
IFDIR := ${SOURCE}/figs
OFDIR := ${SINK}/figs

MKDIRS := ${IFDIR} ${OFDIR} ${IFDIR}/epi ${IFDIR}/intervention_timing ${IFDIR}/adjusted

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds

include ../support.makefile

# overview of all raw data
${IFDIR}/epi/%.png: epi_overview.R ${RAWDATA} | ${IFDIR}/epi
	${Rstar}

# after smoothing on presumed missing data
${IFDIR}/adjusted/%.png: epi_adjusted.R ${EPIDATA} | ${IFDIR}/adjusted
	${Rstar}

${IFDIR}/intervention_timing/%.png: assess_interventions.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SINK}/phylo.rds | ${IFDIR}/intervention_timing
	${Rstar}

${IFDIR}/prefitting/%.png: pre_fitting.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SINK}/phylo.rds | ${IFDIR}/intervention_timing
	${Rstar}

epireview: $(call alliso,${IFDIR}/epi,png)
adjreview: $(call alliso,${IFDIR}/adjusted,png)
timing: $(call alliso,${IFDIR}/intervention_timing,png)
allfigs: epireview adjreview timing

.PRECIOUS: %.rds

# may need to send inputs dependencies back to parent
%.rds:
	$(MAKE) -C .. $@