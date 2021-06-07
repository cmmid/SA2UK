# organization file for making figures

IFDIR := ${SOURCE}/figs
OFDIR := ${SINK}/figs

# overview of all raw data
${IFDIR}/epi/%.png: epi_overview.R ${RAWDATA} | ${IFDIR}/epi
	${Rstar}

# after smoothing on presumed missing data
${IFDIR}/adjusted/%.png: epi_adjusted.R ${EPIDATA} | ${IFDIR}/adjusted
	${Rstar}

${IFDIR}/intervention_timing/%.png: assess_interventions.R ${EPIDATA} ${SINK}/intervention_timing/%.rds ${SINK}/phylo.rds | ${IFDIR}/intervention_timing
	${Rstar}

epireview: $(patsubst %,${IFDIR}/epi/%.png,${ISOS})

adjreview: $(patsubst %,${IFDIR}/adjusted/%.png,${ISOS})

timing: $(patsubst %,${IFDIR}/intervention_timing/%.png,${ISOS})

MKDIRS := ${IFDIR} ${OFDIR} ${IFDIR}/intervention_timing

include ../support.makefile