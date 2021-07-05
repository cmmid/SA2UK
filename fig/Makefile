
TARISO ?= GHA

ALLFIGS := epi adjusted intervention_timing prefitting

MKDIRS := $(patsubst %,${FIGS}/%,${ALLFIGS})

GENITEMS := ${GEND}/intervention_timing/%.rds ${GEND}/pops/%.rds ${GEND}/yuqs/%.rds

include ../support.makefile

ISOTARS := $(patsubst %,${FIGS}/%/${TARISO}.png,${ALLFIGS})

# a directory with symbolic links to all the figure targets for a particular ISO
figures: ${FIGS}/${TARISO}

# overview of all raw data
${FIGS}/epi/%.png: epi_overview.R ${RAWDATA} | ${FIGS}/epi
	${Rstar}

# after smoothing on presumed missing data
${FIGS}/adjusted/%.png: epi_adjusted.R ${EPIDATA} | ${FIGS}/adjusted
	${Rstar}

${FIGS}/intervention_timing/%.png: assess_interventions.R ${EPIDATA} ${GEND}/intervention_timing/%.rds ${SOURCE}/phylo.rds | ${FIGS}/intervention_timing
	${Rstar}

${FIGS}/prefitting/%.png: pre_fitting.R ${EPIDATA} ${GENITEMS} ${GEND}/mobility.rds ${ESTD}/r0/%.rds | ${FIGS}/prefitting
	${Rstar}

define FIGln
ln -s $(1) $(notdir $(subst /${TARISO}.png,.png,$(1)));
endef

${FIGS}/${TARISO}: FORCE | ${ISOTARS}
	mkdir -p $@
	cd $@ && $(foreach fig,$|,$(call FIGln,${fig}))

.PRECIOUS: %.rds

# may need to send inputs dependencies back to parent
${GEND}/%:
	$(MAKE) -C .. $@