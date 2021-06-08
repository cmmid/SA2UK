
MKDIRS := ${SOURCE}/pops ${SOURCE}/yuqs ${SINK}/intervention_timing

include ../support.makefile

${SOURCE}/urbanization.rds: urbanization.R | ${SOURCE}
	${R}

${SOURCE}/matrices.rds: matrices.R | ${SOURCE}
	${R}

${SOURCE}/mortality.rds: mortality.R | ${SOURCE}
	${R}

${SOURCE}/fertility.rds: fertility.R | ${SOURCE}
	${R}

GOOGURL := https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv
OXURLS := $(patsubst %,https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/timeseries/c1_%.csv,school_closing flag)

${SOURCE}/google_mobility.csv:
	${WGET} ${GOOGURL}

${SOURCE}/ox_si_schools.csv:
	${WGET} $(firstword ${OXURLS})

${SOURCE}/ox_si_schools_flag.csv:
	${WGET} $(lastword ${OXURLS})

MOBDATA := $(patsubst %,${SOURCE}/%.csv,google_mobility ox_si_schools ox_si_schools_flag)
rawmob: ${MOBDATA}

${SOURCE}/mobility.rds: mobility.R ${MOBDATA} | ${SOURCE}
	${RSCRIPT} $^ ${MOBURLS} $@

REFDATA := $(patsubst %,${SOURCE}/%.rds,mortality fertility urbanization matrices mobility)

${SINK}/intervention_timing/%.rds: r0_est_timing.R ${PROJRT}/interventions.csv | ${SINK}/intervention_timing
	${Rstar}

${SOURCE}/pops/%.rds: covidm_pop.R ${REFDATA} | ${COVIDM} ${SOURCE}/pops
	${RSCRIPT} $^ $* $(firstword $|) $@

${SOURCE}/yuqs/%.rds: reference_qs.R ${SOURCE}/covidm_fit_yu.qs ${SOURCE}/pops/%.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/mobility.rds | ${COVIDM} ${SOURCE}/yuqs
	${RSCRIPT} $^ $* $(firstword $|) $@

allpops: $(call alliso,${SOURCE}/pops,rds)
allyuqs: $(call alliso,${SOURCE}/yuqs,rds)
allgen: allpops allyuqs

.PRECIOUS: %.rds

# may need to send dependencies back to parent - can't do this, since other dependencies?
# %.rds:
#	$(MAKE) -C .. $@