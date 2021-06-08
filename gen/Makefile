
MKDIRS := ${SOURCE}/pops ${SOURCE}/yuqs

include ../support.makefile

${SOURCE}/urbanization.rds: urbanization.R | ${SOURCE}
	${R}

${SOURCE}/matrices.rds: matrices.R | ${SOURCE}
	${R}

${SOURCE}/mortality.rds: mortality.R | ${SOURCE}
	${R}

${SOURCE}/fertility.rds: fertility.R | ${SOURCE}
	${R}

${SOURCE}/pops/%.rds: covidm_pop.R $(patsubst %,${SOURCE}/%.rds,mortality fertility urbanization matrices) | ${COVIDM} ${SOURCE}/pops
	${RSCRIPT} $^ $* $(firstword $|) $@

${SOURCE}/yuqs/%.rds: reference_qs.R ${SOURCE}/covidm_fit_yu.qs ${SOURCE}/pops/%.rds ${SINK}/intervention_timing/%.rds ${SOURCE}/mobility.rds | ${COVIDM} ${SOURCE}/yuqs
	${RSCRIPT} $^ $* $(firstword $|) $@

allpops: $(call alliso,${SOURCE}/pops,rds)
allyuqs: $(call alliso,${SOURCE}/yuqs,rds)

allgen: allpops allyuqs

.PRECIOUS: %.rds

# may need to send dependencies back to parent
%.rds:
	$(MAKE) -C .. $@