
# can define other variables here as well...
DBROOT := ~/Dropbox
DBPTH := ${DBROOT}/covidLMIC

# ...which can be used when overriding defaults
# n.b.: these paths cannot be identical, but this is *not*
# enforced by the Makefile
SOURCE := ${DBPTH}/inputs
SINK := ${DBPTH}/outputs