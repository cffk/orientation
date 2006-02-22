# $Id$

TARGETS=$(patsubst %.grid,%.quat,$(wildcard *.grid)) $(patsubst %.grid,%.euler,$(wildcard c48??.grid c48???.grid c48?[1-3]??.grid))


all: $(TARGETS) orient.tgz

orient.tgz: $(wildcard *.grid) $(wildcard *.quat) $(wildcard *.euler) $(TARGETS)
	cd .. && tar cfz data/orient.tgz data/*.grid data/*.quat data/*.euler

%.euler : %.grid
	../ExpandSet -e < $< > $@

%.quat : %.grid
	../ExpandSet < $< > $@
