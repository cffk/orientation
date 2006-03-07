# $Id$

TARGETS=$(patsubst %.grid,%.quat,$(wildcard c48??.grid c48???.grid c48????.grid c48?????.grid)) $(patsubst %.grid,%.euler,$(wildcard c48??.grid c48???.grid c48?[1-3]??.grid))


all: $(TARGETS) orientgrid.tgz orientgridall.tgz orientquat.tgz orienteuler.tgz

orientgrid.tgz: $(wildcard c48??.grid c48???.grid c48????.grid c48?????.grid) $(wildcard c600*.quat) $(TARGETS) ../ExpandSet.cpp
	cd .. && tar cfz data/$@ ExpandSet.cpp data/c48?{?,??,???,????}.grid data/c600*.quat

orientgridall.tgz: $(wildcard *.grid) $(wildcard c600*.quat) $(TARGETS) ../ExpandSet.cpp
	cd .. && tar cfz data/$@ ExpandSet.cpp data/*.grid data/c600*.quat

orientquat.tgz: $(wildcard *.quat) $(TARGETS)
	cd .. && tar cfz data/$@ data/*.quat

orienteuler.tgz: $(wildcard *.euler) $(TARGETS)
	cd .. && tar cfz data/$@ data/*.euler

%.euler : %.grid
	../ExpandSet -e < $< > $@

%.quat : %.grid
	../ExpandSet < $< > $@
