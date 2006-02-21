SOURCES = ExpandSet.cpp

# Quaternion.cpp
#HEADERS = Quaternion.h
CONFIG -= qt
CONFIG += warn_on
unix:OBJECTS_DIR = .obj
unix:MOC_DIR = .moc
INCLUDEPATH	+= ../random
LANGUAGE	= C++
QMAKE_CXXFLAGS_RELEASE += -march=i686 -mtune=i686 -O2 -funroll-loops -g
#QMAKE_CXXFLAGS_RELEASE = -g
