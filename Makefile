CXXFLAGS   = -g -std=c++11 -O2 -Wall -fPIC -lm
ROOT = $(shell $(ROOTSYS)/bin/root-config --glibs) $(shell $(ROOTSYS)/bin/root-config --cflags)

all: coulex_ang_corr

coulex_ang_corr: coulex_ang_corr.cpp coulex_ang_corr.h read_config.cpp build_matrices.cpp
	g++ $(CXXFLAGS) coulex_ang_corr.cpp  $(ROOT) -o coulex_ang_corr
clean:
	rm -rf *~ *.o coulex_ang_corr *tmpdatafile*