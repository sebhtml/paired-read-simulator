CXX=g++
CXXFLAGS=-O3 -Wall -std=c++98

all: VirtualNextGenSequencer

VirtualNextGenSequencer: simulatePairedReads.cpp
	@$(CXX) -o $@ $< $(CXXFLAGS)
	@echo "  CXX $<"


