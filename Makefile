CXX = g++

ROOT_DIR=<Your-Path-To-ROOT>
ROOT_CFLAGS=$(shell $(ROOT_DIR)/bin/root-config --cflags)
ROOT_LDFLAGS=$(shell $(ROOT_DIR)/bin/root-config --ldflags)
ROOT_LIBS=$(shell $(ROOT_DIR)/bin/root-config --libs)

HEPMC2_DIR=<Your-Path-To-HepMC2>
HEPMC2_INC=$(HEPMC2_DIR)/include
HEPMC2_LIB=$(HEPMC2_DIR)/lib

FASTJET_DIR=<Your-Path-To-FastJet>
FASTJET_CFLAGS=$(shell $(FASTJET_DIR)/bin/fastjet-config --cxxflags)
FASTJET_LDFLAGS=$(shell $(FASTJET_DIR)/bin/fastjet-config --libs)

CFLAGS=-I$(HEPMC2_INC) $(FASTJET_CFLAGS) $(ROOT_CFLAGS) -Wall
LDFLAGS=-L$(HEPMC2_LIB) -Wl,-rpath $(HEPMC2_LIB) $(FASTJET_LDFLAGS) $(ROOT_LDFALGS)
LIBS=-lHepMC $(ROOT_LIBS)

SRC=main.cpp Isolation.cpp Pdg_check.cpp Btag.cpp
OBJ=main.o Isolation.o Pdg_check.o Btag.o

all: $(SRC) $(OBJ) main

main: $(OBJ)
	$(CXX) $^ -o $@ $(LDFLAGS) $(LIBS)

$(OBJ): $(SRC)
	$(CXX) $(CFLAGS) -c $^

clean:
	rm *.o main
