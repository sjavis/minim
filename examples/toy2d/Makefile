ROOT_DIR = ../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -Wall -fopenmp -march=native -DPARALLEL

all:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) toy2d.cpp $(LIBS) -o run.exe

clean:
	rm run.exe
