ROOT_DIR = ../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -Wall -fopenmp -DPARALLEL

all:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm run.exe *.txt
