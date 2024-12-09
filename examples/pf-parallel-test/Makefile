ROOT_DIR = ../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpicxx
CXXFLAGS = -DPARALLEL -Wall -O2 -g
# CXXFLAGS = -Wall -O2 -g

all:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm run.exe
