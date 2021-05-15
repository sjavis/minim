ROOT_DIR = ../..
INC = $(ROOT_DIR)/include
BIN = $(ROOT_DIR)/bin
LIBS = -lminim

CXX = g++
CXX = mpic++
CXXFLAGS = -I$(INC) -L$(BIN) $(LIBS) -DPARALLEL

all:
	$(CXX) lj3d.cpp -o run.exe $(CXXFLAGS)

clean:
	rm run.exe
