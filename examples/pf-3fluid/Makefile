ROOT_DIR = ../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -DPARALLEL -Wall

all:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) 3fluid.cpp $(LIBS) -o run.exe

clean:
	rm run.exe
