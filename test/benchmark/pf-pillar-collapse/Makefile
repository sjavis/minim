ROOT_DIR = ../../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -DPARALLEL -Wall

all: run.exe
	mpirun -np 1 run.exe > /dev/null
	python3 check.py 0 || exit 1

plot: run.exe
	mpirun -np 1 run.exe > /dev/null
	python3 check.py 1

run.exe: main.cpp
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm -rf run.exe outputs/*
