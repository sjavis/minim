ROOT_DIR = ../../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -Wall -fopenmp -march=native -DPARALLEL

all: run.exe
	mpirun -np 2 run.exe 1 > /dev/null || exit 1

report: run.exe
	mpirun -np 2 run.exe

run.exe:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm -f run.exe
