ROOT_DIR = ../../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -DPARALLEL -Wall

all: run.exe
	mpirun -np 1 run.exe 1 > /dev/null || exit 1

report: run.exe
	mpirun -np 1 run.exe

run.exe:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm -f run.exe
