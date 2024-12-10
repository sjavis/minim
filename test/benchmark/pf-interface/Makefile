ROOT_DIR = ../../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -DPARALLEL -Wall

all: run.exe
	mpirun -np 2 run.exe "1" > /dev/null
	python3 check.py "1" 0 || exit 1

plot: run.exe
	mpirun -np 2 run.exe "0.25 0.5 1 2 4" > /dev/null
	python3 check.py "0.25 0.5 1 2 4"

run.exe:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	rm -f run.exe width-*.txt
