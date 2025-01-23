ROOT_DIR = ../../..
INC = -I$(ROOT_DIR)/include
BIN = -L$(ROOT_DIR)/bin
LIBS = -lminim

CXX = mpic++
CXXFLAGS = -DPARALLEL -Wall

debug ?= 0
ifeq ($(debug), 1)
	CXXFLAGS += -g -fno-omit-frame-pointer
else ifeq ($(debug), 2)
	CXXFLAGS += -g -O0
else ifeq ($(debug), 3)
	CXXFLAGS += -g -O0 -fsanitize=address
endif

all: run.exe
	mpirun -np 1 run.exe
	mpirun -np 2 run.exe
	mpirun -np 4 run.exe
	mpirun -np 5 run.exe
	mpirun -np 8 run.exe
	python3 check.py

run.exe:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	$(RM) run.exe *.txt
