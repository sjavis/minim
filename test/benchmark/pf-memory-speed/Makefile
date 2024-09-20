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

all:
	$(CXX) $(CXXFLAGS) $(INC) $(BIN) main.cpp $(LIBS) -o run.exe

clean:
	$(RM) run.exe *.txt
