TARGET = libminim.a
SRC_DIR = src/core src/minimisers src/potentials src/utils
INC_DIR = include
BUILD_DIR = bin

CXX      = mpic++             # C++ compiler
CXXFLAGS = -Wall -DPARALLEL   # Flags for the C++ compiler

VPATH = $(SRC_DIR)
SRC = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJ = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(SRC)))
INC = $(addprefix -I, $(INC_DIR))
TARGET := $(BUILD_DIR)/$(TARGET)

.PHONY: all clean

all: $(TARGET)

clean:
	rm $(TARGET) $(OBJ)

$(TARGET): $(OBJ)
	ar -rcs $(TARGET) $(OBJ)

$(OBJ): bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@
