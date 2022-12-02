TARGET = libminim.a
SRC_DIR = src/core src/minimisers src/potentials src/utils
INC_DIR = include
BUILD_DIR = bin
LIB_DIR = lib
LIBS =
HLIBS = gsl

CXX      = mpicxx#            C++ compiler
CXXFLAGS = -Wall -DPARALLEL#  Flags for the C++ compiler

TARGET := $(BUILD_DIR)/$(TARGET)
VPATH = $(SRC_DIR)
SRC = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJ = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(SRC)))
INC = $(addprefix -I, $(INC_DIR))
LIB = $(patsubst %,$(BUILD_DIR)/lib%.a, $(LIBS))
HLIB = $(patsubst %,$(INC_DIR)/%, $(HLIBS))
LDLIBS = $(addprefix -l, $(LIBS))
LDFLAGS = $(addprefix -L, $(BUILD_DIR))

.PHONY: all debug deps clean check $(LIB)

all: $(TARGET)

debug: CXXFLAGS+=-g
debug: $(TARGET)

deps: $(LIB) $(HLIB)

clean:
	rm $(TARGET) $(OBJ) $(LIB)

$(TARGET): $(OBJ)
	ar -rcs $(TARGET) $(OBJ)

$(OBJ): $(BUILD_DIR)/%.o: %.cpp $(LIB) $(HLIB)
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LDFLAGS) $(LDLIBS) -o $@

$(LIB): $(BUILD_DIR)/lib%.a:
	git submodule update --init $(LIB_DIR)/$*
	$(MAKE) -C $(LIB_DIR)/$*
	ln -sfn ../$(LIB_DIR)/$*/$(INC_DIR) $(INC_DIR)/$*
	cp $(LIB_DIR)/$*/$@ $@

$(HLIB): $(INC_DIR)/%:
	git submodule update --init $(LIB_DIR)/$*
	ln -sfn ../$(LIB_DIR)/$*/include/$* $@

check:
	$(MAKE) -C test gtest
	$(MAKE) -C test
