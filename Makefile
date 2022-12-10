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

.PHONY: all debug deps clean check $(LIBS)

all: deps $(TARGET)

debug: CXXFLAGS += -g
debug: SUBTARGET = debug
debug: deps $(TARGET)

deps: $(LIBS) $(HLIB)

clean:
	rm $(TARGET) $(OBJ) $(LIB)

$(TARGET): $(OBJ)
	@echo "Making library: $@"
	@ar -rcs $@ $(OBJ)
	@echo "CREATE $@" >> tmp.mri
	@echo "ADDLIB $@" >> tmp.mri
	@for LIBFILE in $(LIB); do echo "ADDLIB $$LIBFILE" >> tmp.mri; done
	@echo "SAVE" >> tmp.mri
	@echo "END" >> tmp.mri
	@ar -M < tmp.mri
	@rm tmp.mri

$(OBJ): $(BUILD_DIR)/%.o: %.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(INC) -c $< $(LDFLAGS) $(LDLIBS) -o $@

$(LIBS): %:
	@git submodule update --init $(LIB_DIR)/$*
	$(MAKE) --no-print-directory -C $(LIB_DIR)/$* $(SUBTARGET)
	@ln -sfn ../$(LIB_DIR)/$*/$(INC_DIR) $(INC_DIR)/$*
	@if [ ! -f $(BUILD_DIR)/lib$*.a ] || [ $(LIB_DIR)/$*/$(BUILD_DIR)/lib$*.a -nt $(BUILD_DIR)/lib$*.a ]; then\
	  cp $(LIB_DIR)/$*/$(BUILD_DIR)/lib$*.a $(BUILD_DIR)/lib$*.a;\
	fi

$(LIB): $(BUILD_DIR)/lib%.a:
	$(MAKE) --no-print-directory $*

$(HLIB): $(INC_DIR)/%:
	git submodule update --init $(LIB_DIR)/$*
	ln -sfn ../$(LIB_DIR)/$*/include/$* $@

check:
	@echo Testing...
	@$(MAKE) --no-print-directory -C test gtest
	@$(MAKE) --no-print-directory -C test -k
