TARGET = libminim.a
SRC_DIR = src
INC_DIR = include
BUILD_DIR = bin
LIB_DIR = lib
LIBS =
HLIBS =

CXX      = mpicxx#            C++ compiler
CXXFLAGS = -O3 -Wall -DPARALLEL -std=c++14#  Flags for the C++ compiler

TARGET := $(BUILD_DIR)/$(TARGET)
SRC_DIRS := $(shell find $(SRC_DIR) -type d)
VPATH = $(SRC_DIRS)
SRC = $(foreach sdir, $(SRC_DIRS), $(wildcard $(sdir)/*.cpp))
OBJ = $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(notdir $(SRC)))
DEPS= $(patsubst %.cpp,$(BUILD_DIR)/%.d,$(notdir $(SRC)))
INC = $(addprefix -I, $(INC_DIR))
LIB = $(patsubst %,$(BUILD_DIR)/lib%.a, $(LIBS))
HLIB= $(patsubst %,$(INC_DIR)/%, $(HLIBS))
LDLIBS = $(addprefix -l, $(LIBS))
LDFLAGS = $(addprefix -L, $(BUILD_DIR))

.PHONY: all libs clean check $(LIBS)

debug ?= 0
MAKE_ARGS = "debug=$(debug)"
ifeq ($(debug), 1)
	CXXFLAGS += -g -fno-omit-frame-pointer
else ifeq ($(debug), 2)
	CXXFLAGS += -g -O2
else ifeq ($(debug), 3)
	CXXFLAGS += -g -O0
else ifeq ($(debug), 4)
	CXXFLAGS += -g -O0 -fsanitize=address
endif


all: libs $(TARGET)

libs: $(LIBS) $(HLIB)

clean:
	$(RM) $(TARGET) $(OBJ) $(DEPS) $(LIB)

$(TARGET): $(OBJ) $(LIB)
	@echo "Making library: $@"
	@ar -rcs $@ $(OBJ)
	@echo "CREATE $@" >> tmp.mri
	@echo "ADDLIB $@" >> tmp.mri
	@for LIBFILE in $(LIB); do echo "ADDLIB $$LIBFILE" >> tmp.mri; done
	@echo "SAVE" >> tmp.mri
	@echo "END" >> tmp.mri
	@ar -M < tmp.mri
	@rm tmp.mri

-include $(DEPS)
$(OBJ): $(BUILD_DIR)/%.o: %.cpp $(LIB)
	$(CXX) $(CXXFLAGS) $(INC) -MMD -MP -c $< $(LDFLAGS) $(LDLIBS) -o $@

$(LIBS): %:
	@git submodule update --init $(LIB_DIR)/$*
	$(MAKE) --no-print-directory -C $(LIB_DIR)/$* $(MAKE_ARGS)
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
	@$(MAKE) --no-print-directory -C test/unit gtest
	@$(MAKE) --no-print-directory -C test/unit
