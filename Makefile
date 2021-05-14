SRC_DIR = src/core src/minimisers src/potentials src/utils
SRC = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJ = $(patsubst %.cpp,bin/%.o,$(notdir $(SRC)))
INC = $(addprefix -I, $(SRC_DIR) include)

VPATH = $(SRC_DIR)

CXX      = mpic++       # C++ compiler
CXXFLAGS = $(INC) -DPARALLEL # Flags for the C++ compiler

all: $(OBJ)
	ar -rcs bin/libminim.a $(OBJ)

bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm bin/*.o bin/libminim.a
