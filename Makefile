SRC_DIR = src/minimiser src/model src/utils
SRC = $(foreach sdir, $(SRC_DIR), $(wildcard $(sdir)/*.cpp))
OBJ = $(patsubst %.cpp,bin/%.o,$(notdir $(SRC)))
INC = $(addprefix -I, $(SRC_DIR) include)

VPATH = $(SRC_DIR)

CXX      = gcc       # C++ compiler
CXXFLAGS = $(INC) # Flags for the C++ compiler

all: $(OBJ)
	ar -rcs bin/libminim.a $(OBJ)

bin/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm bin/*.o bin/libminim.a
