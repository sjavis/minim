OBJ_ = Minimiser.o GradDescent.o Lbfgs.o System.o vec.o
OBJ = $(OBJ_:%=bin/%)

CXX      = gcc       # C++ compiler
CXXFLAGS = -Iinclude # Flags for the C++ compiler

all: $(OBJ)
	ar -rcs bin/libminim.a $(OBJ)

bin/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm bin/*.o bin/libminim.a
