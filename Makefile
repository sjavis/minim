all:
	gcc -c src/Minimiser.cpp -o bin/Minimiser.o
	gcc -c src/GradDescent.cpp -o bin/GradDescent.o
	gcc -c src/Lbfgs.cpp -o bin/Lbfgs.o
	gcc -c src/System.cpp -o bin/System.o
	ar -rcs bin/libminimiser.a bin/Minimiser.o bin/GradDescent.o bin/Lbfgs.o bin/System.o
