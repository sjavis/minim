# Minim

## Overview
A C++ software library containing energy minimisation subroutines to find equilibium states of generic system potentials.
These can utilise MPI parallelisation using built-in frameworks.
This library can be used in conjuction with [ELLib](https://github.com/sjavis/ellib) to perform more complex energy landscape methods.

## Installation and use
Download the repository and call `make` in the root directory to compile the library.

To use in a program, include the `minim.h` header file and compile with the `-lminim` flag.
For example: `mpic++ -I$(MINIM)/include -L$(MINIM)/bin -lminim -DPARALLEL script.cpp -o run.exe`

Refer to the examples folder for simple demonstrations of how to use the library.

## Library structure

This library is split into several core components:
- `Potential` classes provide the interface for calculating the energy and
gradient for a given set of coordinates, and includes any potential specific parameters.
- `Minimiser` classes are used to perform the energy minimisation.
- The `State` class is used to create systems consisting of a
potential and a set of coordinates. These are the objects acted upon by the minimisation methods.
- `Communicator` classes are used to abstract away the MPI communication used in each `State` object.

`Potential` classes:
- `LjNd`: 2D and 3D Lennard-Jones particle potential.
- `PhaseField`: A phase-field potential for multicomponent fluid systems.
- `BarAndHinge`: A triangular mesh bar-and-hinge potential for simulating elastic surfaces.

`Minimiser` classes:
- `Lbfgs`: L-BFGS
- `Fire`: FIRE
- `GradDescent`: Gradient descent
- `Anneal`: Simulated annealing
