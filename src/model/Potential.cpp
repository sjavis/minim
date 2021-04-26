#include "Potential.h"
#include "State.h"
#include <iostream>

typedef std::vector<double> Vector;


double Potential::energy(Vector coords, Args &args) {
  std::cout << "Lj3d" << std::endl;
  return _energy(coords, args);
}


Vector Potential::gradient(Vector coords, Args &args) {
  return _gradient(coords, args);
}


State Potential::newState(Vector coords) {
  Args *args = new Args(coords.size());
  return State(*this, coords, *args);
}

State Potential::newState(int ndof) {
  Vector coords(ndof);
  return newState(coords);
}
