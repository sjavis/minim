#include "Potential.h"
#include "State.h"

typedef std::vector<double> Vector;


double Potential::energy(Vector coords, Args &args) {
  return _energy(coords, args);
}


Vector Potential::gradient(Vector coords, Args &args) {
  return _gradient(coords, args);
}


Args* Potential::newArgs(int ndof) {
  return new Args(ndof);
}


State Potential::newState(Vector coords) {
  Args *args = newArgs(coords.size());
  return State(*this, coords, *args);
}

State Potential::newState(int ndof) {
  Vector coords(ndof);
  return newState(coords);
}
