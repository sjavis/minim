#include <vector>
#include "State.h"
#include "Potential.h"


State::State(Potential &pot, std::vector<double> coords, Args &args)
: _pot(pot), coords(coords), args(args), ndof(coords.size())
{}


double State::energy() {
  return energy(coords);
}

double State::energy(std::vector<double> coords) {
  return _pot.energy(coords, args);
}


std::vector<double> State::gradient() {
  return gradient(coords);
}

std::vector<double> State::gradient(std::vector<double> coords) {
  return _pot.gradient(coords, args);
}


double State::operator[](int i) {
  return coords[i];
}
