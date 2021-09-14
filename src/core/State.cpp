#include <vector>
#include "State.h"
#include "Potential.h"
#include "Communicator.h"
#include "minimMpi.h"


State::State(Potential &pot, const std::vector<double> &coords, Args &args)
  : _pot(pot), args(args), ndof(coords.size()), comm(ndof,args)
{
  setCoords(coords);
}


double State::energy() {
  return energy(_coords);
}

double State::energy(const std::vector<double> &coords) {
  return minim::mpi.sum(_pot.energy(coords, args));
}


std::vector<double> State::gradient() {
  return gradient(_coords);
}

std::vector<double> State::gradient(const std::vector<double> &coords) {
  std::vector<double> grad = _pot.gradient(coords, args);
  comm.communicate(grad);
  return grad;
}


double State::operator[](int i) {
  return comm.get(_coords, i);
}


std::vector<double> State::getCoords() {
  return comm.gather(_coords, -1);
}

void State::setCoords(const std::vector<double> &in) {
  _coords = comm.scatter(in, -1);
}


std::vector<double> State::blockCoords() {
  return _coords;
}

void State::blockCoords(const std::vector<double> &in) {
  _coords = in;
}


void State::communicate() {
  comm.communicate(_coords);
}
