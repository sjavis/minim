#include "State.h"

#include <vector>
#include "Potential.h"
#include "Communicator.h"
#include "utils/mpi.h"

namespace minim {

  typedef std::vector<double> Vector;


  State::State(Potential &pot, const Vector &coords, Args &args)
    : _pot(pot), args(args), ndof(coords.size()), comm(ndof,args)
  {
    setCoords(coords);
  }


  double State::energy() {
    return energy(_coords);
  }

  double State::energy(const Vector &coords) {
    return minim::mpi.sum(_pot.energy(coords, args));
  }


  Vector State::gradient() {
    return gradient(_coords);
  }

  Vector State::gradient(const Vector &coords) {
    Vector grad = _pot.gradient(coords, args);
    comm.communicate(grad);
    return grad;
  }


  double State::operator[](int i) {
    return comm.get(_coords, i);
  }


  Vector State::getCoords() {
    return comm.gather(_coords, -1);
  }

  void State::setCoords(const Vector &in) {
    _coords = comm.scatter(in, -1);
  }


  Vector State::blockCoords() {
    return _coords;
  }

  void State::blockCoords(const Vector &in) {
    _coords = in;
  }


  void State::communicate() {
    comm.communicate(_coords);
  }

}
