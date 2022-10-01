#include "State.h"

#include <vector>
#include "Potential.h"
#include "Communicator.h"
#include "utils/mpi.h"

namespace minim {

  typedef std::vector<double> Vector;


  State::State(Potential& pot, const Vector& coords, std::unique_ptr<Potential::Args>& args)
    : ndof(coords.size()), comm(ndof,*args), args(std::move(args)), _pot(pot)
  {
    setCoords(coords);
  }

  State::State(const State& state)
    : ndof(state.ndof),
      convergence(state.convergence),
      comm(state.comm),
      args(state.args->clone()),
      _istart(state._istart),
      _iend(state._iend),
      _coords(state._coords),
      _pot(state._pot)
  {}

  State& State::operator=(const State& state) {
    ndof = state.ndof;
    convergence = state.convergence;
    comm = state.comm;
    args = state.args->clone();
    _istart = state._istart;
    _iend = state._iend;
    _coords = state._coords;
    _pot = state._pot;
    return *this;
  }


  double State::energy() const {
    return energy(_coords);
  }

  double State::energy(const Vector &coords) const {
    return minim::mpi.sum(_pot.energy(coords, *args));
  }


  Vector State::gradient() const {
    return gradient(_coords);
  }

  Vector State::gradient(const Vector &coords) const {
    Vector grad = _pot.gradient(coords, *args);
    comm.communicate(grad);
    return grad;
  }


  double State::operator[](int i) {
    return comm.get(_coords, i);
  }


  Vector State::getCoords() const {
    return comm.gather(_coords, -1);
  }

  void State::setCoords(const Vector &in) {
    _coords = comm.scatter(in, -1);
  }


  Vector State::blockCoords() const {
    return _coords;
  }

  void State::blockCoords(const Vector &in) {
    _coords = in;
  }


  void State::communicate() {
    comm.communicate(_coords);
  }

}
