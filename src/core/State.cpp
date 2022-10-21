#include "State.h"

#include "utils/mpi.h"

namespace minim {

  class Potential;
  typedef std::vector<double> Vector;


  State::State(const Potential& pot, const Vector& coords)
    : ndof(coords.size()), pot(pot.clone()), comm(ndof,*this->pot)
  {
    setCoords(coords);
  }

  State::State(const State& state)
    : ndof(state.ndof),
      convergence(state.convergence),
      pot(state.pot->clone()),
      comm(state.comm),
      _istart(state._istart),
      _iend(state._iend),
      _coords(state._coords)
  {}

  State& State::operator=(const State& state) {
    ndof = state.ndof;
    convergence = state.convergence;
    pot = state.pot->clone();
    comm = state.comm;
    _istart = state._istart;
    _iend = state._iend;
    _coords = state._coords;
    return *this;
  }


  double State::energy() const {
    return energy(_coords);
  }

  double State::energy(const Vector& coords) const {
    if (pot->totalEnergyDef()) {
      return pot->energy(coords);
    } else if (pot->blockEnergyDef()) {
      return minim::mpi.sum(pot->blockEnergy(coords));
    } else {
      throw std::logic_error("Energy function not defined");
    }
  }


  Vector State::gradient() const {
    return gradient(_coords);
  }

  Vector State::gradient(const Vector& coords) const {
    if (pot->totalEnergyDef()) {
      return pot->gradient(coords);
    } else if (pot->blockEnergyDef()) {
      Vector grad = pot->gradient(coords);
      comm.communicate(grad); // TODO: Check / fix this
      return grad;
    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }


  void State::energyGradient(double* e, Vector* g) const {
    return energyGradient(_coords, e, g);
  }

  void State::energyGradient(const Vector& coords, double* e, Vector* g) const {
    if (pot->totalEnergyDef()) {
      pot->energyGradient(coords, e, g);
    } else if (pot->blockEnergyDef()) {
      pot->blockEnergyGradient(coords, e, g); // TODO: Accumulate e and g
    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  double State::operator[](int i) {
    return comm.get(_coords, i);
  }


  Vector State::getCoords() const {
    return comm.gather(_coords, -1);
  }

  void State::setCoords(const Vector& in) {
    _coords = comm.scatter(in, -1);
  }


  Vector State::blockCoords() const {
    return _coords;
  }

  void State::blockCoords(const Vector& in) {
    _coords = in;
  }


  void State::communicate() {
    comm.communicate(_coords);
  }

}
