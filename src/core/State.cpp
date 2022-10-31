#include "State.h"

#include "utils/mpi.h"
#include <stdexcept>

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


  // Energy and gradient functions using the state coordinates
  // Total energy / gradient
  double State::energy() const {
    if (pot->blockEnergyDef() && !pot->totalEnergyDef()) {
      return minim::mpi.sum(pot->blockEnergy(_coords)); // This is to avoid an unnecessary gather and scatter
    }
    return energy(comm.gather(_coords));
  }

  Vector State::gradient() const {
    if (pot->blockEnergyDef() && !pot->totalEnergyDef()) {
      return comm.gather(pot->blockGradient(_coords)); // This is to avoid an unnecessary gather and scatter
    }
    return gradient(comm.gather(_coords));
  }

  void State::energyGradient(double* e, Vector* g) const {
    if (pot->blockEnergyDef() && !pot->totalEnergyDef()) {
      // This is to avoid an unnecessary gather and scatter
      pot->blockEnergyGradient(_coords, e, g);
      if (e != nullptr) *e = mpi.sum(*e);
      if (g != nullptr) *g = comm.gather(*g);
    } else {
      energyGradient(comm.gather(_coords), e, g);
    }
  }

  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy() const {
    return blockEnergy(_coords);
  }

  Vector State::blockGradient() const {
    return blockGradient(_coords);
  }

  void State::blockEnergyGradient(double* e, Vector* g) const {
    blockEnergyGradient(_coords, e, g);
  }

  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy() const {
    return procEnergy(_coords);
  }

  Vector State::procGradient() const {
    return procGradient(_coords);
  }

  void State::procEnergyGradient(double* e, Vector* g) const {
    procEnergyGradient(_coords, e, g);
  }


  // Energy and gradient functions using given coordinates
  // Total energy / gradient
  double State::energy(const Vector& coords) const {
    if (pot->totalEnergyDef()) {
      return pot->energy(coords);
    } else if (pot->blockEnergyDef()) {
      return minim::mpi.sum(pot->blockEnergy(comm.scatter(coords)));
    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  Vector State::gradient(const Vector& coords) const {
    if (pot->totalEnergyDef()) {
      return pot->gradient(coords);
    } else if (pot->blockEnergyDef()) {
      return comm.gather(pot->blockGradient(comm.scatter(coords)));
    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::energyGradient(const Vector& coords, double* e, Vector* g) const {
    if (pot->totalEnergyDef()) {
      pot->energyGradient(coords, e, g);
    } else if (pot->blockEnergyDef()) {
      pot->blockEnergyGradient(comm.scatter(coords), e, g);
      if (e != nullptr) *e = mpi.sum(*e);
      if (g != nullptr) *g = comm.gather(*g);
    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy(const Vector& coords) const {
    if (pot->blockEnergyDef()) {
      return pot->blockEnergy(coords);
    } else if (pot->totalEnergyDef()) {
      Vector allCoords = comm.gather(coords, 0);
      if (mpi.rank != 0) {
        return 0;
      } else {
        return pot->energy(allCoords);
      }
    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  Vector State::blockGradient(const Vector& coords) const {
    if (pot->blockEnergyDef()) {
      return pot->blockGradient(coords);
    } else if (pot->totalEnergyDef()) {
      return comm.scatter(pot->gradient(comm.gather(coords)));
    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::blockEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (pot->blockEnergyDef()) {
      pot->blockEnergyGradient(coords, e, g);
    } else if (pot->totalEnergyDef()) {
      if (g == nullptr) {
        pot->energyGradient(comm.gather(coords), e, nullptr);
      } else {
        Vector gAll(ndof);
        pot->energyGradient(comm.gather(coords), e, &gAll);
        *g = comm.scatter(gAll);
      }
      if (mpi.rank != 0 && e != nullptr) *e = 0;
    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy(const Vector& coords) const {
    return blockEnergy(coords);
  }

  Vector State::procGradient(const Vector& coords) const {
    Vector g = blockGradient(coords);
    if (pot->blockEnergyDef()) comm.communicate(g);
    return g;
  }

  void State::procEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    blockEnergyGradient(coords, e, g);
    if (pot->blockEnergyDef()) comm.communicate(*g);
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
