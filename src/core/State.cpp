#include "State.h"

#include <stdexcept>

namespace minim {

  class Potential;
  typedef std::vector<double> Vector;


  State::State(const Potential& pot, const Vector& coords, const std::vector<int>& ranks)
    : ndof(coords.size()), pot(pot.clone()), comm(*this->pot,ndof,ranks)
  {
    this->usesThisProc = comm.usesThisProc;
    this->coords(coords);
  }

  State::State(const State& state)
    : ndof(state.ndof),
      convergence(state.convergence),
      pot(state.pot->clone()),
      comm(state.comm),
      usesThisProc(state.usesThisProc),
      _coords(state._coords)
  {}

  State& State::operator=(const State& state) {
    ndof = state.ndof;
    convergence = state.convergence;
    pot = state.pot->clone();
    comm = state.comm;
    usesThisProc = state.usesThisProc;
    _coords = state._coords;
    return *this;
  }


  // Energy and gradient functions using the state coordinates
  // Total energy / gradient
  double State::energy() const {
    if (!usesThisProc) return 0;
    return energy(_coords);
  }

  Vector State::gradient() const {
    if (!usesThisProc) return Vector();
    return gradient(_coords);
  }

  void State::energyGradient(double* e, Vector* g) const {
    if (!usesThisProc) return;
    energyGradient(_coords, e, g);
  }

  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy() const {
    if (!usesThisProc) return 0;
    return blockEnergy(_coords);
  }

  Vector State::blockGradient() const {
    if (!usesThisProc) return Vector();
    return blockGradient(_coords);
  }

  void State::blockEnergyGradient(double* e, Vector* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(_coords, e, g);
  }

  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy() const {
    if (!usesThisProc) return 0;
    return procEnergy(_coords);
  }

  Vector State::procGradient() const {
    if (!usesThisProc) return Vector();
    return procGradient(_coords);
  }

  void State::procEnergyGradient(double* e, Vector* g) const {
    if (!usesThisProc) return;
    procEnergyGradient(_coords, e, g);
  }


  // Energy and gradient functions using given coordinates
  // Total energy / gradient
  double State::energy(const Vector& coords) const {
    if (!usesThisProc) return 0;
    if (pot->totalEnergyDef()) {
      Vector allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      return pot->energy(allCoords);

    } else if (pot->blockEnergyDef()) {
      Vector blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      return comm.sum(pot->blockEnergy(blockCoords, comm));

    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  Vector State::gradient(const Vector& coords) const {
    if (!usesThisProc) return Vector();

    if (pot->totalEnergyDef()) {
      Vector allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      return pot->gradient(allCoords);

    } else if (pot->blockEnergyDef()) {
      Vector blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      return comm.gather(pot->blockGradient(blockCoords, comm));

    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::energyGradient(const Vector& coords, double* e, Vector* g) const {
    if (!usesThisProc) return;

    if (pot->totalEnergyDef()) {
      Vector allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      pot->energyGradient(allCoords, e, g);

    } else if (pot->blockEnergyDef()) {
      Vector blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      pot->blockEnergyGradient(blockCoords, comm, e, g);
      if (e != nullptr) *e = comm.sum(*e);
      if (g != nullptr) *g = comm.gather(*g);

    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy(const Vector& coords) const {
    if (!usesThisProc) return 0;

    if (pot->blockEnergyDef()) {
      return pot->blockEnergy(coords, comm);

    } else if (pot->totalEnergyDef()) {
      Vector allCoords = comm.gather(coords, 0);
      if (comm.rank() != 0) {
        return 0;
      } else {
        return pot->energy(allCoords);
      }

    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  Vector State::blockGradient(const Vector& coords) const {
    if (!usesThisProc) return Vector();

    if (pot->blockEnergyDef()) {
      return pot->blockGradient(coords, comm);

    } else if (pot->totalEnergyDef()) {
      return comm.scatter(pot->gradient(comm.gather(coords)));

    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::blockEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (!usesThisProc) return;

    if (pot->blockEnergyDef()) {
      pot->blockEnergyGradient(coords, comm, e, g);

    } else if (pot->totalEnergyDef()) {
      if (g == nullptr) {
        pot->energyGradient(comm.gather(coords), e, nullptr);
      } else {
        Vector gAll(ndof);
        pot->energyGradient(comm.gather(coords), e, &gAll);
        *g = comm.scatter(gAll);
      }
      if (comm.rank() != 0 && e != nullptr) *e = 0;

    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy(const Vector& coords) const {
    if (!usesThisProc) return 0;
    return blockEnergy(coords);
  }

  Vector State::procGradient(const Vector& coords) const {
    if (!usesThisProc) return Vector();
    Vector g = blockGradient(coords);
    if (pot->blockEnergyDef()) comm.communicate(g);
    return g;
  }

  void State::procEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(coords, e, g);
    if (pot->blockEnergyDef()) comm.communicate(*g);
  }


  double State::operator[](int i) {
    return comm.get(_coords, i);
  }


  Vector State::coords() const {
    return comm.gather(_coords, -1);
  }

  void State::coords(const Vector& in) {
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
