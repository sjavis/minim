#include "State.h"

#include <stdexcept>
#include "utils/mpi.h"

namespace minim {
  using std::vector;
  class Potential;


  State::State(const Potential& pot, const vector<double>& coords, const vector<int>& ranks)
    : ndof(coords.size()), pot(pot.clone())
  {
    if (this->pot->distributed) throw std::invalid_argument("You cannot create a State with a Potential that has already been distributed.");
    this->pot->init(coords);
    this->comm.setup(*this->pot, ndof, ranks);
    this->usesThisProc = comm.usesThisProc;
    this->coords(coords);
    this->pot->distributeParameters(comm);
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

  vector<double> State::gradient() const {
    if (!usesThisProc) return vector<double>();
    return gradient(_coords);
  }

  void State::energyGradient(double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    energyGradient(_coords, e, g);
  }

  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy() const {
    if (!usesThisProc) return 0;
    return blockEnergy(_coords);
  }

  vector<double> State::blockGradient() const {
    if (!usesThisProc) return vector<double>();
    return blockGradient(_coords);
  }

  void State::blockEnergyGradient(double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(_coords, e, g);
  }

  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy() const {
    if (!usesThisProc) return 0;
    return procEnergy(_coords);
  }

  vector<double> State::procGradient() const {
    if (!usesThisProc) return vector<double>();
    return procGradient(_coords);
  }

  void State::procEnergyGradient(double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    procEnergyGradient(_coords, e, g);
  }


  // Energy and gradient functions using given coordinates
  // Total energy / gradient
  double State::energy(const vector<double>& coords) const {
    if (!usesThisProc) return 0;
    if (pot->serialDef()) {
      vector<double> allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      return pot->energy(allCoords);

    } else if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      return comm.sum(blockEnergy(blockCoords));

    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  vector<double> State::gradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();

    if (pot->serialDef()) {
      vector<double> allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      return pot->gradient(allCoords);

    } else if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      return comm.gather(blockGradient(blockCoords));

    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::energyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    if (pot->serialDef()) {
      vector<double> allCoords = (coords.size() == ndof) ? coords : comm.gather(coords);
      pot->energyGradient(allCoords, e, g);

    } else if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm.scatter(coords) : coords;
      blockEnergyGradient(blockCoords, e, g);
      if (e != nullptr) *e = comm.sum(*e);
      if (g != nullptr) *g = comm.gather(*g);

    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy(const vector<double>& coords) const {
    if (!usesThisProc) return 0;

    if (pot->parallelDef()) {
      double e;
      blockEnergyGradient(coords, &e, nullptr);
      return e;

    } else if (pot->serialDef()) {
      vector<double> allCoords = comm.gather(coords, 0);
      if (comm.rank() != 0) {
        return 0;
      } else {
        return pot->energy(allCoords);
      }

    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  vector<double> State::blockGradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();

    if (pot->parallelDef()) {
      vector<double> g;
      blockEnergyGradient(coords, nullptr, &g);
      return g;

    } else if (pot->serialDef()) {
      return comm.scatter(pot->gradient(comm.gather(coords)));

    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    if (pot->parallelDef()) {
      if (e) *e = 0;
      if (g) *g = vector<double>(coords.size());
      // Compute any system-wide contributions
      pot->blockEnergyGradient(coords, comm, e, g);
      // Compute the energy elements
      for (auto el : pot->elements) {
        pot->elementEnergyGradient(coords, el, e, g);
      }
      // Compute the gradient of the halo energy elements TODO: Remove this and use procEnergyGradient instead
      if (g) {
        for (auto el : pot->elements_halo) {
          pot->elementEnergyGradient(coords, el, nullptr, g);
        }
      }

    } else if (pot->serialDef()) {
      if (g == nullptr) {
        pot->energyGradient(comm.gather(coords), e, nullptr);
      } else {
        vector<double> gAll(ndof);
        pot->energyGradient(comm.gather(coords), e, &gAll);
        *g = comm.scatter(gAll);
      }
      if (comm.rank() != 0 && e != nullptr) *e = 0;

    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Processor energy / gradient (the gradient includes the halo)
  double State::procEnergy(const vector<double>& coords) const {
    if (!usesThisProc) return 0;
    return blockEnergy(coords);
  }

  vector<double> State::procGradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();
    vector<double> g = blockGradient(coords);
    if (pot->parallelDef()) comm.communicate(g);
    return g;
  }

  void State::procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(coords, e, g);
    if (pot->parallelDef()) comm.communicate(*g);
  }


  // Get energy / gradient on all procs, including those not used in the state
  double State::allEnergy() const {
    double e = energy();
    mpi.bcast(e, comm.ranks[0]);
    return e;
  }

  vector<double> State::allGradient() const {
    vector<double> g = gradient();
    mpi.bcast(g, comm.ranks[0]);
    return g;
  }

  void State::allEnergyGradient(double* e, vector<double>* g) const {
    energyGradient(e, g);
    mpi.bcast(*e, comm.ranks[0]);
    mpi.bcast(*g, comm.ranks[0]);
  }


  double State::operator[](int i) {
    return comm.get(_coords, i);
  }


  vector<double> State::coords() const {
    return comm.gather(_coords, -1);
  }

  void State::coords(const vector<double>& in) {
    _coords = comm.scatter(in, -1);
  }


  vector<double> State::blockCoords() const {
    return _coords;
  }

  void State::blockCoords(const vector<double>& in) {
    _coords = in;
  }


  vector<double> State::allCoords() const {
    vector<double> coords = comm.gather(_coords, 0);
    mpi.bcast(coords, comm.ranks[0], ndof);
    return coords;
  }


  void State::communicate() {
    comm.communicate(_coords);
  }

}
