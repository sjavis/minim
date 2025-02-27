#include "State.h"

#include <stdexcept>
#include "utils/mpi.h"


namespace minim {
  using std::vector;
  class Potential;


  State::State(const Potential& pot, const vector<double>& coords, const vector<int>& ranks)
    : ndof(coords.size()), pot(pot.clone())
  {
    // Initialise the potential (globally)
    this->pot->init(coords);
    this->convergence = this->pot->convergence;
    // Set-up the communicator
    this->comm = this->pot->newComm();
    this->comm->setup(*this->pot, ndof, ranks);
    this->usesThisProc = comm->usesThisProc;
    // Initialise the potential (locally)
    this->pot->initLocal(coords, *comm);
    // Initialise the coords
    this->coords(coords);
  }

  State::State(const State& state)
    : ndof(state.ndof),
      convergence(state.convergence),
      pot(state.pot->clone()),
      comm(state.comm->clone()),
      usesThisProc(state.usesThisProc),
      _coords(state._coords)
  {}

  State& State::operator=(const State& state) {
    ndof = state.ndof;
    convergence = state.convergence;
    pot = state.pot->clone();
    comm = state.comm->clone();
    usesThisProc = state.usesThisProc;
    _coords = state._coords;
    return *this;
  }


  inline void basicEG(const Potential& pot, const vector<double>& coords, double* e, vector<double>* g, const Communicator& comm) {
    if (e) *e = 0;
    if (g) *g = vector<double>(coords.size());
    pot.energyGradient(coords, comm, e, g);
    if (g) {
      if (comm.size() > 1) comm.communicateAccumulate(*g); // Get correct gradient on the edges
      pot.applyConstraints(coords, comm, *g);
    }
  }

  inline void elementEG(const Potential& pot, const vector<double>& coords, double* e, vector<double>* g, const Communicator& comm) {
    if (e) *e = 0;
    if (g) *g = vector<double>(coords.size());
    // Compute the energy elements
    for (auto el : pot.elements) {
      pot.elementEnergyGradient(coords, el, e, g);
    }
    // Compute any system-wide contributions
    pot.blockEnergyGradient(coords, comm, e, g);
    if (!g) return;
    // Get the correct gradient on the edges (not halo)
    if (comm.size()>1) {
      // A: By communication
      comm.communicateAccumulate(*g);
      // // B: By computing the gradient of the halo energy elements
      // for (auto el : pot.elements_halo) {
      //   pot.elementEnergyGradient(coords, el, nullptr, g);
      // }
    }
    // Constraints
    pot.applyConstraints(coords, comm, *g);
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
    double e;

    // Serial
    if (pot->potentialType() == Potential::SERIAL) {
      basicEG(*pot, coords, &e, nullptr, *comm);
      return e;
    }

    // Parallel
    const vector<double>& blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
    if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, blockCoords, &e, nullptr, *comm);
    } else {
      basicEG(*pot, blockCoords, &e, nullptr, *comm);
    }
    return comm->sum(e);
  }

  vector<double> State::gradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();
    vector<double> g;

    // Serial
    if (pot->potentialType() == Potential::SERIAL) {
      basicEG(*pot, coords, nullptr, &g, *comm);
      return g;
    }

    // Parallel
    const vector<double>& blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
    if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, blockCoords, nullptr, &g, *comm);
    } else {
      basicEG(*pot, blockCoords, nullptr, &g, *comm);
    }
    return comm->gather(g);
  }

  void State::energyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    // Serial
    if (pot->potentialType() == Potential::SERIAL) {
      basicEG(*pot, coords, e, g, *comm);
      return;
    }

    // Parallel
    const vector<double>& blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
    if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, blockCoords, e, g, *comm);
    } else {
      basicEG(*pot, blockCoords, e, g, *comm);
    }
    if (e != nullptr) *e = comm->sum(*e);
    if (g != nullptr) *g = comm->gather(*g);
  }


  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy(const vector<double>& coords) const {
    if (!usesThisProc) return 0;

    double e;
    if (pot->potentialType() == Potential::SERIAL) {
      if (comm->rank() != 0) return 0;
      basicEG(*pot, coords, &e, nullptr, *comm);
    } else if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, coords, &e, nullptr, *comm);
    } else {
      basicEG(*pot, coords, &e, nullptr, *comm);
    }
    return e;
  }

  vector<double> State::blockGradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();

    vector<double> g;
    if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, coords, nullptr, &g, *comm);
    } else {
      basicEG(*pot, coords, nullptr, &g, *comm);
    }
    return g;
  }

  void State::blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    if (pot->potentialType() == Potential::SERIAL) {
      basicEG(*pot, coords, e, g, *comm);
      if (e && comm->rank() != 0) *e = 0;
    } else if (pot->potentialType() == Potential::UNSTRUCTURED) {
      elementEG(*pot, coords, e, g, *comm);
    } else {
      basicEG(*pot, coords, e, g, *comm);
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
    comm->communicate(g);
    return g;
  }

  void State::procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(coords, e, g);
    comm->communicate(*g);
  }


  // Get energy / gradient on all procs, including those not used in the state
  double State::allEnergy() const {
    double e = energy();
    mpi.bcast(e, comm->ranks[0]);
    return e;
  }

  vector<double> State::allGradient() const {
    vector<double> g = gradient();
    mpi.bcast(g, comm->ranks[0]);
    return g;
  }

  void State::allEnergyGradient(double* e, vector<double>* g) const {
    energyGradient(e, g);
    mpi.bcast(*e, comm->ranks[0]);
    mpi.bcast(*g, comm->ranks[0]);
  }


  double State::operator[](int i) {
    return comm->get(_coords, i);
  }


  vector<double> State::coords() const {
    return comm->gather(_coords, -1);
  }

  void State::coords(const vector<double>& in) {
    _coords = comm->scatter(in, -1);
  }


  const vector<double>& State::blockCoords() const {
    return _coords;
  }

  void State::blockCoords(const vector<double>& in) {
    _coords = in;
  }


  vector<double> State::allCoords() const {
    vector<double> coords = comm->gather(_coords, 0);
    mpi.bcast(coords, comm->ranks[0], ndof);
    return coords;
  }


  double State::componentEnergy(int component) const {
    double e = 0;
    for (auto el : pot->elements) {
      if (el.type == component) pot->elementEnergyGradient(_coords, el, &e, nullptr);
    }
    return e;
  }


  void State::communicate() {
    comm->communicate(_coords);
  }


  void State::applyConstraints(vector<double>& data) const {
    pot->applyConstraints(_coords, *comm, data);
  }


  void State::failed() {
    this->isFailed = true;
  }

}
