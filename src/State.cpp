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
    // Initialise the potential
    this->pot->init(coords);
    this->convergence = this->pot->convergence;
    // Set-up the communicator & distribute the potential
    this->comm = this->pot->newComm();
    this->comm->setup(*this->pot, ndof, ranks);
    this->usesThisProc = comm->usesThisProc;
    this->pot->distributeParameters(*comm);
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


  void serialEnergyGradient(const Potential& pot, const vector<double>& coords, double* e, vector<double>* g) {
    pot.energyGradient(coords, e, g);
    if (g) pot.applyConstraints(coords, *g);
  }

  void parallelEnergyGradient(const Potential& pot, const vector<double>& coords, double* e, vector<double>* g, const Communicator& comm) {
    if (e) *e = 0;
    if (g) *g = vector<double>(coords.size());
    // Compute the energy elements
    for (auto el : pot.elements) {
      pot.elementEnergyGradient(coords, el, e, g);
    }
    // Compute the gradient of the halo energy elements TODO: Remove this and use procEnergyGradient instead
    if (g) {
      for (auto el : pot.elements_halo) {
        pot.elementEnergyGradient(coords, el, nullptr, g);
      }
    }
    // Compute any system-wide contributions
    pot.blockEnergyGradient(coords, comm, e, g);
    // Constraints
    if (g) pot.applyConstraints(coords, *g);
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

    if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
      parallelEnergyGradient(*pot, blockCoords, &e, nullptr, *comm);
      return comm->sum(e);

    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, &e, nullptr);
      return e;

    } else {
      throw std::logic_error("Energy function not defined");
    }
  }

  vector<double> State::gradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();
    vector<double> g;

    if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
      parallelEnergyGradient(*pot, blockCoords, nullptr, &g, *comm);
      return comm->gather(g);

    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, nullptr, &g);
      return g;

    } else {
      throw std::logic_error("Gradient function not defined");
    }
  }

  void State::energyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    if (pot->parallelDef()) {
      vector<double> blockCoords = (coords.size() == ndof) ? comm->scatter(coords) : coords;
      parallelEnergyGradient(*pot, blockCoords, e, g, *comm);
      if (e != nullptr) *e = comm->sum(*e);
      if (g != nullptr) *g = comm->gather(*g);

    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, e, g);

    } else {
      throw std::logic_error("Energy and/or gradient function not defined");
    }
  }


  // Block energy / gradient (the gradient includes the halo, but not required to be correct)
  double State::blockEnergy(const vector<double>& coords) const {
    if (!usesThisProc) return 0;

    double e;
    if (pot->parallelDef()) {
      parallelEnergyGradient(*pot, coords, &e, nullptr, *comm);
    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, &e, nullptr);
      if (comm->rank() != 0) e = 0;
    } else {
      throw std::logic_error("Energy function not defined");
    }
    return e;
  }

  vector<double> State::blockGradient(const vector<double>& coords) const {
    if (!usesThisProc) return vector<double>();

    vector<double> g;
    if (pot->parallelDef()) {
      parallelEnergyGradient(*pot, coords, nullptr, &g, *comm);
    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, nullptr, &g);
    } else {
      throw std::logic_error("Gradient function not defined");
    }
    return g;
  }

  void State::blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;

    if (pot->parallelDef()) {
      parallelEnergyGradient(*pot, coords, e, g, *comm);

    } else if (pot->serialDef()) {
      serialEnergyGradient(*pot, coords, e, g);
      if (e && comm->rank() != 0) *e = 0;

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
    if (pot->parallelDef()) comm->communicate(g);
    return g;
  }

  void State::procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const {
    if (!usesThisProc) return;
    blockEnergyGradient(coords, e, g);
    if (pot->parallelDef()) comm->communicate(*g);
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


  vector<double> State::blockCoords() const {
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


  void State::failed() {
    this->isFailed = true;
  }

}
