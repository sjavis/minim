#include "Anneal.h"

#include <math.h>
#include <time.h>
#include "State.h"

namespace minim {

  Anneal& Anneal::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }


  Anneal& Anneal::setDisplacement(double displacement) {
    this->displacement = displacement;
    return *this;
  }

  Anneal& Anneal::setTempInit(double tempInit) {
    this->tempInit = tempInit;
    return *this;
  }

  Anneal& Anneal::setCoolingRate(double coolingRate) {
    this->coolingRate = coolingRate;
    return *this;
  }

  Anneal& Anneal::setCoolingSchedule(std::function<double(int)> coolingSchedule) {
    this->coolingSchedule = coolingSchedule;
    return *this;
  }

  Anneal& Anneal::setMaxRejections(int maxRejections) {
    this->maxRejections = maxRejections;
    return *this;
  }


  void Anneal::init(State& state) {
    _sinceAccepted = 0;
    _currentState = state.blockCoords();
    _currentE = state.energy();
    srand(time(0)+state.comm.rank()); // Set a different random seed on all processors
  }


  void Anneal::iteration(State& state) {
    if (coolingSchedule) {
      _temp = coolingSchedule(iter);
    } else {
      _temp = tempInit / (1 + coolingRate*iter);
    }

    // Randomly perturb state
    std::vector<double> newState(state.comm.nproc);
    for (size_t i=0; i<state.comm.nblock; i++) {
      double random = 2 * ((double) rand() / RAND_MAX) - 1;
      newState[i] = _currentState[i] + random*displacement;
    }
    state.blockCoords(newState);
    state.communicate(); // Communicate to ensure halo regions are correct and not random

    // Accept or reject state
    double energy = state.energy();
    if (acceptMetropolis(energy)) {
      _currentState = state.blockCoords();
      _currentE = energy;
      _sinceAccepted = 0;
    } else {
      _sinceAccepted++;
    }

    // Set final state
    if ((checkConvergence(state)) || (iter == maxIter-1)) state.blockCoords(_currentState);
  }


  bool Anneal::acceptMetropolis(double energy) {
    double random = (double) rand() / RAND_MAX;
    if (energy < _currentE) {
      return true;
    } else if (random < exp((_currentE-energy)/_temp)) {
      return true;
    } else {
      return false;
    }
  }


  bool Anneal::checkConvergence(const State& state) {
    bool isConverged = (maxRejections > 0) && (_sinceAccepted >= maxRejections);
    return isConverged;
  }

}
