#include "Anneal.h"

#include <math.h>
#include <time.h>
#include "State.h"
#include "utils/mpi.h"

namespace minim {

  Anneal& Anneal::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }


  Anneal& Anneal::setTempInit(double temp_init) {
    _temp_init = temp_init;
    return *this;
  }


  Anneal& Anneal::setDisplacement(double displacement) {
    _displacement = displacement;
    return *this;
  }


  void Anneal::init(State& state) {
    _since_accepted = 0;
    _current_state = state.blockCoords();
    _current_e = state.energy();
    srand(time(0)+minim::mpi.rank); // Set a different random seed on all processors
  }


  void Anneal::iteration(State& state) {
    _temp = _temp_init / (1 + iter);

    // Randomly perturb state
    std::vector<double> new_state(state.comm.nproc);
    for (int i=0; i<state.comm.nblock; i++) {
      double random = 2 * ((double) rand() / RAND_MAX) - 1;
      new_state[i] = _current_state[i] + random*_displacement;
    }
    state.blockCoords(new_state);
    state.communicate(); // Communicate to ensure halo regions are correct and not random

    // Accept or reject state
    double energy = state.energy();
    if (acceptMetropolis(energy)) {
      _current_state = state.blockCoords();
      _current_e = energy;
      _since_accepted = 0;
    } else {
      _since_accepted++;
    }

    // Set final state
    if ((checkConvergence(state)) || (iter == maxIter-1)) state.blockCoords(_current_state);
  }


  bool Anneal::acceptMetropolis(double energy) {
    double random = (double) rand() / RAND_MAX;
    if (energy < _current_e) {
      return true;
    } else if (random < exp((_current_e-energy)/_temp)) {
      return true;
    } else {
      return false;
    }
  }


  bool Anneal::checkConvergence(const State& state) {
    bool is_converged = (_max_rejections > 0) && (_since_accepted >= _max_rejections);
    return is_converged;
  }

}
