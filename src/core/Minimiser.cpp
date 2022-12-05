#include "Minimiser.h"

#include "State.h"

namespace minim {

  Minimiser& Minimiser::setMaxIter(int maxIter) {
    this->maxIter = maxIter;
    return *this;
  }


  std::vector<double> Minimiser::minimise(State& state, AdjustFunc adjustState) {
    if (!state.usesThisProc) return state.coords();

    init(state);
    for (iter=0; iter<maxIter; iter++) {
      if (adjustState) adjustState(iter, state);
      iteration(state);
      if (checkConvergence(state)) break;
    }
    return state.coords();
  }

}
