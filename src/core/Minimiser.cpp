#include "Minimiser.h"

#include "State.h"
#include "vec.h"

namespace minim {

  Minimiser& Minimiser::setMaxIter(int maxIter) {
    this->maxIter = maxIter;
    return *this;
  }


  std::vector<double> Minimiser::minimise(State &state, AdjustFunc adjustState) {
    init(state);
    for (iter=0; iter<maxIter; iter++) {
      if (adjustState) adjustState(iter, state);
      iteration(state);
      if (checkConvergence(state)) break;
    }
    return state.getCoords();
  }

}