#include "Minimiser.h"

#include <stdexcept>
#include "State.h"
#include "utils/vec.h"

namespace minim {

  Minimiser& Minimiser::setMaxIter(int maxIter) {
    this->maxIter = maxIter;
    return *this;
  }

  Minimiser& Minimiser::setLinesearch(std::string method) {
    if (!vec::isIn({"backtracking","none"}, method)) {
      throw std::invalid_argument("Invalid line search method.");
    }
    linesearch = method;
    return *this;
  }


  std::vector<double> Minimiser::minimise(State& state, std::function<void(int,State&)> adjustState) {
    if (!state.usesThisProc) return std::vector<double>();

    init(state);
    for (iter=0; iter<=maxIter; iter++) {
      if (adjustState) adjustState(iter, state);
      iteration(state);
      if (checkConvergence(state)) break;
    }
    return state.coords();
  }

}
