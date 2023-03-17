#include "GradDescent.h"

#include <math.h>
#include "State.h"
#include "linesearch.h"
#include "utils/vec.h"

namespace minim {

  GradDescent& GradDescent::setAlpha(double alpha) {
    _alpha = alpha;
    return *this;
  }


  GradDescent& GradDescent::setMaxIter(int maxIter) {
    Minimiser::setMaxIter(maxIter);
    return *this;
  }


  void GradDescent::init(State& state) {
    _g = std::vector<double>(state.ndof);
  }


  void GradDescent::iteration(State& state) {
    // Get step
    _g = state.procGradient();
    auto step = -_alpha * _g;

    // Perform linesearch
    if (linesearch == "backtracking") {
      double gs = state.comm.dotProduct(_g, step);
      backtrackingLinesearch(state, step, gs);
    } else {
      state.blockCoords(state.blockCoords() + step); // The step is correct on the halo, so no need to communicate the coords
    }
  }


  bool GradDescent::checkConvergence(const State& state) {
    double sum = state.comm.dotProduct(_g, _g);
    double rms = sqrt(sum/state.ndof);
    return (rms < state.convergence);
  }

}
