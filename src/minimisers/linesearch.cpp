#include "linesearch.h"

#include "State.h"
#include "utils/vec.h"

namespace minim {

  double backtrackingLinesearch(State& state, std::vector<double>& step, double de0) {
    const double c = 0.5; // Armijo control parameter
    const double tau = 0.5; // Shrink factor

    double t = - c * de0;
    double e0 = state.energy();
    double step_multiplier = 1;
    vector<double> newCoords = state.blockCoords() + step;

    for (int i=0; i<10; i++) {
      double e = state.energy(newCoords);
      if (e0-e >= t) break;

      step = step * tau;
      newCoords = state.blockCoords() + step;
      step_multiplier = step_multiplier * tau;
      t = t * tau;
    }

    state.blockCoords(newCoords);
    return step_multiplier;
  }

}
