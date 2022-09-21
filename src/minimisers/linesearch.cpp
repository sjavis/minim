#include "linesearch.h"

#include "State.h"
#include "utils/vec.h"

namespace minim {

  double backtrackingLinesearch(State &state, std::vector<double> &step, double de0) {
    const double c = 0.5; // Armijo control parameter
    const double tau = 0.5; // Shrink factor

    double t = - c * de0;
    double e0 = state.energy();
    double step_multiplier = 1;

    for (int i=0; i<10; i++) {
      double e = state.energy(state.blockCoords() + step);
      if (e0-e >= t) break;

      step = step * tau;
      step_multiplier = step_multiplier * tau;
      t = t * tau;
    }

    state.blockCoords(state.blockCoords() + step);
    return step_multiplier;
  }

}
