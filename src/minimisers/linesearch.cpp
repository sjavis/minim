#include "State.h"
#include "linesearch.h"

#include "vec.h"


double backtrackingLinesearch(State &state, std::vector<double> &step, double slope) {
  const double c = 0.5; // Armijo control parameter
  const double tau = 0.5; // Shrink factor

  double t = - c * slope;
  double e0 = state.energy();
  double step_multiplier = 1;

  for (int i=0; i<10; i++) {
    double e = state.energy(vec::sum(state.blockCoords(), step));
    if (e0-e >= t) break;

    step = vec::multiply(step, tau);
    step_multiplier = step_multiplier * tau;
    t = t * tau;
  }

  state.blockCoords(vec::sum(state.blockCoords(), step));
  return step_multiplier;
}
