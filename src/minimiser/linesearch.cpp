#include "State.h"
#include "linesearch.h"

#include "vec.h"


void backtrackingLinesearch(State &state, std::vector<double> &step, double slope) {
  const double c = 0.5; // Armijo control parameter
  const double tau = 0.5; // Shrink factor

  double t = - c * slope;
  double e0 = state.energy();

  for (int i=0; i<10; i++) {
    double e = state.energy(vec::sum(state.coords, step));
    if (e0-e >= t) break;

    step = vec::multiply(step, tau);
    t = t * tau;
  }

  state.coords = vec::sum(state.coords, step);
}
