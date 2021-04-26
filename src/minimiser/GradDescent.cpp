#include <cmath>
#include <vector>
#include "State.h"
#include "Minimiser.h"
#include "GradDescent.h"
#include "vec.h"


GradDescent::GradDescent(State &state, double alpha, int maxIter)
: Minimiser(state, maxIter), alpha(alpha), _g(state.ndof)
{}


void GradDescent::iteration() {
  _g = state.gradient();
  for (int i=0; i<state.ndof; i++) {
    state.coords[i] -= alpha * _g[i];
  }
}
