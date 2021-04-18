#include <cmath>
#include <vector>
#include "System.h"
#include "Minimiser.h"
#include "GradDescent.h"
#include "linalg.h"


GradDescent::GradDescent(System &sys, double a, int maxIter)
: Minimiser(sys, maxIter), alpha(a), g(sys.ndof)
{}


void GradDescent::iteration() {
  sys.gradient(g);
  for (int i=0; i<sys.ndof; i++) {
    sys.state[i] -= alpha * g[i];
  }
}


bool GradDescent::checkConvergence() {
  sys.gradient(g);
  double sum = linalg::dotProduct(g, g);
  double rms = sqrt(sum / sys.ndof);

  return (rms < sys.convergence_rms);
}
