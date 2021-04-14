#include <math.h>
#include "System.h"
#include "Minimiser.h"
#include "GradDescent.h"


GradDescent::GradDescent(System &sys, double a, int maxIter)
 : Minimiser(sys, maxIter), alpha(a)
{}


void GradDescent::iteration() {
  double *g = new double[sys.ndof];
  sys.gradient(g);
  
  for (int i=0; i<sys.ndof; i++) {
    sys.state[i] -= alpha * g[i]; 
  }

  delete g;
}


bool GradDescent::checkConvergence() {
  double *g = new double[sys.ndof];
  sys.gradient(g);

  double sum = 0;
  for (int i=0; i<sys.ndof; i++) {
    sum += g[i]*g[i];
  }
  double rms = sqrt(sum/sys.ndof);

  delete g;
  return (rms < 1e-6);
}
