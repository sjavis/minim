#include "System.h"
#include "Minimiser.h"
#include "Lbfgs.h"


Lbfgs::Lbfgs(System &sys, int maxIter)
 : Minimiser(sys, maxIter)
{}


void Lbfgs::iteration() {
}


bool Lbfgs::checkConvergence() {
  return true;
}
