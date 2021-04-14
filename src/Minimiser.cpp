#include "Minimiser.h"
#include "System.h"


Minimiser::Minimiser(System &sys_, int maxIter_)
 : sys(sys_), maxIter(maxIter_)
{}


void Minimiser::minimise() {
  for (iter=0; iter<maxIter; iter++) {
    iteration();
    if (checkConvergence()) break;
  }
}
