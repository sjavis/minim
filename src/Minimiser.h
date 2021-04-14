#ifndef MINIMISER_H
#define MINIMISER_H

#include "System.h"


// Abstract class for minimisation proceedures
class Minimiser {
  public:
    int maxIter;
    int iter;
    System &sys;

    Minimiser(System &sys_, int maxIter_);
    virtual ~Minimiser() {};

    void minimise();
    virtual void iteration() = 0;
    virtual bool checkConvergence() = 0;
};

#endif
