#ifndef MINIMISER_H
#define MINIMISER_H

class State;


// Abstract class for minimisation proceedures
class Minimiser {
  public:
    int maxIter;
    int iter;
    State &state;

    Minimiser(State &state, int maxIter);
    virtual ~Minimiser() {};

    void minimise();
    virtual void iteration() = 0;
    virtual bool checkConvergence();
};

#endif
