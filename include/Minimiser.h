#ifndef MINIMISER_H
#define MINIMISER_H

#include <vector>

class State;


// Abstract class for minimisation proceedures
class Minimiser {
  public:
    int maxIter = 10000;

    typedef void (*AdjustFunc)(int, State&);
    int iter;
    State &state;

    Minimiser(State &state);
    Minimiser(State &state, AdjustFunc adjustModel);
    virtual ~Minimiser() {};

    void setMaxIter(int maxIter);

    std::vector<double> minimise();
    virtual void iteration() = 0;
    virtual bool checkConvergence() { false };
    void adjustModel();

  private:
    AdjustFunc _adjustModel;
};

#endif
