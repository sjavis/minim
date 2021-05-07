#ifndef MINIMISER_H
#define MINIMISER_H

#include <vector>

class State;
class Args;


// Abstract class for minimisation proceedures
class Minimiser {
  public:
    int maxIter = 10000;

    typedef std::vector<double> Vector;
    typedef void (*AdjustFunc)(int, Vector&, Args&);
    int iter;
    State &state;

    Minimiser(State &state);
    Minimiser(State &state, AdjustFunc adjustModel);
    virtual ~Minimiser() {};

    void setMaxIter(int maxIter);

    void minimise();
    virtual void iteration() = 0;
    virtual bool checkConvergence();
    void adjustModel();

  private:
    AdjustFunc _adjustModel;
};

#endif
