#ifndef MINIMISER_H
#define MINIMISER_H

#include <vector>

namespace minim {
  class State;

  // Abstract class for minimisation proceedures
  class Minimiser {
    public:
      int maxIter = 10000;

      typedef void (*AdjustFunc)(int, State&);
      int iter;

      Minimiser() {};
      virtual ~Minimiser() {};

      virtual Minimiser& setMaxIter(int maxIter);

      std::vector<double> minimise(State &state, AdjustFunc adjustState=nullptr);
      virtual void init(State &state) {};
      virtual void iteration(State &state) = 0;
      virtual bool checkConvergence(const State &state) { return false; };
  };

}

#endif
