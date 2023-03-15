#ifndef ANNEAL_H
#define ANNEAL_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Anneal : public NewMinimiser<Anneal> {
    public:
      Anneal(double tempInit, double displacement) : tempInit(tempInit), displacement(displacement) {};
      ~Anneal() {};

      int maxRejections = 0;
      double tempInit;
      double displacement;
      double coolingRate = 1;

      Anneal& setMaxIter(int maxIter);
      Anneal& setTempInit(double tempInit);
      Anneal& setDisplacement(double displacement);

      void init(State& state);
      void iteration(State& state);
      bool checkConvergence(const State& state) override;

    private:
      int _sinceAccepted;
      double _temp;
      double _currentE;
      std::vector<double> _currentState;

      bool acceptMetropolis(double energy);
  };

}

#endif
