#ifndef FIRE_H
#define FIRE_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Fire : public NewMinimiser<Fire> {
    public:
      Fire() : Fire({}) {};
      Fire(double dtMax);
      ~Fire() {};

      double dtMax = 0;

      Fire& setMaxIter(int maxIter);
      Fire& setDtMax(double dtMax);

      void init(State& state);
      void iteration(State& state);
      bool checkConvergence(const State& state) override;

    private:
      int _nMin = 5;
      double _fInc = 1.1;
      double _fDec = 0.5;
      double _fA = 0.99;
      double _aStart = 0.1;

      int _nSteps;
      double _a;
      double _dt;
      double _gNorm;
      std::vector<double> _g;
      std::vector<double> _v;
  };

}

#endif
