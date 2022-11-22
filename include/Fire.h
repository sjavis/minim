#ifndef FIRE_H
#define FIRE_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Fire : public NewMinimiser<Fire> {
    public:
      Fire() {};
      Fire(double dt_max) : _dt_max(dt_max) {};
      ~Fire() {};

      Fire& setMaxIter(int maxIter);
      Fire& setDtMax(double dt_max);

      void init(State& state);
      void iteration(State& state);
      bool checkConvergence(const State& state) override;

    private:
      int _n_min = 5;
      double _f_inc = 1.1;
      double _f_dec = 0.5;
      double _f_a = 0.99;
      double _a_start = 0.1;

      int _n_steps;
      double _a;
      double _dt;
      double _dt_max = 0;
      double _gnorm;
      std::vector<double> _g;
      std::vector<double> _v;
  };

}

#endif
