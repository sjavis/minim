#ifndef FIRE_H
#define FIRE_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Fire : public Minimiser {
    public:
      Fire(State &state, AdjustFunc adjustModel=NULL);
      Fire(State &state, double dt_max, AdjustFunc adjustModel=NULL);
      ~Fire() {};

      Fire& setMaxIter(int maxIter);
      Fire& setDtMax(double dt_max);

      void iteration();
      bool checkConvergence() override;

    private:
      int _n_min = 5;
      double _f_inc = 1.1;
      double _f_dec = 0.5;
      double _f_a = 0.99;
      double _a_start = 0.1;

      int _n_steps;
      double _a;
      double _dt;
      double _dt_max;
      double _gnorm;
      std::vector<double> _g;
      std::vector<double> _v;
  };

}

#endif
