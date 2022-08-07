#ifndef GRADDESCENT_H
#define GRADDESCENT_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  // Gradient desecent minimisation
  class GradDescent : public Minimiser {
    public:
      GradDescent() {};
      ~GradDescent() {};

      GradDescent& setAlpha(double alpha);
      GradDescent& setMaxIter(int maxIter);

      void init(State &state);
      void iteration(State &state);

      bool checkConvergence(const State &state) override;

    private:
      double _alpha = 1e-3;
      std::vector<double> _g;
  };

}

#endif
