#ifndef GRADDESCENT_H
#define GRADDESCENT_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  // Gradient desecent minimisation
  class GradDescent : public Minimiser {
    public:
      GradDescent(State &state, AdjustFunc adjustModel=NULL);
      ~GradDescent() {};

      GradDescent& setAlpha(double alpha);
      GradDescent& setMaxIter(int maxIter);

      void iteration();

      bool checkConvergence() override;

    private:
      double _alpha = 1e-3;
      std::vector<double> _g;
  };

}

#endif
