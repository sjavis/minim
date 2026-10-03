#ifndef GRADDESCENT_H
#define GRADDESCENT_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  /// Gradient descent minimisation.
  ///
  /// Takes steps in the direction of the negative gradient, scaled by the
  /// step size alpha. Converged when the RMS gradient falls below the
  /// state's convergence criterion.
  class GradDescent : public NewMinimiser<GradDescent> {
    public:
      /// Set the step size alpha.
      GradDescent& setAlpha(double alpha);
      /// Set the maximum number of iterations.
      GradDescent& setMaxIter(int maxIter);

      void init(State& state);
      void iteration(State& state);

      bool checkConvergence(const State& state) override;

    private:
      double _alpha = 1e-1;
      std::vector<double> _g;
  };

}

#endif
