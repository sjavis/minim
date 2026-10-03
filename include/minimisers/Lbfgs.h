/**
 * \file Lbfgs.h
 * \author Sam Avis
 *
 * This file contains the class for the LBFGS algorithm.
 */

#ifndef LBFGS_H
#define LBFGS_H

#include <vector>
#include "Minimiser.h"

namespace minim {
  using std::vector;
  template<typename T> using vector2d = vector<vector<T>>;
  class Communicator;


  /// LBFGS minimisation algorithm.
  ///
  /// Limited-memory Broyden-Fletcher-Goldfarb-Shanno. A quasi-Newton
  /// method that uses the last m iterations of position and gradient
  /// changes to approximate the inverse Hessian, giving Newton-like
  /// convergence without computing second derivatives.
  ///
  /// Converged when the RMS gradient falls below the state's convergence
  /// criterion.
  class Lbfgs : public NewMinimiser<Lbfgs> {
    public:
      /// Set the number of stored iterations m used to approximate the
      /// inverse Hessian (default: 5).
      Lbfgs& setM(int m);
      /// Set the maximum number of iterations.
      Lbfgs& setMaxIter(int maxIter);
      /// Set the maximum allowed step size (0: no limit).
      Lbfgs& setMaxStep(double maxStep);
      /// Set the initial step size used in the first iteration, when no
      /// curvature information is yet available (default: 1).
      Lbfgs& setInitStep(double initStep);

      void init(State& state);
      void iteration(State& state);

      bool checkConvergence(const State& state) override;

    private:
      int _m = 5;
      int _i;
      double _maxStep = 0;
      double _initStep = 1;
      vector<double> _g;
      vector<double> _rho;
      vector2d<double> _s;
      vector2d<double> _y;

      vector<double> getDirection(const Communicator& comm);
  };

}

#endif
