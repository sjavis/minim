#ifndef FIRE_H
#define FIRE_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  /// FIRE (Fast Inertial Relaxation Engine) minimisation.
  ///
  /// A damped molecular dynamics method: velocities follow the forces,
  /// with an adaptive time step that grows during downhill motion and
  /// resets when the motion turns uphill. Converged when the RMS gradient
  /// falls below the state's convergence criterion.
  class Fire : public NewMinimiser<Fire> {
    public:
      Fire() = default;
      /// Construct with a maximum time step.
      Fire(double dtMax);

      /// Maximum time step of the dynamics.
      ///
      /// If left as 0, it is estimated from the initial gradient norm.
      double dtMax = 0;

      /// Set the maximum number of iterations.
      Fire& setMaxIter(int maxIter);
      /// Set the maximum time step of the dynamics.
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
