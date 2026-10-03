#ifndef ANNEAL_H
#define ANNEAL_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  /// Simulated annealing minimisation.
  ///
  /// Randomly perturbs the state each iteration, accepting or rejecting
  /// moves by the Metropolis criterion at the current temperature. The
  /// temperature decreases over the minimisation, by a cooling rate or a
  /// custom schedule, so that increasingly only downhill moves are taken.
  ///
  /// Useful for finding global minma unlike other minimiser which typically
  /// look for local minima.
  class Anneal : public NewMinimiser<Anneal> {
    public:
      /// Construct with an initial temperature and displacement.
      Anneal(double tempInit, double displacement) : displacement(displacement), tempInit(tempInit) {};

      // Simulated annealing parameters
      /// The displacement used for random moves.
      double displacement;
      /// The initial temperature.
      double tempInit;
      /// The factor by which the temperature decreases each iteration.
      double coolingRate = 1;
      /// A custom cooling schedule, as a function of the iteration number.
      std::function<double(int)> coolingSchedule = nullptr;

      /// Set the displacement used for random moves.
      Anneal& setDisplacement(double displacement);
      /// Set the initial temperature.
      Anneal& setTempInit(double tempInit);
      /// Set the factor by which the temperature decreases each iteration.
      Anneal& setCoolingRate(double coolingRate);
      /// Set a custom cooling schedule as a function of the iteration number.
      Anneal& setCoolingSchedule(std::function<double(int)> coolingSchedule);

      // Convergence parameters
      /// Maximum number of consecutive rejected moves before terminating
      /// (0: disabled).
      int maxRejections = 0;

      /// Set the maximum number of iterations.
      Anneal& setMaxIter(int maxIter);
      /// Set the maximum number of consecutive rejected moves before terminating.
      Anneal& setMaxRejections(int maxRejections);

      // Other functions
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
