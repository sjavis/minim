#ifndef ANNEAL_H
#define ANNEAL_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Anneal : public NewMinimiser<Anneal> {
    public:
      Anneal(double tempInit, double displacement) : displacement(displacement), tempInit(tempInit) {};

      // Simulated annealing parameters
      double displacement;
      double tempInit;
      double coolingRate = 1;
      std::function<double(int)> coolingSchedule = nullptr;

      Anneal& setDisplacement(double displacement);
      Anneal& setTempInit(double tempInit);
      Anneal& setCoolingRate(double coolingRate);
      Anneal& setCoolingSchedule(std::function<double(int)> coolingSchedule);

      // Convergence parameters
      int maxRejections = 0;

      Anneal& setMaxIter(int maxIter);
      Anneal& setMaxRejections(int maxRejections);

      // Other functions
      ~Anneal() {};
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
