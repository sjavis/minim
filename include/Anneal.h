#ifndef ANNEAL_H
#define ANNEAL_H

#include <vector>
#include "Minimiser.h"

namespace minim {

  class Anneal : public Minimiser {
    public:
      Anneal(double temp_init, double displacement) : _temp_init(temp_init), _displacement(displacement) {};
      ~Anneal() {};

      Anneal& setMaxIter(int maxIter);
      Anneal& setTempInit(double temp_init);
      Anneal& setDisplacement(double displacement);

      void init(State &state);
      void iteration(State &state);
      bool checkConvergence(const State &state) override;

    private:
      int _since_accepted;
      int _max_rejections;
      double _temp;
      double _temp_init;
      double _displacement;
      double _current_e;
      std::vector<double> _current_state;

      bool acceptMetropolis(double energy);
  };

}

#endif
