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

  //! \class Lbfgs
  //! LBFGS minimisation algorithm
  class Lbfgs : public Minimiser {
    typedef std::vector<double> Vector;

    public:
      Lbfgs() {};
      ~Lbfgs() {};

      Lbfgs& setM(int m);
      Lbfgs& setMaxIter(int maxIter);

      void init(State& state);
      void iteration(State& state);

      bool checkConvergence(const State& state) override;

    private:
      int _m = 5;
      int _i_cycle;
      double _init_hessian = 1e-4;
      Vector _g;
      Vector _gNew;
      Vector _step;
      Vector _rho;
      std::vector<Vector> _s;
      std::vector<Vector> _y;

      Vector getDirection();
  };

}

#endif
