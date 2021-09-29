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
    public:
      Lbfgs(State &state, AdjustFunc adjustModel=NULL);
      ~Lbfgs() {};

      Lbfgs& setM(int m);
      Lbfgs& setMaxIter(int maxIter);

      void iteration();

      bool checkConvergence() override;

    private:
      typedef std::vector<double> Vector;

      int _m;
      int _i_cycle;
      double _init_hessian = 1e-4;
      Vector _g0;
      Vector _g1;
      Vector _step;
      Vector _rho;
      std::vector<Vector> _s;
      std::vector<Vector> _y;

      Vector getDirection();
  };

}

#endif
