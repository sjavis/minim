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


  //! \class Lbfgs
  //! LBFGS minimisation algorithm
  class Lbfgs : public NewMinimiser<Lbfgs> {
    public:
      Lbfgs() {};
      ~Lbfgs() {};

      Lbfgs& setM(int m);
      Lbfgs& setMaxIter(int maxIter);
      Lbfgs& setMaxStep(double maxStep);

      void init(State& state);
      void iteration(State& state);

      bool checkConvergence(const State& state) override;

    private:
      int _m = 5;
      int _i;
      double _init_hessian = 1e-4;
      double _maxStep = 0;
      vector<double> _g;
      vector<double> _gNew;
      vector<double> _step;
      vector<double> _rho;
      vector2d<double> _s;
      vector2d<double> _y;

      vector<double> getDirection();
  };

}

#endif
