/**
 * \file Lbfgs.h
 * \author Sam Avis
 *
 * This file contains the class for the LBFGS algorithm.
 */

#ifndef LBFGS_H
#define LBFGS_H


//! \class Lbfgs
//! LBFGS minimisation algorithm
class Lbfgs : public Minimiser {
  public:
    Lbfgs(State &state);
    Lbfgs(State &state, AdjustFunc adjustModel);
    ~Lbfgs() {};

    Lbfgs& setM(int m);
    Lbfgs& setMaxIter(int maxIter);

    void iteration();

    override bool checkConvergence();

  private:
    int _m;
    int _i_cycle;
    double _init_hessian = 1e-4;
    std::vector<double> _g0;
    std::vector<double> _g1;
    std::vector<double> _step;
    std::vector<double> _rho;
    std::vector<std::vector<double>> _s;
    std::vector<std::vector<double>> _y;

    void getDirection();
};

#endif
