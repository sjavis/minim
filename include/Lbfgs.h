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
    Lbfgs(State &state, int m = 5, int maxIter = 10000);
    ~Lbfgs() {};

    void iteration();

  private:
    int _m;
    int _i_cycle;
    double _init_hessian = 1e-4;
    std::vector<double> _g;
    std::vector<double> _rho;
    std::vector<std::vector<double>> _s;
    std::vector<std::vector<double>> _y;

    void getDirection(std::vector<double> &step);
    void linesearch(std::vector<double> &step);
};

#endif
