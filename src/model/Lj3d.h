#ifndef LJ3D_H
#define LJ3D_H

#include <cmath>
#include "Potential.h"


class Lj3dArgs : public Args {
  public:
    int n_particle;
    double sigma;
    double epsilon;
    
    Lj3dArgs(int ndof, double sigma=1, double epsilon=1);
};


class Lj3d : public Potential {
  public:
    Lj3d() {};
    ~Lj3d() {};

    double energy(std::vector<double> coords, Args &args) override;

    std::vector<double> gradient(std::vector<double> coords, Args &args) override;

    using Potential::newState;
    State newState(std::vector<double> coords) override;
};

#endif
