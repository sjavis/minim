#ifndef LJ3D_H
#define LJ3D_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {
  using std::vector;


  class Lj3d : public NewPotential<Lj3d> {
    public:
      int n_particle;
      double sigma = 1;
      double epsilon = 1;

      Lj3d() { _parallelDef = true; };
      Lj3d(double sigma, double epsilon);
      ~Lj3d() {};

      void init(const vector<double>& coords) override;

      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      Lj3d& setSigma(double sigma);
      Lj3d& setEpsilon(double epsilon);
  };

}

#endif
