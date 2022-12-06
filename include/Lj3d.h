#ifndef LJ3D_H
#define LJ3D_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {


  class Lj3d : public NewPotential<Lj3d> {
    typedef std::vector<double> Vector;

    public:
      int n_particle;
      double sigma = 1;
      double epsilon = 1;

      Lj3d() { _parallelDef = true; };
      ~Lj3d() {};

      void elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const override;

      Lj3d& setSigma(double sigma);
      Lj3d& setEpsilon(double epsilon);

      State newState(const Vector& coords, const std::vector<int>& ranks={}) override;
      State newState(const Vector& coords, double sigma, double epsilon, const std::vector<int>& ranks={});
    
    private:
  };

}

#endif
