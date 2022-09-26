#ifndef LJ3D_H
#define LJ3D_H

#include <vector>
#include "Potential.h"

namespace minim {


  class Lj3d : public Potential {
    typedef std::vector<double> Vector;

    public:
      Lj3d() {};
      ~Lj3d() {};

      double energy(const Vector &coords, const Potential::Args &args) const override;
      Vector gradient(const Vector &coords, const Potential::Args &args) const override;

      using Potential::newArgs;
      Potential::Args* newArgs(int ndof) override;

      class Args : public Potential::Args {
        public:
          int n_particle;
          double sigma;
          double epsilon;
          Args(int ndof, double sigma=1, double epsilon=1);
      };
  };

}

#endif
