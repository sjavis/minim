#ifndef LJ3D_H
#define LJ3D_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {


  class Lj3d : public NewPotential<Lj3d> {
    typedef std::vector<double> Vector;

    public:
      Lj3d() {};
      ~Lj3d() {};

      double energy(const Vector& coords, const Potential::Args& args) const override;
      Vector gradient(const Vector& coords, const Potential::Args& args) const override;

      using Potential::newArgs;
      std::unique_ptr<Potential::Args> newArgs(int ndof) override;

      class Args : public Potential::Args {
        public:
          int n_particle;
          double sigma;
          double epsilon;
          Args(int ndof, double sigma=1, double epsilon=1);
          std::unique_ptr<Potential::Args> clone() const;
      };
  };

}

#endif
