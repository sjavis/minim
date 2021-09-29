#ifndef LJ3D_H
#define LJ3D_H

#include <vector>
#include "Potential.h"

namespace minim {

  class Lj3dArgs : public Args {
    public:
      int n_particle;
      double sigma;
      double epsilon;

      Lj3dArgs(int ndof, double sigma=1, double epsilon=1);
  };


  class Lj3d : public Potential {
    private:
      typedef std::vector<double> Vector;

    public:
      Lj3d() {};
      ~Lj3d() {};

      double energy(const Vector &coords, const Args &args) override;

      Vector gradient(const Vector &coords, const Args &args) override;

      using Potential::newArgs;
      Args* newArgs(int ndof) override;
  };

}

#endif
