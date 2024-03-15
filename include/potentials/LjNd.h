#ifndef LJND_H
#define LJND_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {
  using std::vector;


  class LjNd : public NewPotential<LjNd> {
    public:
      int nDim;
      int nParticle;
      double sigma = 1;
      double epsilon = 1;

      LjNd(int nDim) : nDim(nDim) { _parallelDef = true; };
      LjNd(int nDim, double sigma, double epsilon);
      ~LjNd() {};

      void init(const vector<double>& coords) override;

      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      LjNd& setSigma(double sigma);
      LjNd& setEpsilon(double epsilon);
  };


  class Lj2d : public LjNd {
    public:
      Lj2d() : LjNd(2) {};
      Lj2d(double sigma, double epsilon) : LjNd(2, sigma, epsilon) {};
  };


  class Lj3d : public LjNd {
    public:
      Lj3d() : LjNd(3) {};
      Lj3d(double sigma, double epsilon) : LjNd(3, sigma, epsilon) {};
  };

}

#endif
