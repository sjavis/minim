#ifndef LJND_H
#define LJND_H

#include <vector>
#include <memory>
#include "Potential.h"

namespace minim {
  using std::vector;


  /// Lennard-Jones potential in an arbitrary number of dimensions.
  ///
  /// Pairwise interaction of particles, with energy
  /// \f$ 4 \epsilon \left[ (\sigma/r)^{12} - (\sigma/r)^{6} \right] \f$,
  /// where \f$ r \f$ is the particle separation. The repulsive
  /// \f$ r^{-12} \f$ and attractive \f$ r^{-6} \f$ terms give a potential
  /// well of depth \f$ \epsilon \f$ at separation
  /// \f$ \sigma 2^{1/6} \f$.
  ///
  /// The number of particles is inferred from the length of the
  /// coordinates, and all pairs are used as energy elements.
  class LjNd : public NewPotential<LjNd> {
    public:
      int potentialType() const override { return Potential::UNSTRUCTURED; };

      /// The number of spatial dimensions.
      int nDim;
      /// The number of particles.
      int nParticle;
      /// The Lennard-Jones length parameter.
      double sigma = 1;
      /// The Lennard-Jones energy parameter.
      double epsilon = 1;

      /// Construct a potential for the given number of dimensions.
      LjNd(int nDim) : nDim(nDim) {};
      /// Construct a potential with the given dimension and parameters.
      LjNd(int nDim, double sigma, double epsilon) : nDim(nDim), sigma(sigma), epsilon(epsilon) {};

      void init(const vector<double>& coords) override;

      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      /// Set the Lennard-Jones length parameter.
      LjNd& setSigma(double sigma);
      /// Set the Lennard-Jones energy parameter.
      LjNd& setEpsilon(double epsilon);
  };


  /// Lennard-Jones potential in two dimensions.
  class Lj2d : public LjNd {
    public:
      Lj2d() : LjNd(2) {};
      Lj2d(double sigma, double epsilon) : LjNd(2, sigma, epsilon) {};
  };


  /// Lennard-Jones potential in three dimensions.
  class Lj3d : public LjNd {
    public:
      Lj3d() : LjNd(3) {};
      Lj3d(double sigma, double epsilon) : LjNd(3, sigma, epsilon) {};
  };

}

#endif
