#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {
  using std::vector;
  class Communicator;


  /// The bar-and-hinge model for thin sheet buckling.
  ///
  /// Represents an elastic sheet as a triangulated surface: bars (edges)
  /// stretch with energy \f$ k (l - l_0)^2 / 2 \f$, and hinges (pairs of
  /// adjacent triangles) bend with energy
  /// \f$ k (\theta - \theta_0)^2 / 2 \f$. Optionally includes substrate
  /// and external force interactions.
  ///
  /// The bar and hinge lists are generated from a triangulation, or can
  /// be given explicitly. Rest lengths and angles are taken from the
  /// initial coordinates unless set.
  class BarAndHinge : public NewPotential<BarAndHinge> {
    public:
      /// @name Structure
      /// @{

      /// Set the triangulation used to generate the bar and hinge lists.
      BarAndHinge& setTriangulation(const vector2d<int>& triList);
      /// Set the bars (edges) of the model, overriding the generated list.
      BarAndHinge& setBondList(const vector2d<int>& bondList);
      /// Set the hinges of the model, overriding the generated list.
      BarAndHinge& setHingeList(const vector2d<int>& hingeList);

      /// @}

      /// @name Elastic properties
      /// @{

      /// The elastic modulus.
      double modulus = 1;
      /// The Poisson ratio of the material.
      double poissonRatio = 0.3;
      /// Set the elastic modulus.
      BarAndHinge& setModulus(double modulus);
      /// Set the sheet thickness.
      BarAndHinge& setThickness(double thickness);
      /// Set the sheet thickness per bar.
      BarAndHinge& setThickness(const vector<double>& thickness);
      /// Set the stretching and bending rigidities.
      BarAndHinge& setRigidity(double kBond, double kHinge);
      /// Set the stretching and bending rigidities per bar and hinge.
      BarAndHinge& setRigidity(const vector<double>& kBond, const vector<double>& kHinge);
      /// Set the rest length of the bars.
      BarAndHinge& setLength0(double length0);
      /// Set the rest length of each bar.
      BarAndHinge& setLength0(const vector<double>& length0);
      /// Set the rest angle of the hinges.
      BarAndHinge& setTheta0(double theta0);
      /// Set the rest angle of each hinge.
      BarAndHinge& setTheta0(const vector<double>& theta0);

      /// @}

      /// @name Substrate interaction
      /// @{

      /// Whether the substrate (wall) interaction is enabled.
      bool wallOn = false;
      /// Whether adhesion to the substrate is enabled.
      bool wallAdhesion = false;
      /// The energy parameter of the substrate interaction.
      double lj_epsilon = 1e-12;
      /// The length parameter of the substrate interaction.
      double lj_sigma = 1e-5;
      /// Enable or disable the substrate (wall) interaction.
      BarAndHinge& setWall(bool wallOn=true);
      /// Enable or disable adhesion to the substrate.
      BarAndHinge& setWallAdhesion(bool wallAdhesion=true);
      /// Set the substrate interaction parameters.
      BarAndHinge& setWallParams(double epsilon, double sigma);

      /// @}

      /// @name External force
      /// @{

      /// Set an external force on the nodes.
      BarAndHinge& setForce(const vector<double>& force);
      /// Set an external force on the nodes per bar.
      BarAndHinge& setForce(const vector2d<double>& force);

      /// @}

      /// @name Overrides
      /// @{

      int potentialType() const override { return Potential::UNSTRUCTURED; };

      void init(const vector<double>& coords) override;

      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      /// @}

    private:
      vector2d<int> _bondList;
      vector2d<int> _hingeList;
      vector2d<int> _triList;
      vector<double> thickness;
      vector<double> kBond;
      vector<double> kHinge;
      vector<double> length0;
      vector<double> theta0;
      vector<double> force;

      vector2d<int> computeBondList();
      vector2d<int> computeHingeList(const vector2d<int>& bondList);
      void computeRigidities(const vector<double>& coords, const vector2d<int>& bondList, const vector2d<int>& hingeList);
      void computeLength0(const vector<double>& coords, const vector2d<int>& bondList);
      void computeTheta0(const vector<double>& coords, const vector2d<int>& hingeList);

      void stretching(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void bending(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void forceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void substrate(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
  };

}

#endif
