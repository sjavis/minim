#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {

  class BarAndHinge : public NewPotential<BarAndHinge> {
    typedef std::vector<double> Vector;

    public:
      BarAndHinge() { _parallelDef = true; };
      ~BarAndHinge() {};

      void init(const Vector& coords) override;

      void elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const override;

      BarAndHinge& setTriangulation(const std::vector<std::vector<int>>& triList);
      BarAndHinge& setBondList(const std::vector<std::vector<int>>& bondList);
      BarAndHinge& setHingeList(const std::vector<std::vector<int>>& hingeList);

      BarAndHinge& setModulus(double modulus);
      BarAndHinge& setThickness(double thickness);
      BarAndHinge& setThickness(const Vector& thickness);
      BarAndHinge& setRigidity(double kBond, double kHinge);
      BarAndHinge& setRigidity(const Vector& kBond, const Vector& kHinge);
      BarAndHinge& setLength0(double length0);
      BarAndHinge& setLength0(const Vector& length0);
      BarAndHinge& setTheta0(double theta0);
      BarAndHinge& setTheta0(const Vector& theta0);

      BarAndHinge& setWall(bool wallOn=true);
      BarAndHinge& setWallAdhesion(bool wallAdhesion=true);
      BarAndHinge& setWallParams(double epsilon, double sigma);

      BarAndHinge& setForce(const Vector& force);
      BarAndHinge& setForce(const std::vector<Vector>& force);

      BarAndHinge& setFixed(const std::vector<bool>& fixed);

      double modulus = 1;
      double poissonRatio = 0.3;
      bool wallOn = false;
      bool wallAdhesion = false;
      double lj_epsilon = 1e-12;
      double lj_sigma = 1e-5;

    private:
      std::vector<bool> fixed;
      std::vector<std::vector<int>> bondList;
      std::vector<std::vector<int>> hingeList;
      Vector thickness;
      Vector kBond;
      Vector kHinge;
      Vector length0;
      Vector theta0;
      Vector force;

      void stretching(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void bending(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void forceEnergy(const Vector& coords, const Element& el, double* e, Vector* g) const;
      void substrate(const Vector& coords, const Element& el, double* e, Vector* g) const;
  };

}

#endif
