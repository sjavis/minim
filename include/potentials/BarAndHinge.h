#ifndef BARANDHINGE_H
#define BARANDHINGE_H

#include <vector>
#include "Potential.h"

namespace minim {
  using std::vector;
  class Communicator;


  class BarAndHinge : public NewPotential<BarAndHinge> {
    public:
      BarAndHinge() { _parallelDef = true; };

      void init(const vector<double>& coords) override;

      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      BarAndHinge& setTriangulation(const vector2d<int>& triList);
      BarAndHinge& setBondList(const vector2d<int>& bondList);
      BarAndHinge& setHingeList(const vector2d<int>& hingeList);

      BarAndHinge& setModulus(double modulus);
      BarAndHinge& setThickness(double thickness);
      BarAndHinge& setThickness(const vector<double>& thickness);
      BarAndHinge& setRigidity(double kBond, double kHinge);
      BarAndHinge& setRigidity(const vector<double>& kBond, const vector<double>& kHinge);
      BarAndHinge& setLength0(double length0);
      BarAndHinge& setLength0(const vector<double>& length0);
      BarAndHinge& setTheta0(double theta0);
      BarAndHinge& setTheta0(const vector<double>& theta0);

      BarAndHinge& setWall(bool wallOn=true);
      BarAndHinge& setWallAdhesion(bool wallAdhesion=true);
      BarAndHinge& setWallParams(double epsilon, double sigma);

      BarAndHinge& setForce(const vector<double>& force);
      BarAndHinge& setForce(const vector2d<double>& force);

      double modulus = 1;
      double poissonRatio = 0.3;
      bool wallOn = false;
      bool wallAdhesion = false;
      double lj_epsilon = 1e-12;
      double lj_sigma = 1e-5;

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
