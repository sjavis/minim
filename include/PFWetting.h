#ifndef PFWETTING_H
#define PFWETTING_H

#include <vector>
#include <functional>
#include "Potential.h"

namespace minim {

  class PFWetting : public NewPotential<PFWetting> {
    typedef std::vector<double> Vector;

    public:
      PFWetting() { _parallelDef = true; };
      ~PFWetting() {};

      void init();
      void blockEnergyGradient(const Vector& coords, const Communicator& comm, double* e, Vector* g) const override;
      void elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const override;

      State newState(const Vector& coords, const std::vector<int>& ranks={}) override;

      PFWetting& setGridSize(std::array<int,3> gridSize);
      PFWetting& setNFluid(int nFluid);
      PFWetting& setInterfaceSize(double interfaceSize);
      PFWetting& setSurfaceTension(double surfaceTension);
      PFWetting& setResolution(double resolution);
      PFWetting& setPressure(Vector pressure);
      PFWetting& setVolume(Vector volume, double volConst=1e5);
      PFWetting& setSolid(std::vector<bool> solid);
      PFWetting& setSolid(std::function<bool(int,int,int)> solidFn);
      PFWetting& setContactAngle(Vector contactAngle);
      PFWetting& setContactAngle(std::function<double(int,int,int)> contactAngleFn);
      PFWetting& setForce(Vector force, std::vector<int> iFluid={});

      std::array<int,3> gridSize;
      std::vector<bool> solid;
      int nFluid = 1;
      double resolution = 1;
      Vector interfaceSize = {1};
      Vector surfaceTension = {1};
      Vector pressure;
      Vector volume;
      double volConst = 1e5;
      Vector nodeVol;
      Vector kappa;
      Vector kappaP;
      Vector contactAngle;
      std::vector<Vector> force;

    private:
      int nGrid() const;
      std::array<int,3> getCoord(int i) const;
      int getType(int i) const;
      void assignKappa();
  };

}

#endif
