#ifndef PFWETTING_H
#define PFWETTING_H

#include <vector>
#include <functional>
#include "Potential.h"

namespace minim {

  class PFWetting : public NewPotential<PFWetting> {
    typedef std::vector<double> Vector;

    public:
      PFWetting() { _blockEnergyGradientDef = true; };
      ~PFWetting() {};

      void init();
      void blockEnergyGradient(const Vector& coords, const Communicator& comm, double* e, Vector* g) const override;

      State newState(const Vector& coords) override;

      PFWetting& setGridSize(std::array<int,3> gridSize);
      PFWetting& setEpsilon(double epsilon);
      PFWetting& setResolution(double resolution);
      PFWetting& setPressure(double pressure);
      PFWetting& setVolume(double volume, double volConst=1e5);
      PFWetting& setSolid(std::vector<bool> solid);
      PFWetting& setSolid(std::function<bool(int,int,int)> solidFn);
      PFWetting& setContactAngle(Vector contactAngle);
      PFWetting& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      std::array<int,3> gridSize;
      double epsilon = 1;
      double resolution = 1;
      double pressure = 0;
      double volume = 0;
      double volConst = 1e5;
      std::vector<bool> solid;
      Vector nodeVol;
      Vector contactAngle;
      Vector force1;
      Vector force2;

    private:
      void elementEnergyGradient(const Element el, const Vector& coords, double* e, Vector* g) const;
      int nGrid() const;
      std::array<int,3> getCoord(int i) const;
      int getType(int i) const;
  };

}

#endif
