#ifndef PFWETTING_H
#define PFWETTING_H

#include <vector>
#include <functional>
#include "Potential.h"

namespace minim {
  using std::vector;


  class PFWetting : public NewPotential<PFWetting> {
    public:
      PFWetting() { _parallelDef = true; };
      ~PFWetting() {};

      // System size
      int nFluid = 1;
      std::array<int,3> gridSize;
      PFWetting& setNFluid(int nFluid);
      PFWetting& setGridSize(std::array<int,3> gridSize);

      // Fluid interfaces
      double resolution = 1;
      vector<double> interfaceSize;
      vector<double> surfaceTension = {1};
      PFWetting& setResolution(double resolution);
      PFWetting& setInterfaceSize(double interfaceSize);
      PFWetting& setInterfaceSize(vector<double> interfaceSize);
      PFWetting& setSurfaceTension(double surfaceTension);
      PFWetting& setSurfaceTension(vector<double> surfaceTension);

      // Density constraint
      int densityConstraint = 0;
      PFWetting& setDensityConstraint(std::string method);

      // Total volume / pressure constraints
      double volConst = 1e5;
      vector<double> pressure;
      vector<double> volume;
      PFWetting& setPressure(vector<double> pressure);
      PFWetting& setVolume(vector<double> volume, double volConst=1e5);

      // Solid nodes
      vector<bool> solid;
      vector<double> kappa;
      vector<double> kappaP;
      vector<double> contactAngle;
      PFWetting& setSolid(vector<bool> solid);
      PFWetting& setSolid(std::function<bool(int,int,int)> solidFn);
      PFWetting& setContactAngle(vector<double> contactAngle);
      PFWetting& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      // External force
      vector<vector<double>> force;
      PFWetting& setForce(vector<double> force, vector<int> iFluid={});

      // Frozen fluid method
      vector<bool> fixFluid;
      vector<double> confinementStrength;
      PFWetting& setFixFluid(int iFluid, bool fix=true);
      PFWetting& setConfinement(vector<double> strength);

      void init(const vector<double>& coords) override;
      void distributeParameters(const Communicator& comm) override;

      void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;
      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;


    private:
      vector<double> nodeVol;
      vector<int> fluidType;
      int nGrid() const;
      std::array<int,3> getCoord(int i) const;
      int getType(int i) const;
      void assignKappa();
  };

}

#endif
