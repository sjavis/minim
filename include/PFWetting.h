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

      void init(const vector<double>& coords) override;
      void distributeParameters(const Communicator& comm) override;

      void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;
      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;

      State newState(const vector<double>& coords, const vector<int>& ranks={}) override;

      PFWetting& setGridSize(std::array<int,3> gridSize);
      PFWetting& setNFluid(int nFluid);
      PFWetting& setInterfaceSize(double interfaceSize);
      PFWetting& setSurfaceTension(double surfaceTension);
      PFWetting& setResolution(double resolution);
      PFWetting& setPressure(vector<double> pressure);
      PFWetting& setVolume(vector<double> volume, double volConst=1e5);
      PFWetting& setSolid(vector<bool> solid);
      PFWetting& setSolid(std::function<bool(int,int,int)> solidFn);
      PFWetting& setContactAngle(vector<double> contactAngle);
      PFWetting& setContactAngle(std::function<double(int,int,int)> contactAngleFn);
      PFWetting& setForce(vector<double> force, vector<int> iFluid={});
      PFWetting& setFixFluid(int iFluid, bool fix=true);

      std::array<int,3> gridSize;
      vector<bool> solid;
      int nFluid = 1;
      double resolution = 1;
      vector<double> interfaceSize = {1};
      vector<double> surfaceTension = {1};
      vector<double> pressure;
      vector<double> volume;
      double volConst = 1e5;
      vector<double> kappa;
      vector<double> kappaP;
      vector<double> contactAngle;
      vector<vector<double>> force;
      vector<bool> fixFluid = {false};

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
