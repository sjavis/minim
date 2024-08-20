#ifndef PFWETTING_H
#define PFWETTING_H

#include <array>
#include <vector>
#include <functional>
#include "Potential.h"

namespace minim {
  using std::vector;


  class PhaseField : public NewPotential<PhaseField> {
    public:
      PhaseField() { _parallelDef = true; };

      // System size
      int nFluid = 1;
      double resolution = 1;
      std::array<int,3> gridSize;
      PhaseField& setNFluid(int nFluid);
      PhaseField& setGridSize(std::array<int,3> gridSize);
      PhaseField& setResolution(double resolution);

      // Fluid interfaces
      // Order: 1-2, 1-3, ..., 1-N, 2-3 (continue to (N-1)-N for N-comp)
      vector<double> interfaceSize;
      vector<double> surfaceTension = {1};
      PhaseField& setInterfaceSize(double interfaceSize);
      PhaseField& setInterfaceSize(vector<double> interfaceSize);
      PhaseField& setSurfaceTension(double surfaceTension);
      PhaseField& setSurfaceTension(vector<double> surfaceTension);

      // Density constraint
      int densityConstraint = 0; // Default: Hard constraint
      double densityConst = 1;
      PhaseField& setDensityConstraint(std::string method);

      // Total volume / pressure constraints
      bool volumeFixed = false;
      double volConst = 0.01;
      vector<double> pressure;
      vector<double> volume;
      PhaseField& setPressure(vector<double> pressure);
      PhaseField& setVolume(vector<double> volume, double volConst=0.01);
      PhaseField& setVolumeFixed(bool volumeFixed, double volConst=0.01);

      // Solid nodes
      vector<bool> solid;
      vector<double> contactAngle;
      PhaseField& setSolid(vector<bool> solid);
      PhaseField& setSolid(std::function<bool(int,int,int)> solidFn);
      PhaseField& setContactAngle(vector<double> contactAngle);
      PhaseField& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      // External force
      vector<vector<double>> force;
      PhaseField& setForce(vector<double> force, vector<int> iFluid={});

      // Diffuse solid method
      vector<bool> fixFluid;
      vector<double> confinementStrength;
      PhaseField& setFixFluid(int iFluid, bool fix=true);
      PhaseField& setConfinement(vector<double> strength);

      vector<double> diffuseSolid(vector<bool> solid, int iFluid=0, bool twoStep=false);
      static vector<double> diffuseSolid(vector<bool> solid, PhaseField potential, int iFluid=0, bool twoStep=false);
      static vector<double> diffuseSolid(vector<bool> solid, std::array<int,3> gridSize, int nFluid=2, int iFluid=0, bool twoStep=false);

      // Overrides
      void init(const vector<double>& coords) override;
      void distributeParameters(const Communicator& comm) override;

      void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;
      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;


      // Read only
      int nGrid;
      double surfaceTensionMean;
      vector<double> kappa;
      vector<double> kappaP;
      vector<double> nodeVol;
      vector<int> fluidType;

    private:
      std::array<int,3> getCoord(int i) const;
      int getType(int i) const;
      void setDefaults();
      void checkArraySizes();
      void assignFluidCoefficients();

      enum{ MODEL_BASIC=0, MODEL_NCOMP=1 };
      int model = MODEL_BASIC;

      void fluidEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void fluidEnergyAll(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void fluidPairEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void pressureEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void densityConstraintEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void surfaceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void forceEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
      void ffConfinementEnergy(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const;
  };

}

#endif
