#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include <vector>
#include <map>
#include <functional>
#include "Potential.h"

namespace minim {
  using std::vector;


  class PhaseField : public NewPotential<PhaseField> {
    public:
      int potentialType() const override { return Potential::GRID; };

      // System size
      int nFluid = 1;
      double resolution = 1;
      PhaseField& setNFluid(int nFluid);
      PhaseField& setGridSize(vector<int> gridSize);
      PhaseField& setResolution(double resolution);

      // Fluid interfaces
      // Order: 1-2, 1-3, ..., 1-N, 2-3 (continue to (N-1)-N for N-comp)
      vector<double> interfaceSize;
      vector<double> surfaceTension = {1};
      PhaseField& setInterfaceSize(double interfaceSize);
      PhaseField& setInterfaceSize(vector<double> interfaceSize);
      PhaseField& setSurfaceTension(double surfaceTension);
      PhaseField& setSurfaceTension(vector<double> surfaceTension);
      PhaseField& setSurfaceTension(std::function<vector<double>(int,int,int)> surfaceTensionFn);

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
      vector<char> solid;
      vector<double> contactAngle;
      PhaseField& setSolid(vector<char> solid);
      PhaseField& setSolid(std::function<bool(int,int,int)> solidFn);
      PhaseField& setContactAngle(double contactAngle);
      PhaseField& setContactAngle(vector<double> contactAngle);
      PhaseField& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      // External force
      vector<vector<double>> force;
      PhaseField& setForce(vector<double> force, vector<int> iFluid={});

      // Diffuse solid method
      vector<char> fixFluid;
      vector<double> confinementStrength;
      PhaseField& setFixFluid(int iFluid, bool fix=true);
      PhaseField& setConfinement(vector<double> strength);

      vector<double> diffuseSolid(vector<char> solid, int iFluid=0, bool twoStep=false);
      static vector<double> diffuseSolid(vector<char> solid, PhaseField potential, int iFluid=0, bool twoStep=false);
      static vector<double> diffuseSolid(vector<char> solid, vector<int> gridSize, int nFluid=2, int iFluid=0, bool twoStep=false);

      // Overrides
      void init(const vector<double>& coords) override;
      void initLocal(const vector<double>& coords, const Communicator& comm) override;

      void energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;

      std::map<std::string,vector<double>> energyComponents(const vector<double>& coords, const Communicator& comm) const;


      // Read only
      int nParams;
      double surfaceTensionMean;
      vector<double> kappa;
      vector<double> kappaP;
      vector<double> nodeVol;
      vector<double> surfaceArea;
      vector<int> fluidType;
      vector2d<int> neighbours;

      double totalVolume2;
      vector<double> ffInit;
      vector<double> fMag;
      vector2d<double> fNorm;

      int nGrid;
      vector<int> procSizes;
      vector<int> procStart;
      vector<int> haloWidths;

    private:
      int getType(int i) const;
      void setDefaults();
      void checkArraySizes();
      void assignFluidCoefficients();

      enum{ MODEL_BASIC=0, MODEL_NCOMP=1 };
      int model = MODEL_BASIC;

      void phaseGradient(const vector<double>& coords, int iGrid, int iFluid, const vector<int>& xGrid,
                         const vector<int>& neighbours, double factor, double* e, vector<double>* g) const;
      void phasePairGradient(const vector<double>& coords, int iGrid, int iFluid1, int iFluid2, const vector<int>& xGrid,
                             const vector<int>& neighbours, double factor, double* e, vector<double>* g) const;

      void fluidEnergy(const vector<double>& coords, int iNode, const vector<int>& xGrid, double* e, vector<double>* g) const;
      void fluidPairEnergy(const vector<double>& coords, int iNode, const vector<int>& xGrid, double* e, vector<double>* g) const;
      void pressureEnergy(const vector<double>& coords, int iNode, double* e, vector<double>* g) const;
      void densityConstraintEnergy(const vector<double>& coords, int iNode, double* e, vector<double>* g) const;
      void surfaceEnergy(const vector<double>& coords, int iNode, double* e, vector<double>* g) const;
      void forceEnergy(const vector<double>& coords, int iNode, const vector<int>& xGrid, double* e, vector<double>* g) const;
      void ffConfinementEnergy(const vector<double>& coords, int iNode, double* e, vector<double>* g) const;
      vector<bool> volumeConstraintEnergy(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const;
      void specialisedConstraints(const vector<double>& coords, const Communicator& comm, vector<double>& data) const override;
  };

}

#endif
