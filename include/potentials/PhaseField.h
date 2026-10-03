#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include <vector>
#include <map>
#include <functional>
#include "Potential.h"

namespace minim {
  using std::vector;


  /// Multi-component phase field potential on a structured grid.
  class PhaseField : public NewPotential<PhaseField> {
    public:
      /// @name System size
      /// @{

      /// The number of fluid components.
      int nFluid = 1;
      /// The grid resolution.
      double resolution = 1;
      /// Set the number of fluid components.
      PhaseField& setNFluid(int nFluid);
      /// Set the size of the grid along each dimension.
      PhaseField& setGridSize(vector<int> gridSize);
      /// Set the grid resolution.
      PhaseField& setResolution(double resolution);

      /// @}

      /// @name Fluid interfaces
      /// Order: 1-2, 1-3, ..., 1-N, 2-3 (continue to (N-1)-N for N-comp).
      /// @{

      /// The diffuse interface size of each fluid pair.
      vector<double> interfaceSize;
      /// The surface tension of each fluid pair.
      vector<double> surfaceTension = {1};
      /// Set the interface size for all fluid pairs.
      PhaseField& setInterfaceSize(double interfaceSize);
      /// Set the interface size of each fluid pair.
      PhaseField& setInterfaceSize(vector<double> interfaceSize);
      /// Set the surface tension for all fluid pairs.
      PhaseField& setSurfaceTension(double surfaceTension);
      /// Set the surface tension of each fluid pair.
      PhaseField& setSurfaceTension(vector<double> surfaceTension);
      /// Set the surface tension via a function of the two fluid types.
      PhaseField& setSurfaceTension(std::function<vector<double>(int,int,int)> surfaceTensionFn);

      /// @}

      /// @name Density constraint
      /// @{

      /// The density constraint method: 0 for hard constraint (default).
      int densityConstraint = 0;
      /// The strength of the soft density constraint.
      double densityConst = 1;
      /// Set the density constraint method ("hard" or "soft").
      PhaseField& setDensityConstraint(std::string method);

      /// @}

      /// @name Volume and pressure constraints
      /// @{

      /// Whether the total fluid volume is fixed.
      bool volumeFixed = false;
      /// The strength of the volume constraint.
      double volConst = 0.01;
      /// The target pressure of each fluid.
      vector<double> pressure;
      /// The target volume of each fluid.
      vector<double> volume;
      /// Set the target pressure of each fluid.
      PhaseField& setPressure(vector<double> pressure);
      /// Set the target volume of each fluid.
      PhaseField& setVolume(vector<double> volume, double volConst=0.01);
      /// Set whether the total fluid volume is fixed.
      PhaseField& setVolumeFixed(bool volumeFixed, double volConst=0.01);

      /// @}

      /// @name Solid nodes
      /// @{

      /// Which grid nodes are solid.
      vector<char> solid;
      /// The wall contact angle of each fluid pair.
      vector<double> contactAngle;
      /// Set which grid nodes are solid.
      PhaseField& setSolid(vector<char> solid);
      /// Set which grid nodes are solid via a function of the grid indices.
      PhaseField& setSolid(std::function<bool(int,int,int)> solidFn);
      /// Set the wall contact angle for all fluid pairs.
      PhaseField& setContactAngle(double contactAngle);
      /// Set the wall contact angle of each fluid pair.
      PhaseField& setContactAngle(vector<double> contactAngle);
      /// Set the wall contact angle via a function of the grid indices and fluid pair.
      PhaseField& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      /// @}

      /// @name External force
      /// @{

      /// An external force applied to the fluids.
      vector<vector<double>> force;
      /// Set an external force on the given fluid (default: all).
      PhaseField& setForce(vector<double> force, vector<int> iFluid={});

      /// @}

      /// @name Diffuse solid method
      /// @{

      vector<char> fixFluid;
      vector<double> confinementStrength;
      /// Set whether the given fluid is fixed in place.
      PhaseField& setFixFluid(int iFluid, bool fix=true);
      /// Set the confinement strength of each fluid.
      PhaseField& setConfinement(vector<double> strength);

      /// Compute a diffuse representation of the solid, usable as initial coordinates.
      vector<double> diffuseSolid(vector<char> solid, int iFluid=0, bool twoStep=false);
      /// Compute a diffuse representation of the solid using the given potential's grid.
      static vector<double> diffuseSolid(vector<char> solid, PhaseField potential, int iFluid=0, bool twoStep=false);
      /// Compute a diffuse representation of the solid using the given grid.
      static vector<double> diffuseSolid(vector<char> solid, vector<int> gridSize, int nFluid=2, int iFluid=0, bool twoStep=false);

      /// @}

      /// @name Overrides
      /// @{

      int potentialType() const override { return Potential::GRID; };
      std::unique_ptr<Communicator> newComm() const override;

      void init(const vector<double>& coords) override;
      void initLocal(const vector<double>& coords, const Communicator& comm) override;

      void energyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;

      /// Compute the named components of the energy, e.g. for analysis.
      std::map<std::string,vector<double>> energyComponents(const vector<double>& coords, const Communicator& comm) const;


      /// @}

      /// @name Read only
      /// @{

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

      /// @}

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
