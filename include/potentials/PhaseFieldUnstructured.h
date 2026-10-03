#ifndef PHASEFIELDUNSTRUCTURED_H
#define PHASEFIELDUNSTRUCTURED_H

#include <array>
#include <vector>
#include <functional>
#include "Potential.h"

namespace minim {
  using std::vector;


  /// Multi-component phase field potential on an unstructured grid.
  class PhaseFieldUnstructured : public NewPotential<PhaseFieldUnstructured> {
    public:
      int potentialType() const override { return Potential::UNSTRUCTURED; };

      /// @name System size
      /// @{

      /// The number of fluid components.
      int nFluid = 1;
      /// The grid resolution.
      double resolution = 1;
      /// Set the number of fluid components.
      PhaseFieldUnstructured& setNFluid(int nFluid);
      /// Set the size of the grid along each dimension.
      PhaseFieldUnstructured& setGridSize(vector<int> gridSize);
      /// Set the grid resolution.
      PhaseFieldUnstructured& setResolution(double resolution);

      /// @}

      /// @name Fluid interfaces
      /// Order: 1-2, 1-3, ..., 1-N, 2-3 (continue to (N-1)-N for N-comp).
      /// @{

      /// The diffuse interface size of each fluid pair.
      vector<double> interfaceSize;
      /// The surface tension of each fluid pair.
      vector<double> surfaceTension = {1};
      /// Set the interface size for all fluid pairs.
      PhaseFieldUnstructured& setInterfaceSize(double interfaceSize);
      /// Set the interface size of each fluid pair.
      PhaseFieldUnstructured& setInterfaceSize(vector<double> interfaceSize);
      /// Set the surface tension for all fluid pairs.
      PhaseFieldUnstructured& setSurfaceTension(double surfaceTension);
      /// Set the surface tension of each fluid pair.
      PhaseFieldUnstructured& setSurfaceTension(vector<double> surfaceTension);

      /// @}

      /// @name Density constraint
      /// @{

      /// The density constraint method: 0 for hard constraint (default).
      int densityConstraint = 0;
      /// The strength of the soft density constraint.
      double densityConst = 1;
      /// Set the density constraint method ("hard" or "soft").
      PhaseFieldUnstructured& setDensityConstraint(std::string method);

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
      PhaseFieldUnstructured& setPressure(vector<double> pressure);
      /// Set the target volume of each fluid.
      PhaseFieldUnstructured& setVolume(vector<double> volume, double volConst=0.01);
      /// Set whether the total fluid volume is fixed.
      PhaseFieldUnstructured& setVolumeFixed(bool volumeFixed, double volConst=0.01);

      /// @}

      /// @name Solid nodes
      /// @{

      /// Which grid nodes are solid.
      vector<char> solid;
      /// The wall contact angle of each fluid pair.
      vector<double> contactAngle;
      /// Set which grid nodes are solid.
      PhaseFieldUnstructured& setSolid(vector<char> solid);
      /// Set which grid nodes are solid via a function of the grid indices.
      PhaseFieldUnstructured& setSolid(std::function<bool(int,int,int)> solidFn);
      /// Set the wall contact angle of each fluid pair.
      PhaseFieldUnstructured& setContactAngle(vector<double> contactAngle);
      /// Set the wall contact angle via a function of the grid indices and fluid pair.
      PhaseFieldUnstructured& setContactAngle(std::function<double(int,int,int)> contactAngleFn);

      /// @}

      /// @name External force
      /// @{

      /// An external force applied to the fluids.
      vector<vector<double>> force;
      /// Set an external force on the given fluid (default: all).
      PhaseFieldUnstructured& setForce(vector<double> force, vector<int> iFluid={});

      /// @}

      /// @name Diffuse solid method
      /// @{

      vector<char> fixFluid;
      vector<double> confinementStrength;
      /// Set whether the given fluid is fixed in place.
      PhaseFieldUnstructured& setFixFluid(int iFluid, bool fix=true);
      /// Set the confinement strength of each fluid.
      PhaseFieldUnstructured& setConfinement(vector<double> strength);

      /// Compute a diffuse representation of the solid, usable as initial coordinates.
      vector<double> diffuseSolid(vector<char> solid, int iFluid=0, bool twoStep=false);
      /// Compute a diffuse representation of the solid using the given potential's grid.
      static vector<double> diffuseSolid(vector<char> solid, PhaseFieldUnstructured potential, int iFluid=0, bool twoStep=false);
      /// Compute a diffuse representation of the solid using the given grid.
      static vector<double> diffuseSolid(vector<char> solid, vector<int> gridSize, int nFluid=2, int iFluid=0, bool twoStep=false);

      /// @}

      /// @name Overrides
      /// @{

      void init(const vector<double>& coords) override;
      void initLocal(const vector<double>& coords, const Communicator& comm) override;

      void blockEnergyGradient(const vector<double>& coords, const Communicator& comm, double* e, vector<double>* g) const override;
      void elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const override;


      /// @}

      /// @name Read only
      /// @{

      int nGrid;
      double surfaceTensionMean;
      vector<double> kappa;
      vector<double> kappaP;
      vector<double> nodeVol;
      vector<int> fluidType;

    /// @}

    private:
      vector<int> getCoord(int i) const;
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
