#ifndef STATE_H
#define STATE_H

#include <vector>
#include <memory>
#include <cstddef>
#include "Potential.h"
#include "Communicator.h"

namespace minim {
  class Potential;
  using std::vector;

  /// The state of a minimisation.
  ///
  /// Holds the coordinates, potential, and communicator, and provides
  /// methods to compute energies and gradients.
  class State {
    public:
      /// Number of degrees of freedom.
      size_t ndof;
      /// Convergence criterion: the minimisation is converged when the
      /// energy gradient falls below this value.
      double convergence;
      /// The potential to be minimised.
      std::unique_ptr<Potential> pot;
      /// The communicator, handling distribution over MPI processes.
      std::unique_ptr<Communicator> comm;
      /// Whether the current processor is used in the minimisation.
      bool usesThisProc = true;


      /// Construct a state from a potential and initial coordinates.
      ///
      /// @param pot The potential to be minimised.
      /// @param coords The initial coordinates.
      /// @param ranks The MPI ranks to distribute the state over
      ///   (default: all ranks).
      State(const Potential& pot, const vector<double>& coords, const vector<int>& ranks={});

      /// Copy constructor.
      State(const State& state);
      /// Copy assignment operator.
      State& operator=(const State& state);

      // Energy / Gradient

      /// Compute the energy of the current state.
      double energy() const;
      /// Compute the energy of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy of.
      /// @return The total energy.
      double energy(const vector<double>& coords) const;
      /// Compute the energy gradient of the current state.
      vector<double> gradient() const;
      /// Compute the energy gradient of the given coordinates.
      ///
      /// @param coords The coordinates to compute the gradient of.
      /// @return The gradient of the total energy.
      vector<double> gradient(const vector<double>& coords) const;
      /// Compute the energy and gradient of the current state.
      ///
      /// @param e Output: the total energy.
      /// @param g Output: the gradient of the total energy.
      void energyGradient(double* e, vector<double>* g) const;
      /// Compute the energy and gradient of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy and gradient of.
      /// @param e Output: the total energy.
      /// @param g Output: the gradient of the total energy.
      void energyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      // Coordinates

      /// Get a single coordinate by index.
      ///
      /// @param i The index of the coordinate.
      /// @return The value of the coordinate.
      double operator[](int i);
      /// Get the coordinates of the state.
      vector<double> coords() const;
      /// Set the coordinates of the state.
      ///
      /// @param in The new coordinates.
      void coords(const vector<double>& in);

      // Parallel functions

      /// Get the coordinates held by this processor's block.
      const vector<double>& blockCoords() const;
      /// Set this processor's block of coordinates.
      ///
      /// @param in The new block coordinates.
      void blockCoords(const vector<double>& in);

      /// Compute the energy of this processor's block of the current state.
      double blockEnergy() const;
      /// Compute the energy of this processor's block of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy of.
      /// @return The block's contribution to the total energy.
      double blockEnergy(const vector<double>& coords) const;
      /// Compute the gradient of this processor's block of the current state.
      vector<double> blockGradient() const;
      /// Compute the gradient of this processor's block of the given coordinates.
      ///
      /// @param coords The coordinates to compute the gradient of.
      /// @return The block's contribution to the energy gradient.
      vector<double> blockGradient(const vector<double>& coords) const;
      /// Compute the energy and gradient of this processor's block of the current state.
      ///
      /// @param e Output: the block's contribution to the total energy.
      /// @param g Output: the block's contribution to the energy gradient.
      void blockEnergyGradient(double* e, vector<double>* g) const;
      /// Compute the energy and gradient of this processor's block of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy and gradient of.
      /// @param e Output: the block's contribution to the total energy.
      /// @param g Output: the block's contribution to the energy gradient.
      void blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      /// Compute the energy of the coordinates assigned to this processor.
      double procEnergy() const;
      /// Compute the energy of the given coordinates assigned to this processor.
      ///
      /// @param coords The coordinates to compute the energy of.
      /// @return The processor's contribution to the total energy.
      double procEnergy(const vector<double>& coords) const;
      /// Compute the gradient of the coordinates assigned to this processor.
      vector<double> procGradient() const;
      /// Compute the gradient of the given coordinates assigned to this processor.
      ///
      /// @param coords The coordinates to compute the gradient of.
      /// @return The processor's contribution to the energy gradient.
      vector<double> procGradient(const vector<double>& coords) const;
      /// Compute the energy and gradient of the coordinates assigned to this processor.
      ///
      /// @param e Output: the processor's contribution to the total energy.
      /// @param g Output: the processor's contribution to the energy gradient.
      void procEnergyGradient(double* e, vector<double>* g) const;
      /// Compute the energy and gradient of the given coordinates assigned to this processor.
      ///
      /// @param coords The coordinates to compute the energy and gradient of.
      /// @param e Output: the processor's contribution to the total energy.
      /// @param g Output: the processor's contribution to the energy gradient.
      void procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      /// Compute the energy over all coordinates, without communication.
      double allEnergy() const;
      /// Compute the gradient over all coordinates, without communication.
      vector<double> allGradient() const;
      /// Compute the energy and gradient over all coordinates, without communication.
      ///
      /// @param e Output: the total energy.
      /// @param g Output: the gradient of the total energy.
      void allEnergyGradient(double* e, vector<double>* g) const;
      /// Get all coordinates, without communication.
      vector<double> allCoords() const;

      /// Compute the energy of a single component of the potential.
      ///
      /// @param component The index of the component.
      /// @return The energy of the component.
      double componentEnergy(int component) const;

      /// Communicate the coordinates between processors.
      void communicate();

      // Constraints

      /// Apply any constraints from the potential to the given data.
      ///
      /// @param data The data to apply the constraints to.
      void applyConstraints(vector<double>& data) const;

      // Failure early completion checks
      /// Whether the minimisation has failed and should be terminated early.
      bool isFailed = false;
      /// Mark the state as failed, causing the minimiser to terminate early.
      void failed();

      vector<double> _coords;
  };

}

#endif
