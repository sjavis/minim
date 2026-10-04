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
  ///
  /// A state may be distributed over any subset of the available MPI
  /// ranks (all ranks by default). Each rank holds the local data for the
  /// degrees of freedom it requires. The block is the set of degrees of
  /// freedom allocated primarily to that rank, and the proc data is the
  /// block plus the halo of degrees of freedom, owned by neighbouring
  /// ranks, that are needed to compute the energy and gradient of the
  /// block. The local data therefore has size Communicator::nproc.
  ///
  /// The plain energy and gradient functions are intended for end users:
  /// they accept either the full set of coordinates or the local proc
  /// data, and return the total energy and gradient over the whole
  /// system. The remaining parallel functions are intended for writing
  /// new minimisation algorithms that make use of the distributed data.
  class State {
    public:

      /// The total number of degrees of freedom.
      size_t ndof;
      /// Convergence criterion: the minimisation is converged when the
      /// energy gradient falls below this value.
      double convergence;
      /// The potential to be minimised.
      std::unique_ptr<Potential> pot;
      /// The communicator, handling distribution over MPI processes.
      std::unique_ptr<Communicator> comm;
      /// Whether the current MPI rank is used by the state.
      bool usesThisProc = true;


      /// Construct a state from a potential and initial coordinates.
      ///
      /// @param pot The potential to be minimised.
      /// @param coords The initial coordinates.
      /// @param ranks The MPI ranks to distribute the state over
      ///   (default: all ranks).
      State(const Potential& pot, const vector<double>& coords, const vector<int>& ranks={});

      State(const State& state);
      State& operator=(const State& state);

      /// @name Energy and gradient
      /// Compute the total energy and gradient over the whole system.
      /// These are the functions intended for end users; any MPI
      /// communication required is handled internally.
      /// @{

      /// Compute the total energy of the current state.
      ///
      /// @return The total energy (0 on ranks not used by the state).
      double energy() const;
      /// Compute the total energy of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy of: either
      ///   the full set of degrees of freedom, or the local proc data.
      /// @return The total energy, identical on all ranks used by the
      ///   state (0 on unused ranks).
      double energy(const vector<double>& coords) const;
      /// Compute the gradient of the total energy of the current state.
      ///
      /// @return The gradient over all degrees of freedom (empty on
      ///   ranks not used by the state).
      vector<double> gradient() const;
      /// Compute the gradient of the total energy of the given
      /// coordinates.
      ///
      /// @param coords The coordinates to compute the gradient of: either
      ///   the full set of degrees of freedom, or the local proc data.
      /// @return The gradient over all degrees of freedom, identical on
      ///   all ranks used by the state (empty on unused ranks).
      vector<double> gradient(const vector<double>& coords) const;
      /// Compute the total energy and gradient of the current state.
      ///
      /// @param e Output: the total energy (may be null).
      /// @param g Output: the gradient over all degrees of freedom (may
      ///   be null).
      void energyGradient(double* e, vector<double>* g) const;
      /// Compute the total energy and gradient of the given coordinates.
      ///
      /// @param coords The coordinates to compute the energy and gradient
      ///   of: either the full set of degrees of freedom, or the local
      ///   proc data.
      /// @param e Output: the total energy (may be null).
      /// @param g Output: the gradient over all degrees of freedom (may
      ///   be null).
      void energyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      /// Compute the energy of a single component of the potential, e.g.
      /// for analysing the different contributions to the total energy.
      ///
      /// Unlike energy(), the result is not summed over the ranks used by
      /// the state: each rank returns the contribution of the elements it
      /// holds.
      ///
      /// @param component The index of the component.
      /// @return The energy of the component.
      double componentEnergy(int component) const;

      /// @}


      /// @name Coordinates
      /// Read and assign the coordinates of the whole system.
      /// @{

      /// Get a single coordinate by its global index.
      ///
      /// @param i The global index of the degree of freedom.
      /// @return The value of the coordinate.
      double operator[](int i);
      /// Get the full set of coordinates of the state.
      ///
      /// @return All coordinates, identical on all ranks used by the
      ///   state.
      vector<double> coords() const;
      /// Set the coordinates of the state.
      ///
      /// @param in The new coordinates, given as the full set of degrees
      ///   of freedom (distributed over the ranks used by the state).
      void coords(const vector<double>& in);

      /// @}


      /// @name Parallel functions
      /// Compute with the distributed data. These are intended for
      /// writing new minimisation algorithms; see the class description
      /// for the block / proc terminology.
      /// @{

      /// Get the local coordinates held by this processor (its block plus
      /// halo).
      ///
      /// @return The local proc data, of size Communicator::nproc.
      const vector<double>& blockCoords() const;
      /// Set the local coordinates held by this processor.
      ///
      /// The halo region is set exactly as given; call communicate()
      /// afterwards if it needs updating from the neighbouring blocks.
      ///
      /// @param in The new local coordinates, of size
      ///   Communicator::nproc.
      void blockCoords(const vector<double>& in);

      /// Compute the energy contribution of this processor's block of the
      /// current state.
      ///
      /// @return The block's contribution to the total energy: the total
      ///   is the sum of the block energies over all ranks used by the
      ///   state.
      double blockEnergy() const;
      /// Compute the energy contribution of this processor's block of the
      /// given coordinates.
      ///
      /// @param coords The local proc data to compute the energy of.
      /// @return The block's contribution to the total energy: the total
      ///   is the sum of the block energies over all ranks used by the
      ///   state.
      double blockEnergy(const vector<double>& coords) const;
      /// Compute the gradient contribution of this processor's block of
      /// the current state.
      ///
      /// @return The gradient in the local proc layout, of size
      ///   Communicator::nproc. The block entries give the total gradient
      ///   for those degrees of freedom, but the halo entries are not
      ///   guaranteed to be correct; use procGradient() if the halo
      ///   entries are also required.
      vector<double> blockGradient() const;
      /// Compute the gradient contribution of this processor's block of
      /// the given coordinates.
      ///
      /// @param coords The local proc data to compute the gradient of.
      /// @return The gradient in the local proc layout, of size
      ///   Communicator::nproc. The block entries give the total gradient
      ///   for those degrees of freedom, but the halo entries are not
      ///   guaranteed to be correct; use procGradient() if the halo
      ///   entries are also required.
      vector<double> blockGradient(const vector<double>& coords) const;
      /// Compute the energy and gradient contributions of this
      /// processor's block of the current state.
      ///
      /// @param e Output: the block's contribution to the total energy.
      /// @param g Output: the gradient in the local proc layout, of size
      ///   Communicator::nproc, as given by blockGradient().
      void blockEnergyGradient(double* e, vector<double>* g) const;
      /// Compute the energy and gradient contributions of this
      /// processor's block of the given coordinates.
      ///
      /// @param coords The local proc data to compute the energy and
      ///   gradient of.
      /// @param e Output: the block's contribution to the total energy.
      /// @param g Output: the gradient in the local proc layout, of size
      ///   Communicator::nproc, as given by blockGradient().
      void blockEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      /// Compute the energy contribution of the coordinates assigned to
      /// this processor.
      ///
      /// This is identical to blockEnergy(); it is provided for symmetry
      /// with procGradient().
      ///
      /// @return The processor's contribution to the total energy.
      double procEnergy() const;
      /// Compute the energy contribution of the given coordinates
      /// assigned to this processor.
      ///
      /// This is identical to blockEnergy(); it is provided for symmetry
      /// with procGradient().
      ///
      /// @param coords The local proc data to compute the energy of.
      /// @return The processor's contribution to the total energy.
      double procEnergy(const vector<double>& coords) const;
      /// Compute the gradient of the current state over this processor's
      /// data.
      ///
      /// This is blockGradient() followed by a halo communication, so
      /// that the returned gradient is correct over the whole proc region
      /// (block plus halo). This allows minimisation steps to be applied
      /// locally to the proc data without further communication.
      ///
      /// @return The gradient in the local proc layout, of size
      ///   Communicator::nproc.
      vector<double> procGradient() const;
      /// Compute the gradient of the given coordinates over this
      /// processor's data.
      ///
      /// This is blockGradient() followed by a halo communication, so
      /// that the returned gradient is correct over the whole proc region
      /// (block plus halo). This allows minimisation steps to be applied
      /// locally to the proc data without further communication.
      ///
      /// @param coords The local proc data to compute the gradient of.
      /// @return The gradient in the local proc layout, of size
      ///   Communicator::nproc.
      vector<double> procGradient(const vector<double>& coords) const;
      /// Compute the energy and gradient of the current state over this
      /// processor's data.
      ///
      /// This is blockEnergyGradient() followed by the halo communication
      /// of procGradient(): the energy is the processor's contribution to
      /// the total, and the gradient is correct over the whole proc
      /// region (block plus halo).
      ///
      /// @param e Output: the processor's contribution to the total
      ///   energy.
      /// @param g Output: the gradient in the local proc layout, of size
      ///   Communicator::nproc.
      void procEnergyGradient(double* e, vector<double>* g) const;
      /// Compute the energy and gradient of the given coordinates over
      /// this processor's data.
      ///
      /// This is blockEnergyGradient() followed by the halo communication
      /// of procGradient(): the energy is the processor's contribution to
      /// the total, and the gradient is correct over the whole proc
      /// region (block plus halo).
      ///
      /// @param coords The local proc data to compute the energy and
      ///   gradient of.
      /// @param e Output: the processor's contribution to the total
      ///   energy.
      /// @param g Output: the gradient in the local proc layout, of size
      ///   Communicator::nproc.
      void procEnergyGradient(const vector<double>& coords, double* e, vector<double>* g) const;

      /// Compute the total energy of the current state, returning the
      /// result on all MPI ranks.
      ///
      /// Unlike energy(), which returns the result only on the ranks used
      /// by the state, this broadcasts it to every rank. This is needed
      /// when the state is distributed over a subset of the ranks.
      ///
      /// @return The total energy.
      double allEnergy() const;
      /// Compute the gradient of the total energy of the current state,
      /// returning the result on all MPI ranks.
      ///
      /// Unlike gradient(), which returns the result only on the ranks
      /// used by the state, this broadcasts it to every rank. This is
      /// needed when the state is distributed over a subset of the ranks.
      ///
      /// @return The gradient over all degrees of freedom.
      vector<double> allGradient() const;
      /// Compute the total energy and gradient of the current state,
      /// returning the results on all MPI ranks.
      ///
      /// Unlike energyGradient(), which returns the results only on the
      /// ranks used by the state, this broadcasts them to every rank.
      /// This is needed when the state is distributed over a subset of
      /// the ranks.
      ///
      /// @param e Output: the total energy.
      /// @param g Output: the gradient over all degrees of freedom.
      void allEnergyGradient(double* e, vector<double>* g) const;
      /// Get the full set of coordinates, returning them on all MPI
      /// ranks.
      ///
      /// Unlike coords(), which returns the coordinates only on the ranks
      /// used by the state, this broadcasts them to every rank. This is
      /// needed when the state is distributed over a subset of the ranks.
      ///
      /// @return All coordinates.
      vector<double> allCoords() const;

      /// Communicate the local coordinates between processors.
      ///
      /// Updates the halo region of the local proc data using the block
      /// data of the neighbouring ranks.
      void communicate();

      /// @}


      /// @name Constraints
      /// @{

      /// Apply any constraints from the potential to the given data.
      ///
      /// @param data The data to apply the constraints to.
      void applyConstraints(vector<double>& data) const;

      /// @}


      /// @name Failure
      /// Early termination of the minimisation.
      /// @{

      /// Whether the minimisation has failed and should be terminated
      /// early.
      bool isFailed = false;
      /// Mark the state as failed, causing the minimiser to terminate
      /// early.
      void failed();

      /// @}

      vector<double> _coords;
  };

}

#endif
