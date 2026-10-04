#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>
#include <memory>

namespace minim {
  using std::vector;
  class Potential;

  /// Handles the distribution of data over MPI processes.
  ///
  /// An abstract base class. Derived classes define how the degrees of
  /// freedom are divided into blocks and assigned to processors:
  /// CommGrid for data on a structured grid, and CommUnstructured for
  /// element-based data.
  ///
  /// The degrees of freedom are divided into blocks, one for each MPI
  /// rank used by the communicator. Each rank holds its block plus a
  /// halo of degrees of freedom, owned by neighbouring ranks, that are
  /// needed to compute the energy and gradient of the block. This
  /// processor data has size nproc. A communicator is created
  /// automatically when a State is created, and handles all
  /// communication between the ranks used by the state.
  class Communicator {
    public:
      /// Total number of degrees of freedom.
      size_t ndof;
      /// Number of degrees of freedom held by this processor (block
      /// plus halo).
      size_t nproc;
      /// Number of degrees of freedom assigned to this processor
      /// (excluding halo).
      size_t nblock;

      /// Whether the current processor is used by the communicator.
      bool usesThisProc = true;
      /// The global MPI ranks used by the communicator.
      vector<int> ranks = vector<int>();

      /// The MPI rank of this processor within the communicator
      /// (-1 if this processor is not used).
      int rank() const;
      /// The number of MPI ranks used by the communicator.
      int size() const;

      /// @name Assign data
      /// Convert global data to the local block or processor data.
      /// @{

      /// Assign the local block from global or block data.
      ///
      /// The layout of the returned data depends on the derived class:
      /// CommUnstructured returns just the block (size nblock), while
      /// CommGrid returns the processor data (size nproc), with the
      /// halo assigned if global data is given and zeroed otherwise.
      ///
      /// @param in The data to assign: either the full set of degrees
      ///   of freedom, or the local block data.
      /// @return The local block data.
      virtual vector<int> assignBlock(const vector<int>& in) const = 0;
      /// Assign the local block from global or block data.
      /// @param in The data to assign: either the full set of degrees
      ///   of freedom, or the local block data.
      /// @return The local block data.
      virtual vector<char> assignBlock(const vector<char>& in) const = 0;
      /// Assign the local block from global or block data.
      /// @param in The data to assign: either the full set of degrees
      ///   of freedom, or the local block data.
      /// @return The local block data.
      virtual vector<double> assignBlock(const vector<double>& in) const = 0;

      /// Assign the local processor data (block plus halo) from global
      /// data.
      ///
      /// @param in The full set of degrees of freedom (size ndof).
      /// @return The local processor data (size nproc).
      virtual vector<int> assignProc(const vector<int>& in) const = 0;
      /// Assign the local processor data (block plus halo) from global
      /// data.
      /// @param in The full set of degrees of freedom (size ndof).
      /// @return The local processor data (size nproc).
      virtual vector<char> assignProc(const vector<char>& in) const = 0;
      /// Assign the local processor data (block plus halo) from global
      /// data.
      /// @param in The full set of degrees of freedom (size ndof).
      /// @return The local processor data (size nproc).
      virtual vector<double> assignProc(const vector<double>& in) const = 0;

      /// @}

      /// @name Access data
      /// @{

      /// Get the block (processor) that owns a given global index.
      ///
      /// @param loc The global index of the degree of freedom.
      /// @return The index of the owning block.
      virtual int getBlock(int loc) const = 0;
      /// Get the local index of a given global index, or -1 if the
      /// degree of freedom is not held by this processor.
      ///
      /// @param loc The global index of the degree of freedom.
      /// @param block The block to look within (default: the block
      ///   owning the given index).
      /// @return The local index, or -1 if not held by this processor.
      virtual int getLocalIdx(int loc, int block=-1) const = 0;
      /// Get the value at a global index from local (processor) data,
      /// communicating with the owning processor if needed.
      ///
      /// @param vector The local processor data.
      /// @param loc The global index of the degree of freedom.
      /// @return The value at the given index.
      double get(const vector<double>& vector, int loc) const;

      /// @}

      /// @name Communication
      /// @{

      /// Communicate the halo regions of the local data.
      ///
      /// Sends the edges of this processor's block to the halo regions
      /// of the neighbouring processors, and receives this processor's
      /// halo from the neighbours' block edges.
      ///
      /// @param vector The local data whose halo is to be updated.
      void communicate(vector<double>& vector) const;
      /// Communicate the halo regions, accumulating into the existing
      /// values.
      ///
      /// Adds the values of this processor's halo region onto the
      /// corresponding block locations of the neighbouring processors.
      /// Used to accumulate contributions to the gradient from
      /// neighbouring processors.
      ///
      /// @param vector The local data to accumulate from.
      void communicateAccumulate(vector<double>& vector) const;
      /// Gather the local blocks to the given root processor
      /// (default: all).
      ///
      /// @param block The local processor data; only the block region
      ///   is used.
      /// @param root The rank to gather onto (default: all ranks).
      /// @return The full data (size ndof) on the gathering ranks,
      ///   empty on the others.
      vector<double> gather(const vector<double>& block, int root=-1) const;
      /// Scatter global data to the local processor data from the given
      /// root processor.
      ///
      /// If root is -1 (default), the data must already be present on
      /// all ranks, and each rank extracts its own portion without
      /// communication.
      ///
      /// @param data The full data (size ndof).
      /// @param root The rank holding the data (default: all ranks).
      /// @return The local processor data (size nproc).
      vector<double> scatter(const vector<double>& data, int root=-1) const;

      /// Broadcast an integer from the given root processor.
      void bcast(int& value, int root=0) const;
      /// Broadcast a double from the given root processor.
      void bcast(double& value, int root=0) const;
      /// Broadcast a vector of doubles from the given root processor.
      void bcast(vector<double>& value, int root=0) const;

      /// @}

      /// @name MPI reduction functions
      /// Perform operations on local vectors before reducing the result
      /// across all the processors used by the communicator. Every rank
      /// receives the result.
      /// @{

      /// Sum a value over all processors.
      double sum(double a) const;
      /// Sum the elements of a vector over all processors.
      ///
      /// Note that the whole local vector is summed, including any halo
      /// entries.
      double sum(const vector<double>& a) const;
      /// Compute the L2 norm of the vector over all processors.
      double norm(const vector<double>& a) const;
      /// Compute the dot product of two vectors over all processors.
      ///
      /// Only the block entries contribute; halo entries are excluded
      /// (the default implementation uses the whole vector).
      virtual double dotProduct(const vector<double>& a, const vector<double>& b) const;

      /// @}

      // Internal functions
      virtual ~Communicator() = default;
      /// Create a deep copy of the communicator.
      virtual std::unique_ptr<Communicator> clone() const = 0;
      /// Set up the distribution of the given potential over the given
      /// ranks.
      virtual void setup(Potential& pot, size_t ndof, vector<int> ranks) = 0;

    protected:
      int commSize;
      int commRank;
      vector<int> nGather;
      vector<int> iGather;

      #ifdef PARALLEL
      struct CommunicateObj {
        int rank;
        int tag;
        std::shared_ptr<MPI_Datatype> type;
      };
      MPI_Comm comm;
      vector<CommunicateObj> haloTypes;         // Objects containing halo region MPI derived datatypes for each MPI send
      vector<CommunicateObj> edgeTypes;         // Objects containing edge region MPI derived datatypes for each MPI recv
      std::shared_ptr<MPI_Datatype> blockType;  // MPI derived datatype to send the local block
      std::shared_ptr<MPI_Datatype> gatherType; // MPI derived datatype to receive the blocks for gathering
      static void mpiTypeDeleter(MPI_Datatype* type);
      bool mpiTypesCommitted = false;
      #endif

      void defaultSetup(const Potential& pot, size_t ndof, vector<int> ranks);
      void setComm(vector<int> ranks);
      virtual void makeMPITypes() = 0;
  };

}

#endif
