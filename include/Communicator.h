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
  /// freedom are divided into blocks and assigned to processors.
  class Communicator {
    public:
      /// Total number of degrees of freedom.
      size_t ndof;
      /// Number of degrees of freedom on the processor (including halo).
      size_t nproc;
      /// Number of degrees of freedom assigned to the processor (excluding halo).
      size_t nblock;

      /// Whether the current processor is used.
      bool usesThisProc = true;
      /// The MPI ranks used by the communicator.
      vector<int> ranks = vector<int>();

      /// The MPI rank of the current processor.
      int rank() const;
      /// The number of processors in the communicator.
      int size() const;

      // Assign data
      /// Assign the local block from block or global data.
      virtual vector<int> assignBlock(const vector<int>& in) const = 0;
      /// Assign the local block from block or global data.
      virtual vector<char> assignBlock(const vector<char>& in) const = 0;
      /// Assign the local block from block or global data.
      virtual vector<double> assignBlock(const vector<double>& in) const = 0;

      /// Assign the local processor (including halo) from global data.
      virtual vector<int> assignProc(const vector<int>& in) const = 0;
      /// Assign the local processor (including halo) from global data.
      virtual vector<char> assignProc(const vector<char>& in) const = 0;
      /// Assign the local processor (including halo) from global data.
      virtual vector<double> assignProc(const vector<double>& in) const = 0;

      // Access data
      /// Get the processor that owns a given global index.
      virtual int getBlock(int loc) const = 0;
      /// Get the local index of a global index, or -1 if not owned by this processor.
      virtual int getLocalIdx(int loc, int block=-1) const = 0;
      /// Get the value at a global index from local (processor) data.
      double get(const vector<double>& vector, int loc) const;

      // Communication
      /// Communicate the halo regions of the local data.
      void communicate(vector<double>& vector) const;
      /// Communicate the halo regions, accumulating into the existing values.
      void communicateAccumulate(vector<double>& vector) const;
      /// Gather the local blocks to the given root processor (default: all).
      vector<double> gather(const vector<double>& block, int root=-1) const;
      /// Scatter global data to the local blocks from the given root processor.
      vector<double> scatter(const vector<double>& data, int root=-1) const;

      /// Broadcast an integer from the given root processor.
      void bcast(int& value, int root=0) const;
      /// Broadcast a double from the given root processor.
      void bcast(double& value, int root=0) const;
      /// Broadcast a vector of doubles from the given root processor.
      void bcast(vector<double>& value, int root=0) const;

      // MPI reduction functions
      /// Sum a value over all processors.
      double sum(double a) const;
      /// Sum a vector over all processors.
      double sum(const vector<double>& a) const;
      /// Compute the L2 norm of the vector over all processors.
      double norm(const vector<double>& a) const;
      /// Compute the dot product of two vectors over all processors.
      virtual double dotProduct(const vector<double>& a, const vector<double>& b) const;

      // Internal functions
      virtual ~Communicator() = default;
      /// Create a deep copy of the communicator.
      virtual std::unique_ptr<Communicator> clone() const = 0;
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
