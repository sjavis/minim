#ifndef MINIM_MPI_H
#define MINIM_MPI_H

#ifndef PARALLEL
#ifdef MPI_VERSION
#define PARALLEL
#endif
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>

namespace minim {
  using std::vector;

  /// Initialise MPI.
  void mpiInit();
  /// Initialise MPI with command line arguments.
  void mpiInit(int* argc, char*** argv);

  /// A helper for MPI operations without a communicator.
  class Mpi {
    public:
      /// The number of processors.
      int size;
      /// The rank of this processor.
      int rank;

      Mpi();
      ~Mpi();
      /// Initialise MPI, if not already initialised.
      void init();
      /// Initialise MPI with command line arguments.
      void init(int* argc, char*** argv);
#ifdef PARALLEL
      /// Initialise MPI with the given communicator.
      void init(MPI_Comm comm);
      /// Set the size and rank from the given communicator.
      void getSizeRank(MPI_Comm comm);
#endif

      /// Sum a value over all processors.
      double sum(double a) const;
      /// Sum a vector over all processors.
      double sum(const vector<double>& a) const;
      /// Compute the dot product of two vectors over all processors.
      double dotProduct(const vector<double>& a, const vector<double>& b) const;

      /// Broadcast an integer from the given root processor.
      void bcast(int& value, int root=0) const;
      /// Broadcast a double from the given root processor.
      void bcast(double& value, int root=0) const;
      /// Broadcast a vector of doubles from the given root processor.
      void bcast(vector<double>& data, int root=0, int nData=0) const;

      /// Wait for all processors to reach this point.
      void barrier() const;

    private:
      bool _init = false;
  };

  /// The global MPI helper instance.
  extern Mpi mpi;
}

#endif
