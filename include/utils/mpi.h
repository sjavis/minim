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

  /// A helper for MPI operations over all processors, without a
  /// communicator.
  ///
  /// Programs should not create their own instances of this class;
  /// use the global ::minim::mpi object, which is created automatically
  /// when the library is included. It is initialised by calling
  /// mpi.init() (or mpiInit()) at the start of the program, after which
  /// the size and rank members give the number of processors and the
  /// rank of the current processor.
  ///
  /// For operations specific to a state, use its Communicator instead,
  /// which is restricted to the processors used by that state.
  class Mpi {
    public:
      /// The number of processors.
      int size;
      /// The rank of this processor.
      int rank;

      Mpi();
      ~Mpi();

      /// @name Initialisation
      /// @{

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

      /// @}

      /// @name MPI reduction functions
      /// Perform operations on local values before reducing the result
      /// over all processors. Every processor receives the result.
      /// @{

      /// Sum a value over all processors.
      double sum(double a) const;
      /// Sum the elements of a vector over all processors.
      double sum(const vector<double>& a) const;
      /// Compute the dot product of two vectors over all processors.
      double dotProduct(const vector<double>& a, const vector<double>& b) const;

      /// @}

      /// @name Communication
      /// @{

      /// Broadcast an integer from the given root processor.
      void bcast(int& value, int root=0) const;
      /// Broadcast a double from the given root processor.
      void bcast(double& value, int root=0) const;
      /// Broadcast a vector of doubles from the given root processor.
      void bcast(vector<double>& data, int root=0, int nData=0) const;

      /// Wait for all processors to reach this point.
      void barrier() const;

      /// @}

    private:
      bool _init = false;
  };

  /// The global MPI helper instance.
  ///
  /// This is the only instance of Mpi that should be used; it is shared
  /// by the whole library.
  extern Mpi mpi;
}

#endif
