#ifndef MINIM_MPI_H
#define MINIM_MPI_H

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace minim {

  void mpiInit();
  void mpiInit(int* argc, char*** argv);

  class Mpi {
    public:
      int size;
      int rank;

      Mpi();
      ~Mpi();
      void init(int* argc, char*** argv);
#ifdef PARALLEL
      void getSizeRank(MPI_Comm comm);
#endif

      double sum(double summand);

    private:
      bool _init = false;
  };

  extern Mpi mpi;
}

#endif
