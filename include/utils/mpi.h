#ifndef MINIM_MPI_H
#define MINIM_MPI_H

namespace minim {

  void mpiInit();
  void mpiInit(int *argc, char ***argv);

  class Mpi {
    public:
      int size;
      int rank;

      Mpi();
      ~Mpi();
      void init(int *argc, char ***argv);

      double sum(double summand);

    private:
      bool _init = false;
  };

  extern Mpi mpi;
}

#endif
