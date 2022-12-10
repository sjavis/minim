#ifndef MINIM_MPI_H
#define MINIM_MPI_H

#ifdef PARALLEL
#include <mpi.h>
#endif

#include <vector>

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

      double sum(double a) const;
      double sum(const std::vector<double>& a) const;
      double dotProduct(const std::vector<double>& a, const std::vector<double>& b) const;

      void bcast(int& value, int root=0) const;
      void bcast(double& value, int root=0) const;
      void bcast(std::vector<double>& data, int root=0, int nData=0) const;

    private:
      bool _init = false;
  };

  extern Mpi mpi;
}

#endif
