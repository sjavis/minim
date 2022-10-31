#include "utils/mpi.h"

#include <iostream>


namespace minim {

  void mpiInit() {
    minim::mpi.init(0, 0);
  }

  void mpiInit(int* argc, char*** argv) {
    minim::mpi.init(argc, argv);
  }


  Mpi mpi;

  Mpi::Mpi() {
    size = 1;
    rank = 0;
  }


  void Mpi::init(int* argc, char*** argv) {
#ifdef PARALLEL
    MPI_Init(argc, argv);
    getSizeRank(MPI_COMM_WORLD);
    _init = true;
#else
    std::cerr << "Warning: MPI not initialised. -DPARALLEL flag must be passed to the compiler." << std::endl;
#endif
  }


#ifdef PARALLEL
  void Mpi::getSizeRank(MPI_Comm comm) {
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
  }
#endif


  Mpi::~Mpi() {
#ifdef PARALLEL
    if (_init == true) {
      MPI_Finalize();
    }
#endif
  }


  double Mpi::sum(double summand) {
#ifdef PARALLEL
    if (size > 1) {
      double sum;
      MPI_Allreduce(&summand, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      return sum;
    }
#endif
    return summand;
  }

}
