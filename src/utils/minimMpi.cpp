#ifdef PARALLEL
#include <mpi.h>
#endif

#include <iostream>
#include "minimMpi.h"


void mpiInit() {
  minim::mpi.init(0, 0);
}

void mpiInit(int *argc, char ***argv) {
  minim::mpi.init(argc, argv);
}


namespace minim {

  Mpi mpi;
  
  Mpi::Mpi() {
    size = 1;
    rank = 0;
  }
  

  void Mpi::init(int *argc, char ***argv) {
#ifdef PARALLEL
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    _init = true;
#else
    std::cerr << "Warning: MPI not initialised. -DPARALLEL flag must be passed to the compiler." << std::endl;
#endif
  }
  

  Mpi::~Mpi() {
#ifdef PARALLEL
    if (_init = true) {
      MPI_Finalize();
    }
#endif
  }

  
  double Mpi::sum(double summand) {
    if ((_init = true) && (size > 1)) {
      double sum;
#ifdef PARALLEL
      MPI_Allreduce(&summand, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      return sum;
    } else {
      return summand;
    }
  }

}
