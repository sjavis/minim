#include "utils/print.h"

#include <iostream>
#include "utils/mpi.h"


namespace minim {

  void print() {
    if (mpi.rank == 0) {
      std::cout << std::endl;
    }
  }

  void printAll() {
    std::cout << "[Rank " << mpi.rank << "/" << mpi.size << "]" << std::endl;
  }

  void printAllPlain() {
    std::cout << std::endl;
  }

}
