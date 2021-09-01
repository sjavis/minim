#include <iostream>
#include "minimMpi.h"


namespace minim {

  void print() {
    if (mpi.rank == 0) {
      std::cout << std::endl;
    }
  }

  void printAll() {
    std::cout << std::endl;
  }

}
