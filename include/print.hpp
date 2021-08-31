#ifndef MINIM_PRINT_HPP
#define MINIM_PRINT_HPP

#include <iostream>
#include "minimMpi.h"


namespace minim {

  template <typename T, typename ... Args>
  void print(T first, Args ... args) {
    if (mpi.rank == 0) {
      std::cout << first << " ";
      print(args...);
    }
  }

  template <typename T, typename ... Args>
  void print(std::vector<T> first, Args ... args) {
    if (mpi.rank == 0) {
      for (auto elem : first) {
        std::cout << elem << " ";
      }
      print(args...);
    }
  }

}

#endif
