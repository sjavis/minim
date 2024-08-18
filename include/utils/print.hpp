#ifndef MINIM_PRINT_HPP
#define MINIM_PRINT_HPP

#include <vector>
#include <iostream>
#include "mpi.h"

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

  template <typename T, typename ... Args>
  void printAll(T first, Args ... args) {
    std::cout << "[Rank " << mpi.rank << "/" << mpi.size << "] " << first << " ";
    printAllPlain(args...);
  }

  template <typename T, typename ... Args>
  void printAll(std::vector<T> first, Args ... args) {
    std::cout << "[Rank " << mpi.rank << "/" << mpi.size << "] ";
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    printAllPlain(args...);
  }

  template <typename T, typename ... Args>
  void printAllPlain(T first, Args ... args) {
    std::cout << first << " ";
    printAllPlain(args...);
  }

  template <typename T, typename ... Args>
  void printAllPlain(std::vector<T> first, Args ... args) {
    for (auto elem : first) {
      std::cout << elem << " ";
    }
    printAllPlain(args...);
  }

}

#endif
