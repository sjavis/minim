#ifndef MINIM_PRINT_H
#define MINIM_PRINT_H


namespace minim {

  void print();

  template <typename T, typename ... Args>
  void print(T first, Args ... args);

  template <typename T, typename ... Args>
  void print(std::vector<T> first, Args ... args);
}

#include "print.hpp"

#endif
