#ifndef MINIM_PRINT_H
#define MINIM_PRINT_H

#include <vector>


namespace minim {

  void print();
  template <typename T, typename ... Args>
  void print(T first, Args ... args);
  template <typename T, typename ... Args>
  void print(std::vector<T> first, Args ... args);

  void printAll();
  template <typename T, typename ... Args>
  void printAll(T first, Args ... args);
  template <typename T, typename ... Args>
  void printAll(std::vector<T> first, Args ... args);

}

#include "print.hpp"

#endif
