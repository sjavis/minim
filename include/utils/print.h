#ifndef MINIM_PRINT_H
#define MINIM_PRINT_H

#include <vector>

namespace minim {

  /// Print any combination of arguments, from the root processor only.
  void print();
  template <typename T, typename ... Args>
  void print(T first, Args ... args);
  template <typename T, typename ... Args>
  void print(std::vector<T> first, Args ... args);

  /// Print any combination of arguments from all processors, with rank labels.
  void printAll();
  template <typename T, typename ... Args>
  void printAll(T first, Args ... args);
  template <typename T, typename ... Args>
  void printAll(std::vector<T> first, Args ... args);

  /// Print any combination of arguments from all processors, without rank labels.
  void printAllPlain();
  template <typename T, typename ... Args>
  void printAllPlain(T first, Args ... args);
  template <typename T, typename ... Args>
  void printAllPlain(std::vector<T> first, Args ... args);

}

#include "print.hpp"

#endif
