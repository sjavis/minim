#ifndef VEC_H
#define VEC_H

#include <vector>

namespace vec {
  typedef std::vector<double> Vector;

  Vector sum(double a, Vector b);
  Vector sum(Vector a, double b);
  Vector sum(Vector a, Vector b);

  Vector diff(double a, Vector b);
  Vector diff(Vector a, double b);
  Vector diff(Vector a, Vector b);

  Vector multiply(double a, Vector b);
  Vector multiply(Vector a, double b);

  double dotProduct(Vector a, Vector b);
  Vector crossProduct(Vector a, Vector b);

  double norm(Vector a);

  bool any(std::vector<bool> a);
  bool all(std::vector<bool> a);

  void insert_unique(std::vector<int> &vec, int value);
}

#endif
