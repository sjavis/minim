#ifndef VEC_H
#define VEC_H

#include <vector>

namespace minim {
  typedef std::vector<double> Vector;

  Vector operator+(double a, Vector b);
  Vector operator+(const Vector &a, double b);
  Vector operator+(const Vector &a, const Vector &b);
  Vector& operator+=(Vector &a, double b);
  Vector& operator+=(Vector &a, const Vector &b);

  Vector operator-(const Vector &a);
  Vector operator-(double a, const Vector &b);
  Vector operator-(const Vector &a, double b);
  Vector operator-(const Vector &a, const Vector &b);
  Vector& operator-=(Vector &a, double b);
  Vector& operator-=(Vector &a, const Vector &b);

  Vector operator*(double a, const Vector &b);
  Vector operator*(const Vector &a, double b);
  Vector& operator*=(Vector &a, double b);

  namespace vec {
    double dotProduct(const Vector &a, const Vector &b);
    Vector crossProduct(const Vector &a, const Vector &b);

    double norm(const Vector &a);

    bool any(const std::vector<bool> &a);
    bool all(const std::vector<bool> &a);

    void insert_unique(std::vector<int> &vec, int value);
  }
}

#endif
