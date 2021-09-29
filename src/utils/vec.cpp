#include "vec.h"

#include <math.h>
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <functional>

namespace vec {

  typedef std::vector<double> Vector;


  // Sum
  Vector sum(double a, Vector b) {
    return sum(b, a);
  }

  Vector sum(Vector a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return b+x; });
    return c;
  }

  Vector sum(Vector a, Vector b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<double>());
    return c;
  }


  // Diff
  Vector diff(double a, Vector b) {
    Vector c(b.size());
    std::transform(b.begin(), b.end(), c.begin(), [a](double &x){ return a-x; });
    return c;
  }

  Vector diff(Vector a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return x-b; });
    return c;
  }

  Vector diff(Vector a, Vector b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<double>());
    return c;
  }


  // Multiply
  Vector multiply(double a, Vector b) {
    return multiply(b, a);
  }

  Vector multiply(Vector a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return b*x; });
    return c;
  }


  // Dot Product
  double dotProduct(Vector a, Vector b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  }


  // Cross Product
  Vector crossProduct(Vector a, Vector b) {
    assert(a.size()==3);
    assert(b.size()==3);
    Vector c(3);
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
    return c;
  }


  // Norm
  double norm(Vector a) {
    return sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.0));
  }


  // Logical
  bool any(std::vector<bool> a) {
    return std::any_of(a.begin(), a.end(), [](bool i){ return i; });
  }

  bool all(std::vector<bool> a) {
    return std::all_of(a.begin(), a.end(), [](bool i){ return i; });
  }


  // Insert
  void insert_unique(std::vector<int> &vec, int value) {
    for (int i=0; i<vec.size(); i++) {
      if (value <= vec[i]) {
        if (value != vec[i]) vec.insert(vec.begin()+i, value);
        return;
      }
    }
    vec.push_back(value);
  }
}
