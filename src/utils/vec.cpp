#include "vec.h"

#include <math.h>
#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>
#include <functional>

namespace minim {

  typedef std::vector<double> Vector;


  // Sum
  Vector operator+(double a, const Vector &b) {
    Vector c(b.size());
    std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a+x; });
    return c;
  }
  Vector operator+(const Vector &a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x+b; });
    return c;
  }
  Vector operator+(const Vector &a, const Vector &b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<double>());
    return c;
  }

  Vector& operator+=(Vector &a, double b) {
    std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x+b; });
    return a;
  }
  Vector& operator+=(Vector &a, const Vector &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::plus<double>());
    return a;
  }


  // Diff
  Vector operator-(const Vector &a) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [](auto x){ return -x; });
    return c;
  }

  Vector operator-(double a, const Vector &b) {
    Vector c(b.size());
    std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a-x; });
    return c;
  }
  Vector operator-(const Vector &a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x-b; });
    return c;
  }
  Vector operator-(const Vector &a, const Vector &b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<double>());
    return c;
  }

  Vector& operator-=(Vector &a, double b) {
    std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x-b; });
    return a;
  }
  Vector& operator-=(Vector &a, const Vector &b) {
    std::transform(a.begin(), a.end(), b.begin(), a.begin(), std::minus<double>());
    return a;
  }


  // Multiply
  Vector operator*(double a, const Vector &b) {
    Vector c(b.size());
    std::transform(b.begin(), b.end(), c.begin(), [a](auto x){ return a*x; });
    return c;
  }
  Vector operator*(const Vector &a, double b) {
    Vector c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](auto x){ return x*b; });
    return c;
  }

  Vector& operator*=(Vector &a, double b) {
    std::transform(a.begin(), a.end(), a.begin(), [b](auto x){ return x*b; });
    return a;
  }


  namespace vec {

    // Dot Product
    double dotProduct(const Vector &a, const Vector &b) {
      return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
    }


    // Cross Product
    Vector crossProduct(const Vector &a, const Vector &b) {
      assert(a.size()==3);
      assert(b.size()==3);
      Vector c(3);
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = a[2]*b[0] - a[0]*b[2];
      c[2] = a[0]*b[1] - a[1]*b[0];
      return c;
    }


    // Norm
    double norm(const Vector &a) {
      return sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.0));
    }


    // Logical
    bool any(const std::vector<bool> &a) {
      return std::any_of(a.begin(), a.end(), [](bool i){ return i; });
    }

    bool all(const std::vector<bool> &a) {
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

}
