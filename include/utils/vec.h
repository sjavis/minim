#ifndef VEC_H
#define VEC_H

#include <vector>
using std::vector;


/// Element-wise addition of a scalar and a vector.
template<typename T, typename U> auto operator+(T a, const vector<U>& b);
/// Element-wise addition of a vector and a scalar.
template<typename T, typename U> auto operator+(const vector<T>& a, U b);
/// Element-wise addition of two vectors.
template<typename T, typename U> auto operator+(const vector<T>& a, const vector<U>& b);
/// Add a scalar to each element of the vector.
template<typename T, typename U> auto& operator+=(vector<T>& a, U b);
/// Add the elements of another vector to this one.
template<typename T, typename U> auto& operator+=(vector<T>& a, const vector<U>& b);

/// Element-wise negation of a vector.
template<typename T> auto operator-(const vector<T>& a);
/// Element-wise subtraction of a vector from a scalar.
template<typename T, typename U> auto operator-(T a, const vector<U>& b);
/// Element-wise subtraction of a scalar from a vector.
template<typename T, typename U> auto operator-(const vector<T>& a, U b);
/// Element-wise subtraction of two vectors.
template<typename T, typename U> auto operator-(const vector<T>& a, const vector<U>& b);
/// Subtract a scalar from each element of the vector.
template<typename T, typename U> auto& operator-=(vector<T>& a, U b);
/// Subtract the elements of another vector from this one.
template<typename T, typename U> auto& operator-=(vector<T>& a, const vector<U>& b);

/// Element-wise multiplication of a scalar and a vector.
template<typename T, typename U> auto operator*(T a, const vector<U>& b);
/// Element-wise multiplication of a vector and a scalar.
template<typename T, typename U> auto operator*(const vector<T>& a, U b);
/// Element-wise multiplication of two vectors.
template<typename T, typename U> auto operator*(const vector<T>& a, const vector<U>& b);
/// Multiply each element of the vector by a scalar.
template<typename T, typename U> auto& operator*=(vector<T>& a, U b);
/// Multiply the elements of this vector by those of another.
template<typename T, typename U> auto& operator*=(vector<T>& a, const vector<U>& b);

/// Element-wise division of a scalar by a vector.
template<typename T, typename U> auto operator/(T a, const vector<U>& b);
/// Element-wise division of a vector by a scalar.
template<typename T, typename U> auto operator/(const vector<T>& a, U b);
/// Element-wise division of two vectors.
template<typename T, typename U> auto operator/(const vector<T>& a, const vector<U>& b);
/// Divide each element of the vector by a scalar.
template<typename T, typename U> auto& operator/=(vector<T>& a, U b);
/// Divide the elements of this vector by those of another.
template<typename T, typename U> auto& operator/=(vector<T>& a, const vector<U>& b);

/// Vector operations and helpers, in the style of NumPy.
namespace vec {
  /// Compute the dot product of two vectors.
  template<typename T, typename U> auto dotProduct(const vector<T>& a, const vector<U>& b);
  /// Compute the cross product of two vectors.
  template<typename T, typename U> auto crossProduct(const vector<T>& a, const vector<U>& b);
  /// Compute the sum of the elements.
  template<typename T> auto sum(const vector<T>& a);
  /// Compute the product of the elements.
  template<typename T> auto product(const vector<T>& a);
  /// Compute the L2 norm of the vector.
  template<typename T> auto norm(const vector<T>& a);
  /// Compute the root-mean-square of the elements.
  template<typename T> auto rms(const vector<T>& a);

  /// Compute the element-wise absolute value.
  template<typename T> vector<T> abs(const vector<T>& a);
  /// Compute the element-wise square root.
  template<typename T> vector<T> sqrt(const vector<T>& a);
  /// Raise each element to the given power.
  template<typename T, typename U> vector<T> pow(const vector<T>& a, U n);
  /// Raise a scalar to each of the given powers.
  template<typename T, typename U> vector<T> pow(T a, const vector<U>& n);

  /// Whether any element is truthy.
  template<typename T> bool any(const vector<T>& a);
  /// Whether all elements are truthy.
  template<typename T> bool all(const vector<T>& a);

  /// Whether the vector contains the given value.
  template<typename T> bool isIn(const vector<T>& vec, T value);

  /// Element-wise comparison of a vector against a scalar.
  template<typename T> std::vector<char> lessThan(const std::vector<T>& v, T s);
  /// Element-wise comparison of two vectors.
  template<typename T> std::vector<char> lessThan(const std::vector<T>& v1, const std::vector<T>& v2);
  /// Element-wise comparison of a vector against a scalar.
  template<typename T> std::vector<char> greaterThan(const std::vector<T>& v, T s);
  /// Element-wise comparison of two vectors.
  template<typename T> std::vector<char> greaterThan(const std::vector<T>& v1, const std::vector<T>& v2);

  /// Extract the elements at the given indices.
  template<typename T> vector<T> slice(const vector<T>& in, const vector<int>& index);
  /// Sort the elements, optionally outputting the sort indices.
  template<typename T> vector<T> sort(const vector<T>& in, vector<int>* index=nullptr);
  /// Get the unique elements, optionally with the index of the first occurrence.
  template<typename T> vector<T> unique(const vector<T>& in, vector<int>* index=nullptr);

  /// Insert a value into a vector if not already present.
  template<typename T> void insert_unique(vector<T>& vec, T value);

  /// Fill the vector with random values in [0, max).
  void random(vector<double>& vec, double max);
  /// Generate a vector of random values in [0, max).
  vector<double> random(int n, double max);

  /// Generate the integers 0, 1, ..., n-1, with an optional offset.
  vector<int> iota(int n, int start=0);
  /// Generate values from start to stop with the given step, optionally inclusive.
  template<typename T, typename T1, typename T2> vector<T> arange(T1 start, T2 stop, T step, bool inclusive=false);
}

#include "vec.hpp"

#endif
