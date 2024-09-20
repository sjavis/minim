#ifndef VEC_H
#define VEC_H

#include <vector>
using std::vector;


template<typename T, typename U> auto operator+(T a, const vector<U>& b);
template<typename T, typename U> auto operator+(const vector<T>& a, U b);
template<typename T, typename U> auto operator+(const vector<T>& a, const vector<U>& b);
template<typename T, typename U> auto& operator+=(vector<T>& a, U b);
template<typename T, typename U> auto& operator+=(vector<T>& a, const vector<U>& b);

template<typename T> auto operator-(const vector<T>& a);
template<typename T, typename U> auto operator-(T a, const vector<U>& b);
template<typename T, typename U> auto operator-(const vector<T>& a, U b);
template<typename T, typename U> auto operator-(const vector<T>& a, const vector<U>& b);
template<typename T, typename U> auto& operator-=(vector<T>& a, U b);
template<typename T, typename U> auto& operator-=(vector<T>& a, const vector<U>& b);

template<typename T, typename U> auto operator*(T a, const vector<U>& b);
template<typename T, typename U> auto operator*(const vector<T>& a, U b);
template<typename T, typename U> auto operator*(const vector<T>& a, const vector<U>& b);
template<typename T, typename U> auto& operator*=(vector<T>& a, U b);
template<typename T, typename U> auto& operator*=(vector<T>& a, const vector<U>& b);

template<typename T, typename U> auto operator/(T a, const vector<U>& b);
template<typename T, typename U> auto operator/(const vector<T>& a, U b);
template<typename T, typename U> auto operator/(const vector<T>& a, const vector<U>& b);
template<typename T, typename U> auto& operator/=(vector<T>& a, U b);
template<typename T, typename U> auto& operator/=(vector<T>& a, const vector<U>& b);

namespace vec {
  template<typename T, typename U> auto dotProduct(const vector<T>& a, const vector<U>& b);
  template<typename T, typename U> auto crossProduct(const vector<T>& a, const vector<U>& b);
  template<typename T> auto sum(const vector<T>& a);
  template<typename T> auto product(const vector<T>& a);
  template<typename T> auto norm(const vector<T>& a);
  template<typename T> auto rms(const vector<T>& a);

  template<typename T> vector<T> abs(const vector<T>& a);
  template<typename T> vector<T> sqrt(const vector<T>& a);
  template<typename T, typename U> vector<T> pow(const vector<T>& a, U n);
  template<typename T, typename U> vector<T> pow(T a, const vector<U>& n);

  template<typename T> bool any(const vector<T>& a);
  template<typename T> bool all(const vector<T>& a);

  template<typename T> bool isIn(const vector<T>& vec, T value);

  template<typename T> std::vector<char> lessThan(const std::vector<T>& v, T s);
  template<typename T> std::vector<char> lessThan(const std::vector<T>& v1, const std::vector<T>& v2);
  template<typename T> std::vector<char> greaterThan(const std::vector<T>& v, T s);
  template<typename T> std::vector<char> greaterThan(const std::vector<T>& v1, const std::vector<T>& v2);

  template<typename T> vector<T> slice(const vector<T>& in, const vector<int>& index);
  template<typename T> vector<T> sort(const vector<T>& in, vector<int>* index=nullptr);
  template<typename T> vector<T> unique(const vector<T>& in, vector<int>* index=nullptr);

  template<typename T> void insert_unique(vector<T>& vec, T value);

  void random(vector<double>& vec, double max);
  vector<double> random(int n, double max);

  vector<int> iota(int n, int start=0);
  template<typename T, typename T1, typename T2> vector<T> arange(T1 start, T2 stop, T step, bool inclusive=false);
}

#include "vec.hpp"

#endif
