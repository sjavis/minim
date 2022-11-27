#ifndef VEC_H
#define VEC_H

#include <vector>

#if __has_include(<span>)
#include <span>
#endif
#ifdef __cpp_lib_span
using std::span;
#else
#include <gsl/span>
using gsl::span;
#endif


std::vector<double> operator+(double a, span<const double> b);
std::vector<double> operator+(span<const double> a, double b);
std::vector<double> operator+(span<const double> a, span<const double> b);
std::vector<double>& operator+=(std::vector<double>& a, double b);
std::vector<double>& operator+=(std::vector<double>& a, span<const double> b);

std::vector<double> operator-(span<const double> a);
std::vector<double> operator-(double a, span<const double> b);
std::vector<double> operator-(span<const double> a, double b);
std::vector<double> operator-(span<const double> a, span<const double> b);
std::vector<double>& operator-=(std::vector<double>& a, double b);
std::vector<double>& operator-=(std::vector<double>& a, span<const double> b);

std::vector<double> operator*(double a, span<const double> b);
std::vector<double> operator*(span<const double> a, double b);
std::vector<double> operator*(span<const double> a, span<const double> b);
std::vector<double>& operator*=(std::vector<double>& a, double b);
std::vector<double>& operator*=(std::vector<double>& a, span<const double> b);

std::vector<double> operator/(double a, span<const double> b);
std::vector<double> operator/(span<const double> a, double b);
std::vector<double> operator/(span<const double> a, span<const double> b);
std::vector<double>& operator/=(std::vector<double>& a, double b);
std::vector<double>& operator/=(std::vector<double>& a, span<const double> b);

namespace vec {
  double dotProduct(span<const double> a, span<const double> b);
  std::vector<double> crossProduct(span<const double> a, span<const double> b);
  double norm(span<const double> a);

  std::vector<double> abs(span<const double> a);
  std::vector<double> sqrt(span<const double> a);
  std::vector<double> pow(span<const double> a, double n);

  bool any(const std::vector<bool>& a);
  bool all(const std::vector<bool>& a);

  void insert_unique(std::vector<int>& vec, int value);

  void random(span<double> vec, double max);
  std::vector<double> random(int n, double max);
}

#endif
