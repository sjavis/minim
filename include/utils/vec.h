#ifndef VEC_H
#define VEC_H

#include <vector>


std::vector<double> operator+(double a, std::vector<double> b);
std::vector<double> operator+(const std::vector<double>& a, double b);
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double>& operator+=(std::vector<double>& a, double b);
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b);

std::vector<double> operator-(const std::vector<double>& a);
std::vector<double> operator-(double a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, double b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double>& operator-=(std::vector<double>& a, double b);
std::vector<double>& operator-=(std::vector<double>& a, const std::vector<double>& b);

std::vector<double> operator*(double a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& a, double b);
std::vector<double>& operator*=(std::vector<double>& a, double b);

namespace vec {
  double dotProduct(const std::vector<double>& a, const std::vector<double>& b);
  std::vector<double> crossProduct(const std::vector<double>& a, const std::vector<double>& b);

  double norm(const std::vector<double>& a);

  bool any(const std::vector<bool>& a);
  bool all(const std::vector<bool>& a);

  void insert_unique(std::vector<int>& vec, int value);
}

#endif
