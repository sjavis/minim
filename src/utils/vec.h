#ifndef VEC_H
#define VEC_H

#include <vector>

namespace vec {

  std::vector<double> sum(double a, std::vector<double> b);
  std::vector<double> sum(std::vector<double> a, double b);
  std::vector<double> sum(std::vector<double> a, std::vector<double> b);

  std::vector<double> diff(double a, std::vector<double> b);
  std::vector<double> diff(std::vector<double> a, double b);
  std::vector<double> diff(std::vector<double> a, std::vector<double> b);

  std::vector<double> multiply(double a, std::vector<double> b);
  std::vector<double> multiply(std::vector<double> a, double b);

  double dotProduct(std::vector<double> a, std::vector<double> b);

  double norm(std::vector<double> a);
}

#endif
