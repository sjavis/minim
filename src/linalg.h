#ifndef LINALG_H
#define LINALG_H

#include <vector>

namespace linalg {

  std::vector<double> sum(double a, std::vector<double> b);
  std::vector<double> sum(std::vector<double> a, double b);
  std::vector<double> sum(std::vector<double> a, std::vector<double> b);

  std::vector<double> diff(double a, std::vector<double> b);
  std::vector<double> diff(std::vector<double> a, double b);
  std::vector<double> diff(std::vector<double> a, std::vector<double> b);

  std::vector<double> multiply(double a, std::vector<double> b);
  std::vector<double> multiply(std::vector<double> a, double b);

  double dotProduct(std::vector<double> a, std::vector<double> b);

}

#endif
