#include <cmath>
#include <vector>
#include <numeric>
#include <algorithm>
#include <functional>
#include "vec.h"

namespace vec {


  // Sum
  std::vector<double> sum(double a, std::vector<double> b) {
    return sum(b, a);
  }

  std::vector<double> sum(std::vector<double> a, double b) {
    std::vector<double> c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return b+x; });
    return c;
  }

  std::vector<double> sum(std::vector<double> a, std::vector<double> b) {
    std::vector<double> c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::plus<double>());
    return c;
  }


  // Diff
  std::vector<double> diff(double a, std::vector<double> b) {
    std::vector<double> c(b.size());
    std::transform(b.begin(), b.end(), c.begin(), [a](double &x){ return a-x; });
    return c;
  }

  std::vector<double> diff(std::vector<double> a, double b) {
    std::vector<double> c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return x-b; });
    return c;
  }

  std::vector<double> diff(std::vector<double> a, std::vector<double> b) {
    std::vector<double> c(a.size());
    std::transform(a.begin(), a.end(), b.begin(), c.begin(), std::minus<double>());
    return c;
  }


  // Multiply
  std::vector<double> multiply(double a, std::vector<double> b) {
    return multiply(b, a);
  }

  std::vector<double> multiply(std::vector<double> a, double b) {
    std::vector<double> c(a.size());
    std::transform(a.begin(), a.end(), c.begin(), [b](double &x){ return b*x; });
    return c;
  }


  // Dot Product
  double dotProduct(std::vector<double> a, std::vector<double> b) {
    return  std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  }


  // Norm
  double norm(std::vector<double> a) {
    return  sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.0));
  }


}
