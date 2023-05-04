#include "utils/vec.h"

#include <random>
#include <numeric>

namespace vec {

  // Generate random vector
  std::random_device rd; 
  std::mt19937 gen(rd());

  void random(vector<double>& vec, double max) {
    if (max < 0) max = - max;
    std::uniform_real_distribution<double> distr(-max, max);
    std::generate(vec.begin(), vec.end(), [&](){ return distr(gen); });
  }

  vector<double> random(int n, double max) {
    vector<double> vec(n);
    random(vec, max);
    return vec;
  }

  vector<int> iota(int n, int start) {
    vector<int> vec(n);
    std::iota(vec.begin(), vec.end(), start);
    return vec;
  }

}
