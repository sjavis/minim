#ifndef SVPARAM_H
#define SVPARAM_H

#include <vector>

namespace minim {

  class SVParam {
    private:
      bool isScalar;
      double scalar;
      std::vector<double> vector;

    public:
      SVParam(double scalar) : isScalar(true), scalar(scalar) {}

      SVParam(std::vector<double> vector) : isScalar(true), vector(vector) {}

      double& operator[](std::size_t index) {
        if (isScalar) {
          return scalar;
        } else {
          return vector[index];
        }
      }

      const double& operator[](std::size_t index) const {
        if (isScalar) {
          return scalar;
        } else {
          return vector[index];
        }
      }

      std::size_t size() const {
        return isScalar ? 1 : vector.size();
      }

  }
}

#endif
