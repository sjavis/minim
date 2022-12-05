#include "Lj3d.h"

#include <math.h>
#include <stdexcept>
#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;


  void Lj3d::blockEnergyGradient(const Vector& coords, const Communicator& comm, double* e, Vector* g) const {
    if (e != nullptr) *e = 0;
    if (g != nullptr) *g = Vector(coords.size());

    for (auto el : elements) {
      elementEnergyGradient(el, coords, e, g);
    }

    // Compute the gradient of the halo energy elements
    if (g != nullptr) {
      for (auto el : elements_halo) {
        elementEnergyGradient(el, coords, nullptr, g);
      }
    }
  }


  void Lj3d::elementEnergyGradient(const Element el, const Vector& coords, double* e, Vector* g) const {
    double dx = coords[el.idof[0]] - coords[el.idof[3]];
    double dy = coords[el.idof[1]] - coords[el.idof[4]];
    double dz = coords[el.idof[2]] - coords[el.idof[5]];
    double r2 = dx*dx + dy*dy + dz*dz;
    double lj6 = pow(sigma, 6) / pow(r2, 3);

    if (e != nullptr) *e += 4*epsilon * (lj6*lj6 - lj6);

    if (g != nullptr) {
      double r = sqrt(r2);
      double dedr = -24 * epsilon * (2*lj6*lj6 - lj6) / r;
      (*g)[el.idof[0]] += dx/r * dedr;
      (*g)[el.idof[1]] += dy/r * dedr;
      (*g)[el.idof[2]] += dz/r * dedr;
      (*g)[el.idof[3]] -= dx/r * dedr;
      (*g)[el.idof[4]] -= dy/r * dedr;
      (*g)[el.idof[5]] -= dz/r * dedr;
    }
  }


  Lj3d& Lj3d::setSigma(double sigma) {
    this->sigma = sigma;
    return *this;
  }


  Lj3d& Lj3d::setEpsilon(double epsilon) {
    this->epsilon = epsilon;
    return *this;
  }


  State Lj3d::newState(const Vector& coords, const std::vector<int>& ranks) {
    int ndof = coords.size();
    if (ndof % 3 != 0) throw std::invalid_argument("Length of coords must be a multiple of 3.");
    n_particle = ndof / 3;
    // Generate energy elements
    elements = {};
    for (int i=0; i<n_particle; i++) {
      for (int j=i+1; j<n_particle; j++) {
        elements.push_back({0, {3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2}});
      }
    }
    return State(*this, coords, ranks);
  }

  State Lj3d::newState(const Vector& coords, double sigma, double epsilon, const std::vector<int>& ranks) {
    this->sigma = sigma;
    this->epsilon = epsilon;
    return newState(coords, ranks);
  }

}
