#include "Lj3d.h"

#include <math.h>
#include <stdexcept>
#include "State.h"

namespace minim {

  typedef std::vector<double> Vector;


  double Lj3d::energy(const Vector& coords) const {
    double energy = 0;

    for (auto el : elements) {
      double dx = coords[el.idof[0]] - coords[el.idof[3]];
      double dy = coords[el.idof[1]] - coords[el.idof[4]];
      double dz = coords[el.idof[2]] - coords[el.idof[5]];

      double r2 = dx*dx + dy*dy + dz*dz;
      double lj6 = pow(sigma, 6) / pow(r2, 3);
      energy += 4*epsilon * (lj6*lj6 - lj6);
    }
    return energy;
  }


  Vector Lj3d::gradient(const Vector& coords) const {
    Vector g(coords.size());

    int ne1 = elements.size();
    int ne2 = elements_halo.size();
    for (int ie=0; ie<(ne1+ne2); ie++) {
      auto el = (ie<ne1) ? elements[ie] : elements_halo[ie-ne1];

      double dx = coords[el.idof[0]] - coords[el.idof[3]];
      double dy = coords[el.idof[1]] - coords[el.idof[4]];
      double dz = coords[el.idof[2]] - coords[el.idof[5]];
      double r2 = dx*dx + dy*dy + dz*dz;
      double r = sqrt(r2);

      double lj6 = pow(sigma, 6) / pow(r2, 3);
      double dedr = -24 * epsilon * (2*lj6*lj6 - lj6) / r;

      g[el.idof[0]] += dx/r * dedr;
      g[el.idof[1]] += dy/r * dedr;
      g[el.idof[2]] += dz/r * dedr;
      g[el.idof[3]] -= dx/r * dedr;
      g[el.idof[4]] -= dy/r * dedr;
      g[el.idof[5]] -= dz/r * dedr;
    }
    return g;
  }


  Lj3d& Lj3d::setSigma(double sigma) {
    this->sigma = sigma;
    return *this;
  }


  Lj3d& Lj3d::setEpsilon(double epsilon) {
    this->epsilon = epsilon;
    return *this;
  }


  State Lj3d::newState(const Vector& coords) {
    int ndof = coords.size();
    if (ndof % 3 != 0) throw std::invalid_argument("Length of coords must be a multiple of 3.");
    n_particle = ndof / 3;
    // Generate energy elements
    elements = {};
    int id = 0;
    for (int i=0; i<n_particle; i++) {
      for (int j=i+1; j<n_particle; j++) {
        Potential::Element el = {id, 0, {3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2}};
        elements.push_back(el);
        id++;
      }
    }
    return State(*this, coords);
  }

  State Lj3d::newState(const Vector& coords, double sigma, double epsilon) {
    this->sigma = sigma;
    this->epsilon = epsilon;
    return newState(coords);
  }

}
