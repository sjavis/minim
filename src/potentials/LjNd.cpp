#include "LjNd.h"

#include <math.h>
#include <stdexcept>
#include "State.h"
#include "utils/vec.h"

namespace minim {
  using std::vector;


  LjNd::LjNd(int nDim, double sigma, double epsilon) : LjNd(nDim) {
    this->nDim = nDim;
    this->sigma = sigma;
    this->epsilon = epsilon;
  }


  void LjNd::init(const vector<double>& coords) {
    int nDof = coords.size();
    if (nDof % nDim != 0) throw std::invalid_argument("LjNd: Length of coords must be a multiple of the number of dimensions.");
    int nParticle = nDof / nDim;
    // Generate energy elements
    elements = {};
    for (int i=0; i<nParticle; i++) {
      for (int j=i+1; j<nParticle; j++) {
        vector<int> iDof(2*nDim);
        std::iota(iDof.begin(), iDof.begin()+nDim, nDim*i);
        std::iota(iDof.begin()+nDim, iDof.end(), nDim*j);
        elements.push_back({0, iDof});
      }
    }
  }


  void LjNd::elementEnergyGradient(const vector<double>& coords, const Element& el, double* e, vector<double>* g) const {
    vector<double> dx(nDim);
    for (int iDim=0; iDim<nDim; iDim++) {
      dx[iDim] = coords[el.idof[iDim]] - coords[el.idof[nDim+iDim]];
    }
    double r2 = vec::dotProduct(dx, dx);
    // double dx = coords[el.idof[0]] - coords[el.idof[3]];
    // double dy = coords[el.idof[1]] - coords[el.idof[4]];
    // double dz = coords[el.idof[2]] - coords[el.idof[5]];
    // double r2 = dx*dx + dy*dy + dz*dz;
    double lj6 = pow(sigma, 6) / pow(r2, 3);
    double lj12 = lj6 * lj6;

    if (e != nullptr) *e += 4*epsilon * (lj12 - lj6);

    if (g != nullptr) {
      double factor = -24 * epsilon * (2*lj12 - lj6) / r2;
      for (int iDim=0; iDim<nDim; iDim++) {
        int iDof1 = el.idof[iDim];
        int iDof2 = el.idof[nDim+iDim];
        (*g)[iDof1] += dx[iDim] * factor;
        (*g)[iDof2] -= dx[iDim] * factor;
      }
    }
  }


  LjNd& LjNd::setSigma(double sigma) {
    this->sigma = sigma;
    return *this;
  }


  LjNd& LjNd::setEpsilon(double epsilon) {
    this->epsilon = epsilon;
    return *this;
  }

}
