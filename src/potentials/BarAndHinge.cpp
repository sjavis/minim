#include "BarAndHinge.h"

#include <math.h>
#include "utils/vec.h"
#include "utils/print.h"

namespace minim {

  typedef std::vector<double> Vector;

  constexpr double pi() { return 4*atan(1); }


  void BarAndHinge::elementEnergyGradient(const Vector& coords, const Element& el, double* e, Vector* g) const {
    switch (el.type) {
      case 0:
        stretching(coords, el, e, g);
        break;
      case 1:
        bending(coords, el, e, g);
        break;
      case 2:
        forceEnergy(coords, el, e, g);
        break;
    }
  }


  void BarAndHinge::stretching(const Vector& coords, const Element& el, double* e, Vector* g) const {
    Vector x1(coords.cbegin()+el.idof[0], coords.cbegin()+el.idof[2]+1);
    Vector x2(coords.cbegin()+el.idof[3], coords.cbegin()+el.idof[5]+1);

    // Compute distance
    auto dx = x1 - x2;
    auto l = vec::norm(dx);
    auto dl = l - el.parameters[1];

    if (e != nullptr) {
      *e += el.parameters[0] * pow(dl, 2);
    }
    if (g != nullptr) {
      auto g_factor = 2 * el.parameters[0] * dl / l;
      if (fixed.empty() || !fixed[el.idof[0]]) (*g)[el.idof[0]] += g_factor * dx[0];
      if (fixed.empty() || !fixed[el.idof[1]]) (*g)[el.idof[1]] += g_factor * dx[1];
      if (fixed.empty() || !fixed[el.idof[2]]) (*g)[el.idof[2]] += g_factor * dx[2];
      if (fixed.empty() || !fixed[el.idof[3]]) (*g)[el.idof[3]] -= g_factor * dx[0];
      if (fixed.empty() || !fixed[el.idof[4]]) (*g)[el.idof[4]] -= g_factor * dx[1];
      if (fixed.empty() || !fixed[el.idof[5]]) (*g)[el.idof[5]] -= g_factor * dx[2];
    }
  }


  void BarAndHinge::bending(const Vector& coords, const Element& el, double* e, Vector* g) const {
    Vector x1(coords.cbegin()+el.idof[0], coords.cbegin()+el.idof[2]+1);
    Vector x2(coords.cbegin()+el.idof[3], coords.cbegin()+el.idof[5]+1);
    Vector x3(coords.cbegin()+el.idof[6], coords.cbegin()+el.idof[8]+1);
    Vector x4(coords.cbegin()+el.idof[9], coords.cbegin()+el.idof[11]+1);

    // Compute bond vectors
    auto b1 = x2 - x1;
    auto b2 = x3 - x2;
    auto b3 = x4 - x3;
    double b2m = vec::norm(b2);

    // Compute normal vectors
    auto n1 = vec::crossProduct(b1, b2);
    auto n2 = vec::crossProduct(b2, b3);
    double n1sq = vec::dotProduct(n1, n1);
    double n2sq = vec::dotProduct(n2, n2);
    double n12m = sqrt(n1sq * n2sq);

    // Compute cos
    double c = vec::dotProduct(n1, n2) / n12m;
    // E = k (1 - cos(\theta - \theta_0))
    // double s = b2m / n12m * vec::dotProduct(n1, b3);
    // double c0 = cos(el.parameters[1]);
    // double s0 = sin(el.parameters[1]);
    // E = k/2 (\theta - \theta_0)^2
    double n1b3 = vec::dotProduct(n1, b3);
    double theta = (n1b3>=0) ? acos(c) : 2*pi()-acos(c);
    double dtheta = theta - el.parameters[1];
    dtheta = fmod(dtheta + pi(), 2*pi()) - pi();

    if (e != nullptr) {
      // E = k (1 - cos(\theta - \theta_0))
      // *e += el.parameters[0] * (1 - c*c0 + s*s0); // Double angle formula for cos(t - t0)
      // E = k/2 (\theta - \theta_0)^2
      *e += el.parameters[0] * 0.5 * pow(theta-el.parameters[1], 2);
    }
    if (g != nullptr) {
      // E = k (1 - cos(\theta - \theta_0))
      // double gfactor = el.parameters[0] * (s*c0 + c*s0);
      // E = k/2 (\theta - \theta_0)^2
      double gfactor = el.parameters[0] * (theta - el.parameters[1]);
      // Normal vectors with magnitude 1 / triangle height
      auto n1h = b2m/n1sq * n1;
      auto n2h = b2m/n2sq * n2;
      // Quantify triangular skew, 0.5 if symmetric
      double skew1 = -vec::dotProduct(b1, b2) / (b2m*b2m);
      double skew2 = -vec::dotProduct(b3, b2) / (b2m*b2m);
      if (fixed.empty() || !fixed[el.idof[0]])  (*g)[el.idof[0]] += -gfactor * n1h[0];
      if (fixed.empty() || !fixed[el.idof[1]])  (*g)[el.idof[1]] += -gfactor * n1h[1];
      if (fixed.empty() || !fixed[el.idof[2]])  (*g)[el.idof[2]] += -gfactor * n1h[2];
      if (fixed.empty() || !fixed[el.idof[3]])  (*g)[el.idof[3]] +=  gfactor * ((1-skew1)*n1h[0] - skew2*n2h[0]);
      if (fixed.empty() || !fixed[el.idof[4]])  (*g)[el.idof[4]] +=  gfactor * ((1-skew1)*n1h[1] - skew2*n2h[1]);
      if (fixed.empty() || !fixed[el.idof[5]])  (*g)[el.idof[5]] +=  gfactor * ((1-skew1)*n1h[2] - skew2*n2h[2]);
      if (fixed.empty() || !fixed[el.idof[6]])  (*g)[el.idof[6]] += -gfactor * ((1-skew2)*n2h[0] - skew1*n1h[0]);
      if (fixed.empty() || !fixed[el.idof[7]])  (*g)[el.idof[7]] += -gfactor * ((1-skew2)*n2h[1] - skew1*n1h[1]);
      if (fixed.empty() || !fixed[el.idof[8]])  (*g)[el.idof[8]] += -gfactor * ((1-skew2)*n2h[2] - skew1*n1h[2]);
      if (fixed.empty() || !fixed[el.idof[9]])  (*g)[el.idof[9]]  += gfactor * n2h[0];
      if (fixed.empty() || !fixed[el.idof[10]]) (*g)[el.idof[10]] += gfactor * n2h[1];
      if (fixed.empty() || !fixed[el.idof[11]]) (*g)[el.idof[11]] += gfactor * n2h[2];
    }
  }


  void BarAndHinge::forceEnergy(const Vector& coords, const Element& el, double* e, Vector* g) const {
    if (fixed.empty() || !fixed[el.idof[0]]) {
      if (e != nullptr) (*e) -= force[2] * coords[el.idof[0]];
      if (g != nullptr) (*g)[el.idof[0]] -= force[2];
    }
  }


  BarAndHinge& BarAndHinge::setFixed(const std::vector<bool>& fixed) {
    this->fixed = fixed;
    return *this;
  }

}
