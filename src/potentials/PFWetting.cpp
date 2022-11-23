#include "PFWetting.h"

#include <stdexcept>
#include "State.h"


namespace minim {
  typedef std::vector<double> Vector;


  PFWetting::blockEnergyGradient(const Vector& coords, double* e, Vector* g) const {
    if (e) *e = 0;
    if (g) *g = Vector(coords.size());

    for (auto el : elements) {
      elementEnergyGradient(el, coords, e, g);
    }

    // Compute the gradient of the halo energy elements
    if (g) {
      for (auto el : elements_halo) {
        elementEnergyGradient(el, coords, nullptr, g);
      }
    }
  }


  PFWetting::elementEnergyGradient(const Element el, const Vector& coords, double* e, Vector* g) const {
    switch (el.type) {
      case 0:
        vol = el.parameters[0];
        // Bulk energy
        phi = coords[el.idof[0]];
        if (e) *e += 0.25/epsilon * (pow(phi,4) - 2*pow(phi,2) + 1) * vol;
        if (g) (*g)[el.idof[0]] += 1/epsilon * (pow(phi,3) - phi) * vol;
        // Gradient energy
        double diffzm = phi - coords[el.idof[1]];
        double diffzp = phi - coords[el.idof[2]];
        double diffym = phi - coords[el.idof[3]];
        double diffyp = phi - coords[el.idof[4]];
        double diffxm = phi - coords[el.idof[5]];
        double diffxp = phi - coords[el.idof[6]];
        double gradx2 = pow(diffxm, 2) + pow(diffxp, 2);
        double grady2 = pow(diffym, 2) + pow(diffyp, 2);
        double gradz2 = pow(diffzm, 2) + pow(diffzp, 2);
        double grad2 = gradx2 + grady2 + gradz2;
        double factor = 0.5 * epsilon * vol;
        if (e) *e += factor * grad2;
        if (g) {
          (*g)[el.idof[0]] += factor * (diffxm + diffxp + diffym + diffyp + diffzm + diffzp);
          (*g)[el.idof[1]] -= factor * diffzm;
          (*g)[el.idof[2]] -= factor * diffzp;
          (*g)[el.idof[3]] -= factor * diffym;
          (*g)[el.idof[4]] -= factor * diffyp;
          (*g)[el.idof[5]] -= factor * diffxm;
          (*g)[el.idof[6]] -= factor * diffxp;
        }
        break;

      case 1:
        // Surface energy
        // parameter[0]: Surface area
        // parameter[1]: Wetting parameter sqrt(2)cos(theta)
        phi = coords[el.idof[0]];
        if (e): *e += el.parameters[1]/6 * (pow(phi,3) - 3*phi - 2) * el.parameters[0];
        if (g): *g[el.idof[0]] += el.parameters[0]*0.5 * (pow(phi,2) - 1) * el.parameters[0];
        break;

      case 2:
        // External force
        // Parameters:
        //   0: Volume
        //   1: Magnitude of force on component 1
        //   2: Magnitude of force on component 2
        //   3-5: Direction force on component 1
        //   6-8: Direction force on component 2
        phi = coords[el.idof[0]];
        double f1 = el.parameters[1];
        double f2 = el.parameters[2];
        Vector f1Norm = {el.parameters[3], el.parameters[4], el.parameters[5]};
        Vector f2Norm = {el.parameters[6], el.parameters[7], el.parameters[8]};
        Vector coord = getCoord(el.idof[0]);
        double h1 = - vec::dot_product(coord, f1Norm);
        double h2 = - vec::dot_product(coord, f2Norm);
        if (e) *e += 0.5 * ((1+phi)*f1*h1 + (1-phi)*f2*h2) * el.parameters[0];
        if (g) (*g)[el.idof[0]] += 0.5 * (f1*h1 - f2*h2) * el.parameters[0];
        break;

      default:
        std::invalid_argument("Unknown energy element type.");
    }
  }


  Vector PFWetting::getCoord(int i) const {
    int z = i % gridSize[2];
    int y = (i - z) / gridSize[2] % ly;
    int x = i / (gridSize[1] * gridSize[2]);
    return {x, y, z};
  }


  PFWetting::newState(const Vector& coords) {
    return State(*this, coords);
  }

}
