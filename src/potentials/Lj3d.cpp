#include "Lj3d.h"

#include <math.h>
#include "State.h"
#include "Potential.h"

typedef std::vector<double> Vector;


Lj3dArgs::Lj3dArgs(int ndof, double sigma, double epsilon)
  : Args(ndof), n_particle(ndof/3), sigma(sigma), epsilon(epsilon) {

  // Generate energy elements
  int id = 0;
  for (int i=0; i<n_particle; i++) {
    for (int j=i+1; j<n_particle; j++) {
      Args::Element el = {id, 0, {3*i, 3*i+1, 3*i+2, 3*j, 3*j+1, 3*j+2}};
      elements.push_back(el);
      id++;
    }
  }
}


double Lj3d::energy(const Vector &coords, const Args &args_tmp) {
  const Lj3dArgs &args = static_cast<const Lj3dArgs&> (args_tmp);
  double energy = 0;

  for (Args::Element el : args.elements) {
    double dx = coords[el.idof[0]] - coords[el.idof[3]];
    double dy = coords[el.idof[1]] - coords[el.idof[4]];
    double dz = coords[el.idof[2]] - coords[el.idof[5]];

    double r2 = dx*dx + dy*dy + dz*dz;
    double lj6 = pow(args.sigma, 6) / pow(r2, 3);
    energy += 4*args.epsilon * (lj6*lj6 - lj6);
  }
  return energy;
}


Vector Lj3d::gradient(const Vector &coords, const Args &args_tmp) {
  const Lj3dArgs &args = static_cast<const Lj3dArgs&> (args_tmp);
  Vector g(coords.size());

  int ne1 = args.elements.size();
  int ne2 = args.elements_halo.size();
  for (int ie=0; ie<(ne1+ne2); ie++) {
    auto el = (ie<ne1) ? args.elements[ie] : args.elements_halo[ie-ne1];

    double dx = coords[el.idof[0]] - coords[el.idof[3]];
    double dy = coords[el.idof[1]] - coords[el.idof[4]];
    double dz = coords[el.idof[2]] - coords[el.idof[5]];
    double r2 = dx*dx + dy*dy + dz*dz;
    double r = sqrt(r2);

    double lj6 = pow(args.sigma, 6) / pow(r2, 3);
    double dedr = -24 * args.epsilon * (2*lj6*lj6 - lj6) / r;

    g[el.idof[0]] += dx/r * dedr;
    g[el.idof[1]] += dy/r * dedr;
    g[el.idof[2]] += dz/r * dedr;
    g[el.idof[3]] -= dx/r * dedr;
    g[el.idof[4]] -= dy/r * dedr;
    g[el.idof[5]] -= dz/r * dedr;
  }
  return g;
}


Args* Lj3d::newArgs(int ndof) {
  return new Lj3dArgs(ndof);
}
