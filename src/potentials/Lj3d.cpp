#include <math.h>
#include "Lj3d.h"
#include "State.h"


Lj3dArgs::Lj3dArgs(int ndof, double sigma, double epsilon)
  : Args(ndof), n_particle(ndof/3), sigma(sigma), epsilon(epsilon)
{}


double Lj3d::energy(std::vector<double> coords, Args &args_tmp) {
  Lj3dArgs &args = static_cast <Lj3dArgs&> (args_tmp); 

  double energy = 0;
  for (int i=0; i<args.n_particle; i++) {
    for (int j=i+1; j<args.n_particle; j++) {
      double dx = coords[3*i] - coords[3*j];
      double dy = coords[3*i+1] - coords[3*j+1];
      double dz = coords[3*i+2] - coords[3*j+2];

      double r2 = dx*dx + dy*dy + dz*dz;
      double lj6 = pow(args.sigma, 6) / pow(r2, 3);
      energy += 4*args.epsilon * (lj6*lj6 - lj6);
    }
  }
  return energy;
}


std::vector<double> Lj3d::gradient(std::vector<double> coords, Args &args_tmp) {
  Lj3dArgs &args = static_cast <Lj3dArgs&> (args_tmp); 

  std::vector<double> g(args.ndof);
  for (int i=0; i<args.n_particle; i++) {
    for (int j=i+1; j<args.n_particle; j++) {
      double dx = coords[3*i] - coords[3*j];
      double dy = coords[3*i+1] - coords[3*j+1];
      double dz = coords[3*i+2] - coords[3*j+2];
      double r2 = dx*dx + dy*dy + dz*dz;

      double lj6 = pow(args.sigma, 6) / pow(r2, 3);
      double dedr = 12 * args.epsilon * (-lj6*lj6 + lj6) / r2;

      g[3*i] += 2 * dx * dedr;
      g[3*i+1] += 2 * dy * dedr;
      g[3*i+2] += 2 * dz * dedr;
      g[3*j] -= 2 * dx * dedr;
      g[3*j+1] -= 2 * dy * dedr;
      g[3*j+2] -= 2 * dz * dedr;
    }
  }
  return g;
}


State Lj3d::newState(std::vector<double> coords) {
  Args *args = new Lj3dArgs(coords.size());
  return State(*this, coords, *args);
}
