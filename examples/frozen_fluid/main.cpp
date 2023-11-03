#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

using namespace minim;
using std::vector;

double PI = acos(-1);

int nx = 100;
int ny = 100;
int nz = 1;

// Solid parameters
double solidHeight = 9.5;
double confinementStrength = 10; // units: energy per unit area

// Liquid parameters
double surfaceTension = 1; // units: energy per unit area
double contactAngle = 150;
double dropRad = 30;


double gammaDiff = surfaceTension * cos(contactAngle*PI/180);
vector<double> surfaceTensions = {surfaceTension-gammaDiff/2, surfaceTension+gammaDiff/2, surfaceTension};


vector<double> initialiseSolid() {
  vector<double> data(3*nx*ny*nz);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
  for (int z=0; z<nz; z++) {
    int i = x*ny*nz + y*nz + z;

    if (y < solidHeight) {
      data[3*i+0] = 1;
    }

  }
  }
  }
  return data;
}


vector<double> initialiseFluid(vector<double> solid) {
  vector<double> data(3*nx*ny*nz);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
  for (int z=0; z<nz; z++) {
    int i = x*ny*nz + y*nz + z;

    // Initialise a spherical cap
    double y0 = solidHeight - dropRad * cos(contactAngle*PI/180);
    double r2 = pow(x - 0.5*(nx-1), 2) + pow(y - y0, 2);
    if (r2 < pow(dropRad,2)) {
      // Liquid
      data[3*i+1] = 1;
    } else {
      // Gas
      data[3*i+2] = 1;
    }

    // Incorporate the solid
    data[3*i+0] = solid[3*i];
    data[3*i+1] *= (1 - solid[3*i]);
    data[3*i+2] *= (1 - solid[3*i]);
  }
  }
  }
  return data;
}


vector<double> getTotVolumes(vector<double> coords) {
  vector<double> volumes(3);
  for (int i=0; i<nx*ny*nz; i++) {
    volumes[0] += coords[3*i+0];
    volumes[1] += coords[3*i+1];
    volumes[2] += coords[3*i+2];
  }
  return volumes;
}


void output(vector<double> data, std::string filename) {
  if (mpi.rank == 0) {
    std::ofstream f(filename);
    for (auto x: data) {
      f << x << std::endl;
    }
    f.close();
  }
}


void energyComponents(State& state) {
  double e0 = state.componentEnergy(0);
  double e2 = state.componentEnergy(2);
  double e3 = state.componentEnergy(3);
  double e4 = state.componentEnergy(4);
  double e_pv = 0;
  state.pot->blockEnergyGradient(state.blockCoords(), state.comm, &e_pv, nullptr);
  e_pv = state.comm.sum(e_pv);
  print("Interfaces:", e0, "Surfaces:", e2, "Force:", e3, "Confinement:", e4, "Pressure/volume:", e_pv);
}



int main(int argc, char** argv) {
  mpi.init(&argc, &argv);

  auto log = [](int iter, State& state) {
    if (iter%100 == 0) {
      print("ITER:", iter, "ENERGY:", state.energy(), "GRAD:", vec::norm(state.gradient()));
    }
  };

  // Set up minimiser and potential
  Lbfgs min;
  min.setLinesearch("none");
  PFWetting pot;
  pot.setNFluid(3);
  pot.setGridSize({nx, ny, nz});
  pot.setSurfaceTension(surfaceTensions);
  // Set a simulation boundary (optional)
  pot.setSolid([](int x, int y, int z) { return (x==0 || x==nx-1 || y==0 || y==ny-1); });

  print("Step 1: Evolving solid");
  pot.setConfinement({confinementStrength, 0, 0});
  pot.setDensityConstraint("none");
  State state1(pot, initialiseSolid());
  min.setMaxIter(100);
  auto solid = min.minimise(state1, log);

  print("Step 2: Evolving fluid");
  auto initFluid = initialiseFluid(solid);
  pot.setFixFluid(0);
  pot.setConfinement({});
  pot.setDensityConstraint("hard");
  pot.setVolume(getTotVolumes(initFluid), 1e-2);
  State state2(pot, initFluid);
  min.setMaxIter(10000);
  auto minimum = min.minimise(state2, log);
  output(minimum, "minimum.txt");

  return 0;
}
