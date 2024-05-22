#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace minim;
using std::vector;

double radius = 20;
int nx = 4*radius;
int ny = 4*radius;


vector<double> initialiseFluid() {
  vector<double> data(2*nx*ny);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
    int i = x*ny + y;

    double x0 = (nx - 1) / 2.0;
    double y0 = (ny - 1) / 2.0;
    double r = sqrt(pow(x-x0, 2) + pow(y-y0, 2));
    if (r < radius) {
      data[2*i+0] = 1; // Liquid
    } else {
      data[2*i+1] = 1; // Gas
    }
  }
  }
  return data;
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

double run(double volConst) {
  double surfaceTension = 1e4;
  double resolution = 1e-8;

  PhaseField pf;
  pf.setNFluid(2);
  pf.setGridSize({nx, ny, 1});
  pf.setResolution(resolution);
  pf.setSurfaceTension(surfaceTension);

  pf.setDensityConstraint("hard");
  pf.setVolumeFixed(true, volConst);

  auto init = initialiseFluid();
  State state(pf, init);
  state.convergence = 1e-8*surfaceTension*resolution*resolution;

  Lbfgs min;
  min.setLinesearch("none");
  min.setMaxIter(1000);
  auto minimum = min.minimise(state);

  double v1 = 0;
  double v2 = 0;
  double v01 = 0;
  double v02 = 0;
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
    int i = x*ny + y;
    v1 += minimum[2*i+0];
    v2 += minimum[2*i+1];
    v01 += init[2*i+0];
    v02 += init[2*i+1];
  }
  }
  double volError = std::max(abs(v1-v01)/v01, abs(v2-v02)/v02);
  print("Volume constraint", volConst, ": Error =", volError, ", Iterations =", min.iter);
  return volError;
}

int main(int argc, char** argv) {
  mpi.init(&argc, &argv);

  if (argc>1 && atoi(argv[1])==1) {
    if (run(0.01) > 0.003) return 1;

  } else {
    for (auto volConst : {0.1, 0.01, 0.001}) {
      run(volConst);
    }
  }


  return 0;
}
