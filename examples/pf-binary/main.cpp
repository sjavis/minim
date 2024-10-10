#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

using namespace minim;
using std::vector;

int nx = 120;
int ny = 60;


vector<double> initCoords() {
  // Make a semi-circular droplet on the solid
  double dropRad = 30;
  vector<double> coords(nx*ny, -1);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
    int i = ny*x+y;
    double r = sqrt(pow(x-0.5*(nx-1), 2) + pow(y-1, 2));
    if (r < dropRad) coords[i] = 1;
  }
  }
  return coords;
}


void output(vector<double> data) {
  if (mpi.rank == 0) {
    std::ofstream f("output.txt");
    for (auto x: data) {
      f << x << std::endl;
    }
    f.close();
  }
}


int main(int argc, char** argv) {
  mpi.init(&argc, &argv);

  PhaseField pot;
  pot.setGridSize({nx, ny, 1});
  pot.setSolid([](int x, int y, int z){ return (y==0); });
  pot.setContactAngle(60);
  pot.setVolumeFixed(true);
  State state(pot, initCoords());

  Lbfgs min;
  min.setMaxIter(5000);
  min.minimise(state, "eg-100");

  output(state.coords());

  return 0;
}
