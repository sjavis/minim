#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

using namespace minim;
using std::vector;

int nx = 100;
int ny = 100;


vector<double> initCoords() {
  // Make a semi-circular droplet on the solid
  double dropRad = 30;
  vector<double> coords(nx*ny);
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

  auto log = [](int iter, State& state) {
    if (iter%100 == 0) {
      print("ITER:", iter, "ENERGY:", state.energy(), "GRAD:", vec::norm(state.gradient()));
    }
  };

  PFWetting pot;
  pot.setGridSize({nx, ny, 1});
  pot.setSolid([](int x, int y, int z){ return (y==0); });
  pot.setVolumeFixed(true, 1e-4);
  State state(pot, initCoords());

  Lbfgs min;
  min.setMaxIter(5000);
  min.minimise(state, log);

  output(state.coords());

  return 0;
}
