#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

using namespace minim;
using std::vector;

int nx = 100;
int ny = 100;


vector<double> initCoords() {
  vector<double> data(3*nx*ny);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
    int i = ny*x+y;

    if (y<=ny/3) {
      data[3*i+0] = 1;
      data[3*i+1] = 0;
    } else if (x>=nx/2) {
      data[3*i+0] = 0;
      data[3*i+1] = 1;
    } else {
      data[3*i+0] = 0;
      data[3*i+1] = 0;
    }

    data[3*i+2] = 1 - data[3*i] - data[3*i+1];
  }
  }
  return data;
}

vector<char> getSolid() {
  vector<char> solid(nx*ny, false);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
    int i = ny*x+y;
    if (x==0 || x==nx-1 || y==0 || y==ny-1) solid[i] = true;
  }
  }
  return solid;
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

  PhaseField pot;
  pot.setNFluid(3);
  pot.setGridSize({nx, ny, 1});
  pot.setSolid(getSolid());
  pot.setDensityConstraint("hard");
  pot.setVolumeFixed(true);
  State state(pot, initCoords());

  Lbfgs min;
  min.setMaxIter(5000);
  min.minimise(state, log);

  output(state.coords());

  return 0;
}
