#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace minim;
using std::vector;

int nx = 60;


vector<double> initialiseFluid() {
  vector<double> data(2*nx);
  for (int x=0; x<nx; x++) {
    if (x > (nx-1)/2.0) {
      data[2*x+0] = 1;
    } else {
      data[2*x+1] = 1;
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


void run(double width) {
  print("Width:", width);
  double surfaceTension = 1e4;
  double resolution = 1e-8;

  PhaseField pf;
  pf.setNFluid(2);
  pf.setGridSize({nx, 1, 1});
  pf.setResolution(resolution);
  pf.setSurfaceTension(surfaceTension);
  pf.setInterfaceSize(width*resolution);
  pf.setSolid([](int x, int y, int z){ return (x==0 || x==nx-1);});

  State state(pf, initialiseFluid());

  Lbfgs min;
  min.setLinesearch("none");
  min.setMaxIter(1000);
  auto minimum = min.minimise(state);

  std::string filename = (std::stringstream() << "width-" << width <<".txt").str();
  output(minimum, filename);
}


int main(int argc, char** argv) {
  mpi.init(&argc, &argv);

  if (argc != 2) return 1;
  std::istringstream ss(argv[1]);

  double width;
  while (ss >> width) {
    run(width);
  }

  return 0;
}
