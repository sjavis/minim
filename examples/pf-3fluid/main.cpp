#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

using namespace minim;
using std::vector;

int nx = 100;
int ny = 100;
double c1 = 0.33;
double c2 = 0.33;


vector<double> initCoords(int n, double c1, double c2) {
  double randMax = 0.1;
  vector<double> rand(n);
  std::generate(rand.begin(), rand.end(), [](){ return (std::rand()%1000)/1000.0; });
  auto c1v = c1 + randMax * rand;
  std::generate(rand.begin(), rand.end(), [](){ return (std::rand()%1000)/1000.0; });
  auto c2v = c2 + randMax * rand;
  for (int x=0; x<100; x++) {
  for (int y=0; y<100; y++) {
    int i = 100*x+y;
    if (y<=33) {
      c1v[i] = 1;
      c2v[i] = 0;
    } else if (x>=50) {
      c1v[i] = 0;
      c2v[i] = 1;
    } else {
      c1v[i] = 0;
      c2v[i] = 0;
    }
  }
  }
  auto c3v = 1 - c1v - c2v;
  vector<double> data(3*n);
  for (int i=0; i<n; i++) {
    data[3*i] = c1v[i];
    data[3*i+1] = c2v[i];
    data[3*i+2] = c3v[i];
  }
  return data;
}

vector<bool> getSolid() {
  vector<bool> solid(10000, false);
  for (int x=0; x<100; x++) {
  for (int y=0; y<100; y++) {
    int i = 100*x+y;
    if (x==0 || x==99 || y==0 || y==99) solid[i] = true;
  }
  }
  return solid;
}


void output(vector<double> data) {
  char *filename;
  asprintf(&filename, "%.2f-%.2f.txt", c1, c2);

  if (mpi.rank == 0) {
    std::ofstream f(filename);
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

  auto pot = PhaseField();
  pot.setNFluid(3);
  pot.setGridSize({nx, ny, 1});
  pot.setSolid(getSolid());
  pot.setDensityConstraint("hard");
  pot.setVolumeFixed(true, 1e-4);
  State state(pot, initCoords(nx*ny, c1, c2));

  auto min = Lbfgs();
  min.setMaxIter(5000);
  min.minimise(state, log);

  output(state.coords());

  return 0;
}
