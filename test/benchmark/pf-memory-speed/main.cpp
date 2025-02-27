#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>

#include <sys/resource.h>
#include <malloc.h>
#include <chrono>

using namespace minim;
using std::vector;

int nx = 200;
int ny = 200;

vector<float> getMemoryUsage() {
  vector<float> memory(4); // MaxRSS, allocated heap, free heap, total heap
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  memory[0] = mpi.sum(usage.ru_maxrss / 1024.0) / mpi.size;
  struct mallinfo info = mallinfo();
  memory[1] = mpi.sum(info.uordblks / pow(1024.0, 2)) / mpi.size;
  memory[2] = mpi.sum(info.fordblks / pow(1024.0, 2)) / mpi.size;
  memory[3] = mpi.sum(info.arena / pow(1024.0, 2)) / mpi.size;
  return memory;
}


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

  auto time0 = std::chrono::high_resolution_clock::now();

  PhaseField pot;
  pot.setNFluid(3);
  pot.setGridSize({nx, ny, 1});
  pot.setSolid(getSolid());
  pot.setDensityConstraint("hard");
  pot.setVolumeFixed(true);
  State state(pot, initCoords());

  Lbfgs min;
  min.setMaxIter(5000);
  min.minimise(state, "eg-100");

  output(state.coords());

  auto time1 = std::chrono::high_resolution_clock::now();
  float timeElapsed = std::chrono::duration_cast<std::chrono::milliseconds>(time1 - time0).count() * 1e-3;
  auto memory = getMemoryUsage();
  print("TIME (s):", timeElapsed);
  print("MEMORY: maximum RSS (MB):", memory[0]);
  print("Allocated heap size (MB):", memory[1]);
  print("Free heap size      (MB):", memory[2]);
  print("Total heap size     (MB):", memory[3]);

  // Write the results to a file
  std::string filename = "scaling.txt";
  bool fileExists = std::ifstream(filename).good();
  std::ofstream f(filename, std::ios::app);
  if (mpi.rank == 0) {
    if (!fileExists) {
      f << "N TIME(S) MEM(MB)" << std::endl;
    }
    f << mpi.size << " " << timeElapsed << " " << memory[0] << std::endl;
  }
  f.close();

  return 0;
}
