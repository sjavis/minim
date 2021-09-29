#include <vector>
#include <math.h>
#include "minim.h"

#include <chrono>
#include <algorithm>

using namespace std::chrono;

int main(int argc, char **argv) {
  mpiInit(&argc, &argv);
  minim::print();

  Lj3d potential = Lj3d();
  State state = potential.newState({0,0,0, 2,0,0, 1,1,0, 5,1,1});
  state.convergence = 1e-4;

  Lbfgs min = Lbfgs(state);
  //Fire min = Fire(state);
  //GradDescent min = GradDescent(state).setMaxIter(100000);
  //Anneal min = Anneal(state, 1, 0.0001).setMaxIter(1000000);
  auto result = min.minimise();

  minim::print("Complete after", min.iter, "iterations.");
  minim::print(result[0], result[1], result[2]);
  minim::print(result[3], result[4], result[5]);
  minim::print(result[6], result[7], result[8]);
  minim::print(result[9], result[10], result[11]);
  minim::print("Separations:");
  minim::print(sqrt(pow(result[0]-result[3], 2) + pow(result[1]-result[4] , 2) + pow(result[2]-result[5] , 2)));
  minim::print(sqrt(pow(result[0]-result[6], 2) + pow(result[1]-result[7] , 2) + pow(result[2]-result[8] , 2)));
  minim::print(sqrt(pow(result[0]-result[9], 2) + pow(result[1]-result[10], 2) + pow(result[2]-result[11], 2)));
  minim::print(sqrt(pow(result[3]-result[6], 2) + pow(result[4]-result[7] , 2) + pow(result[5]-result[8] , 2)));
  minim::print(sqrt(pow(result[3]-result[9], 2) + pow(result[4]-result[10], 2) + pow(result[5]-result[11], 2)));
  minim::print(sqrt(pow(result[6]-result[9], 2) + pow(result[7]-result[10], 2) + pow(result[8]-result[11], 2)));


  minim::print();
  minim::print("Measuring minimiser speed...");
  srand(2);//time(NULL));
  int nruns = 1000;
  std::vector<double> init(12);
  auto start = high_resolution_clock::now();
  for (int i=0; i<nruns; i++) {
    std::generate(init.begin(), init.end(), [](){ return rand()%1000/100.-5; });
    state.comm.bcast(init);
    min.state.setCoords(init);
    auto result = min.minimise();
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  minim::print("Time taken by minimiser:", duration.count()/nruns, "microseconds");
}
