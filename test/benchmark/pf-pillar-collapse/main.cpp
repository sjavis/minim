#include "minim.h"
#include "utils/vec.h"
#include <algorithm>
#include <fstream>
#include <sstream>

using namespace minim;
using std::vector;

double PI = acos(-1);


int logIter = -100; // <=0 for no log
bool saveCoords = false;

// Geometry
int nx = 80;
int ny = 100;
int nz = 1;
double res = 1e-3;
double yTop = 20;
double yBottom = 9.5;

// Liquid parameters
double surfaceTension = 100;
double maxPressure = 0.1*surfaceTension/res;

// Run parameters
vector<double> spacings = {40};
vector<double> contactAngles = {30, 50, 70, 90, 110, 130, 150};
std::string directory;


// Initialise

vector<char> initialiseSolid(double spacing) {
  vector<char> data(nx*ny*nz, false);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
  for (int z=0; z<nz; z++) {
    int i = x*ny*nz + y*nz + z;

    double pillarWidth = nx - spacing;
    bool isBase = (y < yBottom);
    bool isPillar = (y < ny-1-yTop) && ((x <= pillarWidth/2) || (x > nx-1-pillarWidth/2));
    if (isBase || isPillar) {
      data[i] = true;
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

    if (y > ny-1-yTop) {
      data[3*i+1] = 1; // Liquid
    } else {
      data[3*i+2] = 1; // Gas
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


// Utilities

void output(vector<double> data, std::string filename) {
  if (mpi.rank == 0) {
    std::ofstream f(filename);
    for (auto x: data) {
      f << x << std::endl;
    }
    f.close();
  }
}

bool checkFailed(const vector<double>& coords) {
  // Check for liquid near bottom
  int y = yBottom + 5;
  int x = (nx-1) / 2;
  int k = x*ny + y;
  return (coords[3*k+1]>0.5);
}


void logFn(int iter, State& state) {
  if (checkFailed(state.coords())) state.failed();
  if (iter%logIter == 0 && logIter>0) {
    print("ITER:", iter, "ENERGY:", state.energy(), "GRAD:", vec::norm(state.gradient()));
  }
};


// Run

bool evolveFluid(PhaseField pot, double pressure, std::vector<double> init) {
  pot.setPressure({0, pressure, 0});
  pot.setFixFluid(0);
  pot.setDensityConstraint("hard");
  Lbfgs min;
  min.setLinesearch("none");
  min.setMaxIter(5000);

  State state(pot, init);
  auto minimum = min.minimise(state, logFn);
  if (saveCoords) output(minimum, directory+"/p-"+std::to_string(pressure)+".txt");
  return state.isFailed;
}


bool checkComplete(double pressure[2]) {
  double pressureMean = (pressure[0] + pressure[1]) / 2;
  double pressureErr = abs(pressure[1] - pressure[0]) / 2;
  if (pressureErr > 1e-3*maxPressure) return false;

  std::string filename = directory + "/critical_pressure.txt";
  std::ofstream f(filename);
  f << pressureMean << " ± " << pressureErr << std::endl;
  f << pressure[0] << " " << pressure[1] << std::endl;
  f.close();
  print("Critical Pressure:", pressureMean, "±", pressureErr);
  return true;
}


void run(double spacing, double contactAngle) {
  print("\nRes:", res, "Surface Tension:", surfaceTension, "Spacing:", spacing, "Contact Angle:", contactAngle);

  // Make potential
  directory = (std::stringstream() << "outputs/res-"<<res<<"/st-"<<surfaceTension<<"/s-"<<spacing<<"/ca-"<<contactAngle).str();
  if (mpi.rank==0) {
    int result = system(("mkdir -p " + directory).c_str());
    if (result > 0) throw std::logic_error("Error making output directory");
  }

  double gammaDiff = surfaceTension * cos(contactAngle*PI/180);
  vector<double> surfaceTensions = {surfaceTension-gammaDiff/2, surfaceTension+gammaDiff/2, surfaceTension};

  // Set up potential
  PhaseField pot;
  pot.setNFluid(3);
  pot.setGridSize({nx, ny, nz});
  pot.setResolution(res);
  pot.setSurfaceTension(surfaceTensions);
  pot.setSolid([](int x, int y, int z) { return (y==0 || y==ny-1); });

  // Initialise
  print("Initialising");
  auto solid = pot.diffuseSolid(initialiseSolid(spacing));
  auto init = initialiseFluid(solid);

  // Critical potential binary search
  double pressureBounds[2] = {-maxPressure, maxPressure};
  for (int iPressure=0; iPressure<50; iPressure++) {
    // Get new pressure
    double pressure = 0.5 * (pressureBounds[0] + pressureBounds[1]);
    print("Pressure:", pressure, "\tBounds:", pressureBounds[0], pressureBounds[1]);
    // Run
    int status = evolveFluid(pot, pressure, init);
    // Update bounds
    pressureBounds[status] = pressure;
    if (checkComplete(pressureBounds)) break;
  }
}


int main(int argc, char** argv) {
  mpi.init(&argc, &argv);

  for (auto spacing: spacings) {
    for (auto contactAngle: contactAngles) {
      run(spacing, contactAngle);
    }
  }

  return 0;
}
