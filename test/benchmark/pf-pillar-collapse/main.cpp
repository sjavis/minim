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
int nx = 120;
int ny = 100;
int nz = 120;
double res = 1e-3;
double yTop = 30;
double yBottom = 9.5;

// Reentrant Geometry
int reentrant_height = 30;
double cone_angle = 80 * PI/180; 

// Pillar Geometry
int pillar_radius = 50;
int pillar_height = 21;
int pillar_spacing = 30;
int hole_radius = 30;
int total_pillars = 6;
int pillar_gaps[4][2] = {
  {1,1},
  {3,3},
  {5,2},
  {3,5}
};
int total_gaps = sizeof pillar_gaps / sizeof pillar_gaps[0];

// Torus Geometry
double torus_r = 2;
int torus_R = pillar_radius;
int total_curves = floor((2*yTop-yTop/2)/(2*torus_r));

// Liquid parameters
double surfaceTension = 100;
double maxPressure = 0.1*surfaceTension/res;

// Run parameters
vector<double> spacings = {2}; 
vector<double> contactAngles = {110};
std::string directory;
bool fixed_pressure = true;
vector<double> defined_pressure = {0};



// Initialise

vector<char> initialiseSolid(double spacing) {
  vector<char> data(nx*ny*nz, false);
  for (int x=0; x<nx; x++) {
  for (int y=0; y<ny; y++) {
  for (int z=0; z<nz; z++) {
    int i = x*ny*nz + y*nz + z;

    double pillarWidth = nx - spacing;

    // bool isBase = (y < yBottom);

    bool isBase = (y < ny-1-2*yTop && y > ny-1-2*yTop-yBottom);

    //bool isBase = ((y < ny-1-yTop) && (y > ny-1-yTop-2*yBottom) && !(sqrt(pow(x,2)+pow(z,2)) <= hole_radius) && !(sqrt(pow(x,2)+pow(nz-z,2)) <= hole_radius) && !(sqrt(pow(nx-x,2)+pow(nz/2-z,2)) <= hole_radius));

    // bool isPillar = (y < ny-1-yTop) && ((x <= pillarWidth/2) || (x > nx-1-pillarWidth/2)); //simple pillar

    //bool isPillar = (y < ny-1-yTop && y > ny-1-yTop-20) && ((x <= pillarWidth/2) || (x > nx-1-pillarWidth/2)) && ((z <= pillarWidth/2) || (z > nz-1-pillarWidth/2)); //posts

    //bool isPillar = (y < ny-1-yTop && y > ny-1-yTop-20) && ((sqrt(pow(x,2)+pow(z,2)) < pillar_radius) || (sqrt(pow(nx-x,2)+pow(z,2)) < pillar_radius) || (sqrt(pow(x,2)+pow(nz-z,2)) < pillar_radius) || (sqrt(pow(nx-x,2)+pow(nz-z,2)) < pillar_radius));


    //bool isPillar = (y < ny-1-yTop) && ((x <= pillarWidth/2 && z <= nz) || (x > nx-1-pillarWidth/2 && z <= nz) || (z <= pillarWidth/2 && x <= nx) || (z > nz-1-pillarWidth/2 && x <= nx));

    //bool isPillar = ((y < ny-1-yTop) && (y > ny-1-yTop-2*yBottom)) && ((x <= pillarWidth/2 && z <= nz) || (x > nx-1-pillarWidth/2 && z <= nz) || (z > nz-1-pillarWidth/2 && x <= nx) || (z <= pillarWidth/2 && x <= nx)); //half hole


    //Torus shape

    bool isPillar = (y <= ny-1-yTop/2 && y > ny-1-2*yTop) && ((sqrt(pow(x,2) + pow(z,2)) <= pillar_radius) || (sqrt(pow(x-nx,2) + pow(z-nz,2)) <= pillar_radius));


    // if (pow(sqrt(pow(x,2)+pow(z,2)) - torus_R,2) + pow(y-(ny-1-yTop-torus_r),2) <= pow(torus_r,2) || pow(sqrt(pow(x-nx,2)+pow(z-nz,2)) - torus_R,2) + pow(y-(ny-1-yTop-torus_r),2) <= pow(torus_r,2)){

    //   isPillar = false;
    // }

    for (int curve_count = 0; curve_count < total_curves; curve_count++){
    if (pow(sqrt(pow(x,2)+pow(z,2)) - torus_R,2) + pow(y-(ny-1-yTop/2-curve_count*2*torus_r),2) <= pow(torus_r,2) || pow(sqrt(pow(x-nx,2)+pow(z-nz,2)) - torus_R,2) + pow(y-(ny-1-yTop/2-curve_count*2*torus_r),2) <= pow(torus_r,2)){

      isPillar = false;
    }
    }

    // for (int curves = 0; curves < total_curves; curves++){

    // if (pow(sqrt(pow(x,2)+pow(z,2)) - torus_R,2) + pow(y-(ny-1-yTop-torus_r*curves),2) <= pow(torus_r,2) || pow(sqrt(pow(x-nx,2)+pow(z-nz,2)) - torus_R,2) + pow(y-(ny-1-yTop-torus_r*curves),2) <= pow(torus_r,2)){

    //   isPillar = false;

    // }
    // }

    // bool isPillar = (y <= ny-1-yTop && y >ny-1-yTop-50) && ((sqrt(pow(x,2) + pow(z,2)) <= pillar_radius) || (sqrt(pow(x-nx,2) + pow(z-nz,2)) <= pillar_radius));

    // if (pow(sqrt(pow(x,2)+pow(z,2)) - torus_R,2) + pow(y-ny-1-yTop-25,2) <= pow(torus_r,2)){

    //   isPillar = false;

    // }


    //Circular pillar array for n number of pillars per axis
    // bool isPillar = false;
    // for (int pillar_x = 0; pillar_x <= total_pillars; pillar_x++){
    // for (int pillar_z = 0; pillar_z <= total_pillars; pillar_z++){

    // if (sqrt(pow(x + nx/(2*total_pillars) - pillar_x*nx/total_pillars,2)+pow(z + nz/(2*total_pillars)- pillar_z*nz/total_pillars,2)) <= pillar_radius && (y <= ny-1-yTop+pillar_height && y >= ny-1-yTop)){

    //   isPillar = true;
    //   break;
    //     }

    // }
    // }

    // for (int no_pillar = 0; no_pillar < total_gaps; no_pillar++){
    //   if (sqrt(pow(x + nx/(2*total_pillars) - pillar_gaps[no_pillar][0]*nx/total_pillars,2)+pow(z + nz/(2*total_pillars)- pillar_gaps[no_pillar][1]*nz/total_pillars,2)) <= pillar_radius && (y <= ny-1-yTop+pillar_height && y >= ny-1-yTop)){

    //   isPillar = false;
    //   break;
    //     }

    // }


    //for (int pillar_array; pillar_array < pillar_number; pillar_array++)
      
      //!((x <= pillar_radius && z <= pillar_radius && sqrt(pow(x,2) + pow(z,2)) <= pillar_radius) &&
      //((x >= nx-pillar_radius && x <= nx) && z <= pillar_radius && sqrt(pow(nx-x,2) + pow(z,2)) <= pillar_radius) &&
      //((x >= nx/2-pillar_radius && x <= nx/2 + pillar_radius) && z >= nz-pillar_radius && sqrt(pow(nx/2-x,2) + pow(nz-z,2)) <= pillar_radius));

    //bool isPillar = ((y < ny-1-yTop) && ((x <= pillarWidth/2) || (x > nx-1-pillarWidth/2))) ||
    //                ((y < -tan(cone_angle)*(x-pillarWidth/2) + (ny-1-yTop)) && (x >= pillarWidth/2 && x <= nx/2)) ||
    //                ((y < tan(cone_angle)*(x-nx+pillarWidth/2) + (ny-1-yTop)) && (x >= nx/2 && x <= nx-1-pillarWidth/2)); //cone

                    
    //bool isPillar = ((((y > reentrant_height) && (y < ny-1-yTop)) && ((x <= pillarWidth/2) || (x > nx-1-pillarWidth/2)))
     //                 || (((y <= reentrant_height) && (y > yBottom)) && ((x <= pillarWidth/4) || (x > nx-1-pillarWidth/4))) 
     //                 || (y <= yBottom)); //reentrant geometry



    if (isPillar || isBase) {
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
  int y = ny-1-2*yTop+4;
  int x = nx/2;
  int z = nx/2;
  // int k = x*ny + y;
  int k = x*ny*nz +y*nz + z;
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
  pot.setSolid([](int x, int y, int z) { return (y==0 || y==ny-1 || x==0 || x==nx-1 || z==0 || z==nz-1); });
  // pot.setSolid([](int x, int y, int z) { return (y==0 || y==ny-1); });

  // Initialise
  print("Initialising");
  auto solid = pot.diffuseSolid(initialiseSolid(spacing));
  auto init = initialiseFluid(solid);

  // Critical potential binary search
  if (fixed_pressure == false) {
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
    } } 
  
  else {
      for (int  i = 0; i < defined_pressure.size(); i++) {
        int pressure = defined_pressure[i];
        
        print("Pressure:", pressure);
        int status = evolveFluid(pot, pressure, init);

      }





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
