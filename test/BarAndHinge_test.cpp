#include "test_main.cpp"
#include "BarAndHinge.h"

#include <math.h>
#include "State.h"
#include "Lbfgs.h"
#include "utils/mpi.h"
#include "utils/print.h"

using namespace minim;
typedef std::vector<double> Vector;

TEST(BarAndHingeTest, TestMin) {
  double pi = 4*atan(1);

  BarAndHinge pot;
  pot.setTriangulation({{0,1,2}, {1,3,2}});
  pot.setRigidity(100, 1);
  State initState(pot, {-1,0,1, 0,-sqrt(2),0, 0,sqrt(2),0, 1,0,1});
  int iBend;
  for (int i=0; i<(int)initState.pot->elements.size(); i++) {
    if (initState.pot->elements[i].type==1) iBend = i;
  }
  Lbfgs min;

  State s1 = initState;
  if (mpi.rank==0) s1.pot->elements[iBend].parameters[1] = 1e-6;
  Vector x1 = min.minimise(s1);
  double dz = sqrt(2)/2;
  Vector x1_expected = {0,0,0.5+dz, 0,-sqrt(2),0.5-dz, 0,sqrt(2),0.5-dz, 0,0,0.5+dz};
  EXPECT_TRUE(ArraysNear(x1, x1_expected, 1e-4));

  State s2 = initState;
  if (mpi.rank==0) s2.pot->elements[iBend].parameters[1] = pi/2;
  Vector x2 = min.minimise(s2);
  Vector x2_expected = {-1,0,1, 0,-sqrt(2),0, 0,sqrt(2),0, 1,0,1};
  EXPECT_TRUE(ArraysNear(x2, x2_expected, 1e-4));

  State s3 = initState;
  if (mpi.rank==0) s3.pot->elements[iBend].parameters[1] = pi;
  Vector x3 = min.minimise(s3);
  Vector x3_expected = {-sqrt(2),0,0.5, 0,-sqrt(2),0.5, 0,sqrt(2),0.5, sqrt(2),0,0.5};
  EXPECT_TRUE(ArraysNear(x3, x3_expected, 1e-4));

  State s4 = initState;
  if (mpi.rank==0) s4.pot->elements[iBend].parameters[1] = 3*pi/2;
  Vector x4 = min.minimise(s4);
  Vector x4_expected = {-1,0,0, 0,-sqrt(2),1, 0,sqrt(2),1, 1,0,0};
  EXPECT_TRUE(ArraysNear(x4, x4_expected, 1e-4));

  State s5 = initState;
  if (mpi.rank==0) s5.pot->elements[iBend].parameters[1] = 2*pi-1e-6;
  Vector x5 = min.minimise(s5);
  Vector x5_expected = {0,0,0.5-dz, 0,-sqrt(2),0.5+dz, 0,sqrt(2),0.5+dz, 0,0,0.5-dz};
  EXPECT_TRUE(ArraysNear(x5, x5_expected, 1e-4));
}
