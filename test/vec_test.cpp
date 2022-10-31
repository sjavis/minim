#include "utils/vec.h"

#include "gtest/gtest.h"
#include "gtest-mpi-listener.hpp"
#include "ArraysMatch.h"

#include "utils/mpi.h"

typedef std::vector<double> Vector;


TEST(VecTest, TestMultiply) {
  Vector v1 = {5, 1, -0.2};
  Vector v2 = {1, 2, 4};
  double s = -2;
  int i = -2;
  Vector v1ms = v1;
  Vector v1mi = v1;
  Vector v1mv2 = v1;
  v1ms *= s;
  v1mi *= i;
  v1mv2 *= v2;
  EXPECT_TRUE(ArraysNear(s*v1, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1*s, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1ms, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(i*v1, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1*i, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1mi, {-10, -2, 0.4}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1*v2, {5, 2, -0.8}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1mv2, {5, 2, -0.8}, 1e-6));
}


TEST(VecTest, TestDivide) {
  Vector v1 = {5, 1, -0.2};
  Vector v2 = {1, 2, 4};
  double s = -2;
  int i = -2;
  Vector v1ds = v1;
  Vector v1di = v1;
  Vector v1dv2 = v1;
  v1ds /= s;
  v1di /= i;
  v1dv2 /= v2;
  EXPECT_TRUE(ArraysNear(s/v1, {-0.4, -2, 10}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1/s, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1ds, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(i/v1, {-0.4, -2, 10}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1/i, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1di, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1/v2, {5, 0.5, -0.05}, 1e-6));
  EXPECT_TRUE(ArraysNear(v1dv2, {5, 0.5, -0.05}, 1e-6));
}


TEST(VecTest, TestSquareRoot) {
  Vector v = {25, 1, 0.36};
  EXPECT_TRUE(ArraysNear(vec::sqrt(v), {5, 1, 0.6}, 1e-6));
}


TEST(VecTest, TestPow) {
  Vector v1 = {4, 1, 0.04};
  Vector v2 = {-1};
  EXPECT_TRUE(ArraysNear(vec::pow(v1,1), v1, 1e-6));
  EXPECT_TRUE(ArraysNear(vec::pow(v1,2), {16, 1, 0.0016}, 1e-6));
  EXPECT_TRUE(ArraysNear(vec::pow(v1,-1), {0.25, 1, 25}, 1e-6));
  EXPECT_TRUE(ArraysNear(vec::pow(v1,0.5), {2, 1, 0.2}, 1e-6));
  EXPECT_TRUE(ArraysNear(vec::pow(v2,-1), {-1}, 1e-6));
  EXPECT_TRUE(ArraysNear(vec::pow(v2,2), {1}, 1e-6));
}


TEST(VecTest, TestAbs) {
  EXPECT_TRUE(ArraysNear(vec::abs({5, 0, -0.2}), {5, 0, 0.2}, 1e-6));
}


int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  minim::mpiInit(&argc, &argv);

  // Add an MPI listener (https://github.com/LLNL/gtest-mpi-listener)
  ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
  ::testing::TestEventListener *l = listeners.Release(listeners.default_result_printer());
  listeners.Append(new GTestMPIListener::MPIWrapperPrinter(l, MPI_COMM_WORLD));

  return RUN_ALL_TESTS();
}
