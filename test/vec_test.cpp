#include "utils/vec.h"

#include "gtest/gtest.h"
#include "ArraysMatch.h"

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


TEST(VecTest, TestAbs) {
  EXPECT_TRUE(ArraysNear(vec::abs({5, 0, -0.2}), {5, 0, 0.2}, 1e-6));
}
