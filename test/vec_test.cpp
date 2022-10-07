#include "utils/vec.h"

#include "gtest/gtest.h"
#include "ArraysMatch.h"

typedef std::vector<double> Vector;


TEST(VecTest, TestDivide) {
  Vector v = {5, 1, -0.2};
  double s = -2;
  int i = -2;
  Vector vds = v;
  Vector vdi = v;
  vds /= s;
  vdi /= i;
  EXPECT_TRUE(ArraysNear(s/v, {-0.4, -2, 10}, 1e-6));
  EXPECT_TRUE(ArraysNear(v/s, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(vds, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(i/v, {-0.4, -2, 10}, 1e-6));
  EXPECT_TRUE(ArraysNear(v/i, {-2.5, -0.5, 0.1}, 1e-6));
  EXPECT_TRUE(ArraysNear(vdi, {-2.5, -0.5, 0.1}, 1e-6));
}


TEST(VecTest, TestAbs) {
  EXPECT_TRUE(ArraysNear(vec::abs({5, 0, -0.2}), {5, 0, 0.2}, 1e-6));
}
