#include "test_main.cpp"
#include "utils/range.h"

using namespace minim;


TEST(RangeTest, RangeXSize) {
  std::vector<int> xs, ys, zs;
  for (auto xyz : RangeX({2,2,3})) {
    xs.push_back(xyz[0]);
    ys.push_back(xyz[1]);
    zs.push_back(xyz[2]);
  }
  EXPECT_TRUE(ArraysMatch(xs, {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1}));
  EXPECT_TRUE(ArraysMatch(ys, {0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1}));
  EXPECT_TRUE(ArraysMatch(zs, {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2}));
}


TEST(RangeTest, RangeXSizeHalo) {
  std::vector<int> xs, ys, zs;
  for (auto xyz : RangeX({4,4,5}, 1)) {
    xs.push_back(xyz[0]);
    ys.push_back(xyz[1]);
    zs.push_back(xyz[2]);
  }
  EXPECT_TRUE(ArraysMatch(xs, { 1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2}));
  EXPECT_TRUE(ArraysMatch(ys, { 1,  1,  1,  2,  2,  2,  1,  1,  1,  2,  2,  2}));
  EXPECT_TRUE(ArraysMatch(zs, { 1,  2,  3,  1,  2,  3,  1,  2,  3,  1,  2,  3}));
}


TEST(RangeTest, RangeISize) {
  std::vector<int> is;
  for (int i : RangeI({2,2,3})) {
    is.push_back(i);
  }
  EXPECT_TRUE(ArraysMatch(is, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
}


TEST(RangeTest, RangeISizeHalo) {
  std::vector<int> is;
  for (int i : RangeI({4,4,5}, 1)) {
    is.push_back(i);
  }
  EXPECT_TRUE(ArraysMatch(is, {26, 27, 28, 31, 32, 33, 46, 47, 48, 51, 52, 53}));
}
