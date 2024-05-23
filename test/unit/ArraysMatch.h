#include "gtest/gtest.h"
#include <vector>


template<typename T>
::testing::AssertionResult ArraysMatch(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size()) {
    return ::testing::AssertionFailure() << "array 1 has size " << a.size() << " != expected size " << b.size();
  }

  if (a != b) {
    return ::testing::AssertionFailure() << ::testing::PrintToString(a) << " != " << ::testing::PrintToString(b);
  }

  return ::testing::AssertionSuccess();
}


template<typename T>
::testing::AssertionResult ArraysNear(const std::vector<T>& a, const std::vector<T>& b, float delta=1e-6) {
  if (a.size() != b.size()) {
    return ::testing::AssertionFailure() << "array 1 has size " << a.size() << " != expected size " << b.size();
  }

  for (size_t i(0); i < a.size(); ++i) {
    if (std::abs(a[i] - b[i]) > delta) {
      return ::testing::AssertionFailure() << ::testing::PrintToString(a) << " != " << ::testing::PrintToString(b);
    }
  }

  return ::testing::AssertionSuccess();
}
