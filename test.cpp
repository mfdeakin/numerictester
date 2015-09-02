
#include <gtest/gtest.h>

#include "numerictester.hpp"

template <typename fptype>
class NTest;

template <typename fptype>
class NTestCase : public NumericTester::TestCase {
 public:
  NTestCase(fptype correctVal, fptype estimateVal)
      : NumericTester::TestCase(), estimate(estimateVal) {
    correct = correctVal;
  }
  mpfr::mpreal estimate;
};

TEST(Statistics, average) {
  float knownAvg[] = {2, 4, 8, 16};
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
