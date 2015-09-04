
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

template <typename fptype>
class NTest : public NumericTester::NumericTest {
 public:
  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    assert(typeid(testCase) ==
           typeid(const NTestCase<fptype>));
    const NTestCase<fptype> *test =
        static_cast<const NTestCase<fptype> *>(&testCase);
    startTimer();
    stopTimer();
    addStatistic(test->estimate, test->correctValue());
  }
};

TEST(Statistics, average) {
  const float epsilon = 1.0 - std::nextafter(1.0, 0.0);
  constexpr const float knownAvg[] = {2.0,  4.0,  8.0,
                                      16.0, 32.0, 64.0};
  constexpr const unsigned numTests =
      sizeof(knownAvg) / sizeof(knownAvg[0]);
  NTest<float> test;
  float accum = 0.0;
  for(unsigned i = 0; i < numTests; i++) {
    constexpr const float correctVal = 1.0;
    NTestCase<float> testcase(correctVal,
                              knownAvg[i] + correctVal);
    test.updateStats(testcase);
    accum += knownAvg[i];
  }
  const float avg = accum / numTests;
  const float ulp = epsilon * avg;
  mpfr::mpreal result = test.calcRelErrorAvg();
  EXPECT_NEAR(static_cast<double>(result), avg, 2 * ulp);
}

TEST(Statistics, variance) {
  const float epsilon = 1.0 - std::nextafter(1.0, 0.0);
  constexpr const float known[] = {2.0,  4.0,  8.0,
                                   16.0, 32.0, 64.0};
  constexpr const unsigned numTests =
      sizeof(known) / sizeof(known[0]);
  NTest<float> test;
  double sum = 0.0, sumSq = 0.0;
  for(unsigned i = 0; i < numTests; i++) {
    constexpr const float correctVal = 1.0;
    NTestCase<float> testcase(correctVal,
                              known[i] + correctVal);
    test.updateStats(testcase);
    sum += known[i];
    sumSq += known[i] * known[i];
  }
  const double avg = sum / numTests;
  const double variance =
      (-avg * sum + sumSq) / (numTests - 1);
  const float ulp = variance * epsilon;
  mpfr::mpreal result = test.calcRelErrorVar();
  EXPECT_NEAR(static_cast<double>(result), variance, ulp);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
