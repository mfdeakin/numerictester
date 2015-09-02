
#include "numerictester.hpp"
#include "mpreal.h"

#include <random>
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <assert.h>
#include <time.h>

template <typename fptype>
class DotProdTest;

template <typename fptype>
class DotProdCase : public NumericTester::TestCase {
 public:
  DotProdCase(std::mt19937_64 &rgen,
              std::uniform_real_distribution<fptype> &dist,
              unsigned dim)
      : NumericTester::TestCase(),
        v1(new fptype[dim]),
        v2(new fptype[dim]),
        dim(dim) {
    mpfr::mpreal val1, val2;
    for(unsigned i = 0; i < dim; i++) {
      v1[i] = dist(rgen);
      val1 = v1[i];
      v2[i] = dist(rgen);
      val2 = v2[i];
      correct = correct + val1 * val2;
    }
  }

  virtual ~DotProdCase() {
    delete[] v1;
    delete[] v2;
  }

  template <typename>
  friend class DPNaiveTest;

 private:
  fptype *v1;
  fptype *v2;
  const unsigned dim;
};

template <typename fptype>
class DPNaiveTest : public NumericTester::NumericTest {
 public:
  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    assert(typeid(testCase) ==
           typeid(const DotProdCase<fptype>));
    const DotProdCase<fptype> *dpCase =
        static_cast<const DotProdCase<fptype> *>(&testCase);
    startTimer();
    fptype accumulator = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++)
      accumulator += dpCase->v1[i] * dpCase->v2[i];
    stopTimer();
    mpfr::mpreal estimate(accumulator);
    addStatistic(estimate, testCase.correctValue());
  }
};

int main(int argc, char **argv) {
  mpfr::mpreal::set_default_prec(1024);
  typedef double fptype;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag,
                                              maxMag);
  DPNaiveTest<float> naive;
  for(int i = 0; i < 6e6; i++) {
    DotProdCase<float> test(engine, rgenf, 4);
    naive.updateStats(test);
  }
  struct timespec t = naive.totalRunTime();
  constexpr const int nsDigits = 9;
  std::cout << "Running Time: " << t.tv_sec << "."
            << std::setw(nsDigits) << std::setfill('0')
            << t.tv_nsec << "\n";
  std::cout << "Relative Error Mean: "
            << naive.calcRelErrorAvg() << "\n";
  std::cout << "Relative Error Variance: "
            << naive.calcRelErrorVar() << "\n";
  std::cout << "Relative Error Skew: "
            << naive.calcRelErrorSkew() << "\n";
  return 0;
}
