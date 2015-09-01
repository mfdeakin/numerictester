
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
              unsigned dim, unsigned precision)
      : NumericTester::TestCase(precision),
        v1(new fptype[dim]),
        v2(new fptype[dim]),
        dim(dim) {
    mpfr::mpreal val1(precision), val2(precision);
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
    for(int i = 0; i < dpCase->dim; i++)
      accumulator += dpCase->v1[i] * dpCase->v2[i];
    stopTimer();
    mpfr::mpreal absErr =
        abs(dpCase->correctValue() - accumulator);
    absErrors.push_back(absErr);
    mpfr::mpreal relErr =
        abs(absErr / dpCase->correctValue());
    relErrors.push_back(relErr);
    std::cout << "Absolute Error: " << absErr << "\n"
              << "Relative Error: " << relErr << "\n"
              << "Correct Bits: " << log(relErr) << "\n"
              << "Run Time: " << runningTime.tv_sec << "."
              << std::setfill('0') << std::setw(9)
              << runningTime.tv_nsec << " s\n";
  }
};

int main(int argc, char **argv) {
  typedef double fptype;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag,
                                              maxMag);
  DotProdCase<float> t1(engine, rgenf, 4, 512);
  DPNaiveTest<float> naive;
  naive.updateStats(t1);
  return 0;
}
