
#include "numerictester.hpp"
#include "mpreal.h"

#include <random>

#include <assert.h>
#include <time.h>

template <typename fptype>
class DotProdTest;

template <typename fptype, unsigned dim, unsigned precision>
class DotProdCase : public NumericTester::TestCase {
 public:
  DotProdCase(std::mt19937_64 &rgen,
              std::uniform_real_distribution<fptype> &dist)
      : NumericTester::TestCase(precision),
        v1(new fptype[dim]),
        v2(new fptype[dim]) {
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
  friend class DotProdTest;

 private:
  fptype *v1;
  fptype *v2;

  constexpr static const int _dim = dim;
};

template <typename fptype>
class DotProdTest : public NumericTester::NumericTest {
 public:
  virtual void updateStats(const DotProdCase &dpCase) {
		startTimer();
		fptype accumulator = 0.0;
		for(int i = 0; i < dpCase._dim; i++)
			accumulator += dpCase.v1[i] * dpCase.v2[i];
		stopTimer();
		mpfr::mpreal err = dpCase.correctValue() - accumulator;
	}
};

int main(int argc, char **argv) {
  typedef double fptype;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag,
                                              maxMag);
  DotProdCase<float, 4, 256> t1(engine, rgenf);
  return 0;
}
