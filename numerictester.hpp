
#ifndef _NUMERICTESTER_HPP_
#define _NUMERICTESTER_HPP_

#include <vector>
#include <string>
#include <array>

#include <iostream>
#include <time.h>

#include "mpreal.h"

namespace NumericTester {

class TestCase {
 public:
  TestCase() : correct() {}
  virtual ~TestCase() {}
  virtual const mpfr::mpreal correctValue() const {
    return correct;
  }

 protected:
  mpfr::mpreal correct;
};

class NumericTest {
 public:
  NumericTest()
      : runningTime({0, 0}),
        startTime({0, 0}),
        endTime({0, 0}),
        absErrors(),
        relErrors(),
        accumRelErr(0),
        maxRelErr(NAN),
        minRelErr(NAN){};

  virtual ~NumericTest(){};

  virtual void updateStats(const TestCase &) = 0;

  virtual struct timespec totalRunTime() const;

  virtual std::string testName() = 0;
  virtual void printStats(std::ostream &out = std::cout);
  virtual void dumpData(std::ostream &out = std::cout);

  /* These methods either return the specified statistic,
   * or they throw a NoElementsError
   */
  mpfr::mpreal calcRelErrorAvg();
  mpfr::mpreal calcRelErrorMed();
  mpfr::mpreal calcRelErrorVar();
  /* Calculates the upper and lower frac percentile,
   * with 0 < frac < 1.
   * ie. frac=0.99 calculates the 99th percentile,
   * and returns {{lowerPercentile, upperPercentile}}
   */
  std::array<mpfr::mpreal, 2> calcRelErrorPercentile(
      double frac);
  /* Calculates the univariate data skew */
  mpfr::mpreal calcRelErrorSkew();
  mpfr::mpreal calcRelErrorKurtosis();
  mpfr::mpreal calcRelErrorMax();
  mpfr::mpreal calcRelErrorMin();

  /* Calculates the central moment of the data */
  template <unsigned moment>
  mpfr::mpreal calcRelErrorMoment() {
    if(moment == 1)
      return mpfr::mpreal(0);
    else if(moment > 1 && relErrors.size() <= 1)
      throw NoElementsError();
    mpfr::mpreal accumulator;
    accumulator = 0;
    for(auto &err : relErrors) {
      mpfr::mpreal delta = err - calcRelErrorAvg();
      mpfr::mpreal power = delta;
      if(moment == 0) {
        if(err > 0) accumulator += 1;
      } else {
        for(unsigned i = 1; i < moment; i++) power *= delta;
        accumulator += power;
      }
    }
    mpfr::mpreal result;
    if(moment == 0)
      result = accumulator;
    else
      result = accumulator / (relErrors.size() - 1);
    return result;
  }

  class TimerError {};
  class NoElementsError {};
  class BadPercentileError {};

 protected:
  /* startTimer and stopTimer are timing critical;
   * don't waste time on function calls
   */
  void startTimer() {
    int err =
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startTime);
  }
  __attribute__((always_inline));

  void stopTimer() {
    int err =
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endTime);
    struct timespec elapsed = calcDeltaTime();
    addTime(elapsed);
  }
  __attribute__((always_inline));

  void addStatistic(mpfr::mpreal estimate,
                    mpfr::mpreal correct);

  struct timespec calcDeltaTime();
  void addTime(struct timespec len);

  struct timespec runningTime;
  struct timespec startTime, endTime;
  std::vector<mpfr::mpreal> absErrors;
  std::vector<mpfr::mpreal> relErrors;
  mpfr::mpreal accumRelErr, maxRelErr, minRelErr;
};
};

#endif
