
#ifndef _NUMERICTESTER_HPP_
#define _NUMERICTESTER_HPP_

#include <vector>

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
        avgRelErr(0){};

  virtual ~NumericTest(){};

  virtual void updateStats(const TestCase &) = 0;

  virtual struct timespec totalRunTime() const;

  /* These methods either return the specified statistic,
   * or they throw a NoElementsError
   */
  mpfr::mpreal calcRelErrorAvg();
  mpfr::mpreal calcRelErrorMed();
  mpfr::mpreal calcRelErrorVar();
  mpfr::mpreal calcRelErrorSkew();

  /* Calculates the central moment of the data */
  template <unsigned moment>
  mpfr::mpreal calcRelErrorMoment() {
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
      result = 0;
    else if(relErrors.size() > 1)
      result = accumulator / (relErrors.size() - 1);
    else
      throw NoElementsError();
    return result;
  }

  class TimerError {};
  class NoElementsError {};

 protected:
  /* startTimer and stopTimer are timing critical;
   * don't waste time on function calls
   */
  void startTimer() {
    int err =
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &startTime);
    if(err) throw TimerError();
  }
  __attribute__((always_inline));

  void stopTimer() {
    int err =
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endTime);
    if(err) throw TimerError();
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
  mpfr::mpreal avgRelErr;
};
};

#endif
