
#ifndef _NUMERICTESTER_HPP_
#define _NUMERICTESTER_HPP_

#include <vector>

#include <time.h>

#include "mpreal.h"

namespace NumericTester {

class TestCase {
 public:
  TestCase(unsigned precision) : correct(precision) {}
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
        relErrors(){};

  virtual ~NumericTest(){};

  virtual void updateStats(const TestCase &) = 0;

  virtual struct timespec totalRunTime() const;

  mpfr::mpreal calcRelErrorAvg() { return avgRelErr; }
  mpfr::mpreal calcRelErrorMed();
  mpfr::mpreal calcRelErrorVar(unsigned precision);
  mpfr::mpreal calcRelErrorSkew(unsigned precision);

  template <unsigned moment>
  mpfr::mpreal calcRelErrorMoment(unsigned precision) {
    mpfr::mpreal accumulator(precision);
    accumulator = 0;
    for(auto &err : relErrors) {
      mpfr::mpreal delta(precision);
      delta = err - avgRelErr;
      if(moment == 0) {
        if(delta > 0) accumulator += 1;
      } else {
        for(unsigned i = 1; i < moment; i++) delta *= delta;
      }
      accumulator += delta;
    }
    return accumulator;
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
