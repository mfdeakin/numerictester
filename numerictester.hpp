
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
        timerRunning(false),
        absErrors(),
        relErrors(){};

  virtual ~NumericTest(){};

  virtual void updateStats(const TestCase &) = 0;

  virtual struct timespec totalRunTime() const {
    return runningTime;
  };

  class TimerError {};

 protected:
  void startTimer() {
    if(!timerRunning) {
      int err = clock_gettime(CLOCK_PROCESS_CPUTIME_ID,
                              &startTime);
      if(err) {
        throw new TimerError;
      }
      timerRunning = true;
    }
  }

  void stopTimer() {
    if(timerRunning) {
      int err =
          clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &endTime);
      if(err) {
        throw new TimerError;
      }
      timerRunning = false;
      struct timespec elapsed = calcDeltaTime();
      addTime(elapsed);
    }
  }

  struct timespec calcDeltaTime() {
    struct timespec delta;
    delta.tv_sec = endTime.tv_sec - startTime.tv_sec;
    delta.tv_nsec = endTime.tv_nsec - startTime.tv_nsec;
    if(delta.tv_nsec < 0) {
      delta.tv_sec--;
      constexpr const int nsPerS = 1e9;
      delta.tv_nsec += nsPerS;
    }
    return delta;
  }

  void addTime(struct timespec len) {
    runningTime.tv_nsec += len.tv_nsec;
    runningTime.tv_sec += len.tv_sec;
    constexpr const int maxNSec = 1e9;
    if(runningTime.tv_nsec > maxNSec) {
      runningTime.tv_nsec -= maxNSec;
      runningTime.tv_sec++;
    }
  }

  struct timespec runningTime;
  struct timespec startTime, endTime;
  bool timerRunning;
  std::vector<mpfr::mpreal> absErrors;
  std::vector<mpfr::mpreal> relErrors;
};
};

#endif
