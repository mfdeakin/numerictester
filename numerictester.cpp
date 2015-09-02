
#include "numerictester.hpp"

namespace NumericTester {

struct timespec NumericTest::totalRunTime() const {
  return runningTime;
};

mpfr::mpreal NumericTest::calcRelErrorAvg() {
  if(relErrors.size() == 0) throw NoElementsError();
  return avgRelErr / relErrors.size();
}

mpfr::mpreal NumericTest::calcRelErrorVar() {
  if(relErrors.size() <= 1) throw NoElementsError();
  return calcRelErrorMoment<2>() / (relErrors.size() - 1);
}

mpfr::mpreal NumericTest::calcRelErrorMed() {
  if(relErrors.size() <= 0) throw NoElementsError();
  std::sort(relErrors.begin(), relErrors.end());
  int medianPos = relErrors.size() / 2,
      isOdd = relErrors.size() % 2;
  mpfr::mpreal median = relErrors[medianPos];
  if(!isOdd) {
    median += relErrors[medianPos - 1];
    median /= 2;
  }
  return median;
}

mpfr::mpreal NumericTest::calcRelErrorSkew() {
  mpfr::mpreal stddev = sqrt(calcRelErrorVar());
  mpfr::mpreal denominator = stddev * stddev * stddev;
  mpfr::mpreal moment3 = calcRelErrorMoment<3>();
  mpfr::mpreal numerator = moment3 / relErrors.size();
  return numerator / denominator;
}

struct timespec NumericTest::calcDeltaTime() {
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

void NumericTest::addTime(struct timespec len) {
  runningTime.tv_nsec += len.tv_nsec;
  runningTime.tv_sec += len.tv_sec;
  constexpr const int maxNSec = 1e9;
  if(runningTime.tv_nsec > maxNSec) {
    runningTime.tv_nsec -= maxNSec;
    runningTime.tv_sec++;
  }
}

void NumericTest::addStatistic(mpfr::mpreal estimate,
                               mpfr::mpreal correct) {
  mpfr::mpreal absErr(512);
  absErr = abs(estimate - correct);
  absErrors.push_back(absErr);
  mpfr::mpreal relErr(512);
  relErr = abs(absErr / correct);
  relErrors.push_back(relErr);
  avgRelErr += relErr;
}
};
