
#include "numerictester.hpp"

#include <iomanip>
#include <assert.h>

namespace NumericTester {

struct timespec NumericTest::totalRunTime() const {
  return runningTime;
};

mpfr::mpreal NumericTest::calcRelErrorAvg() {
  if(relErrors.size() == 0) throw NoElementsError();
  return accumRelErr / relErrors.size();
}

mpfr::mpreal NumericTest::calcRelErrorVar() {
  if(relErrors.size() <= 1) throw NoElementsError();
  return calcRelErrorMoment<2>();
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

mpfr::mpreal NumericTest::calcRelErrorMax() {
  return maxRelErr;
}

mpfr::mpreal NumericTest::calcRelErrorMin() {
  return minRelErr;
}

std::array<mpfr::mpreal, 2>
NumericTest::calcRelErrorPercentile(double frac) {
  if(relErrors.size() <= 0) throw NoElementsError();
  if(frac < 0.0 || frac > 1.0) throw BadPercentileError();
  std::sort(relErrors.begin(), relErrors.end());
  int botPos =
      (int)std::floor((1.0 - frac) * relErrors.size());
  assert(botPos >= 0);
  assert(botPos < relErrors.size());
  int topPos = (int)std::ceil(frac * relErrors.size());
  assert(topPos >= 0);
  assert(topPos < relErrors.size());
  std::array<mpfr::mpreal, 2> ret;
  ret[0] = relErrors[botPos];
  ret[1] = relErrors[topPos];
  return ret;
}

mpfr::mpreal NumericTest::calcRelErrorSkew() {
  mpfr::mpreal stddev = sqrt(calcRelErrorVar());
  mpfr::mpreal denominator = stddev * stddev * stddev;
  mpfr::mpreal moment3 = calcRelErrorMoment<3>();
  mpfr::mpreal numerator = moment3 / relErrors.size();
  return numerator / denominator;
}

mpfr::mpreal NumericTest::calcRelErrorKurtosis() {
  mpfr::mpreal var = calcRelErrorVar();
  mpfr::mpreal denominator = var * var;
  mpfr::mpreal moment4 = calcRelErrorMoment<4>();
  mpfr::mpreal kurtosis = moment4 / denominator - 3.0;
  return kurtosis;
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
  mpfr::mpreal absErr = abs(estimate - correct);
  absErrors.push_back(absErr);
  mpfr::mpreal relErr = abs(absErr / correct);
  relErrors.push_back(relErr);
  accumRelErr += relErr;
  if(isnan(maxRelErr) || relErr > maxRelErr)
    maxRelErr = relErr;
  if(isnan(minRelErr) || relErr < minRelErr)
    minRelErr = relErr;
}

void NumericTest::dumpData(std::ostream &out) {
  printStats(out);
  out << "Absolute Error, Relative Error\n";
  for(unsigned i = 0; i < relErrors.size(); i++) {
    mpfr::mpreal absErr = absErrors[i];
    out << absErr << ", ";
    mpfr::mpreal relErr = relErrors[i];
    out << relErr << "\n";
  }
}

void NumericTest::printStats(std::ostream &out) {
  constexpr const int nsDigits = 9;
  std::array<mpfr::mpreal, 2> percent =
      calcRelErrorPercentile(0.99);
  out << testName() << "\n"
      << "Running Time: " << runningTime.tv_sec << "."
      << std::setw(nsDigits) << std::setfill('0')
      << runningTime.tv_nsec << "\n"
      << "Relative Error Average: " << calcRelErrorAvg()
      << "\n"
      << "Relative Error Variance: " << calcRelErrorVar()
      << "\n"
      << "Relative Error Skew: " << calcRelErrorSkew()
      << "\n"
      << "Relative Error Kurtosis: "
      << calcRelErrorKurtosis() << "\n"
      << "Relative Error 99th Percentile: " << percent[0]
      << ", " << percent[1] << "\n";
}
};
