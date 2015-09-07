
#include "numerictester.hpp"
#include "genericfp.hpp"
#include "mpreal.h"

#include <random>
#include <typeinfo>
#include <cmath>
#include <fstream>

#include <assert.h>
#include <time.h>

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
		correctRounded = correct.toLDouble();
  }

  virtual ~DotProdCase() {
    delete[] v1;
    delete[] v2;
  }

  template <typename>
  friend class DPNaiveTest;
  template <typename>
  friend class DPFMATest;
  template <typename>
  friend class DPFMAKahanTest;

 private:
  fptype *v1;
  fptype *v2;
	fptype correctRounded;
  const unsigned dim;
};

template <typename fptype>
class DPNullTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Null Dot Product");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    fptype result = 0.0;
    if(typeid(testCase) ==
       typeid(const DotProdCase<float>)) {
      const DotProdCase<float> *dpCase =
          static_cast<const DotProdCase<float> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<double>)) {
      const DotProdCase<double> *dpCase =
          static_cast<const DotProdCase<double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<long double>)) {
      const DotProdCase<long double> *dpCase =
          static_cast<const DotProdCase<long double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    }
    mpfr::mpreal estimate(result);
    addStatistic(estimate, testCase.correctValue());
  }

  template <typename intype>
  fptype runTest(const DotProdCase<intype> *dpCase) {
    return 0.0;
  }
};

template <typename fptype>
class DPNaiveTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Naive Dot Product with ") +
           fpconvert<fptype>::fpname;
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    fptype result = 0.0;
    if(typeid(testCase) ==
       typeid(const DotProdCase<float>)) {
      const DotProdCase<float> *dpCase =
          static_cast<const DotProdCase<float> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<double>)) {
      const DotProdCase<double> *dpCase =
          static_cast<const DotProdCase<double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<long double>)) {
      const DotProdCase<long double> *dpCase =
          static_cast<const DotProdCase<long double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    }
    mpfr::mpreal estimate(result);
    addStatistic(estimate, testCase.correctValue());
  }

  template <typename intype>
  fptype runTest(const DotProdCase<intype> *dpCase) {
    fptype accumulator = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++)
      accumulator += dpCase->v1[i] * dpCase->v2[i];
    return accumulator;
  }
};

template <typename fptype>
class DPFMATest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("FMA Dot Product with ") +
           fpconvert<fptype>::fpname;
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    fptype result = 0.0;
    if(typeid(testCase) ==
       typeid(const DotProdCase<float>)) {
      const DotProdCase<float> *dpCase =
          static_cast<const DotProdCase<float> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<double>)) {
      const DotProdCase<double> *dpCase =
          static_cast<const DotProdCase<double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<long double>)) {
      const DotProdCase<long double> *dpCase =
          static_cast<const DotProdCase<long double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    }
    mpfr::mpreal estimate(result);
    addStatistic(estimate, testCase.correctValue());
  }

  template <typename intype>
  fptype runTest(const DotProdCase<intype> *dpCase) {
    fptype accumulator = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++)
      accumulator = std::fma(dpCase->v1[i], dpCase->v2[i],
                             accumulator);
    return accumulator;
  }
};

template <typename fptype>
class DPFMAKahanTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Kahan FMA Dot Product with ") +
           fpconvert<fptype>::fpname;
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    fptype result = 0.0;
    if(typeid(testCase) ==
       typeid(const DotProdCase<float>)) {
      const DotProdCase<float> *dpCase =
          static_cast<const DotProdCase<float> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<double>)) {
      const DotProdCase<double> *dpCase =
          static_cast<const DotProdCase<double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<long double>)) {
      const DotProdCase<long double> *dpCase =
          static_cast<const DotProdCase<long double> *>(
              &testCase);
      startTimer();
      result = runTest(dpCase);
      stopTimer();
    }
    mpfr::mpreal estimate(result);
    addStatistic(estimate, testCase.correctValue());
  }

  template <typename intype>
  fptype runTest(const DotProdCase<intype> *dpCase) {
    fptype accumulator = 0.0;
    fptype c = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++) {
      fptype mod =
          std::fma(dpCase->v1[i], dpCase->v2[i], -c);
      fptype tmp = accumulator + mod;
      c = (tmp - accumulator) - mod;
      accumulator = tmp;
    }
    return accumulator;
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
  DPNullTest<float> null;
  NumericTester::NumericTest *tests[] = {
      new DPNaiveTest<float>(),
      new DPFMATest<float>(),
      new DPFMAKahanTest<float>(),
      new DPNaiveTest<double>(),
      new DPFMATest<double>(),
      new DPFMAKahanTest<double>(),
      new DPNaiveTest<long double>(),
      new DPFMATest<long double>(),
      new DPFMAKahanTest<long double>()};
  for(int i = 0; i < 1e3; i++) {
    DotProdCase<float> testcase(engine, rgenf, 4);
    null.updateStats(testcase);
    for(auto t : tests) t->updateStats(testcase);
  }
  null.printStats();
  for(auto t : tests) {
    t->printStats();
    std::string fname = t->testName().append(".csv");
    std::ofstream results(fname, std::ios::out);
    t->dumpData(results);
    delete t;
  }
  return 0;
}
