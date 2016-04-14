
#include "numerictester.hpp"
#include "genericfp.hpp"
#include "kobbelt.hpp"
#include "mpreal.h"
#include "accurate_math.hpp"

#include <random>
#include <typeinfo>
#include <cmath>
#include <fstream>

#include <assert.h>
#include <time.h>

template <typename fptype>
class DotProdCase : public NumericTester::TestCase {
 public:
  DotProdCase(std::mt19937_64 &rgen, auto &signDist,
              auto &expDist, auto &manDist, unsigned dim)
      : NumericTester::TestCase(),
        v1(new fptype[dim]),
        v2(new fptype[dim]),
        dim(dim) {
    mpfr::mpreal val1, val2;
    for(unsigned i = 0; i < dim; i++) {
      v1[i] =
          generateFPVal(rgen, signDist, expDist, manDist);
      val1 = v1[i];
      v2[i] =
          generateFPVal(rgen, signDist, expDist, manDist);
      val2 = v2[i];
      correct = correct + val1 * val2;
    }
    correctRounded = correct.toLDouble();
  }

  virtual ~DotProdCase() {
    delete[] v1;
    delete[] v2;
  }

  static fptype generateFPVal(std::mt19937_64 &rgen,
                              auto &signDist, auto &expDist,
                              auto &manDist) {
    for(;;) {
      union {
        GenericFP::fpconvert<fptype> fpBits;
        fptype fpVal;
      } fpBuf;
      fpBuf.fpBits.sign = signDist(rgen);
      fpBuf.fpBits.exponent = expDist(rgen);
      fpBuf.fpBits.mantissa = manDist(rgen);
      if(!std::isnan(fpBuf.fpVal) &&
         !std::isinf(fpBuf.fpVal))
        return fpBuf.fpVal;
    }
  }

  template <typename>
  friend class DPNaiveTest;
  template <typename>
  friend class DPFMATest;
  template <typename>
  friend class DPKahanTest;
  template <typename>
  friend class DPFMAKahanTest;
  template <typename>
  friend class DPExactFMACompTest;
  template <typename>
  friend class DPKobbeltTest;

 private:
  fptype *v1;
  fptype *v2;
  fptype correctRounded;
  const unsigned dim;
};

/* Use the Curiously Recurring Template Pattern (CRTP)
 * to implement static polymorphism here
 */
template <typename fptype, typename derived>
class DPTestInterface : public NumericTester::NumericTest {
 public:
  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    fptype result = 0.0;
    if(typeid(testCase) ==
       typeid(const DotProdCase<float>)) {
      const DotProdCase<float> *dpCase =
          static_cast<const DotProdCase<float> *>(
              &testCase);
      startTimer();
      result =
          static_cast<derived *>(this)->runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<double>)) {
      const DotProdCase<double> *dpCase =
          static_cast<const DotProdCase<double> *>(
              &testCase);
      startTimer();
      result =
          static_cast<derived *>(this)->runTest(dpCase);
      stopTimer();
    } else if(typeid(testCase) ==
              typeid(const DotProdCase<long double>)) {
      const DotProdCase<long double> *dpCase =
          static_cast<const DotProdCase<long double> *>(
              &testCase);
      startTimer();
      result =
          static_cast<derived *>(this)->runTest(dpCase);
      stopTimer();
    }
    mpfr::mpreal estimate(result);
    addStatistic(estimate, testCase.correctValue());
  }
};

template <typename fptype>
class DPNaiveTest
    : public DPTestInterface<fptype, DPNaiveTest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string("Naive Dot Product with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
    fptype accumulator = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++)
      accumulator += dpCase->v1[i] * dpCase->v2[i];
    return accumulator;
  }
};

template <typename fptype>
class DPFMATest
    : public DPTestInterface<fptype, DPFMATest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string("FMA Dot Product with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
    intype accumulator = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++) {
      accumulator = std::fma(dpCase->v1[i], dpCase->v2[i],
                             accumulator);
    }
    return accumulator;
  }
};

template <typename fptype>
class DPKahanTest
    : public DPTestInterface<fptype, DPKahanTest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string("Kahan Dot Product with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
    fptype accumulator = 0.0;
    fptype c = 0.0;
    for(unsigned i = 0; i < dpCase->dim; i++) {
      fptype mod = dpCase->v1[i] * dpCase->v2[i] - c;
      fptype tmp = accumulator + mod;
      c = (tmp - accumulator) - mod;
      accumulator = tmp;
    }
    return accumulator;
  }
};

template <typename fptype>
class DPFMAKahanTest
    : public DPTestInterface<fptype,
                             DPFMAKahanTest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string("Kahan FMA Dot Product with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
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

template <typename fptype>
class DPExactFMACompTest
    : public DPTestInterface<fptype,
                             DPExactFMACompTest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string(
               "Exact FMA Compensated Dot Product "
               "with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
    std::array<intype, 2> terms(
        twoProd(dpCase->v1[0], dpCase->v2[0]));
    for(unsigned i = 1; i < dpCase->dim; i++) {
      std::array<intype, 3> accumulated(
          threeFMA(dpCase->v1[i], dpCase->v2[i], terms[0]));
      terms[0] = accumulated[0];
      terms[1] += (accumulated[1] + accumulated[2]);
    }
    return fptype(terms[0] + terms[1]);
  }
};

template <typename fptype>
class DPKobbeltTest
    : public DPTestInterface<fptype,
                             DPKobbeltTest<fptype>> {
 public:
  virtual std::string testName() {
    return std::string("Kobbelt Dot Product with ") +
           GenericFP::fpconvert<fptype>::fpname;
  }

  template <typename intype>
  fptype __attribute__((noinline))
  runTest(const DotProdCase<intype> *dpCase) {
    return kobbeltDotProd<intype, fptype>(
        dpCase->v1, dpCase->v2, dpCase->dim);
  }
};

void runTests(const int numTests, const int vecSize);

int main(int argc, char **argv) {
  int numTests = 1e5;
  int vecSize = 4;
  if(argc > 1) {
    numTests = atoi(argv[1]);
    if(numTests < 1) {
      printf("Number of tests must be greater than 0\n");
      return -1;
    }
    if(argc > 2) {
      vecSize = atoi(argv[2]);
      if(vecSize < 1) {
        printf("Vector size must be greater than 0\n");
        return -1;
      }
    }
  }
  runTests(numTests, vecSize);
  return 0;
}

void runTests(const int numTests, const int vecSize) {
  mpfr::mpreal::set_default_prec(1024);
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_int_distribution<int> rgenExp(
      0, GenericFP::fpconvert<float>::centralExp + 20);
  std::uniform_int_distribution<int> rgenMan(
      0, GenericFP::fpconvert<float>::maxMantissa);
  std::uniform_int_distribution<int> rgenSign(0, 1);
  NumericTester::NumericTest *tests[] = {
      new DPNaiveTest<float>(),
      new DPFMATest<float>(),
      new DPKahanTest<float>(),
      new DPFMAKahanTest<float>(),
      new DPExactFMACompTest<float>(),
      new DPKobbeltTest<float>(),

      new DPNaiveTest<double>(),
      new DPFMATest<double>(),
      new DPKahanTest<double>(),
      new DPFMAKahanTest<double>(),
      new DPExactFMACompTest<double>(),
      new DPKobbeltTest<double>(),

      new DPNaiveTest<long double>(),
      new DPFMATest<long double>(),
      new DPKahanTest<long double>(),
      new DPFMAKahanTest<long double>(),
      new DPExactFMACompTest<long double>(),
      new DPKobbeltTest<long double>()};
  for(int i = 0; i < numTests; i++) {
    DotProdCase<float> testcase(engine, rgenSign, rgenExp,
                                rgenMan, vecSize);
    for(auto t : tests) t->updateStats(testcase);
  }
  for(auto t : tests) {
    t->printStats();
    std::cout << "\n\n";
    std::string fname = t->testName().append(".csv");
    std::ofstream results(fname, std::ios::out);
    t->dumpData(results);
    delete t;
  }
}
