
#include "numerictester.hpp"
#include "mpreal.h"

#include <random>
#include <typeinfo>
#include <cmath>
#include <fstream>

#include <assert.h>
#include <time.h>

template <typename fptype>
class QuadricTestCase : public NumericTester::TestCase {
 public:
  QuadricTestCase() : NumericTester::TestCase() {
    for(unsigned i = 0; i < dim; i++) {
      pos[i] = 0.0;
      trans[i] = 0.0;
    }
    radius = 0.0;
  }
  static constexpr const unsigned dim = 3;
  fptype pos[dim];
  fptype trans[dim];
  fptype radius;
};

template <typename fptype>
class SphereTransCase : public QuadricTestCase<fptype> {
 public:
  SphereTransCase(
      std::mt19937_64 &rgen,
      std::uniform_real_distribution<fptype> &dist)
      : QuadricTestCase<fptype>() {
    this->radius = std::fabs(dist(rgen));
    this->correct = -this->radius;
    this->correct *= this->radius;
    for(unsigned i = 0; i < this->dim; i++) {
      this->pos[i] = dist(rgen);
      this->trans[i] = dist(rgen);
      mpfr::mpreal tmp(this->pos[i]);
      tmp += this->trans[i];
      tmp *= tmp;
      this->correct += tmp;
    }
  }
};

template <typename fptype>
class AxisCylinderTransCase
    : public QuadricTestCase<fptype> {
 public:
  AxisCylinderTransCase(
      std::mt19937_64 &rgen,
      std::uniform_real_distribution<fptype> &dist)
      : QuadricTestCase<fptype>() {
    unsigned axis =
        ((unsigned)std::floor(dist(rgen))) % this->dim;
    this->radius = std::fabs(dist(rgen));
    this->correct = -this->radius;
    this->correct *= this->radius;
    for(unsigned i = 0; i < this->dim; i++) {
      if(i == axis)
        this->trans[i] = 0.0;
      else
        this->trans[i] = dist(rgen);
      this->pos[i] = dist(rgen);
      mpfr::mpreal tmp(this->pos[i]);
      tmp += this->trans[i];
      tmp *= tmp;
      this->correct += tmp;
    }
  }
};

template <typename fptype>
class QuadNullTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Null Quadric Evaluation");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    const QuadricTestCase<fptype> *stCase =
        dynamic_cast<const QuadricTestCase<fptype> *>(
            &testCase);
    assert(stCase != NULL);
    startTimer();
    fptype accumulator = NAN;
    stopTimer();
    mpfr::mpreal estimate(accumulator);
    addStatistic(estimate, testCase.correctValue());
  }
};

template <typename fptype>
class QuadNaiveTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Naive Quadric Evaluation");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    const QuadricTestCase<fptype> *stCase =
        dynamic_cast<const QuadricTestCase<fptype> *>(
            &testCase);
    assert(stCase != NULL);
    startTimer();
    fptype moddedPt[stCase->dim + 1];
    fptype transSum = -stCase->radius * stCase->radius;
    for(unsigned i = 0; i < stCase->dim; i++) {
      moddedPt[i] = stCase->pos[i] + stCase->trans[i];
      transSum += stCase->pos[i] * stCase->trans[i] +
                  stCase->trans[i] * stCase->trans[i];
    }
    fptype accumulator = transSum;
    for(unsigned i = 0; i < stCase->dim; i++) {
      accumulator += stCase->pos[i] * moddedPt[i];
    }
    stopTimer();
    mpfr::mpreal estimate(accumulator);
    addStatistic(estimate, testCase.correctValue());
  }
};

template <typename fptype>
class QuadFMATest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("FMA Quadric Evaluation");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    const QuadricTestCase<fptype> *stCase =
        dynamic_cast<const QuadricTestCase<fptype> *>(
            &testCase);
    assert(stCase != NULL);
    startTimer();
    fptype moddedPt[stCase->dim + 1];
    fptype transSum = -stCase->radius * stCase->radius;
    for(unsigned i = 0; i < stCase->dim; i++) {
      moddedPt[i] = stCase->pos[i] + stCase->trans[i];
      transSum = std::fma(stCase->pos[i], stCase->trans[i],
                          transSum);
      transSum = std::fma(stCase->trans[i],
                          stCase->trans[i], transSum);
    }
    fptype accumulator = transSum;
    for(unsigned i = 0; i < stCase->dim; i++) {
      accumulator = std::fma(stCase->pos[i], moddedPt[i],
                             accumulator);
    }
    stopTimer();
    mpfr::mpreal estimate(accumulator);
    addStatistic(estimate, testCase.correctValue());
  }
};

template <typename fptype>
class QuadKahanFMATest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Kahan FMA Quadric Evaluation");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    const QuadricTestCase<fptype> *stCase =
        dynamic_cast<const QuadricTestCase<fptype> *>(
            &testCase);
    assert(stCase != NULL);
    startTimer();
    fptype moddedPt[stCase->dim + 1];
    fptype transSum = -stCase->radius * stCase->radius;
    fptype c1 = 0.0, c2 = 0.0;
    for(unsigned i = 0; i < stCase->dim; i++) {
      moddedPt[i] = stCase->pos[i] + stCase->trans[i];
      fptype mod1 =
          std::fma(stCase->pos[i], stCase->trans[i], -c1);
      fptype tmp = transSum + mod1;
      c1 = (tmp - transSum) - mod1;
      transSum = tmp;
      fptype mod2 =
          std::fma(stCase->trans[i], stCase->trans[i], -c2);
      tmp = transSum + mod2;
      c2 = (tmp - transSum) - mod2;
      transSum = tmp;
    }
    fptype accumulator = transSum;
    for(unsigned i = 0; i < stCase->dim; i++) {
      accumulator = std::fma(stCase->pos[i], moddedPt[i],
                             accumulator);
    }
    stopTimer();
    mpfr::mpreal estimate(accumulator);
    addStatistic(estimate, testCase.correctValue());
  }
};

template <typename testtype, typename fptype>
void runQuadricTests(
    std::mt19937_64 engine,
    std::uniform_real_distribution<fptype> rgenf,
    const int n, const std::string testclass) {
  NumericTester::NumericTest *tests[] = {
      new QuadNullTest<fptype>(),
      new QuadNaiveTest<fptype>(),
      new QuadFMATest<fptype>(),
      new QuadKahanFMATest<fptype>(),
      new QuadNullTest<fptype>(),
      new QuadNaiveTest<fptype>(),
      new QuadFMATest<fptype>(),
      new QuadKahanFMATest<fptype>()};
  constexpr const int numTests =
      sizeof(tests) / sizeof(tests[0]);
  for(int i = 0; i < n; i++) {
    testtype testcase(engine, rgenf);
    for(auto t : tests) {
      t->updateStats(testcase);
    }
  }
  std::cout << testclass << "\n";
  for(auto t : tests) t->printStats();
  for(int i = numTests / 2; i < numTests; i++) {
    auto t = tests[i];
    std::string fname =
        t->testName().append(testclass).append(".csv");
    std::ofstream results(fname, std::ios::out);
    t->dumpData(results);
  }
  for(auto t : tests) delete t;
}

int main(int argc, char **argv) {
  mpfr::mpreal::set_default_prec(128);
  using fptype = float;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<fptype> rgenf(-maxMag,
                                               maxMag);
  runQuadricTests<SphereTransCase<fptype>, fptype>(
      engine, rgenf, 1e4, std::string("Sphere Tests"));
  runQuadricTests<AxisCylinderTransCase<fptype>, fptype>(
      engine, rgenf, 1e4,
      std::string("Axis Aligned Cylinder Tests"));
  return 0;
}
