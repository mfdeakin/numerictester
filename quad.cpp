
#include "numerictester.hpp"
#include "mpreal.h"

#include <random>
#include <typeinfo>
#include <cmath>
#include <fstream>

#include <assert.h>
#include <time.h>

template <typename fptype>
class SphereTransCase : public NumericTester::TestCase {
 public:
  SphereTransCase(
      std::mt19937_64 &rgen,
      std::uniform_real_distribution<fptype> &dist)
      : NumericTester::TestCase() {
    radius = std::fabs(dist(rgen));
    correct = -radius;
    correct *= radius;
    for(unsigned i = 0; i < dim; i++) {
      pos[i] = dist(rgen);
      trans[i] = dist(rgen);
      mpfr::mpreal tmp(pos[i]);
      tmp += trans[i];
      tmp *= tmp;
      correct += tmp;
    }
  }

  static constexpr const unsigned dim = 3;
  fptype pos[dim];
  fptype trans[dim];
  fptype radius;
};

template <typename fptype>
class QuadNaiveTest : public NumericTester::NumericTest {
 public:
  virtual std::string testName() {
    return std::string("Naive Quadric Evaluation");
  }

  virtual void updateStats(
      const NumericTester::TestCase &testCase) {
    assert(typeid(testCase) ==
           typeid(const SphereTransCase<fptype>));
    const SphereTransCase<fptype> *stCase =
        static_cast<const SphereTransCase<fptype> *>(
            &testCase);
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
    assert(typeid(testCase) ==
           typeid(const SphereTransCase<fptype>));
    const SphereTransCase<fptype> *stCase =
        static_cast<const SphereTransCase<fptype> *>(
            &testCase);
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
    assert(typeid(testCase) ==
           typeid(const SphereTransCase<fptype>));
    const SphereTransCase<fptype> *stCase =
        static_cast<const SphereTransCase<fptype> *>(
            &testCase);
    startTimer();
    fptype moddedPt[stCase->dim + 1];
    fptype transSum = -stCase->radius * stCase->radius;
		fptype c1 = 0.0, c2 = 0.0;
    for(unsigned i = 0; i < stCase->dim; i++) {
      moddedPt[i] = stCase->pos[i] + stCase->trans[i];
      fptype mod1 = std::fma(stCase->pos[i], stCase->trans[i], -c1);
			fptype tmp = transSum + mod1;
			c1 = (tmp - transSum) - mod1;
			transSum = tmp;
			fptype mod2 = std::fma(stCase->trans[i], stCase->trans[i], -c2);
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

int main(int argc, char **argv) {
  mpfr::mpreal::set_default_prec(1024);
  typedef double fptype;
  constexpr const fptype maxMag = 1024.0 * 1024.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<float> rgenf(-maxMag,
                                              maxMag);
  NumericTester::NumericTest *tests[] = {
      new QuadNaiveTest<float>(), new QuadFMATest<float>(),
      new QuadKahanFMATest<float>()};
  for(int i = 0; i < 1e4; i++) {
    SphereTransCase<float> testcase(engine, rgenf);
    for(auto t : tests) {
      t->updateStats(testcase);
    }
  }
  for(auto t : tests) {
    std::string fname = t->testName().append(".csv");
    t->printStats();
    std::ofstream results(fname, std::ios::out);
    t->dumpData(results);
  }
  return 0;
}
