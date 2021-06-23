/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitCore                                                       *
 * @(#)root/roofitcore:$Id$
 * Authors:                                                                  *
 *   PB, Patrick Bos,     NL eScience Center, p.bos@esciencecenter.nl        *
 *   IP, Inti Pelupessy,  NL eScience Center, i.pelupessy@esciencecenter.nl  *
 *   VC, Vince Croft,     DIANA / NYU,        vincent.croft@cern.ch          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

// MultiProcess back-end
//#include <MultiProcess/BidirMMapPipe.h>
#include <MultiProcess/messages.h>
#include <MultiProcess/TaskManager.h>
#include <MultiProcess/Vector.h>

#include <cstdlib>  // std::_Exit
#include <cmath>
#include <vector>
#include <map>
//#include <exception>
#include <numeric> // accumulate
#include <tuple>   // for google test Combine in parameterized test

#include <RooRealVar.h>
#include <RooRandom.h>

#include "gtest/gtest.h"
#include "test_lib.h"

class xSquaredPlusBVectorSerial {
 public:
  xSquaredPlusBVectorSerial(double b, std::vector<double> x_init) :
      _b("b", "b", b),
      x(std::move(x_init)),
      result(x.size()) {}

  virtual void evaluate() {
    // call evaluate_task for each task
    for (std::size_t ix = 0; ix < x.size(); ++ix) {
      result[ix] = std::pow(x[ix], 2) + _b.getVal();
    }
  }

  std::vector<double> get_result() {
    evaluate();
    return result;
  }

 protected:
  RooRealVar _b;
  std::vector<double> x;
  std::vector<double> result;
};


using RooFit::MultiProcess::JobTask;

class xSquaredPlusBVectorParallel : public RooFit::MultiProcess::Vector<xSquaredPlusBVectorSerial> {
 public:
  xSquaredPlusBVectorParallel(std::size_t NumCPU, double b_init, std::vector<double> x_init) :
      RooFit::MultiProcess::Vector<xSquaredPlusBVectorSerial>(NumCPU, b_init,
           x_init) // NumCPU stands for everything that defines the parallelization behaviour (number of cpu, strategy, affinity etc)
  {}

  void evaluate() override {
    if (get_manager()->is_master()) {
      // master fills queue with tasks
      for (std::size_t task_id = 0; task_id < x.size(); ++task_id) {
        JobTask job_task(id, task_id);
        get_manager()->to_queue(job_task);
      }
      waiting_for_queued_tasks = true;

      // wait for task results back from workers to master
      gather_worker_results();
      
      // put task results in desired container (same as used in serial class)
      for (std::size_t task_id = 0; task_id < x.size(); ++task_id) {
        result[task_id] = results[task_id];
      }
    }
  }


 private:
  void evaluate_task(std::size_t task) override {
    assert(get_manager()->is_worker());
    result[task] = std::pow(x[task], 2) + _b.getVal();
  }

  double get_task_result(std::size_t task) override {
    assert(get_manager()->is_worker());
    return result[task];
  }

};

class MultiProcessVectorSingleJob : public ::testing::TestWithParam<std::size_t> {
  // You can implement all the usual fixture class members here.
  // To access the test parameter, call GetParam() from class
  // TestWithParam<T>.
};




TEST_P(MultiProcessVectorSingleJob, getResult) {
  // Simple test case: calculate x^2 + b, where x is a vector. This case does
  // both a simple calculation (squaring the input vector x) and represents
  // handling of state updates in b.
  std::vector<double> x{0, 1, 2, 3};
  double b_initial = 3.;

  // start serial test

  xSquaredPlusBVectorSerial x_sq_plus_b(b_initial, x);

  auto y = x_sq_plus_b.get_result();
  std::vector<double> y_expected{3, 4, 7, 12};

  EXPECT_EQ(Hex(y[0]), Hex(y_expected[0]));
  EXPECT_EQ(Hex(y[1]), Hex(y_expected[1]));
  EXPECT_EQ(Hex(y[2]), Hex(y_expected[2]));
  EXPECT_EQ(Hex(y[3]), Hex(y_expected[3]));

  std::size_t NumCPU = GetParam();

  // start parallel test

  xSquaredPlusBVectorParallel x_sq_plus_b_parallel(NumCPU, b_initial, x);

  auto y_parallel = x_sq_plus_b_parallel.get_result();

  EXPECT_EQ(Hex(y_parallel[0]), Hex(y_expected[0]));
  EXPECT_EQ(Hex(y_parallel[1]), Hex(y_expected[1]));
  EXPECT_EQ(Hex(y_parallel[2]), Hex(y_expected[2]));
  EXPECT_EQ(Hex(y_parallel[3]), Hex(y_expected[3]));
}


INSTANTIATE_TEST_SUITE_P(NumberOfWorkerProcesses,
                        MultiProcessVectorSingleJob,
                        ::testing::Values(1,2,3));


class MultiProcessVectorMultiJob : public ::testing::TestWithParam<std::size_t> {};

TEST_P(MultiProcessVectorMultiJob, getResult) {
  // Simple test case: calculate x^2 + b, where x is a vector. This case does
  // both a simple calculation (squaring the input vector x) and represents
  // handling of state updates in b.
  std::vector<double> x{0, 1, 2, 3};
  double b_initial = 3.;

  std::vector<double> y_expected{3, 4, 7, 12};

  std::size_t NumCPU = GetParam();

  // define jobs
  xSquaredPlusBVectorParallel x_sq_plus_b_parallel(NumCPU, b_initial, x);
  xSquaredPlusBVectorParallel x_sq_plus_b_parallel2(NumCPU, b_initial + 1, x);

  // do stuff
  auto y_parallel = x_sq_plus_b_parallel.get_result();
  auto y_parallel2 = x_sq_plus_b_parallel2.get_result();

  EXPECT_EQ(Hex(y_parallel[0]), Hex(y_expected[0]));
  EXPECT_EQ(Hex(y_parallel[1]), Hex(y_expected[1]));
  EXPECT_EQ(Hex(y_parallel[2]), Hex(y_expected[2]));
  EXPECT_EQ(Hex(y_parallel[3]), Hex(y_expected[3]));

  EXPECT_EQ(Hex(y_parallel2[0]), Hex(y_expected[0] + 1));
  EXPECT_EQ(Hex(y_parallel2[1]), Hex(y_expected[1] + 1));
  EXPECT_EQ(Hex(y_parallel2[2]), Hex(y_expected[2] + 1));
  EXPECT_EQ(Hex(y_parallel2[3]), Hex(y_expected[3] + 1));
}


INSTANTIATE_TEST_SUITE_P(NumberOfWorkerProcesses,
                        MultiProcessVectorMultiJob,
                        ::testing::Values(2,1,3));
