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

#include <cmath>
#include <vector>

#include <RooFit/MultiProcess/Job.h>
#include <RooFit/MultiProcess/types.h>  // JobTask
// needed to complete type returned from...
#include <RooFit/MultiProcess/JobManager.h>      // ... Job::get_manager()
#include <RooFit/MultiProcess/ProcessManager.h>  // ... JobManager::process_manager()
#include <RooFit/MultiProcess/Queue.h>           // ... JobManager::queue()

#include "gtest/gtest.h"
#include "utils.h"

class xSquaredPlusBVectorSerial {
public:
   xSquaredPlusBVectorSerial(double b, std::vector<double> x_init) : _b(b), x(std::move(x_init)), result(x.size()) {}

   void evaluate()
   {
      // call evaluate_task for each task
      for (std::size_t ix = 0; ix < x.size(); ++ix) {
         result[ix] = std::pow(x[ix], 2) + _b;
      }
   }

   std::vector<double> get_result()
   {
      evaluate();
      return result;
   }

   // for simplicity of the example (avoiding getters/setters) we make data members public as well
   double _b;
   std::vector<double> x;
   std::vector<double> result;
};

class xSquaredPlusBVectorParallel : public RooFit::MultiProcess::Job {
public:
   xSquaredPlusBVectorParallel(xSquaredPlusBVectorSerial* serial)
      : serial(serial)
   {
   }

   void evaluate()
   {
      if (get_manager()->process_manager().is_master()) {
         // master fills queue with tasks
         for (std::size_t task_id = 0; task_id < serial->x.size(); ++task_id) {
            RooFit::MultiProcess::JobTask job_task(id, task_id);
            get_manager()->queue().add(job_task);
         }

         // wait for task results back from workers to master
         gather_worker_results();
      }
   }

   std::vector<double> get_result()
   {
      evaluate();
      return serial->result;
   }

   void update_real(std::size_t /*ix*/, double val, bool /*is_constant*/) override {
      if (get_manager()->process_manager().is_worker()) {
         serial->_b = val;
      }
   }

   void update_bool(std::size_t /*ix*/, bool /*value*/) override {
      // pass
   }

   // -- BEGIN plumbing --

   void receive_task_result_on_queue(std::size_t task, std::size_t worker_id) override {
      double result = get_manager()->messenger().template receive_from_worker_on_queue<double>(worker_id);
      serial->result[task] = result;
   }

   void send_back_task_result_from_worker(std::size_t task) override {
      get_manager()->messenger().template send_from_worker_to_queue(serial->result[task]);
   }

   void send_back_results_from_queue_to_master() override {
//      get_manager()->messenger().send_from_queue_to_master(serial->result.size());
      for (std::size_t task_ix = 0ul; task_ix < serial->result.size(); ++task_ix) {
         get_manager()->messenger().send_from_queue_to_master(task_ix, serial->result[task_ix]);
      }
   }

   void clear_results() override {
      // no need to clear any results cache since we just reuse the result vector on the queue
   }

   void receive_results_on_master() override {
//      std::size_t N_job_tasks = get_manager()->messenger().template receive_from_queue_on_master<std::size_t>();
//      for (std::size_t task_ix = 0ul; task_ix < N_job_tasks; ++task_ix) {
      for (std::size_t task_ix = 0ul; task_ix < serial->result.size(); ++task_ix) {
         std::size_t task_id = get_manager()->messenger().template receive_from_queue_on_master<std::size_t>();
         serial->result[task_id] = get_manager()->messenger().template receive_from_queue_on_master<double>();
      }
   }

   bool receive_task_result_on_master(const zmq::message_t & /*message*/) override {
      // TODO: implement; this no-op placeholder is just to make everything compile first so I can check whether merge was successful
      return true;
   }

   // -- END plumbing --


private:
   void evaluate_task(std::size_t task) override
   {
      assert(get_manager()->process_manager().is_worker());
      serial->result[task] = std::pow(serial->x[task], 2) + serial->_b;
   }

   xSquaredPlusBVectorSerial* serial;
};

class TestMPJob : public ::testing::TestWithParam<std::size_t> {
   // You can implement all the usual fixture class members here.
   // To access the test parameter, call GetParam() from class
   // TestWithParam<T>.
};

TEST_P(TestMPJob, singleJobGetResult)
{
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

   xSquaredPlusBVectorParallel x_sq_plus_b_parallel(&x_sq_plus_b);
   RooFit::MultiProcess::JobManager::default_N_workers = NumCPU;

   auto y_parallel = x_sq_plus_b_parallel.get_result();

   EXPECT_EQ(Hex(y_parallel[0]), Hex(y_expected[0]));
   EXPECT_EQ(Hex(y_parallel[1]), Hex(y_expected[1]));
   EXPECT_EQ(Hex(y_parallel[2]), Hex(y_expected[2]));
   EXPECT_EQ(Hex(y_parallel[3]), Hex(y_expected[3]));
}


TEST_P(TestMPJob, multiJobGetResult)
{
   // Simple test case: calculate x^2 + b, where x is a vector. This case does
   // both a simple calculation (squaring the input vector x) and represents
   // handling of state updates in b.
   std::vector<double> x{0, 1, 2, 3};
   double b_initial = 3.;

   std::vector<double> y_expected{3, 4, 7, 12};

   std::size_t NumCPU = GetParam();

   // define jobs
   xSquaredPlusBVectorSerial x_sq_plus_b(b_initial, x);
   xSquaredPlusBVectorSerial x_sq_plus_b2(b_initial + 1, x);
   xSquaredPlusBVectorParallel x_sq_plus_b_parallel(&x_sq_plus_b);
   xSquaredPlusBVectorParallel x_sq_plus_b_parallel2(&x_sq_plus_b2);
   RooFit::MultiProcess::JobManager::default_N_workers = NumCPU;

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


// TODO: implement update_real test!

//TEST_P(TestMPJob, singleJobUpdateReal)
//{
//// Simple test case: calculate x^2 + b, where x is a vector. This case does
//// both a simple calculation (squaring the input vector x) and represents
//// handling of state updates in b.
//std::vector<double> x{0, 1, 2, 3};
//double b_initial = 3.;
//
//// start serial test
//
//xSquaredPlusBVectorSerial x_sq_plus_b(b_initial, x);
//
//auto y = x_sq_plus_b.get_result();
//std::vector<double> y_expected{3, 4, 7, 12};
//
//EXPECT_EQ(Hex(y[0]), Hex(y_expected[0]));
//EXPECT_EQ(Hex(y[1]), Hex(y_expected[1]));
//EXPECT_EQ(Hex(y[2]), Hex(y_expected[2]));
//EXPECT_EQ(Hex(y[3]), Hex(y_expected[3]));
//
//std::size_t NumCPU = GetParam();
//
//// start parallel test
//
//xSquaredPlusBVectorParallel x_sq_plus_b_parallel(NumCPU, &x_sq_plus_b);
//
//auto y_parallel = x_sq_plus_b_parallel.get_result();
//
//EXPECT_EQ(Hex(y_parallel[0]), Hex(y_expected[0]));
//EXPECT_EQ(Hex(y_parallel[1]), Hex(y_expected[1]));
//EXPECT_EQ(Hex(y_parallel[2]), Hex(y_expected[2]));
//EXPECT_EQ(Hex(y_parallel[3]), Hex(y_expected[3]));
//}


INSTANTIATE_TEST_SUITE_P(NumberOfWorkerProcesses, TestMPJob, ::testing::Values(1, 2, 3));
