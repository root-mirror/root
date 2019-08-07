#include <benchmark/benchmark.h>
#include <iostream>
#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <map>
#include <vector>
#include <array>
#include <utility>

#include "TInterpreter.h" // for gInterpreter

#include "bdt_helpers.h"

#include "unique_bdt.h"
#include "array_bdt.h"
#include "forest.h"

#include <xgboost/c_api.h> // for xgboost
#include "generated_files/evaluate_forest2.h"

using json = nlohmann::json;

#define safe_xgboost(call)                                                                          \
   {                                                                                                \
      int err = (call);                                                                             \
      if (err != 0) {                                                                               \
         fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError()); \
         exit(1);                                                                                   \
      }                                                                                             \
   }

/// Benchmark eval unique_bdts
static void BM_EvalUniqueBdt(benchmark::State &state)
{
   Forest<unique_bdt::Tree> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions(events_vector, preds);
   }
   write_csv(preds_file, preds);
}
// /*
BENCHMARK(BM_EvalUniqueBdt)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });
// */

/// Benchmark eval unique_bdts
static void BM_EvalUniqueBdt_batch_32(benchmark::State &state)
{
   Forest<unique_bdt::Tree> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions_batch2(events_vector, preds, 32);
   }
   write_csv(preds_file, preds);
}
BENCHMARK(BM_EvalUniqueBdt_batch_32)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });

/// Benchmark eval unique_bdts
static void BM_EvalUniqueBdt_batch_128(benchmark::State &state)
{
   Forest<unique_bdt::Tree> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions_batch2(events_vector, preds, 128);
   }
   write_csv(preds_file, preds);
}
BENCHMARK(BM_EvalUniqueBdt_batch_128)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });

/// Benchmark eval unique_bdts
static void BM_EvalUniqueBdt_batch_256(benchmark::State &state)
{
   Forest<unique_bdt::Tree> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions_batch2(events_vector, preds, 256);
   }
   write_csv(preds_file, preds);
}
BENCHMARK(BM_EvalUniqueBdt_batch_256)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });

/// Benchmark eval array_bdts
static void BM_EvalArrayBdt(benchmark::State &state)
{

   Forest<array_bdt::Tree> Forest;

   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string events_file = "./data_files/events.csv";
   std::string preds_file  = "./data_files/test.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);
   for (auto _ : state) { // only bench what is inside the loop
      // for (int i = 0; i < 1000; i++)
      Forest.do_predictions(events_vector, preds);
   }
   write_csv(preds_file, preds);
}
// /*
BENCHMARK(BM_EvalArrayBdt)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });
// */
/// Benchmark eval Jitted_bdts
static void BM_EvalJittedBdt(benchmark::State &state)
{
   Forest<std::function<float(std::vector<float>)>> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);
   preds.reserve(events_vector.size());
   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions(events_vector, preds);
   }
   write_csv(preds_file, preds);
}
/*
BENCHMARK(BM_EvalJittedBdt)
  ->Unit(benchmark::kMillisecond)
  ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
     return *(std::min_element(std::begin(v), std::end(v)));
  });
//   */

/// Benchmark eval Jitted_bdts
static void BM_EvalJitForestBdt(benchmark::State &state)
{
   Forest<std::function<bool(std::vector<float>)>> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   std::function<bool(std::vector<float>)> my_func = Forest.trees[0];
   preds.reserve(events_vector.size());
   /*
   for (auto _ : state) { // only bench what is inside the loop
      for (size_t i = 0; i < events_vector.size(); i++) {
         preds[i] = my_func(events_vector[i]);
      }
   }
   */
   // /*
   for (auto _ : state) { // only bench what is inside the loop
      Forest.do_predictions(events_vector, preds);
   }
   // */
   write_csv(preds_file, preds);
}
// /*
BENCHMARK(BM_EvalJitForestBdt)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });
// */

/// Benchmark eval Jitted_bdts
static void BM_EvalJitForestWholeBdt(benchmark::State &state)
{
   Forest<std::function<std::vector<bool>(std::vector<std::vector<float>>)>> Forest;
   Forest.get_Forest("model.json");
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   std::function<std::vector<bool>(std::vector<std::vector<float>>)> my_func = Forest.trees[0];
   preds.reserve(events_vector.size());
   for (auto _ : state) { // only bench what is inside the loop
      preds = my_func(events_vector);
   }
   write_csv(preds_file, preds);
}
BENCHMARK(BM_EvalJitForestWholeBdt)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });

/// Benchmark eval Jitted_bdts
static void BM_StaticForestWholeBdt_event(benchmark::State &state)
{
   std::vector<bool> preds;
   preds.clear();
   std::string preds_file  = "./data_files/test.csv";
   std::string events_file = "./data_files/events.csv";

   std::vector<std::vector<float>> events_vector = read_csv(events_file);
   preds.reserve(events_vector.size());

   std::function<std::vector<bool>(std::vector<std::vector<float>>)> my_func = s_f_event_31564752128::evaluate_forest;

   for (auto _ : state) { // only bench what is inside the loop
      preds = my_func(events_vector);
   }
   write_csv(preds_file, preds);
}
/*
BENCHMARK(BM_StaticForestWholeBdt_event)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });
   // */

/// Benchmark eval xgboost_bdt
static void BM_EvalXgboostBdt(benchmark::State &state)
{
   std::string events_fname = "data_files/events.csv";
   std::string preds_fname  = "data_files/python_predictions.csv";
   const char *model_fname  = "./data/model.rabbit";

   std::vector<std::vector<float>> events;
   std::vector<std::vector<float>> labels;
   events = read_csv(events_fname);
   labels = read_csv(preds_fname);

   int cols = events[0].size();
   int rows = events.size();
   // std::cout << "r" << rows << std::endl;
   // std::cout << cols << std::endl;

   std::vector<std::vector<float>> events2 = events;

   // float train[rows][cols];
   // std::cout << "Trains: " << train << std::endl;
   // std::cout << &train << std::endl;
   // float train3[rows * cols];

   std::vector<float> train2;
   train2.reserve(rows * cols);

   float tmp = 0;
   for (int i = 0; i < rows; i++) {
      // std::cout << i << std::endl;
      for (int j = 0; j < cols; j++) {
         // std::cout << j << std::endl;
         train2[i * cols + j] = events.at(i).at(j);
         // tmp         = events.at(i).at(j);
         // train[i][j] = tmp;
      }
   }

   // for (size_t i = 0; i < rows; i++)
   //  for (size_t j = 0; j < cols; j++) train3[i * cols + j] = events.at(i).at(j);

   // for (size_t i = 0; i < rows; i++)
   //  for (size_t j = 0; j < cols; j++) train[i][j] = events.at(i).at(j);

   float m_labels[rows];
   // for (int i = 0; i < rows; i++) m_labels[i] = labels[i][0];
   // Transform to single vector and pass vector.data();
   DMatrixHandle h_train;
   safe_xgboost(XGDMatrixCreateFromMat((float *)train2.data(), rows, cols, -1, &h_train));
   // safe_xgboost(XGDMatrixCreateFromMat((float *)&train, rows, cols, -1, &h_train));

   BoosterHandle boosterHandle;
   safe_xgboost(XGBoosterCreate(0, 0, &boosterHandle));
   // std::cout << "Loading model \n";
   safe_xgboost(XGBoosterLoadModel(boosterHandle, model_fname));
   XGBoosterSetParam(boosterHandle, "objective", "binary:logistic");

   // std::cout << "***** Predicts ***** \n";
   bst_ulong    out_len;
   const float *f;

   for (auto _ : state) { // only bench what is inside the loop
      // for (int i = 0; i < 1000; i++)
      XGBoosterPredict(boosterHandle, h_train, 0, 0, &out_len, &f);
   }

   std::vector<float> preds;
   for (int i = 0; i < out_len; i++) preds.push_back(f[i]);
   std::string preds_file = "data_files/test.csv";
   write_csv(preds_file, preds);

   // free xgboost internal structures
   safe_xgboost(XGBoosterFree(boosterHandle));
}
// /*
BENCHMARK(BM_EvalXgboostBdt)
   ->Unit(benchmark::kMillisecond)
   ->ComputeStatistics("min", [](const std::vector<double> &v) -> double {
      return *(std::min_element(std::begin(v), std::end(v)));
   });
// */

BENCHMARK_MAIN();
