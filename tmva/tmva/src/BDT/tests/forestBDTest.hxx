#include "gtest/gtest.h"
#include "RForestInference.hxx"
#include "TreeHelpers.hxx"

//#include <xgboost/c_api.h> // for xgboost
//#include "generated_files/evaluate_forest2.h"
std::string events_file     = "./data/events.csv";
std::string preds_file      = "./data/python_predictions.csv";
std::string json_model_file = "./data/model.json";
std::string tmp_file        = "./data/tmp.csv";
int         loop_size       = 1002;

int tree_number = 1;

template <typename T, typename ForestType>
void test_predictions(int loop_size = 1, std::string my_tmp_file = "")
{
   DataStruct<T> _data(events_file, preds_file);

   ForestType Forest;
   Forest.LoadFromJson("my_key", json_model_file);
   if (loop_size > 1) {
      Forest.inference(_data.events_pointer, _data.rows, _data.cols, _data.scores.data(), loop_size);
   } else {
      Forest.inference(_data.events_pointer, _data.rows, _data.cols, _data.scores.data());
   }
   _predict<T>(_data.scores.data(), _data.preds.size(), _data.preds);

   if ((!my_tmp_file.empty())) {
      write_csv<bool>(my_tmp_file, _data.preds);
   }

   for (size_t i = 0; i < _data.preds.size(); i++) {
      ASSERT_EQ(_data.preds[i], _data.groundtruth[i][0]);
   }
}

TEST(forestBDT, BranchedPredictionsSingleEvent)
{
   test_predictions<float, ForestBranched<float>>(1, "./data/tmp2.csv");
   test_predictions<double, ForestBranched<double>>();
   test_predictions<long double, ForestBranched<long double>>();
   //
}
TEST(forestBDT, BranchedPredictionsBatch)
{
   test_predictions<float, ForestBranched<float>>(loop_size, "./data/tmp4.csv");
}
TEST(forestBDT, BranchedPredictionsBatchDoubles)
{
   test_predictions<double, ForestBranched<double>>(loop_size);
   test_predictions<long double, ForestBranched<long double>>(loop_size);
}

TEST(forestBDT, BranchlessPredictions)
{
   test_predictions<float, ForestBranchless<float>>();
   test_predictions<double, ForestBranchless<double>>();
   test_predictions<long double, ForestBranchless<long double>>();
}
TEST(forestBDT, BranchlessPredictionsBatch)
{
   test_predictions<float, ForestBranchless<float>>(loop_size);
   // test_predictions<double, ForestBranchless<double>>(loop_size);
   // test_predictions<long double, ForestBranchless<long double>>(loop_size);
}

TEST(forestBDT, JitForestPredictions)
{
   test_predictions<float, ForestBranchedJIT<float>>();
   test_predictions<double, ForestBranchedJIT<double>>();
   test_predictions<long double, ForestBranchedJIT<long double>>();

   // test_predictions<float, ForestBranchedJIT<float>>(loop_size);
   // test_predictions<double, ForestBranchedJIT<double>>(loop_size);
   // test_predictions<long double, ForestBranchedJIT<long double>>(loop_size);
}
TEST(forestBDT, JitForestPredictionsBatch)
{
   test_predictions<float, ForestBranchedJIT<float>>(loop_size);
   test_predictions<double, ForestBranchedJIT<double>>(loop_size);
   test_predictions<long double, ForestBranchedJIT<long double>>(loop_size);
}

TEST(forestBDT, JitForestBranchless)
{
   test_predictions<float, ForestBranchlessJIT<float>>();
   test_predictions<double, ForestBranchlessJIT<double>>();
   test_predictions<long double, ForestBranchlessJIT<long double>>();
}
TEST(forestBDT, JitForestBranchlessBatch)
{
   test_predictions<float, ForestBranchlessJIT<float>>(loop_size);
   test_predictions<double, ForestBranchlessJIT<double>>(loop_size);
   test_predictions<long double, ForestBranchlessJIT<long double>>(loop_size);
}
