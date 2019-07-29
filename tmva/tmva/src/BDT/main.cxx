#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <map>
#include <vector>
#include <array>
#include <chrono>
#include <ctime>      // for date
#include <functional> // for std::fucntion

#include "bdt_helpers.h"
#include "TInterpreter.h" // for gInterpreter
//#include "TMVA/RTensor.hxx"
#include <xgboost/c_api.h>

#define BDT_KIND 2

#include "unique_bdt.h"
#include "array_bdt.h"
#if BDT_KIND == 1
#include "bdt.h"
using namespace shared;
#elif BDT_KIND == 2
// using namespace unique_bdt;
#elif BDT_KIND == 3
// using namespace array_bdt;
#endif

//#include "array_bdt.h"
#include "jitted_bdt.h"
#include "forest.h"

using json = nlohmann::json;
// ------------------- helper functions -------------------
json read_file(const std::string &filename)
{
   std::ifstream i(filename);
   json          j;
   i >> j;
   return j;
}

void print(const std::string message)
{
   std::cout << message << std::endl;
}

void print_json_type(json j)
{
   std::string my_type = "unknown";
   if (j.type() == json::value_t::null) {
      my_type = "null";
   } else if (j.type() == json::value_t::boolean) {
      my_type = "boolean";
   } else if (j.type() == json::value_t::number_integer) {
      my_type = "number_integer";
   } else if (j.type() == json::value_t::number_unsigned) {
      my_type = "number_unsigned";
   } else if (j.type() == json::value_t::number_float) {
      my_type = "number_float";
   } else if (j.type() == json::value_t::object) {
      my_type = "object";
   } else if (j.type() == json::value_t::array) {
      my_type = "array";
   } else if (j.type() == json::value_t::string) {
      my_type = "string";
   }
   std::cout << "Type: " << my_type << '\n';
}

/* // params of the children
std::string params[7] = {"depth",
                          "no",
                          "split_condition",
                          "nodeid",
                          "split",
                          "yes",
                          "children"};
*/
void check_json(json &jTree)
{
   bool node_has_children = (jTree.find("children") != jTree.end());
   if (!node_has_children) {
      std::cerr << "Warning: node has no childrens!!!\n";
   }
   if (jTree["yes"] != jTree["children"][0]["nodeid"]) {
      std::cerr << "Implementation error for reading nodes.";
   }
}

// ---------------------------  read json  ------------------------------
// TODO
int main()
{
   std::cout << "\n\n\n ########## READING MAIN.CXX ##########\n\n";

   std::cout << "\n ***** READ JSON *****\n";
   // json config = read_file("model.json"); // get the model as a json object
   std::string my_config = read_file_string("model.json"); // get model as string
   // Parse the string with all the model
   auto json_model = json::parse(my_config);
   std::cout << "Json read, there are " << json_model.size() << " trees in the forest.\n";
   int                number_of_trees = json_model.size();
   std::vector<float> event_sample{1., 115., 70., 30.}; // event to test trees

   std::cout << "\n\n ***** Create unique_ptr representation ***** \n";
   unique_bdt::Tree trees[number_of_trees];

   for (int i = 0; i < number_of_trees; i++) {
      unique_bdt::read_nodes_from_tree(json_model[i], trees[i]);
      std::cout << trees[i].nodes->split_variable << "  " << trees[i].nodes->split_threshold << "\n";
   }
   for (auto &tree : trees) {
      std::cout << "unique_ptr pred: " << tree.inference(event_sample) << std::endl;
   }

   std::cout << "\n\n ***** Create array representation ***** \n";
   array_bdt::Tree trees_array[number_of_trees];
   for (int i = 0; i < number_of_trees; i++) {
      array_bdt::read_nodes_from_tree(json_model[i], trees_array[i]);
   }
   for (auto &tree : trees_array) {
      std::cout << "array pred: " << tree.inference(event_sample) << std::endl;
   }

   std::cout << "\n\n ***** Create Jitted representation ***** \n";
   time_t      my_time          = time(0);
   std::string s_namespace_name = std::to_string(my_time);
   std::cout << "current time used as namespace: " << s_namespace_name << std::endl;
   std::string s_trees[number_of_trees];

   for (int i = 0; i < number_of_trees; i++) {
      std::filebuf fb;
      std::string  filename = "./generated_files/generated_tree_" + std::to_string(i) + ".h";
      fb.open(filename, std::ios::out);
      std::ostream os(&fb);
      generate_code_bdt(os, trees[i], i);
      fb.close();

      std::stringstream ss;
      generate_code_bdt(ss, trees[i], i, s_namespace_name);
      s_trees[i] = ss.str();
   }
   // Read functions
   /*
   std::function<float(std::vector<float>)>              func;
   std::vector<std::function<float(std::vector<float>)>> function_vector;
   for (int i = 0; i < number_of_trees; i++) {
      // func = jit_function_reader_file(i); //, s_trees[i]);
      func = jit_function_reader_string(i, s_trees[i], s_namespace_name);

      function_vector.push_back(func);
   }

   for (auto &tree : function_vector) {
      std::cout << "jitted pred: " << tree(event_sample) << std::endl;
   }
  */

   std::cout << "\n\n ***** Entering benchmarking section ***** \n";
   std::string                     data_folder   = "./data_files/";
   std::string                     events_file   = data_folder + "events.csv";
   std::vector<std::vector<float>> events_vector = read_csv(events_file);

   float              prediction = 0; // define used variables
   std::vector<float> preds_tmp;
   std::vector<bool>  preds;
   float              preds_sum;

   std::cout << "\n\n ***** Benchmarking unique ***** \n";
   preds.clear();
   for (auto &event : events_vector) {
      preds_tmp.clear();
      for (auto &tree : trees) {
         prediction = tree.inference(event);
         preds_tmp.push_back(prediction);
      }
      preds_sum = vec_sum(preds_tmp);
      preds.push_back(binary_logistic(preds_sum));
   }
   std::string preds_unique_file = data_folder + "preds_unique_file.csv";
   write_csv(preds_unique_file, preds); // write predictions

   std::cout << "\n\n ***** Benchmarking array ***** \n";
   preds.clear();
   for (auto &event : events_vector) {
      preds_tmp.clear();
      for (auto &tree : trees_array) {
         prediction = tree.inference(event);
         preds_tmp.push_back(prediction);
      }
      preds_sum = vec_sum(preds_tmp);
      preds.push_back(binary_logistic(preds_sum));
   }
   std::string preds_array_file = data_folder + "preds_array_file.csv";
   write_csv(preds_unique_file, preds); // write predictions

   std::cout << "\n\n ***** tests ***** \n";
   Forest<int> test1;
   test1.test();
   test1.get_Forest();

   Forest<unique_bdt::Tree> test2;
   test2.test();
   test2.get_Forest("model.json");
   preds.clear();
   std::string preds_file = "./data_files/test2.csv";
   preds                  = test2.do_predictions(events_vector);
   write_csv(preds_file, preds);

   Forest<array_bdt::Tree> test3;
   test3.test();
   test3.get_Forest("model.json");
   preds.clear();
   preds      = test3.do_predictions(events_vector);
   preds_file = "./data_files/test3.csv";
   write_csv(preds_file, preds);

   Forest<std::function<float(std::vector<float>)>> test4;
   test4.test();
   test4.get_Forest("model.json");
   preds.clear();
   preds = test4.do_predictions(events_vector);
   preds.clear();
   preds = test4.do_predictions(events_vector);
   preds.clear();
   preds      = test4.do_predictions(events_vector);
   preds_file = "./data_files/test4.csv";
   write_csv(preds_file, preds);

   std::cout << std::to_string(get_time()) << std::endl;
   std::cout << get_time_string() << std::endl;
   std::cout << test4.counter << std::endl;

   //   xgboost ------------------------------------------------
   // create the train data

   int   cols = 3, rows = 5;
   float train[rows][cols];
   for (int i = 0; i < rows; i++)
      for (int j = 0; j < cols; j++) train[i][j] = (i + 1) * (j + 1);

   float train_labels[rows];
   for (int i = 0; i < rows; i++) train_labels[i] = 1 + i * i * i;

   // convert to DMatrix
   DMatrixHandle h_train[1];
   XGDMatrixCreateFromMat((float *)train, rows, cols, -1, &h_train[0]);

   // load the labels
   XGDMatrixSetFloatInfo(h_train[0], "label", train_labels, rows);

   // read back the labels, just a sanity check
   bst_ulong    bst_result;
   const float *out_floats;
   XGDMatrixGetFloatInfo(h_train[0], "label", &bst_result, &out_floats);
   for (unsigned int i = 0; i < bst_result; i++) std::cout << "label[" << i << "]=" << out_floats[i] << std::endl;

   // create the booster and load some parameters
   BoosterHandle h_booster;
   XGBoosterCreate(h_train, 1, &h_booster);
   XGBoosterSetParam(h_booster, "booster", "gbtree");
   XGBoosterSetParam(h_booster, "objective", "reg:linear");
   XGBoosterSetParam(h_booster, "max_depth", "3");
   XGBoosterSetParam(h_booster, "eta", "0.1");
   XGBoosterSetParam(h_booster, "min_child_weight", "1");
   XGBoosterSetParam(h_booster, "subsample", "0.5");
   XGBoosterSetParam(h_booster, "colsample_bytree", "1");
   XGBoosterSetParam(h_booster, "num_parallel_tree", "1");

   // perform 200 learning iterations
   for (int iter = 0; iter < 200; iter++) XGBoosterUpdateOneIter(h_booster, iter, h_train[0]);

   // predict
   const int sample_rows = 5;
   float     test[sample_rows][cols];
   for (int i = 0; i < sample_rows; i++)
      for (int j = 0; j < cols; j++) test[i][j] = (i + 1) * (j + 1);
   DMatrixHandle h_test;
   XGDMatrixCreateFromMat((float *)test, sample_rows, cols, -1, &h_test);
   bst_ulong    out_len;
   const float *f;
   XGBoosterPredict(h_booster, h_test, 0, 0, &out_len, &f);

   for (unsigned int i = 0; i < out_len; i++) std::cout << "prediction[" << i << "]=" << f[i] << std::endl;

   XGBoosterDumpModel(h_booster, "", int 0, &out_len, const char ***out_dump_array);

   // free xgboost internal structures
   XGDMatrixFree(h_train[0]);
   XGDMatrixFree(h_test);
   XGBoosterFree(h_booster);

   // */
   std::cout << "\n ########## END MAIN.CXX ##########\n\n\n";
   return 0;
} // End main
