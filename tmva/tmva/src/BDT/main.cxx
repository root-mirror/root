#include "json.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <streambuf>
#include <map>
#include <vector>

#include "bdt.h"

using json = nlohmann::json;

json read_file(const std::string &filename) {
  std::ifstream i(filename);
  json j;
  i >> j;
  //std::cout << "Read file: " << j.type();
  return j;
}

std::string read_file_string(const std::string &filename){
  std::ifstream t(filename);
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

void print(const std::string message){
  std::cout << message << std::endl;
}

void print_json_type(json j){
  std::string my_type = "unknown";
  if (j.type() == json::value_t::null){my_type = "null";}
  else if (j.type() == json::value_t::boolean){my_type = "boolean";}
  else if (j.type() == json::value_t::number_integer){my_type = "number_integer";}
  else if (j.type() == json::value_t::number_unsigned){my_type = "number_unsigned";}
  else if (j.type() == json::value_t::number_float){my_type = "number_float";}
  else if (j.type() == json::value_t::object){my_type = "object";}
  else if (j.type() == json::value_t::array){my_type = "array";}
  else if (j.type() == json::value_t::string){my_type = "string";}
  std::cout << "Type: " << my_type << '\n';
}



void check_params(json j, int max_counter, int counter=0){

  std::string params[7] = {"depth",
                            "no",
                            "split_condition",
                            "nodeid",
                            "split",
                            "yes",
                            "children"};

  for (size_t i = 0; i<j.size(); i++){
    for (auto &mess : params){
      if ((mess == "children") && (counter<max_counter)){
        counter++;
        check_params(j[i][mess],max_counter, counter);
      }
      else{
        std::cout << mess <<": " <<j[0]["children"][i][mess] << std::endl;
      }
    }
  }
}

void check_params2(json j, int max_counter, int counter=0){

  std::string params[7] = {"depth",
                            "no",
                            "split_condition",
                            "nodeid",
                            "split",
                            "yes",
                            "children"};

  for (size_t i = 0; i<j.size(); i++){
    for (auto &mess : params){
    //for (auto &mess : j[i]){
      if (j[i].count("leaf") > 0){
        std::cout<< "FOUND\n";
      }
      else if ((mess == "children") && (counter<max_counter)){
        counter++;
        check_params2(j[i][mess],max_counter, counter);
      }
      else{
        std::cout << mess <<": " <<j[i]["children"][i][mess] << std::endl;
      }
    }
  }
}

void write_node_members(json &jTree, Node &tmp_node){
  tmp_node.split_value = jTree["split_condition"];
  tmp_node.node_id = jTree["nodeid"];
  tmp_node.child_id_true = jTree["yes"];
  tmp_node.child_id_false = jTree["no"];

  std::string tmp_str = jTree["split"].get<std::string>();
  tmp_str.erase(tmp_str.begin(), tmp_str.begin()+1); // remove initial "f"
  tmp_node.split_variable = std::stoi(tmp_str);
}
/// Need a nlohmann::json object from an xgboost saved format
//AbstractNode* _read_nodes(json jTree, std::vector<AbstractNode> *nodes_vector){
AbstractNode* _read_nodes(json jTree, Nodes &nodes_struct){
  bool node_has_children = (jTree.find("children") != jTree.end());
  bool is_leaf_node = (
    (jTree["children"][0].find("leaf") != jTree["children"][0].end())
    && (jTree["children"][0].find("nodeid") != jTree["children"][0].end())
  );

  if (node_has_children){
    if (jTree["yes"] != jTree["children"][0]["nodeid"]){
      std::cerr << "Implementation error for reading nodes.";
    }
    if (is_leaf_node){
      LeafNode tmp_node;
      write_node_members(jTree, tmp_node);
      tmp_node.leaf_true = jTree["children"][0]["leaf"];
      tmp_node.leaf_false = jTree["children"][1]["leaf"];
      nodes_struct.leaf_nodes.push_back(tmp_node);
      return &tmp_node;
    }
    else {
      Node tmp_node;
      write_node_members(jTree, tmp_node);
      tmp_node.child_true = _read_nodes(jTree["children"][0], nodes_struct);
      tmp_node.child_false = _read_nodes(jTree["children"][1], nodes_struct);
      nodes_struct.normal_nodes.push_back(tmp_node);
      return &tmp_node;
    }
  }
  else {
    std::cerr << "Warning: node has no childrens!!!\n";
  }
  //return 1;
}

//std::vector<AbstractNode>
void read_nodes_from_tree(json jTree,Nodes &Nodes){
  //std::vector<AbstractNode> nodes;
  _read_nodes(jTree, Nodes);
  //return nodes;
}

// */


int main() {
  json config = read_file("model.json"); // get the model as a json object
  std::string my_config = read_file_string("model.json"); // get model as string
  //std::cout << "String: " << my_config << std::endl;

  // Parse the string with all the model
  auto json_model = json::parse(my_config);
  //std::cout << "json stringed" << j3.dump() << std::endl;
  //std::cout << "Size: " <<  j3.size() << std::endl;
  //auto my_type = json_model.type();
  //std::cout << my_type.get<std::string>();

  /*
  int count = 0;
  for (auto &tree : j3){
    std::cout << tree << "\n\n";
    count++;
  }
  */
  std::cout << "\n *************************** \n\n";

  int number_of_trees = json_model.size();
  Tree trees[number_of_trees];
  for (int i=0; i<number_of_trees; i++){
    //trees[i].nodes =
    read_nodes_from_tree(json_model[i], trees[i].nodes);
    std::cout << "Number of nodes : "
              << trees[i].nodes.normal_nodes.size()+trees[i].nodes.leaf_nodes.size()
              << std::endl;
  }

  Tree test = trees[0];
  double event[4] = {0,1,2,3};
  /*
  for (auto &node : test.nodes.leaf_nodes){
    std::cout << node.split_value << std::endl;
    std::cout << node.inference(event) << std::endl;
  }*/
/*
  for (auto &node : test.nodes.normal_nodes){
    std::cout << node.split_value << std::endl;
    std::cout << node.inference(event) << std::endl;
  }
*/


/*
  Forest my_forest;
  for (auto &tree : trees){
    my_forest.trees.push_back(tree);
  }
*/
  //my_forest.trees = trees;
  //std::vector<AbstractNode> nodes = read_nodes_from_tree(json_model[0]);
  //std::cout << "Number of nodes : " << nodes.size() << std::endl;
  //tree_1.nodes = nodes;

  // template faillure
  //Node2<double> leaf;
  //Node2<Node2<double>> normal_node;
  Node3 ttt;




  std::cout << "\n***** END *****" << std::endl;
} // End main
