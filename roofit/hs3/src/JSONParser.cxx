#include "JSONParser.h"
#include <sstream>

#include "json.hpp"

// TJSONTree interface

TJSONTree::Node& TJSONTree::rootnode() {
  return root;
}

// TJSONTree methods

TJSONTree::TJSONTree() : root(this) {};

TJSONTree::TJSONTree(std::istream& is) : root(this,is) {};

TJSONTree::~TJSONTree(){
  TJSONTree::_nodecache.clear();  
};

TJSONTree::Node& TJSONTree::incache(const TJSONTree::Node& n){
  _nodecache.push_back(n);
  return _nodecache.back();
}

void TJSONTree::clearcache(){
  TJSONTree::_nodecache.clear();
}

// TJSONTree::Node implementation

class TJSONTree::Node::Impl {
public:
  std::string _key;
  std::string key() const { return _key; }
  virtual nlohmann::json& get() = 0;
  virtual const nlohmann::json& get() const = 0;  
  class BaseNode;
  class NodeRef;
  Impl(const std::string& k) : _key(k) {}
  static TJSONTree::Node& mkNode(TJSONTree* t, const std::string& k, nlohmann::json& n);
};
  
class TJSONTree::Node::Impl::BaseNode : public TJSONTree::Node::Impl {
  nlohmann::json node;
public:
  virtual nlohmann::json& get() override {
    return node;
  }
  virtual const nlohmann::json& get() const override {
    return node;
  }  
  BaseNode(std::istream& is) : Impl(""), node(nlohmann::json::parse(is)) {}
  BaseNode() : Impl("") {}    
};

class TJSONTree::Node::Impl::NodeRef : public TJSONTree::Node::Impl {
  nlohmann::json& node;
public:
  virtual nlohmann::json& get(){
    return node;
  }
  virtual const nlohmann::json& get() const override {
    return node;
  }  
  NodeRef(const std::string& k,nlohmann::json& n) : Impl(k), node(n) {}
  NodeRef(const NodeRef& other) : Impl(other.key()), node(other.node) {}      
};  

TJSONTree::Node& TJSONTree::Node::Impl::mkNode(TJSONTree* t, const std::string& k, nlohmann::json& n){
  Node::Impl::NodeRef ref(k,n);
  return t->incache(Node(t,ref));
}

TJSONTree::Node::Node(TJSONTree* t,std::istream& is) : tree(t), node(new Impl::BaseNode(is)) {}

TJSONTree::Node::Node(TJSONTree* t) : tree(t), node(new Impl::BaseNode()) {}

TJSONTree::Node::Node(TJSONTree* t,Impl& other) : tree(t), node(new Impl::NodeRef(other.key(),other.get())) {}

TJSONTree::Node::Node(const Node& other) : Node(other.tree,*other.node) {}

TJSONTree::Node::~Node(){}

// TJSONNode interface

void TJSONTree::Node::writeJSON(std::ostream& os) const {
  os << node->get();
}

TJSONTree::Node& TJSONTree::Node::operator<< (std::string const& s) {
  node->get() = s;
  return *this;
}

TJSONTree::Node& TJSONTree::Node::operator<< (int i) {
  node->get() = i;
  return *this;
}

TJSONTree::Node& TJSONTree::Node::operator<< (double d) {
  node->get() = d;
  return *this;
}

const TJSONTree::Node& TJSONTree::Node::operator>> (std::string &v) const {
  v = node->get().get<std::string>();
  return *this;
}

TJSONTree::Node& TJSONTree::Node::operator[] (std::string const& k) {
  return Impl::mkNode(tree,k,node->get()[k]);
}

TJSONTree::Node& TJSONTree::Node::operator[] (size_t pos) {
  return Impl::mkNode(tree,"",node->get()[pos]);
}

const TJSONTree::Node& TJSONTree::Node::operator[] (std::string const& k) const {
  return Impl::mkNode(tree,k,node->get()[k]);
}

const TJSONTree::Node& TJSONTree::Node::operator[] (size_t pos) const {
  return Impl::mkNode(tree,"",node->get()[pos]);
}

bool TJSONTree::Node::is_container() const {
  return node->get().is_array() || node->get().is_object();
}

bool TJSONTree::Node::is_map() const {
  return node->get().is_object();
}

bool TJSONTree::Node::is_seq() const {
  return node->get().is_array();
}

void TJSONTree::Node::set_map() {
  if(node->get().type() == nlohmann::json::value_t::null){
    node->get() = nlohmann::json::object();
  } else if(node->get().type() != nlohmann::json::value_t::object){
    throw std::runtime_error("cannot declare "+this->key()+" to be of map-type, already of type "+node->get().type_name());
  }
}

void TJSONTree::Node::set_seq() {
  if(node->get().type() == nlohmann::json::value_t::null){
    node->get() = nlohmann::json::array();
  } else if(node->get().type() != nlohmann::json::value_t::array){
    throw std::runtime_error("cannot declare "+this->key()+" to be of seq-type, already of type "+node->get().type_name());
  }  
}

std::string TJSONTree::Node::key() const {
  return node->key();
}

namespace {
  std::string itoa(int i) {
    std::stringstream ss;
    ss << i;
    return ss.str();
  }
  std::string ftoa(float f) {
    std::stringstream ss;
    ss << f;
    return ss.str();
  }
}

std::string TJSONTree::Node::val() const {
  switch (node->get().type()){
  case nlohmann::json::value_t::string:
    return node->get().get<std::string>();
  case nlohmann::json::value_t::boolean:
    return node->get().get<bool>() ? "true" : "false";
  case nlohmann::json::value_t::number_integer:
    return ::itoa(node->get().get<int>());
  case nlohmann::json::value_t::number_unsigned:
    return ::itoa(node->get().get<unsigned int>());
  case nlohmann::json::value_t::number_float:
    return ::ftoa(node->get().get<float>());
  default:
    throw std::runtime_error(std::string("implicit string conversion for type ")+node->get().type_name()+std::string(" not supported!"));
  }
}

int TJSONTree::Node::val_int() const {
  return node->get().get<int>();
}
float TJSONTree::Node::val_float() const {
  return node->get().get<float>();
}
bool TJSONTree::Node::val_bool() const {
  return node->get().get<bool>();
}    

bool TJSONTree::Node::has_key() const {
  return node->key().size()>0;
}

bool TJSONTree::Node::has_val() const {
  return node->get().is_primitive();
}

bool TJSONTree::Node::has_child(std::string const& c) const {
  return node->get().find(c) != node->get().end();
}

TJSONTree::Node& TJSONTree::Node::append_child() {
  node->get().push_back("");  
  return Impl::mkNode(tree,"",node->get().back());
}

size_t TJSONTree::Node::num_children() const {
  return node->get().size();
}

TJSONTree::Node& TJSONTree::Node::child(size_t pos) {
  auto it=node->get().begin(); 
  for(size_t i=0; i<pos; ++i) ++it;
  return Impl::mkNode(tree,this->is_map() ? it.key() : "",*it);
}

const TJSONTree::Node& TJSONTree::Node::child(size_t pos) const {
  auto it=node->get().begin(); 
  for(size_t i=0; i<pos; ++i) ++it;
  return Impl::mkNode(tree,this->is_map() ? it.key() : "",*it);  
}
