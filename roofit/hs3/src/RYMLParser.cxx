#include "RYMLParser.h"
#include <sstream>
#include <stdexcept>

#include <ryml.hpp>
#include <c4/yml/std/map.hpp>
#include <c4/yml/std/string.hpp>
#include <c4/yml/common.hpp>

namespace {
  size_t count_nlines(c4::csubstr src) {
    // helper function to count the lines
    size_t n = (src.len > 0);
    while(src.len > 0)
      {
        n += (src.begins_with('\n') || src.begins_with('\r'));
        src = src.sub(1);
      }
    return n;
  }

  c4::yml::Tree makeTree(std::istream& is){
    if(!is.good()) throw std::runtime_error("invalid input!");
    std::string s(std::istreambuf_iterator<char>(is), {});
    
    auto src = c4::to_csubstr(s.c_str());
    size_t nlines = count_nlines(src);
    c4::yml::Tree tree(nlines,s.size());
    c4::yml::Parser np;
    np.parse({}, tree.copy_to_arena(src), &tree);
    return tree;
  }
}

// TRYMLTree Implementation
class TRYMLTree::Impl {
  ryml::Tree tree;
public:
  Impl(){};
  Impl(std::istream& is) : tree(makeTree(is)) {
    // constructor parsing the RYML tree
  } 
  inline ryml::Tree& get(){
    // getter for the RYML tree reference
    return this->tree;
  }
  inline const ryml::Tree& get() const {
    // const getter for the RYML tree reference    
    return this->tree;
  }
};

// TRYMLNode Implementation
class TRYMLTree::Node::Impl {
  c4::yml::NodeRef node;
public:
  Impl(const c4::yml::NodeRef& n) : node(n) {
    // constructor taking the RYML node reference
  };
  inline static TRYMLTree::Node& mkNode(TRYMLTree* t,c4::yml::NodeRef node){
    // function for creating a new node based on a RYML node reference
    return t->incache(TRYMLTree::Node(t,TRYMLTree::Node::Impl(node)));
  }  
  inline c4::yml::NodeRef& get(){
    // getter for the RYML node reference        
    return node;
  }
  inline const c4::yml::NodeRef& get() const {
    // getter for the RYML node reference            
    return node;
  }
};

// JSONTree interface implementation

void TRYMLTree::Node::writeJSON(std::ostream& os) const {
  // write the tree as JSON to an ostream
  os << c4::yml::as_json(node->get());
}

void TRYMLTree::Node::writeYML(std::ostream& os) const {
  // write the tree as YML to an ostream  
  os << node->get();
}

void TRYMLTree::Node::set_map() {
  // assign this node to be a map (JSON object)
  node->get() |= c4::yml::MAP;
}

void TRYMLTree::Node::set_seq() {
  // assign this node to be a sequence (JSON array)
  node->get() |= c4::yml::SEQ;
}

TRYMLTree::TRYMLTree(std::istream& is) : tree(new Impl(is)) {
  // constructor taking an istream (for reading)
};

TRYMLTree::TRYMLTree() : tree(new Impl()) {
  // default constructor (for writing)
};

TRYMLTree::~TRYMLTree(){
  // destructor. clears the cache.
  clearcache();
};

// JSONNode interface implementation

TRYMLTree::Node::Node(TRYMLTree* t, const TRYMLTree::Node::Impl& imp) : tree(t), node(new Impl(imp)) {
  // construct a new node from scratch
}

TRYMLTree::Node::Node(const Node& other) : Node(other.tree,other.node->get()) {
  // copy constructor
}   

TRYMLTree::Node& TRYMLTree::rootnode(){
  // obtain the root node of a tree
  return Node::Impl::mkNode(this,tree->get().rootref());
}

TRYMLTree::Node& TRYMLTree::Node::operator<< (std::string const& s) {
  // write a string to this node
  node->get() << s;
  return *this;
}

TRYMLTree::Node& TRYMLTree::Node::operator<< (int i) {
  // write an int to this node  
  node->get() << i;
  return *this;
}

TRYMLTree::Node& TRYMLTree::Node::operator<< (double d) {
  // write an float to this node    
  node->get() << d;
  return *this;
}

const TRYMLTree::Node& TRYMLTree::Node::operator>> (std::string &v) const {
  // read a string from this node
  node->get() >> v;
  return *this;
}

TRYMLTree::Node& TRYMLTree::Node::operator[] (std::string const& k) {
  // get a child node with the given key
  return Impl::mkNode(tree,node->get()[c4::to_csubstr(tree->incache(k))]);
}

TRYMLTree::Node& TRYMLTree::Node::operator[] (size_t pos) {
  // get a child node with the given index
  return Impl::mkNode(tree,node->get()[pos]);     
}

const TRYMLTree::Node& TRYMLTree::Node::operator[] (std::string const& k) const {
  // get a child node with the given key (const version)  
  return Impl::mkNode(tree,node->get()[c4::to_csubstr(tree->incache(k))]);     
}

const TRYMLTree::Node& TRYMLTree::Node::operator[] (size_t pos) const {
  // get a child node with the given index (const version)    
  return Impl::mkNode(tree,node->get()[pos]);          
}

bool TRYMLTree::Node::is_container() const {
  // return true if this node can have child nodes
  return node->get().is_container();
}

bool TRYMLTree::Node::is_map() const {
  // return true if this node is a map (JSON object)
  return node->get().is_map();
}

bool TRYMLTree::Node::is_seq() const {
  // return true if this node is a sequence (JSON array)  
  return node->get().is_seq();
}

std::string TRYMLTree::Node::key() const {
  // obtain the key of this node
  std::stringstream ss;    
  ss << node->get().key();
  return ss.str();
}

std::string TRYMLTree::Node::val() const {;
  // obtain the value of this node (as a string)
  std::stringstream ss;    
  ss << node->get().val();
  return ss.str();
}

TRYMLTree::Node& TRYMLTree::Node::append_child(){
  // append a new child to this node
  return Impl::mkNode(tree,node->get().append_child());
}

bool TRYMLTree::Node::has_key() const {
  // return true if this node has a key
  return node->get().has_key();
}

bool TRYMLTree::Node::has_val() const {
  // return true if this node has a value  
  return node->get().has_val();
}

bool TRYMLTree::Node::has_child(std::string const&s) const {
  // return true if this node has a child with the given key
  return node->get().has_child(c4::to_csubstr(s.c_str()));
}

size_t TRYMLTree::Node::num_children() const {
  // return the number of child nodes of this node
  return node->get().num_children();
}

TRYMLTree::Node& TRYMLTree::Node::child(size_t pos){
  // return the child with the given index
  return Impl::mkNode(tree,node->get().child(pos));
}

const TRYMLTree::Node& TRYMLTree::Node::child(size_t pos) const {
  // return the child with the given index (const version)
  return Impl::mkNode(tree,node->get().child(pos));
}

// specific functions

namespace {
  void error_cb(const char* msg, size_t msg_len, void* /*user_data*/){
    // error callback using std::runtime_error
    if(msg && msg_len > 0){
      throw std::runtime_error(msg);
    } else {
      throw std::runtime_error("error handler invoked without error message");
    }
  }
  
  bool setcallbacks(){
    // set the custom callback functions
    c4::yml::set_callbacks(c4::yml::Callbacks(c4::yml::get_callbacks().m_user_data,
                                              c4::yml::get_callbacks().m_allocate,
                                              c4::yml::get_callbacks().m_free,
                                              &::error_cb));
    return true;
  }
  bool ok = setcallbacks();
}

const char* TRYMLTree::incache(const std::string& str){
  // obtain a string from the string cache
  _strcache.push_back(str);
  return _strcache.back().c_str();
}

TRYMLTree::Node& TRYMLTree::incache(const TRYMLTree::Node& n){
  // obtain a node from the node cache  
  _nodecache.push_back(n);
  return _nodecache.back();
}

void TRYMLTree::clearcache(){
  // clear all caches
  TRYMLTree::_strcache.clear();
  TRYMLTree::_nodecache.clear();  
}

