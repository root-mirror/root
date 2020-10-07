#ifndef JSON_INTERFACE_H
#define JSON_INTERFACE_H

#include <string>
#include <stdexcept>
#include <iostream>
#include <vector>

class JSONNode {
 protected:
  template<class Nd> class child_iterator_t {
    Nd& node;
    size_t pos;
  public:
    child_iterator_t(Nd& n, size_t p) : node(n), pos(p) {};

    child_iterator_t& operator++ () { ++pos; return *this; }
    child_iterator_t& operator-- () { --pos; return *this; }
    
    Nd& operator*  () const { return node.child(pos); }
    Nd& operator-> () const { return node.child(pos); }
    
    bool operator!= (const child_iterator_t& that) const { return this->pos != that.pos; };
    bool operator== (const child_iterator_t& that) const { return this->pos == that.pos; };
  };
  template<class Nd> class children_view_t {
    child_iterator_t<Nd> b, e;
  public:
    inline children_view_t(child_iterator_t<Nd> const& b_, child_iterator_t<Nd> const& e_) : b(b_), e(e_) {}
    
    inline child_iterator_t<Nd> begin() const { return b; }
    inline child_iterator_t<Nd> end  () const { return e; }
  };  

 public:
  virtual void writeJSON(std::ostream& os) const = 0;
  virtual void writeYML(std::ostream&) const { throw std::runtime_error("YML not supported"); }  
    
 public:
  virtual JSONNode& operator<< (std::string const& s) = 0;
  virtual JSONNode& operator<< (int i) = 0;
  virtual JSONNode& operator<< (double d) = 0;
  template<class T> JSONNode& operator<< (const std::vector<T>& v){
    this->set_seq(); for(const auto&e:v){ this->append_child() << e; }; return *this;
  }
  virtual const JSONNode& operator>> (std::string &v) const = 0;
  virtual JSONNode& operator[] (std::string const& k) = 0;
  virtual JSONNode& operator[] (size_t pos) = 0;
  virtual const JSONNode& operator[] (std::string const& k) const = 0;
  virtual const JSONNode& operator[] (size_t pos) const = 0;
  virtual bool is_container() const = 0;
  virtual bool is_map() const = 0;
  virtual bool is_seq() const = 0;
  virtual void set_map() = 0;
  virtual void set_seq() = 0;  

  virtual std::string key() const = 0;
  virtual std::string val() const = 0;
  virtual int val_int() const {
    return atoi(this->val().c_str());
  }
  virtual float val_float() const {
    return atof(this->val().c_str());
  }
  virtual bool val_bool() const {
    return atoi(this->val().c_str());
  }    
  virtual bool has_key() const = 0;
  virtual bool has_val() const = 0;  
  virtual bool has_child(std::string const&)const = 0;
  virtual JSONNode& append_child() = 0;  
  virtual size_t num_children() const = 0;

  using       children_view = children_view_t<      JSONNode>;
  using const_children_view = children_view_t<const JSONNode>;
  
  children_view children() {
    return children_view(child_iterator_t<JSONNode>(*this,0),child_iterator_t<JSONNode>(*this,this->num_children()));
  }
  const_children_view children() const {
    return const_children_view(child_iterator_t<const JSONNode>(*this,0),child_iterator_t<const JSONNode>(*this,this->num_children()));
  }    
  virtual JSONNode& child(size_t pos) = 0;
  virtual const JSONNode& child(size_t pos) const = 0;    
};

class JSONTree {
  virtual JSONNode& rootnode() = 0;
};



#endif
