#ifndef ROO_JSON_LIST_FILE
#define ROO_JSON_LIST_FILE

#include <string>
#include <fstream>
#include <vector>


class RooJsonListFile {
public:
  // ctors
  RooJsonListFile(): _member_index(0) {}
  RooJsonListFile(const std::string & filename);
  // default move ctors
  RooJsonListFile(RooJsonListFile&& other) = default;
  RooJsonListFile& operator=(RooJsonListFile&& other) = default;
  // dtor
  ~RooJsonListFile();

  void open(const std::string & filename);

  template <class Iter>
  void set_member_names(Iter begin, Iter end, bool reset_index = true);

  RooJsonListFile& add_member_name(const std::string &name);

  template <typename T>
  RooJsonListFile& operator<< (const T& obj);

private:
  std::ofstream _out;
  std::vector<std::string> _member_names;
  unsigned long _next_member_index();
  unsigned long _member_index;
};


template <class Iter>
void RooJsonListFile::set_member_names(Iter begin, Iter end, bool reset_index) {
  _member_names.clear();
  for(Iter it = begin; it != end; ++it) {
    _member_names.push_back(*it);
  }
  if (reset_index) {
    _member_index = 0;
  }
}

////////////////////////////////////////////////////////////////////////////////
/// This method assumes that std::ofstream::operator<<(T) exists.

template <typename T>
RooJsonListFile& RooJsonListFile::operator<< (const T& obj)
{
  auto ix = _next_member_index();
  if (ix == 0) {
    _out << "{";
  }

  // `"member name": `
  _out << "\"" << _member_names[ix] << "\": ";
  // `"value"` (comma added below, if not last value in list element)
  _out << "\"" << obj << "\"";

  if (ix == _member_names.size() - 1) {
    _out << "},\n";
  } else {
    _out << ", ";
  }

  return *this;
}

#endif