#include "RooJsonListFile.h"

RooJsonListFile::RooJsonListFile(const std::string & filename) :
    _member_index(0) {
  open(filename);
}

void RooJsonListFile::open(const std::string & filename) {
  // do not use ios::app for opening out!
  // app moves put pointer to end of file before each write, which makes seekp useless.
  // See http://en.cppreference.com/w/cpp/io/basic_filebuf/open
  _out.open(filename, std::ios_base::in | std::ios_base::out);  // "mode r+"
  if (!_out.is_open()) {
    _out.clear();
    // new file
    _out.open(filename, std::ios_base::out);  // "mode w"
    _out << "[\n";
  } else {
    // existing file that, presumably, has been closed with close_json_list() and thus ends with "\n]".
    _out.seekp(-2, std::ios_base::end);
    _out << ",\n";
  }
}

RooJsonListFile::~RooJsonListFile() {
  _out.seekp(-2, std::ios_base::end);
  _out << "\n]";
}

unsigned long RooJsonListFile::_next_member_index() {
  auto current_index = _member_index;
  _member_index = (_member_index + 1) % _member_names.size();
  return current_index;
}

RooJsonListFile& RooJsonListFile::add_member_name(const std::string &name) {
  _member_names.push_back(name);

  return *this;
}