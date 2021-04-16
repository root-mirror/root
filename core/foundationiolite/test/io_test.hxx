#include "RConfigure.h"
#include "ROOT/RRawFile.hxx"
#include "ROOT/RMakeUnique.hxx"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "gtest/gtest.h"

using RRawFile = ROOT::Internal::RRawFile;

namespace {

/**
 * An RAII wrapper around an open temporary file on disk. It cleans up the guarded file when the wrapper object
 * goes out of scope.
 */
class FileRaii {
private:
   std::string fPath;
public:
   FileRaii(const std::string &path, const std::string &content) : fPath(path)
   {
      std::ofstream ostrm(path, std::ios::binary | std::ios::out | std::ios::trunc);
      ostrm << content;
   }
   FileRaii(const FileRaii&) = delete;
   FileRaii& operator=(const FileRaii&) = delete;
   ~FileRaii() {
      std::remove(fPath.c_str());
   }
};

} // anonymous namespace
