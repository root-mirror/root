#include "Minuit2/MnPrint.h"

using ROOT::Minuit2::MnPrint;

#ifndef USE_ROOT_ERROR

#include <iostream>

#ifndef MN_OS
#define MN_OS std::cerr
#endif

void MnPrint::Impl(MnPrint::Verbosity level, const std::string &s)
{
   const char *label[4] = {"[Error]", "[Warn]", "[Info]", "[Debug]"};
   const int ilevel = static_cast<int>(level);
   MN_OS << label[ilevel] << " " << s << std::endl;
}

#else // use ROOT error reporting system

#include "TError.h"
#include <sstream>

void MnPrint::Impl(MnPrint::Verbosity level, const std::string &s)
{
   switch (level) {
   case MnPrint::eError: ::Error("Minuit2", "%s", s.c_str()); break;
   case MnPrint::eWarn: ::Warning("Minuit2", "%s", s.c_str()); break;
   case MnPrint::eInfo:
   case MnPrint::eDebug: ::Info("Minuit2", "%s", s.c_str()); break;
   }
}

#endif // USE_ROOT_ERROR
