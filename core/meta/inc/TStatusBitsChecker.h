// @(#)root/meta:$Id$
// Author: Philippe Canal, 2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TStatusBitsChecker
#define ROOT_TStatusBitsChecker

#include <map>
#include <list>

#include "RtypesCore.h"

class TClass;

namespace ROOT {
namespace Detail {

class TStatusBitsChecker {
protected:
   static UChar_t ConvertToBit(Long64_t constant, TClass &classRef, const char *constantName);

public:
   class Registry {
   protected:
      struct Info;

      std::map<UChar_t, std::list<Info>> fRegister; ///<! Register of bits seen so far.

   public:

      void RegisterBits(TClass &classRef);

      bool Check(TClass &classRef, bool verbose = false);

      Registry();  // Implemented in source file to allow hiding of the Info struct.
      ~Registry(); // Implemented in source file to allow hiding of the Info struct.
   };

   static bool Check(TClass &classRef, bool verbose = false);
   static bool Check(const char *classname, bool verbose = false);
   static bool CheckAllClasses(bool verbosity = false);
};

} // Details
} // ROOT

#endif // ROOT__TStatusBitsChecker
