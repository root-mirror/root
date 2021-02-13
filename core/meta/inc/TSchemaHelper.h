// @(#)root/core:$Id$
// author: Lukasz Janyst <ljanyst@cern.ch>

/*************************************************************************
 * Copyright (C) 1995-2010, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSchemaHelper
#define ROOT_TSchemaHelper

#include "RtypesCore.h"

#include <string>

namespace ROOT {
namespace Internal {
   struct TSchemaHelper
   {
      TSchemaHelper(): fTarget(), fSourceClass(),
       fSource(), fCode(), fVersion(), fChecksum(),
       fInclude(), fEmbed(true), fFunctionPtr(nullptr),
       fAttributes() {}
      std::string fTarget;
      std::string fSourceClass;
      std::string fSource;
      std::string fCode;
      std::string fVersion;
      std::string fChecksum;
      std::string fInclude;
      Bool_t      fEmbed;
      void*       fFunctionPtr;
      std::string fAttributes;

      TSchemaHelper(const TSchemaHelper &tsh) :
       fTarget(tsh.fTarget), fSourceClass(tsh.fSourceClass),
       fSource(tsh.fSource), fCode(tsh.fCode), fVersion(tsh.fVersion),fChecksum(tsh.fChecksum),
       fInclude(tsh.fInclude), fEmbed(tsh.fEmbed), fFunctionPtr(tsh.fFunctionPtr),
       fAttributes(tsh.fAttributes) {}

      TSchemaHelper& operator=(const TSchemaHelper &) {return *this;} // Not implemented
   };
}
}

#endif // ROOT_TSchemaHelper
