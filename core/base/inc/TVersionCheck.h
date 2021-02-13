// @(#)root/base:$Id$
// Author: Fons Rademakers   9/5/2007

/*************************************************************************
 * Copyright (C) 1995-2007, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVersionCheck
#define ROOT_TVersionCheck

#include "RVersion.h"

/** \class TVersionCheck
\ingroup Base

Used to check if the shared library or plugin is compatible with
the current version of ROOT.
*/

class TVersionCheck {
public:
   TVersionCheck(int versionCode);  // implemented in TSystem.cxx
};

namespace ROOT {
namespace Internal {
static TVersionCheck gVersionCheck(ROOT_VERSION_CODE);
}
}
#endif
