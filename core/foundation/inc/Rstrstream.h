// @(#)root/base:$Id$
// Author: Fons Rademakers   19/7/02

/*************************************************************************
 * Copyright (C) 1995-2002, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_Rstrstream
#define ROOT_Rstrstream

#include <ROOT/RConfig.hxx>

#if defined(R__ANSISTREAM)
#  if defined(R__SSTREAM)
#    include <sstream>
#  else
#    include <strstream>
#  endif
#else
#  ifndef R__WIN32
#    include <strstream.h>
#  else
#    include <strstrea.h>
#  endif
#endif

#endif
