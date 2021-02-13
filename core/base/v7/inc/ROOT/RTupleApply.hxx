/// \file ROOT/RTupleApply.hxx
/// \ingroup Base StdExt ROOT7
/// \author Axel Naumann <axel@cern.ch>
/// \date 2015-09-06
/// \warning This is part of the ROOT 7 prototype! It will change without notice. It might trigger earthquakes. Feedback is welcome!

/*************************************************************************
 * Copyright (C) 1995-2015, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT7_RTupleApply
#define ROOT7_RTupleApply

#include "RConfigure.h"

#ifdef R__HAS_STD_APPLY
# include <tuple>
#else
# include "ROOT/impl_tuple_apply.hxx"
#endif


#endif
