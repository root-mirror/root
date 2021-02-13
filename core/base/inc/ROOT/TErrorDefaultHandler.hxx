// @(#)root/thread:$Id$
// Author: Jakob Blomer, June 2020

/*************************************************************************
 * Copyright (C) 1995-2006, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TErrorDefaultHandler
#define ROOT_TErrorDefaultHandler

namespace ROOT {
namespace Internal {

/// Destructs resources that are taken by using the default error handler.
/// This function is called during the destruction of gROOT.
void ReleaseDefaultErrorHandler();

}
}

#endif
