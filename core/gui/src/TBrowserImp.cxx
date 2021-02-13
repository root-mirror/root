// @(#)root/base:$Id$
// Author: Fons Rademakers   15/11/95

/*************************************************************************
 * Copyright (C) 1995-2021, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

/** \class TBrowserImp
\ingroup Base

ABC describing GUI independent browser implementation protocol.
*/

#include "TBrowserImp.h"

ClassImp(TBrowserImp);

///////////////////////////////////////////////////////////////////
/// Default constructor

TBrowserImp::TBrowserImp(TBrowser *b) :
   fBrowser(b), fShowCycles(kFALSE)
{
}

///////////////////////////////////////////////////////////////////
/// Constructor with browser width and height

TBrowserImp::TBrowserImp(TBrowser *, const char *, UInt_t, UInt_t, Option_t *) :
   fBrowser(nullptr), fShowCycles(kFALSE)
{
}

///////////////////////////////////////////////////////////////////
/// Constructor with browser x, y, width and height

TBrowserImp::TBrowserImp(TBrowser *, const char *, Int_t, Int_t, UInt_t, UInt_t, Option_t *)
   : fBrowser(nullptr), fShowCycles(kFALSE)
{
}
