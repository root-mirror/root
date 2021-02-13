// @(#)root/hist:$Id$
// Author: Rene Brun   31/08/99

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TVirtualPaveStats
#define ROOT_TVirtualPaveStats

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TVirtualPaveStats                                                    //
//                                                                      //
// Abstract base class for PaveStats                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "Rtypes.h"

class TObject;

class TVirtualPaveStats {

public:
   virtual ~TVirtualPaveStats() = default;

   virtual TObject *GetParent() const = 0;
   virtual void SetParent(TObject *) = 0;

   ClassDef(TVirtualPaveStats, 0)  //Abstract interface for TPaveStats
};

#endif
