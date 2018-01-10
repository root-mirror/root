// $Id$
// Author: Sergey Linev   20/10/2017

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THttpWSEngine
#define ROOT_THttpWSEngine

#include "Rtypes.h"

class THttpCallArg;

class THttpWSEngine {

protected:
   THttpWSEngine() = default;

public:
   virtual ~THttpWSEngine() {}

   void AttachTo(THttpCallArg &);

   virtual UInt_t GetId() const = 0;

   virtual void ClearHandle() = 0;

   virtual void Send(const void *buf, int len) = 0;

   virtual void SendCharStar(const char *str);

   virtual Bool_t PreviewData(THttpCallArg &);

   virtual void PostProcess(THttpCallArg &);
};

#endif
