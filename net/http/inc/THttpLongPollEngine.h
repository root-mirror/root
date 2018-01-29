// $Id$
// Author: Sergey Linev   8/01/2018

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_THttpLongPollEngine
#define ROOT_THttpLongPollEngine

#include "THttpWSEngine.h"

#include <string>
#include <vector>

class THttpLongPollEngine : public THttpWSEngine {
protected:
   THttpCallArg *fPoll;              ///!< polling request, which can be used for the next sending
   std::vector<std::string> fQueue;  ///!< entries submitted to client
   static const char *gLongPollNope; ///!< default reply on the longpoll request
public:
   THttpLongPollEngine() : THttpWSEngine(), fPoll(nullptr), fQueue() {}

   virtual UInt_t GetId() const;

   virtual void ClearHandle();

   virtual void Send(const void *buf, int len);

   virtual void SendCharStar(const char *buf);

   virtual Bool_t PreviewData(THttpCallArg &arg);

   virtual void PostProcess(THttpCallArg &arg);
};

#endif
