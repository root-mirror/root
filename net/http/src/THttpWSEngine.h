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

#include "RtypesCore.h"

#include "THttpCallArg.h"

#include <mutex>
#include <thread>
#include <string>
#include <condition_variable>

class THttpWSHandler;

class THttpWSEngine {

private:
   friend class THttpWSHandler;

   bool fMTSend{false};     ///<!  true when send operation runs, set under locked fMutex from WSHandler
   bool fDisabled{false};   ///<!  true shortly before cleanup, set under locked fMutex from WSHandler

   std::thread fSendThrd;    ///<! dedicated thread for all send operations
   bool fHasSendThrd{false}; ///<! if thread was started one have to call join method for it

   std::mutex fMutex;                                  ///<! protects all data behind
   std::condition_variable fCond;                      ///<! condition used to sync with sending thread
   bool fWaiting{false};                               ///<! if condition wait is called
   bool fSending{false};                               ///<! performing send operation in other thread
   enum { kNone, kData, kHeader, kText } fKind{kNone}; ///<! kind of operation
   std::string fData;                                  ///<! data (binary or text)
   std::string fHdr;                                   ///<! header

protected:
   THttpWSEngine() = default;

   /// Indicate if engine require extra thread to complete postponed thread operation
   virtual Bool_t SupportSendThrd() const { return kFALSE; }

   /// One always can send data to websocket - as long as previous send operation completed
   virtual Bool_t CanSendDirectly() { return kTRUE; }

public:
   virtual ~THttpWSEngine() = default;

   virtual UInt_t GetId() const = 0;

   virtual void ClearHandle(Bool_t) = 0;

   virtual void Send(const void *buf, int len) = 0;

   virtual void SendHeader(const char *hdr, const void *buf, int len) = 0;

   virtual void SendCharStar(const char *str);

   virtual Bool_t PreProcess(std::shared_ptr<THttpCallArg> &arg);

   virtual void PostProcess(std::shared_ptr<THttpCallArg> &arg);
};

#endif
