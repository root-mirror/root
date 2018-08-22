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

#include "THttpCallArg.h"

#include <mutex>
#include <string>

class THttpWSHandler;

class THttpWSEngine {

private:
   friend class THttpWSHandler;

   bool fMTSend{false};     ///<!  true when send operation runs, set under locked fMutex from WSHandler
   bool fDisabled{false};   ///<!  set shortly before cleanup

   std::mutex fDataMutex;                              ///<! protects data submited for send operation
   enum { kNone, kData, kHeader, kText } fKind{kNone}; ///<! kind of operation
   bool fDoingSend{false};                             ///<! doing send operation in other thread
   std::string fData;                                  ///<! data (binary or text)
   std::string fHdr;                                   ///<! header

protected:
   THttpWSEngine() = default;

   /// Indicate if engine require extra thread to complete postponed thread operation
   virtual Bool_t SupportSendThrd() const { return kFALSE; }

   /// One always can send data to websocket - as long as previous send operation completed
   virtual Bool_t CanSendDirectly() { return kTRUE; }

public:
   virtual ~THttpWSEngine() {}

   /// Returns kTRUE when handle is deactivated and need to be destroyed
   Bool_t IsDisabled() const { return fDisabled; }

   virtual UInt_t GetId() const = 0;

   virtual void ClearHandle() = 0;

   virtual void Send(const void *buf, int len) = 0;

   virtual void SendHeader(const char *hdr, const void *buf, int len) = 0;

   virtual void SendCharStar(const char *str);

   virtual Bool_t PreProcess(std::shared_ptr<THttpCallArg> &arg);

   virtual void PostProcess(std::shared_ptr<THttpCallArg> &arg);
};

#endif
