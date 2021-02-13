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
#include <queue>
#include <mutex>

class THttpServer;

class THttpLongPollEngine : public THttpWSEngine {
   friend class THttpServer;

protected:

   enum EBufKind { kNoBuf, kTxtBuf, kBinBuf };

   bool fRaw{false};                    ///!< if true, only content can be used for data transfer
   std::mutex fMutex;                   ///!< protect polling request to use it from different threads
   std::shared_ptr<THttpCallArg> fPoll; ///!< hold polling request, which can be immediately used for the next sending
   EBufKind fBufKind{kNoBuf};           ///!< if buffered data available
   std::string fBuf;                    ///!< buffered data
   std::string fBufHeader;              ///!< buffered header
   static const std::string gLongPollNope;    ///!< default reply on the longpoll request

   std::string MakeBuffer(const void *buf, int len, const char *hdr = nullptr);

   virtual Bool_t CanSendDirectly() override;

public:
   THttpLongPollEngine(bool raw = false);
   virtual ~THttpLongPollEngine() = default;

   UInt_t GetId() const override;

   void ClearHandle(Bool_t) override;

   void Send(const void *buf, int len) override;

   void SendCharStar(const char *buf) override;

   void SendHeader(const char *hdr, const void *buf, int len) override;

   Bool_t PreProcess(std::shared_ptr<THttpCallArg> &arg) override;

   void PostProcess(std::shared_ptr<THttpCallArg> &arg) override;
};

#endif
