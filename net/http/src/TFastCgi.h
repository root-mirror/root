// $Id$
// Author: Sergey Linev   28/12/2013

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFastCgi
#define ROOT_TFastCgi

#include "THttpEngine.h"

#include <thread>
#include <memory>

class TFastCgi : public THttpEngine {
protected:
   Int_t fSocket{0};            ///<! socket used by fastcgi
   Bool_t fDebugMode{kFALSE};   ///<! debug mode, may required for fastcgi debugging in other servers
   TString fTopName;            ///<! name of top item
   std::unique_ptr<std::thread> fThrd;  ///<! thread which takes requests, can be many later
   Bool_t fTerminating{kFALSE};     ///<! set when http server wants to terminate all engines

   void Terminate() override { fTerminating = kTRUE; }

public:
   TFastCgi();
   virtual ~TFastCgi();

   Bool_t Create(const char *args) override;

   Int_t GetSocket() const { return fSocket; }

   Bool_t IsTerminating() const { return fTerminating; }

   Bool_t IsDebugMode() const { return fDebugMode; }

   const char *GetTopName() const { return fTopName.Length() > 0 ? fTopName.Data() : nullptr; }
};

#endif
