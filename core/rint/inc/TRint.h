// @(#)root/rint:$Id$
// Author: Rene Brun   17/02/95

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TRint
#define ROOT_TRint

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// Rint                                                                 //
//                                                                      //
// Rint is the ROOT Interactive Interface. It allows interactive access //
// to the ROOT system via a C++ interpreter.                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TApplication.h"
#include "TString.h"

class TFileHandler;


class TRint : public TApplication {

private:
   Int_t         fNcmd;               // command history number
   TString       fDefaultPrompt;      // default prompt: "root [%d] "
   TString       fNonContinuePrompt;  // default prompt before continue prompt was set
   char          fPrompt[64];         // interpreter prompt
   Bool_t        fInterrupt;          // if true macro execution will be stopped
   Int_t         fCaughtSignal;       // TRint just caught a signal
   TFileHandler *fInputHandler;       // terminal input handler

   TRint(const TRint&) = delete;
   TRint& operator=(const TRint&) = delete;

   void    ExecLogon();
   Long_t  ProcessRemote(const char *line, Int_t *error = nullptr);
   Long_t  ProcessLineNr(const char* filestem, const char *line, Int_t *error = nullptr);

public:
   TRint(const char *appClassName, int *argc, char **argv,
         void *options = 0, int numOptions = 0, Bool_t noLogo = kFALSE);
   virtual             ~TRint();
   virtual char       *GetPrompt();
   virtual const char *SetPrompt(const char *newPrompt);
   virtual void        SetEchoMode(Bool_t mode);
   virtual void        HandleException(Int_t sig);
   virtual Bool_t      HandleTermInput();
   virtual void        PrintLogo(Bool_t lite = kFALSE);
   virtual void        Run(Bool_t retrn = kFALSE);
   virtual void        Terminate(int status);
           void        Interrupt() { fInterrupt = kTRUE; }
   virtual Int_t       TabCompletionHook(char *buf, int *pLoc, std::ostream& out);

   TFileHandler       *GetInputHandler() { return fInputHandler; }

   ClassDef(TRint,0);  //ROOT Interactive Application Interface
};

#endif
