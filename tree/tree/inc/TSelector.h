// @(#)root/tree:$Id$
// Author: Rene Brun   05/02/97

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TSelector
#define ROOT_TSelector


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TSelector                                                            //
//                                                                      //
// A utility class for Trees selections.                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TObject.h"
#include "TString.h"
#include "TSelectorList.h"

class TTree;

class TSelector : public TObject {

public:
   enum EAbort { kContinue, kAbortProcess, kAbortFile };

protected:
   Long64_t       fStatus;  ///<  Selector status
   EAbort         fAbort;   ///<  Abort status
   TString        fOption;  ///<  Option given to TTree::Process
   TObject       *fObject;  ///<! Current object if processing object (vs. TTree)
   TList         *fInput;   ///<  List of objects available during processing
   TSelectorList *fOutput;  ///<! List of objects created during processing

private:
   TSelector(const TSelector&);            // not implemented
   TSelector& operator=(const TSelector&); // not implemented

public:
   TSelector();
   virtual            ~TSelector();

   virtual int         Version() const { return 0; }
   virtual void        Init(TTree *) { }
   virtual void        Begin(TTree *) { }
   virtual void        SlaveBegin(TTree *) { }
   virtual Bool_t      Notify() { return kTRUE; }
   virtual const char *GetOption() const { return fOption; }
   virtual Long64_t    GetStatus() const { return fStatus; }
   virtual Int_t       GetEntry(Long64_t /*entry*/, Int_t /*getall*/ = 0) { return 0; }
   virtual Bool_t      ProcessCut(Long64_t /*entry*/);
   virtual void        ProcessFill(Long64_t /*entry*/);
   virtual Bool_t      Process(Long64_t /*entry*/);
   virtual void        ImportOutput(TList *output);
   virtual void        SetOption(const char *option) { fOption = option; }
   virtual void        SetObject(TObject *obj) { fObject = obj; }
   virtual void        SetInputList(TList *input) { fInput = input; }
   virtual void        SetStatus(Long64_t status) { fStatus = status; }
   virtual TList      *GetInputList() const { return fInput; }
   virtual TList      *GetOutputList() const { return fOutput; }
   virtual void        SlaveTerminate() { }
   virtual void        Terminate() { }
   virtual void        Abort(const char *why, EAbort what = kAbortProcess);
   virtual EAbort      GetAbort() const { return fAbort; }
   virtual void        ResetAbort() { fAbort = kContinue; }

   static  TSelector  *GetSelector(const char *filename);
   static  Bool_t      IsStandardDraw(const char *selec);

   ClassDef(TSelector,2)  //A utility class for tree and object processing
};

#endif

