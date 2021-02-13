// Author: Wim Lavrijsen   March 2008

#ifndef ROOT_TPySelector
#define ROOT_TPySelector

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPySelector                                                              //
//                                                                          //
// Python base class equivalent of PROOF TSelector.                         //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


//- ROOT
#include "TSelector.h"

#include "TTree.h"

// Python
struct _object;
typedef _object PyObject;


class TPySelector : public TSelector {
public:
   using TSelector::fStatus;
// using TSelector::fAbort;
   using TSelector::fOption;
   using TSelector::fObject;
   using TSelector::fInput;
   using TSelector::fOutput;

public:
   TTree* fChain;

public:
// ctor/dtor ... cctor and assignment are private in base class
   TPySelector( TTree* /* tree */ = 0, PyObject* self = 0 );
   virtual ~TPySelector();

// TSelector set of forwarded (overridden) methods
   virtual Int_t  Version() const;
   virtual Int_t  GetEntry( Long64_t entry, Int_t getall = 0 );
   virtual Bool_t Notify();

   virtual void   Init( TTree* tree );
   virtual void   Begin( TTree* tree = 0 /* not used */ );
   virtual void   SlaveBegin( TTree* tree );
   virtual Bool_t Process( Long64_t entry );
   virtual void   SlaveTerminate();
   virtual void   Terminate();

   virtual void Abort( const char* why, EAbort what = kAbortProcess );

   ClassDef( TPySelector, 1 );   //Python equivalent base class for PROOF

private:
// private helpers for forwarding to python
   void SetupPySelf();
   PyObject* CallSelf( const char* method, PyObject* pyobject = 0 );

private:
   PyObject* fPySelf;              //! actual python object
};

#endif
