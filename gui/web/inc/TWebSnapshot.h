// Author: Sergey Linev   6/04/2017

/*************************************************************************
 * Copyright (C) 2017, Sergey Linev                                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/


#ifndef ROOT_TWebSnapshot
#define ROOT_TWebSnapshot

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TWebSnapshot                                                         //
//                                                                      //
// Paint state of object to transfer to JavaScript side                 //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif

#include <vector>

class TWebSnapshot : public TObject {

protected:

   TString    fObjectID;   ///<   object identifier
   TString    fOption;     ///<   object draw option
   Int_t      fKind;       ///<   kind of snapshots
   TObject*   fSnapshot;   ///<   snapshot data

   void SetKind(Int_t kind) { fKind = kind; }

public:

   enum {
     kNone = 0,         // dummy
     kObject = 1,       // object itself
     kSVG = 2,          // list of SVG primitives
     kSubPad = 3        // subpad
   };

   TWebSnapshot();
   virtual ~TWebSnapshot();

   void SetObjectIDAsPtr(void* ptr);
   void SetObjectID(const char* id) { fObjectID = id; }
   const char* GetObjectID() const { return fObjectID.Data(); }

   void SetOption(const char* opt) { fOption = opt; }

   void SetSnapshot(Int_t kind, TObject* shot);
   Int_t GetKind() const { return fKind; }
   TObject* GetSnapshot() const { return fSnapshot; }

   ClassDef(TWebSnapshot,1)  // Object painting snapshot, used for JSROOT
};

class TPadWebSnapshot : public TWebSnapshot {
protected:
   std::vector<TWebSnapshot*> fPrimitives;
public:
   TPadWebSnapshot() : TWebSnapshot(), fPrimitives() { SetKind(kSubPad); }
   virtual ~TPadWebSnapshot();

   void Add(TWebSnapshot *snap) { fPrimitives.push_back(snap); }

   ClassDef(TPadWebSnapshot,1)  // Object painting snapshot, used for JSROOT
};
#endif
