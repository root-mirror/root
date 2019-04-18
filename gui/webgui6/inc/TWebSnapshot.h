// Author:  Sergey Linev, GSI,  6/04/2017

/*************************************************************************
 * Copyright (C) 1995-2018, Rene Brun and Fons Rademakers.               *
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

#include "TObject.h"

#include <vector>
#include <memory>
#include <string>

class TWebSnapshot : public TObject {

protected:
   std::string fObjectID;       ///<   object identifier
   std::string fOption;         ///<   object draw option
   Int_t fKind{0};              ///<   kind of snapshots
   TObject *fSnapshot{nullptr}; ///<   snapshot data
   Bool_t fOwner{kFALSE};       ///<!  if objected owned

   void SetKind(Int_t kind) { fKind = kind; }

public:

   enum {
     kNone = 0,        // dummy
     kObject = 1,      // object itself
     kSVG = 2,         // list of SVG primitives
     kSubPad = 3,      // subpad
     kColors = 4,      // list of ROOT colors + palette
     kStyle = 5        // gStyle object
   };

   virtual ~TWebSnapshot();

   void SetObjectIDAsPtr(void *ptr);
   void SetObjectID(const std::string &id) { fObjectID = id; }
   const char* GetObjectID() const { return fObjectID.c_str(); }

   void SetOption(const std::string &opt) { fOption = opt; }

   void SetSnapshot(Int_t kind, TObject *snapshot, Bool_t owner = kFALSE);
   Int_t GetKind() const { return fKind; }
   TObject *GetSnapshot() const { return fSnapshot; }

   ClassDef(TWebSnapshot,1)  // Object painting snapshot, used for JSROOT
};

// =================================================================================

class TPadWebSnapshot : public TWebSnapshot {
protected:
   bool fActive{false};                                    ///< true when pad is active
   bool fReadOnly{true};                                   ///< when canvas or pad are in readonly mode
   std::vector<std::unique_ptr<TWebSnapshot>> fPrimitives; ///< list of all primitives, drawn in the pad
public:
   TPadWebSnapshot() { SetKind(kSubPad); }

   void SetActive(bool on = true) { fActive = on; }

   void SetReadOnly(bool on = true) { fReadOnly = on; }

   TWebSnapshot &NewPrimitive(TObject *obj = nullptr, const std::string &opt = "");

   TPadWebSnapshot &NewSubPad();

   TWebSnapshot &NewSpecials();

   ClassDef(TPadWebSnapshot,1)  // Pad painting snapshot, used for JSROOT
};

#endif
