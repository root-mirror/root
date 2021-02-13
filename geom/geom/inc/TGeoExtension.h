/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

// Author: Andrei.Gheata@cern.ch  29/05/2013
// Following proposal by Markus Frank

#ifndef ROOT_TGeoExtension
#define ROOT_TGeoExtension

#include <cassert>

#include "TObject.h"

class TGeoExtension : public TObject
{
protected:
   TGeoExtension() : TObject() {}
   virtual ~TGeoExtension() {}

public:
   // Method called whenever requiring a pointer to the extension
   // Equivalent to new()
   virtual TGeoExtension *Grab() = 0;
   // Method called always when the pointer to the extension is not needed
   // Equivalent to delete()
   virtual void           Release() const = 0;

   ClassDef(TGeoExtension, 1)       // User extension for volumes and nodes
};

class TGeoRCExtension : public TGeoExtension
{
protected:
   virtual ~TGeoRCExtension() {delete fUserObject;}
public:
   TGeoRCExtension() : TGeoExtension(), fRC(0), fUserObject(0) { fRC++; }
   TGeoRCExtension(TObject *obj) : TGeoExtension(), fRC(0), fUserObject(obj) { fRC++; }

   TGeoExtension       *Grab() override             { fRC++; return this; }
   void                 Release() const override    { assert(fRC > 0); fRC--; if (fRC==0) delete this; }

   void                 SetUserObject(TObject *obj) { fUserObject = obj; }
   TObject             *GetUserObject() const       { return fUserObject; }


private:
   // Copy constructor and assignment not allowed
   TGeoRCExtension(const TGeoRCExtension &) = delete;
   TGeoRCExtension &operator =(const TGeoRCExtension &) = delete;
   mutable Int_t        fRC{0};           // Reference counter
   TObject             *fUserObject{nullptr};   // Attached user object

   ClassDefOverride(TGeoRCExtension, 1)       // Reference counted extension for volumes and nodes
};

#endif
