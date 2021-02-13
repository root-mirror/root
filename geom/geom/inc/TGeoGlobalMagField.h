// @(#)root/geom:$Id$

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TGeoGlobalMagField
#define ROOT_TGeoGlobalMagField

#include "TObject.h"

#include "TVirtualMagField.h"

class TGeoGlobalMagField : public TObject
{
private:
   static TGeoGlobalMagField *fgInstance;     // Static pointer to the field manager;
   TVirtualMagField       *fField{nullptr};   // Magnetic field
   Bool_t                  fLock{kFALSE};     // Lock flag for global field.

protected:
   TGeoGlobalMagField(const TGeoGlobalMagField&) = delete;
   TGeoGlobalMagField& operator=(const TGeoGlobalMagField&) = delete;
   void                    Unlock() {fLock = kFALSE;}

public:
   TGeoGlobalMagField();
   virtual ~TGeoGlobalMagField();

   // Using SetField() makes a given field global. The field manager owns it from now on.
   TVirtualMagField       *GetField() const {return fField;}
   void                    SetField(TVirtualMagField *field);
   Bool_t                  IsLocked() {return fLock;}
   void                    Lock();

   // The field manager should be accessed via TGeoGlobalMagField::Instance()
   static TGeoGlobalMagField *Instance();
   static TGeoGlobalMagField *GetInstance();

   // Inline access to Field() method
   void                    Field(const Double_t *x, Double_t *B) {if (fField) fField->Field(x,B);}

   ClassDef(TGeoGlobalMagField, 0)              // Global field manager
};

#endif
