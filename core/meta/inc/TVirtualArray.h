// @(#)root/io:$Id$
// Author: Philippe Canal July, 2008

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
  *************************************************************************/

#ifndef ROOT_TVirtualArray
#define ROOT_TVirtualArray



/**
\class TVirtualArray
\ingroup IO
Wrapper around an object and giving indirect access to its content
even if the object is not of a class in the Cint/Reflex dictionary.
*/

#include "TClass.h"
#include "TClassRef.h"

class TVirtualArray {
public:
   using ObjectPtr = TClass::ObjectPtr;

   TClassRef  fClass;
   UInt_t     fCapacity;
   UInt_t     fSize;
   ObjectPtr  fArray; ///< fSize elements

   TVirtualArray( TClass *cl, UInt_t size ) : fClass(cl), fCapacity(size), fSize(size), fArray( ( cl ? cl->NewObjectArray(size) : ObjectPtr{nullptr, nullptr}) ) {};
   ~TVirtualArray() { if (fClass) fClass->DeleteArray( fArray ); }

   TClass *GetClass() { return fClass; }
   char *operator[](UInt_t ind) const { return GetObjectAt(ind); }
   char *GetObjectAt(UInt_t ind) const { return ((char*)fArray.GetPtr())+fClass->Size()*ind; }

   void SetSize(UInt_t size) {
      // Set the used size of this array to 'size'.   If size is greater than the existing
      // capacity, reallocate the array BUT no data is preserved.
      fSize = size;
      if (fSize > fCapacity && fClass) {
         fClass->DeleteArray( fArray );
         fArray = fClass->NewObjectArray(fSize);
         fCapacity = fSize;
      }
   }


};

#endif // ROOT_TVirtualArray
