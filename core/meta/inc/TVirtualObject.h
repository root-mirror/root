// @(#)root/io:$Id$
// Author: Philippe Canal July, 2008

/*************************************************************************
 * Copyright (C) 1995-2000, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
  *************************************************************************/

#ifndef ROOT_TVirtualObject
#define ROOT_TVirtualObject

/**
\class TVirtualObject
\ingroup IO

Wrapper around an object and giving indirect access to its content
even if the object is not of a class in the Cint/Reflex dictionary.
*/

#include "TClassRef.h"


class TVirtualObject {
private:

   TVirtualObject(const TVirtualObject&) = delete;
   TVirtualObject &operator=(const TVirtualObject&) = delete;

public:
   TClassRef  fClass;
   void      *fObject;

   TVirtualObject(TClass *cl) : fClass(cl), fObject(cl ? cl->New() : nullptr) { }
   ~TVirtualObject() { if (fClass) fClass->Destructor(fObject); }


   TClass *GetClass() const { return fClass; }
   void   *GetObject() const { return fObject; }

};

#endif // ROOT_TVirtualObject
