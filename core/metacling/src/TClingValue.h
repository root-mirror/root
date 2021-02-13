// @(#)root/core/meta:$Id$e
// Author: Vassil Vassilev   13/03/2013

/*******************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.                     *
 * All rights reserved.                                                        *
 *                                                                             *
 * For the licensing terms see $ROOTSYS/LICENSE.                               *
 * For the list of contributors see $ROOTSYS/README/CREDITS.                   *
 ******************************************************************************/

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Class representing a value coming from cling. Its main use case           //
//  is to TCallFunc. When TCallFunc returns by-value, i.e. a temporary        //
//  variable, its lifetime has to be extended. TClingValue provides a         //
//  way to extend the temporaries lifetime and gives the user to control it.  //
//                                                                            //
//  The class is used to hide the implementation details of                   //
//  cling::Value.                                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TClingValue
#define ROOT_TClingValue

#include "RtypesCore.h"
#include "TInterpreterValue.h"
#include <string>
#include <utility>

namespace cling {
   class Value;
}

class TClingValue : public TInterpreterValue {
private:
   struct HasTheSameSizeAsClingValue {
      long double fBiggestElementOfUnion;
      int   fStorageType;
      void* fType;
      void* fInterpreter;
   } fValue;

   cling::Value& ToCV() {
      return reinterpret_cast<cling::Value&>(fValue); }
   const cling::Value& ToCV() const {
      return reinterpret_cast<const cling::Value&>(fValue); }

public:
   TClingValue();
   TClingValue(const TClingValue& Other);
   TClingValue& operator=(TClingValue &Other);
   ~TClingValue();

   const void* GetValAddr() const { return &fValue; }
   void* GetValAddr() { return &fValue; }

   std::pair<std::string, std::string> ToTypeAndValueString() const;
   Bool_t      IsValid() const;
   Double_t    GetAsDouble() const;
   Long_t      GetAsLong() const;
   ULong_t     GetAsUnsignedLong() const;
   void*       GetAsPointer() const;
   std::string ToString() const;
};

#endif // ROOT_TClingValue
