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
//  Class representing a value coming from the interpreter. Its main use case //
//  is to TCallFunc. When TCallFunc returns by-value, i.e. a temporary        //
//  variable, its lifetime has to be extended. TInterpreterValue provides a   //
//  way to extend the temporaries lifetime and gives the user to control it.  //
//                                                                            //
//  The class needs to be derived from for the actual interpreter,            //
//  see TClingValue.                                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TInterpreterValue
#define ROOT_TInterpreterValue

#include "RtypesCore.h"

#include <string>
#include <utility>

class TInterpreterValue {
private:
   TInterpreterValue(const TInterpreterValue &) = delete;
   TInterpreterValue& operator=(TInterpreterValue &) = delete;
public:
   TInterpreterValue() { }
   virtual ~TInterpreterValue() { }

   virtual const void* GetValAddr() const = 0;
   virtual void* GetValAddr() = 0;
   virtual std::pair<std::string, std::string> ToTypeAndValueString() const = 0;

   virtual Bool_t      IsValid() const = 0;
   virtual Double_t    GetAsDouble() const = 0;
   virtual Long_t      GetAsLong() const = 0;
   virtual ULong_t     GetAsUnsignedLong() const = 0;
   virtual void*       GetAsPointer() const = 0;
   virtual std::string ToString() const = 0;
};

#endif // ROOT_TInterpreterValue
