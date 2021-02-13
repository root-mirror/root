// @(#)root/base:$Id$
// Author: Victor Perev and Philippe Canal   08/05/02

/*************************************************************************
 * Copyright (C) 1995-2003, Rene Brun, Fons Rademakers and al.           *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TClassStreamer_h
#define ROOT_TClassStreamer_h

#include "TClassRef.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TClassStreamer is used to stream an object of a specific class.      //
//                                                                      //
// The address passed to operator() will be the address of the start    //
// of the  object.                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

class TClassStreamer {
protected:
   TClassStreamer() : fStreamer(0) {};
   TClassStreamer(const TClassStreamer &rhs) : fStreamer(rhs.fStreamer), fOnFileClass() {};
   TClassStreamer &operator=(const TClassStreamer &rhs) {   fOnFileClass = rhs.fOnFileClass; fStreamer = rhs.fStreamer; return *this; }

public:
   TClassStreamer(ClassStreamerFunc_t pointer) : fStreamer(pointer), fOnFileClass() {};

   virtual void SetOnFileClass( const TClass* cl ) { fOnFileClass = const_cast<TClass*>(cl); }
   virtual const TClass* GetOnFileClass() const { return fOnFileClass; }

   virtual TClassStreamer *Generate() const {
      // Virtual copy constructor.
      return new TClassStreamer(*this);
   }

   virtual  ~TClassStreamer(){};
   virtual void operator()(TBuffer &b, void *objp)
   {
      // The address passed to operator() will be the address of the start of the
      // object.

      (*fStreamer)(b,objp);
   }
   virtual void Stream(TBuffer &b, void *objp, const TClass *onfileClass)
   {
      // The address passed to operator() will be the address of the start of the
      // object.   Overload this routine, if your derived class can optimize
      // the handling of the onfileClass (rather than storing and restoring from the
      // fOnFileClass member.

      // Note we can not name this routine 'operator()' has it would be slightly
      // backward incompatible and lead to the following warning/error from the
      // compiler in the derived class overloading the other operator():
//      include/TClassStreamer.h:51: error: ‘virtual void TClassStreamer::operator()(TBuffer&, void*, const TClass*)’ was hidden
//      include/TCollectionProxyFactory.h:180: error:   by ‘virtual void TCollectionClassStreamer::operator()(TBuffer&, void*)’
//   cc1plus: warnings being treated as errors

      SetOnFileClass(onfileClass);
      (*this)(b,objp);
   }

private:
   ClassStreamerFunc_t fStreamer;
protected:
   TClassRef           fOnFileClass;
};

#endif
