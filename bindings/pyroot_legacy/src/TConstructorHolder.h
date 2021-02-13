// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

#ifndef PYROOT_TCONSTRUCTORHOLDER_H
#define PYROOT_TCONSTRUCTORHOLDER_H

// Bindings
#include "TMethodHolder.h"


namespace PyROOT {

   class TConstructorHolder : public TMethodHolder {
   public:
      using TMethodHolder::TMethodHolder;

   public:
      virtual PyObject* GetDocString();
      virtual PyCallable* Clone() { return new TConstructorHolder( *this ); }

   public:
      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt = 0 );

   protected:
      virtual Bool_t InitExecutor_( TExecutor*&, TCallContext* ctxt = 0 );
   };

} // namespace PyROOT

#endif // !PYROOT_TCONSTRUCTORHOLDER_H
