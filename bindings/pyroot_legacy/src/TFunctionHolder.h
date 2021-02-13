// Author: Wim Lavrijsen, Apr 2005

#ifndef PYROOT_TFUNCTIONHOLDER_H
#define PYROOT_TFUNCTIONHOLDER_H

// Bindings
#include "TMethodHolder.h"


namespace PyROOT {

   class TFunctionHolder : public TMethodHolder {
   public:
      using TMethodHolder::TMethodHolder;

      virtual PyCallable* Clone() { return new TFunctionHolder( *this ); }

      virtual PyObject* PreProcessArgs( ObjectProxy*& self, PyObject* args, PyObject* kwds );
      virtual PyObject* Call(
         ObjectProxy*&, PyObject* args, PyObject* kwds, TCallContext* ctx = 0 );
   };

} // namespace PyROOT

#endif // !PYROOT_TFUNCTIONHOLDER_H
