// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Aug 2004

#ifndef PYROOT_TCLASSMETHODHOLDER_H
#define PYROOT_TCLASSMETHODHOLDER_H

// Bindings
#include "TMethodHolder.h"


namespace PyROOT {

   class TClassMethodHolder : public TMethodHolder {
   public:
      using TMethodHolder::TMethodHolder;

      virtual PyCallable* Clone() { return new TClassMethodHolder( *this ); }
      virtual PyObject* Call(
         ObjectProxy*&, PyObject* args, PyObject* kwds, TCallContext* ctxt = 0 );
   };

} // namespace PyROOT

#endif // !PYROOT_TCLASSMETHODHOLDER_H
