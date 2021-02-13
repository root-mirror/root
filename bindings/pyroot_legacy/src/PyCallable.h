// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Aug 2004

#ifndef PYROOT_PYCALLABLE_H
#define PYROOT_PYCALLABLE_H

// Bindings
#include "TCallContext.h"


namespace PyROOT {

   class ObjectProxy;

   class PyCallable {
   public:
      virtual ~PyCallable() {}

   public:
      virtual PyObject* GetSignature() = 0;
      virtual PyObject* GetPrototype() = 0;
      virtual PyObject* GetDocString() { return GetPrototype(); }

      virtual Int_t GetPriority() = 0;

      virtual Int_t GetMaxArgs() = 0;
      virtual PyObject* GetCoVarNames() = 0;
      virtual PyObject* GetArgDefault( Int_t /* iarg */ )  = 0;

      virtual PyObject* GetScopeProxy()  = 0;

      virtual PyCallable* Clone() = 0;

   public:
      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt = 0 ) = 0;
   };

} // namespace PyROOT

#endif // !PYROOT_PYCALLABLE_H
