// Author: Wim Lavrijsen, Oct 2005

#ifndef PYROOT_TSETITEMHOLDER_H
#define PYROOT_TSETITEMHOLDER_H

// Bindings
#include "TMethodHolder.h"


namespace PyROOT {

   class TExecutor;
   class TMemberAdapter;
   class TScopeAdapter;

   class TSetItemHolder : public TMethodHolder {
   public:
      using TMethodHolder::TMethodHolder;

   public:
      virtual PyCallable* Clone() { return new TSetItemHolder( *this ); }
      virtual PyObject* PreProcessArgs( ObjectProxy*& self, PyObject* args, PyObject* kwds );

   protected:
      virtual Bool_t InitExecutor_( TExecutor*&, TCallContext* ctxt = 0 );
   };

} // namespace PyROOT

#endif // !PYROOT_TSETITEMHOLDER_H
