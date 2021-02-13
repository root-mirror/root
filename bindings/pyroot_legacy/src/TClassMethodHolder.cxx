// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Aug 2004

// Bindings
#include "PyROOT.h"
#include "TClassMethodHolder.h"


//- public members --------------------------------------------------------------
PyObject* PyROOT::TClassMethodHolder::Call(
      ObjectProxy*&, PyObject* args, PyObject* kwds, TCallContext* ctxt )
{
// preliminary check in case keywords are accidently used (they are ignored otherwise)
   if ( kwds != 0 && PyDict_Size( kwds ) ) {
      PyErr_SetString( PyExc_TypeError, "keyword arguments are not yet supported" );
      return 0;
   }

// setup as necessary
   if ( ! this->Initialize( ctxt ) )
      return 0;                              // important: 0, not Py_None

// translate the arguments
   if ( ! this->ConvertAndSetArgs( args, ctxt ) )
      return 0;                              // important: 0, not Py_None

// execute function
   return this->Execute( 0, 0, ctxt );
}
