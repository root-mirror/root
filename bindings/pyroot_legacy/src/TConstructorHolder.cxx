// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

// Bindings
#include "PyROOT.h"
#include "TConstructorHolder.h"
#include "Executors.h"
#include "ObjectProxy.h"
#include "TMemoryRegulator.h"

// Standard
#include <string>


//- protected members --------------------------------------------------------
Bool_t PyROOT::TConstructorHolder::InitExecutor_( TExecutor*& executor, TCallContext* )
{
// pick up special case new object executor
   executor = CreateExecutor( "__init__" );
   return kTRUE;
}

//- public members -----------------------------------------------------------
PyObject* PyROOT::TConstructorHolder::GetDocString()
{
// GetMethod() may return an empty function if this is just a special case place holder
   const std::string& clName = Cppyy::GetFinalName( this->GetScope() );
   return PyROOT_PyUnicode_FromFormat( "%s::%s%s",
      clName.c_str(), clName.c_str(), this->GetMethod() ? this->GetSignatureString().c_str() : "()" );
}

////////////////////////////////////////////////////////////////////////////////
/// preliminary check in case keywords are accidently used (they are ignored otherwise)

PyObject* PyROOT::TConstructorHolder::Call(
      ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt )
{
   if ( kwds != 0 && PyDict_Size( kwds ) ) {
      PyErr_SetString( PyExc_TypeError, "keyword arguments are not yet supported" );
      return 0;
   }

// do not allow instantiation of abstract classes
   if ( Cppyy::IsAbstract( this->GetScope() ) ) {
      PyErr_Format( PyExc_TypeError,
         "%s is abstract and can not be instantiated", Cppyy::GetFinalName( this->GetScope() ).c_str() );
      return 0;
   }

// setup as necessary
   if ( ! this->Initialize( ctxt ) )
      return 0;                              // important: 0, not Py_None

// fetch self, verify, and put the arguments in usable order
   if ( ! ( args = this->PreProcessArgs( self, args, kwds ) ) )
      return 0;

// translate the arguments
   if ( ! this->ConvertAndSetArgs( args, ctxt ) ) {
      Py_DECREF( args );
      return 0;
   }

// perform the call, 0 makes the other side allocate the memory
   Long_t address = (Long_t)this->Execute( 0, 0, ctxt );

// done with filtered args
   Py_DECREF( args );

// return object if successful, lament if not
   if ( address != 0 ) {
      Py_INCREF( self );

   // note: constructors are no longer set to take ownership by default; instead that is
   // decided by the method proxy (which carries a creator flag) upon return
      self->Set( (void*)address );

   // allow lookup upon destruction on the ROOT/CINT side for TObjects
      static Cppyy::TCppType_t sTObjectType = (Cppyy::TCppType_t)Cppyy::GetScope( "TObject" );
   // TODO: cache IsSubtype and offset results ...
      if ( Cppyy::IsSubtype( GetScope(), sTObjectType ) ) {
         TObject* object = (TObject*)(address + \
            Cppyy::GetBaseOffset( GetScope(), sTObjectType, (void*)address, 1 /* up-cast */ ) );
         TMemoryRegulator::RegisterObject( self, object );
      }

   // done with self
      Py_DECREF( self );

      Py_INCREF( Py_None );
      return Py_None;                        // by definition
   }

   if ( ! PyErr_Occurred() )   // should be set, otherwise write a generic error msg
      PyErr_SetString( PyExc_TypeError, const_cast< char* >(
         ( Cppyy::GetFinalName( GetScope() ) + " constructor failed" ).c_str() ) );

// do not throw an exception, '0' might trigger the overload handler to choose a
// different constructor, which if all fails will throw an exception
   return 0;
}
