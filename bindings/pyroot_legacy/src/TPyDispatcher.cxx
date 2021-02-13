// Author: Wim Lavrijsen, Aug 2007

// Bindings
#include "PyROOT.h"
#include "TPyDispatcher.h"
#include "RootWrapper.h"

// ROOT
#include "TClass.h"
#include "TObject.h"

// Standard
#include <stdarg.h>


//______________________________________________________________________________
//                         Python callback dispatcher
//                         ==========================
//
// The TPyDispatcher class acts as a functor that can be used for TFn's and GUIs
// to install callbacks from CINT.


//- data ---------------------------------------------------------------------
ClassImp(TPyDispatcher);


//- constructors/destructor --------------------------------------------------
TPyDispatcher::TPyDispatcher( PyObject* callable ) : fCallable( 0 )
{
// Construct a TPyDispatcher from a callable python object. Applies python
// object reference counting.
   Py_XINCREF( callable );
   fCallable = callable;
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Applies python object reference counting.

TPyDispatcher::TPyDispatcher( const TPyDispatcher& other ) : TObject ( other )
{
   Py_XINCREF( other.fCallable );
   fCallable = other.fCallable;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator. Applies python object reference counting.

TPyDispatcher& TPyDispatcher::operator=( const TPyDispatcher& other )
{
   if ( this != &other ) {
      this->TObject::operator=( other );

      Py_XDECREF( fCallable );
      Py_XINCREF( other.fCallable );
      fCallable = other.fCallable;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Reference counting for the held python object is in effect.

TPyDispatcher::~TPyDispatcher() {
   Py_XDECREF( fCallable );
}


//- public members -----------------------------------------------------------
PyObject* TPyDispatcher::DispatchVA( const char* format, ... )
{
// Dispatch the arguments to the held callable python object, using format to
// interpret the types of the arguments. Note that format is in python style,
// not in C printf style. See: https://docs.python.org/2/c-api/arg.html .
   PyObject* args = 0;

   if ( format ) {
      va_list va;
      va_start( va, format );

      args = Py_VaBuildValue( (char*)format, va );

      va_end( va );

      if ( ! args ) {
         PyErr_Print();
         return 0;
      }

      if ( ! PyTuple_Check( args ) ) {    // if only one arg ...
         PyObject* t = PyTuple_New( 1 );
         PyTuple_SET_ITEM( t, 0, args );
         args = t;
      }

   }

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::DispatchVA1( const char* clname, void* obj, const char* format, ... )
{
   PyObject* pyobj = PyROOT::BindCppObject( obj, Cppyy::GetScope( clname ), kFALSE /* isRef */ );
   if ( ! pyobj ) {
      PyErr_Print();
      return 0;
   }

   PyObject* args = 0;

   if ( format ) {
      va_list va;
      va_start( va, format );

      args = Py_VaBuildValue( (char*)format, va );

      va_end( va );

      if ( ! args ) {
         PyErr_Print();
         return 0;
      }

      if ( ! PyTuple_Check( args ) ) {    // if only one arg ...
         PyObject* t = PyTuple_New( 2 );
         PyTuple_SET_ITEM( t, 0, pyobj );
         PyTuple_SET_ITEM( t, 1, args );
         args = t;
      } else {
         PyObject* t = PyTuple_New( PyTuple_GET_SIZE( args ) + 1 );
         PyTuple_SET_ITEM( t, 0, pyobj );
         for ( int i = 0; i < PyTuple_GET_SIZE( args ); i++ ) {
            PyObject* item = PyTuple_GET_ITEM( args, i );
            Py_INCREF( item );
            PyTuple_SET_ITEM( t, i+1, item );
         }
         Py_DECREF( args );
         args = t;
      }
   } else {
      args = PyTuple_New( 1 );
      PyTuple_SET_ITEM( args, 0, pyobj );
   }

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( TPad* selpad, TObject* selected, Int_t event )
{
   PyObject* args = PyTuple_New( 3 );
   PyTuple_SET_ITEM( args, 0, PyROOT::BindCppObject( selpad, Cppyy::GetScope( "TPad" ) ) );
   PyTuple_SET_ITEM( args, 1, PyROOT::BindCppObject( selected, Cppyy::GetScope( "TObject" ) ) );
   PyTuple_SET_ITEM( args, 2, PyInt_FromLong( event ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( Int_t event, Int_t x, Int_t y, TObject* selected )
{
   PyObject* args = PyTuple_New( 4 );
   PyTuple_SET_ITEM( args, 0, PyInt_FromLong( event ) );
   PyTuple_SET_ITEM( args, 1, PyInt_FromLong( x ) );
   PyTuple_SET_ITEM( args, 2, PyInt_FromLong( y ) );
   PyTuple_SET_ITEM( args, 3, PyROOT::BindCppObject( selected, Cppyy::GetScope( "TObject" ) ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( TVirtualPad* pad, TObject* obj, Int_t event )
{
   PyObject* args = PyTuple_New( 3 );
   PyTuple_SET_ITEM( args, 0, PyROOT::BindCppObject( pad, Cppyy::GetScope( "TVirtualPad" ) ) );
   PyTuple_SET_ITEM( args, 1, PyROOT::BindCppObject( obj, Cppyy::GetScope( "TObject" ) ) );
   PyTuple_SET_ITEM( args, 2, PyInt_FromLong( event ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( TGListTreeItem* item, TDNDData* data )
{
   PyObject* args = PyTuple_New( 2 );
   PyTuple_SET_ITEM( args, 0, PyROOT::BindCppObject( item, Cppyy::GetScope( "TGListTreeItem" ) ) );
   PyTuple_SET_ITEM( args, 1, PyROOT::BindCppObject( data, Cppyy::GetScope( "TDNDData" ) ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( const char* name, const TList* attr )
{
   PyObject* args = PyTuple_New( 2 );
   PyTuple_SET_ITEM( args, 0, PyBytes_FromString( name ) );
   PyTuple_SET_ITEM( args, 1, PyROOT::BindCppObject( (void*)attr, Cppyy::GetScope( "TList" ) ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

PyObject* TPyDispatcher::Dispatch( TSlave* slave, TProofProgressInfo* pi )
{
   PyObject* args = PyTuple_New( 2 );
   PyTuple_SET_ITEM( args, 0, PyROOT::BindCppObject( slave, Cppyy::GetScope( "TSlave" ) ) );
   PyTuple_SET_ITEM( args, 1, PyROOT::BindCppObject( pi, Cppyy::GetScope( "TProofProgressInfo" ) ) );

   PyObject* result = PyObject_CallObject( fCallable, args );
   Py_XDECREF( args );

   if ( ! result ) {
      PyErr_Print();
      return 0;
   }

   return result;
}
