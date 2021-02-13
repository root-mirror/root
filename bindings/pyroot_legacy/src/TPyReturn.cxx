// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, May 2004

// Bindings
#include "PyROOT.h"
#include "TPyReturn.h"
#include "ObjectProxy.h"

// ROOT
#include "TObject.h"
#include "TInterpreter.h"

// Standard
#include <stdexcept>


//______________________________________________________________________________
//                        Python expression eval result
//                        =============================
//
// Transport class for bringing objects from python (dynamically typed) to Cling
// (statically typed). It is best to immediately cast a TPyReturn to the real
// type, either implicitly (for builtin types) or explicitly (through a void*
// cast for pointers to ROOT objects).
//
// Examples:
//
//  root [0] TBrowser* b = (void*)TPython::Eval( "ROOT.TBrowser()" );
//  root [1] int i = TPython::Eval( "1+1" );
//  root [2] i
//  (int)2
//  root [3] double d = TPython::Eval( "1+3.1415" );
//  root [4] d
//  (double)4.14150000000000063e+00


//- data ---------------------------------------------------------------------
ClassImp(TPyReturn);


//- constructors/destructor --------------------------------------------------
TPyReturn::TPyReturn()
{
// Construct a TPyReturn object from Py_None.
   Py_INCREF( Py_None );
   fPyObject = Py_None;
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TPyReturn from a python object. The python object may represent
/// a ROOT object. Steals reference to given python object.

TPyReturn::TPyReturn( PyObject* pyobject )
{
   if ( ! pyobject ) {
      Py_INCREF( Py_None );
      fPyObject = Py_None;
   } else
      fPyObject = pyobject;             // steals reference
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor. Applies python object reference counting.

TPyReturn::TPyReturn( const TPyReturn& other )
{
   Py_INCREF( other.fPyObject );
   fPyObject = other.fPyObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator. Applies python object reference counting.

TPyReturn& TPyReturn::operator=( const TPyReturn& other )
{
   if ( this != &other ) {
      Py_INCREF( other.fPyObject );
      Py_DECREF( fPyObject );
      fPyObject = other.fPyObject;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Reference counting for the held python object is in effect.

TPyReturn::~TPyReturn()
{
   Py_DECREF( fPyObject );
}


//- public members -----------------------------------------------------------
TPyReturn::operator char*() const
{
// Cast python return value to C-style string (may fail).
   return (char*)((const char*)*this);
}

////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to C-style string (may fail).

TPyReturn::operator const char*() const
{
   if ( fPyObject == Py_None )     // for void returns
      return 0;

   const char* s = PyROOT_PyUnicode_AsString( fPyObject );
   if ( PyErr_Occurred() ) {
      PyErr_Print();
      return 0;
   }

   return s;
}

////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to C++ char (may fail).

TPyReturn::operator Char_t() const
{
   std::string s = operator const char*();
   if ( s.size() )
      return s[0];

   return '\0';
}

////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to C++ long (may fail).

TPyReturn::operator Long_t() const
{
   Long_t l = PyLong_AsLong( fPyObject );

   if ( PyErr_Occurred() )
      PyErr_Print();

   return l;
}

////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to C++ unsigned long (may fail).

TPyReturn::operator ULong_t() const
{
   ULong_t ul = PyLong_AsUnsignedLong( fPyObject );

   if ( PyErr_Occurred() )
      PyErr_Print();

   return ul;
}


////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to C++ double (may fail).

TPyReturn::operator Double_t() const
{
   Double_t d = PyFloat_AsDouble( fPyObject );

   if ( PyErr_Occurred() )
      PyErr_Print();

   return d;
}

////////////////////////////////////////////////////////////////////////////////
/// Cast python return value to ROOT object with dictionary (may fail; note that
/// you have to use the void* converter, as CINT will not call any other).

TPyReturn::operator void*() const
{
   if ( fPyObject == Py_None )
      return 0;

   if ( PyROOT::ObjectProxy_Check( fPyObject ) ) {
      ((PyROOT::ObjectProxy*)fPyObject)->Release();
      return ((PyROOT::ObjectProxy*)fPyObject)->GetObject();
   } else
      return fPyObject;                 // borrows reference
}

////////////////////////////////////////////////////////////////////////////////
/// Direct return of the held PyObject; note the new reference.

TPyReturn::operator PyObject*() const
{
   if ( fPyObject == Py_None )
      return 0;

   Py_INCREF( fPyObject );
   return fPyObject;
}
