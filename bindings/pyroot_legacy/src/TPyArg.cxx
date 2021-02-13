// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Aug 2013

// Bindings
#include "PyROOT.h"
#include "TPyArg.h"


//______________________________________________________________________________
//                        Generic wrapper for arguments
//                        =============================
//
// Transport class for bringing C++ values and objects from Cling to Python. It
// provides, from the selected constructor, the proper conversion to a PyObject.
// In principle, there should be no need to use this class directly: it relies
// on implicit conversions.


//- data ---------------------------------------------------------------------
ClassImp(TPyArg);

//- constructor dispatcher ---------------------------------------------------
void TPyArg::CallConstructor( PyObject*& pyself, PyObject* pyclass, const std::vector<TPyArg>& args )
{
   int nArgs = args.size();
   PyObject* pyargs = PyTuple_New( nArgs );
   for ( int i = 0; i < nArgs; ++i )
      PyTuple_SET_ITEM( pyargs, i, (PyObject*)args[i] );
   pyself = PyObject_Call( pyclass, pyargs, NULL );
   Py_DECREF( pyargs );
}

////////////////////////////////////////////////////////////////////////////////
void CallConstructor( PyObject*& pyself, PyObject* pyclass )
{
   PyObject* pyargs = PyTuple_New( 0 );
   pyself = PyObject_Call( pyclass, pyargs, NULL );
    Py_DECREF( pyargs );
}

//- generic dispatcher -------------------------------------------------------
PyObject* TPyArg::CallMethod( PyObject* pymeth, const std::vector<TPyArg>& args )
{
   int nArgs = args.size();
   PyObject* pyargs = PyTuple_New( nArgs );
   for ( int i = 0; i < nArgs; ++i )
      PyTuple_SET_ITEM( pyargs, i, (PyObject*)args[i] );
   PyObject* result = PyObject_Call( pymeth, pyargs, NULL );
   Py_DECREF( pyargs );
   return result;
}

//- denstructor dispatcher ----------------------------------------------------
void TPyArg::CallDestructor( PyObject*& pyself, PyObject*, const std::vector<TPyArg>& )
{
   Py_XDECREF( pyself );      // calls actual dtor if ref-count down to 0
}

////////////////////////////////////////////////////////////////////////////////
void TPyArg::CallDestructor( PyObject*& pyself )
{
   Py_XDECREF( pyself );
}

//- constructors/destructor --------------------------------------------------
TPyArg::TPyArg( PyObject* pyobject )
{
// Construct a TPyArg from a python object.
   Py_XINCREF( pyobject );
   fPyObject = pyobject;
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TPyArg from an integer value.

TPyArg::TPyArg( Int_t value )
{
   fPyObject = PyInt_FromLong( value );
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TPyArg from an integer value.

TPyArg::TPyArg( Long_t value )
{
   fPyObject = PyLong_FromLong( value );
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TPyArg from a double value.

TPyArg::TPyArg( Double_t value )
{
   fPyObject = PyFloat_FromDouble( value );
}

////////////////////////////////////////////////////////////////////////////////
/// Construct a TPyArg from a C-string.

TPyArg::TPyArg( const char* value )
{
   fPyObject = PyROOT_PyUnicode_FromString( value );
}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor.

TPyArg::TPyArg( const TPyArg& s )
{
   Py_XINCREF( s.fPyObject );
   fPyObject = s.fPyObject;
}

////////////////////////////////////////////////////////////////////////////////
/// Assignment operator.

TPyArg& TPyArg::operator=( const TPyArg& s )
{
   if ( &s != this ) {
      Py_XINCREF( s.fPyObject );
      fPyObject = s.fPyObject;
   }
   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// Done with held PyObject.

TPyArg::~TPyArg()
{
   Py_XDECREF( fPyObject );
   fPyObject = NULL;
}

//- public members -----------------------------------------------------------
TPyArg::operator PyObject*() const
{
// Extract the python object.
   Py_XINCREF( fPyObject );
   return fPyObject;
}
