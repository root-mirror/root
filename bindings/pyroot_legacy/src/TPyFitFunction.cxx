// Author: Wim Lavrijsen   November 2010

// Bindings
#include "PyROOT.h"
#include "TPyFitFunction.h"
#include "ObjectProxy.h"
#include "MethodProxy.h"
#include "TPyBufferFactory.h"

// Standard
#include <stdexcept>

//______________________________________________________________________________
//                       Python wrapper for Fit functions
//                       ================================
//


//- data ---------------------------------------------------------------------
ClassImp(TPyMultiGenFunction);
ClassImp(TPyMultiGradFunction);


//- helper functions ---------------------------------------------------------
static PyObject* GetOverriddenPyMethod( PyObject* pyself, const char* method )
{
// Retrieve an overriden method on pyself
   PyObject* pymethod = 0;

   if ( pyself && pyself != Py_None ) {
      pymethod = PyObject_GetAttrString( (PyObject*)pyself, const_cast< char* >( method ) );
      if ( ! PyROOT::MethodProxy_CheckExact( pymethod ) )
         return pymethod;

      Py_XDECREF( pymethod );
      pymethod = 0;
   }

   return pymethod;
}

static PyObject* DispatchCall( PyObject* pyself, const char* method, PyObject* pymethod = NULL,
   PyObject* arg1 = NULL, PyObject* arg2 = NULL, PyObject* arg3 = NULL )
{
// Forward <method> to python (need to refactor this with TPySelector).
   PyObject* result = 0;

// get the named method and check for python side overload by not accepting the
// binding's methodproxy
   if ( ! pymethod )
      pymethod = GetOverriddenPyMethod( pyself, method );

   if ( pymethod ) {
      result = PyObject_CallFunctionObjArgs( pymethod, arg1, arg2, arg3, NULL );
   } else {
   // means the method has not been overridden ... simply accept its not there
      result = 0;
      PyErr_Format( PyExc_AttributeError,
         "method %s needs implementing in derived class", const_cast< char* >( method ) );
   }

   Py_XDECREF( pymethod );

   return result;
}


//- constructors/destructor --------------------------------------------------
TPyMultiGenFunction::TPyMultiGenFunction( PyObject* self ) : fPySelf( 0 )
{
// Construct a TPyMultiGenFunction derived with <self> as the underlying
   if ( self ) {
   // steal reference as this is us, as seen from python
      fPySelf = self;
   } else {
      Py_INCREF( Py_None );        // using None allows clearer diagnostics
      fPySelf = Py_None;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Only deref if still holding on to Py_None (circular otherwise).

TPyMultiGenFunction::~TPyMultiGenFunction()
{
   if ( fPySelf == Py_None ) {
      Py_DECREF( fPySelf );
   }
}


//- public functions ---------------------------------------------------------
unsigned int TPyMultiGenFunction::NDim() const
{
// Simply forward the call to python self.
   PyObject* pyresult = DispatchCall( fPySelf, "NDim" );

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyMultiGenFunction::NDim" );
   }

   unsigned int cppresult = (unsigned int)PyLong_AsLong( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

double TPyMultiGenFunction::DoEval( const double* x ) const
{
   PyObject* xbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)x );
   PyObject* pyresult = DispatchCall( fPySelf, "DoEval", NULL, xbuf );
   Py_DECREF( xbuf );

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyMultiGenFunction::DoEval" );
   }

   double cppresult = (double)PyFloat_AsDouble( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}



//- constructors/destructor --------------------------------------------------
TPyMultiGradFunction::TPyMultiGradFunction( PyObject* self )
{
// Construct a TPyMultiGradFunction derived with <self> as the underlying
   if ( self ) {
   // steal reference as this is us, as seen from python
      fPySelf = self;
   } else {
      Py_INCREF( Py_None );        // using None allows clearer diagnostics
      fPySelf = Py_None;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Only deref if still holding on to Py_None (circular otherwise).

TPyMultiGradFunction::~TPyMultiGradFunction()
{
   if ( fPySelf == Py_None ) {
      Py_DECREF( fPySelf );
   }
}


//- public functions ---------------------------------------------------------
unsigned int TPyMultiGradFunction::NDim() const
{
// Simply forward the call to python self.
   PyObject* pyresult = DispatchCall( fPySelf, "NDim" );

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyMultiGradFunction::NDim" );
   }

   unsigned int cppresult = (unsigned int)PyLong_AsLong( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

double TPyMultiGradFunction::DoEval( const double* x ) const
{
   PyObject* xbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)x );
   PyObject* pyresult = DispatchCall( fPySelf, "DoEval", NULL, xbuf );
   Py_DECREF( xbuf );

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyMultiGradFunction::DoEval" );
   }

   double cppresult = (double)PyFloat_AsDouble( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

void TPyMultiGradFunction::Gradient( const double* x, double* grad ) const {
   PyObject* pymethod = GetOverriddenPyMethod( fPySelf, "Gradient" );

   if ( pymethod ) {
      PyObject* xbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)x );
      PyObject* gbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)grad );
      PyObject* pyresult = DispatchCall( fPySelf, "Gradient", pymethod, xbuf, gbuf );
      Py_DECREF( gbuf );
      Py_DECREF( xbuf );

      if ( ! pyresult ) {
         PyErr_Print();
         throw std::runtime_error( "Failure in TPyMultiGradFunction::Gradient" );
      }

      Py_DECREF( pyresult );

   } else
      return ROOT::Math::IMultiGradFunction::Gradient( x, grad );
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

void TPyMultiGradFunction::FdF( const double* x, double& f, double* df ) const
{
   PyObject* pymethod = GetOverriddenPyMethod( fPySelf, "FdF" );

   if ( pymethod ) {
      PyObject* xbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)x );
      PyObject* pyf = PyList_New( 1 );
      PyList_SetItem( pyf, 0, PyFloat_FromDouble( f ) );
      PyObject* dfbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)df );

      PyObject* pyresult = DispatchCall( fPySelf, "FdF", pymethod, xbuf, pyf, dfbuf );
      f = PyFloat_AsDouble( PyList_GetItem( pyf, 0 ) );

      Py_DECREF( dfbuf );
      Py_DECREF( pyf );
      Py_DECREF( xbuf );

      if ( ! pyresult ) {
         PyErr_Print();
         throw std::runtime_error( "Failure in TPyMultiGradFunction::FdF" );
      }

      Py_DECREF( pyresult );

   } else
      return ROOT::Math::IMultiGradFunction::FdF( x, f, df );
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

double TPyMultiGradFunction::DoDerivative( const double * x, unsigned int icoord ) const
{
   PyObject* xbuf = PyROOT::TPyBufferFactory::Instance()->PyBuffer_FromMemory( (Double_t*)x );
   PyObject* pycoord = PyLong_FromLong( icoord );

   PyObject* pyresult = DispatchCall( fPySelf, "DoDerivative", NULL, xbuf, pycoord );
   Py_DECREF( pycoord );
   Py_DECREF( xbuf );

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyMultiGradFunction::DoDerivative" );
   }

   double cppresult = (double)PyFloat_AsDouble( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}

