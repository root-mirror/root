// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

// Bindings
#include "PyROOT.h"
#include "TMethodHolder.h"
#include "Converters.h"
#include "Executors.h"
#include "ObjectProxy.h"
#include "RootWrapper.h"
#include "TPyException.h"
#include "Utility.h"

// ROOT
#include "TClass.h"           // for exception types (to move to Cppyy.cxx)
#include "TException.h"       // for TRY ... CATCH
#include "TVirtualMutex.h"    // for R__LOCKGUARD2
#include "TClassEdit.h"       // demangler
#include "TInterpreter.h"     // for Interpreter exceptions

// Standard
#include <assert.h>
#include <string.h>
#include <exception>
#include <sstream>
#include <string>
#include <typeinfo>
#include <memory>

//- data and local helpers ---------------------------------------------------
namespace PyROOT {
   R__EXTERN PyObject* gRootModule;
}


//- private helpers ----------------------------------------------------------
inline void PyROOT::TMethodHolder::Copy_( const TMethodHolder& /* other */ )
{
// fScope and fMethod handled separately

// do not copy caches
   fExecutor   = 0;
   fArgsRequired = -1;

// being uninitialized will trigger setting up caches as appropriate
   fIsInitialized  = kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// destroy executor and argument converters

inline void PyROOT::TMethodHolder::Destroy_() const
{
   delete fExecutor;

   for ( int i = 0; i < (int)fConverters.size(); ++i )
      delete fConverters[ i ];
}

////////////////////////////////////////////////////////////////////////////////
/// Helper code to prevent some duplication; this is called from CallSafe() as well
/// as directly from TMethodHolder::Execute in fast mode.

inline PyObject* PyROOT::TMethodHolder::CallFast( void* self, ptrdiff_t offset, TCallContext* ctxt )
{
   PyObject* result = nullptr;

   try {       // C++ try block
      result = fExecutor->Execute( fMethod, (Cppyy::TCppObject_t)((Long_t)self + offset), ctxt );
   } catch ( TPyException& ) {
      result = nullptr;           // error already set
   } catch ( std::exception& e ) {
      if (gInterpreter->DiagnoseIfInterpreterException(e)) {
         return result;
      }
   // map user exceptions .. this needs to move to Cppyy.cxx
      TClass* cl = TClass::GetClass( typeid(e) );

      PyObject* pyUserExcepts = PyObject_GetAttrString( gRootModule, "UserExceptions" );
      std::string exception_type;
      if (cl) exception_type = cl->GetName();
      else {
         int errorCode;
         std::unique_ptr<char[]> demangled(TClassEdit::DemangleTypeIdName(typeid(e),errorCode));
         if (errorCode) exception_type = typeid(e).name();
         else exception_type = demangled.get();
      }
      PyObject* pyexc = PyDict_GetItemString( pyUserExcepts, exception_type.c_str() );
      if ( !pyexc ) {
         PyErr_Clear();
         pyexc = PyDict_GetItemString( pyUserExcepts, ("std::"+exception_type).c_str() );
      }
      if ( !pyexc ) {
         PyErr_Clear();
         pyexc = PyDict_GetItemString( pyUserExcepts, ("ROOT::"+exception_type).c_str() );
      }
      Py_DECREF( pyUserExcepts );

      if ( pyexc ) {
         PyErr_Format( pyexc, "%s", e.what() );
      } else {
         PyErr_Format( PyExc_Exception, "%s (C++ exception of type %s)", e.what(), exception_type.c_str() );
      }
      result = nullptr;
   } catch ( ... ) {
      PyErr_SetString( PyExc_Exception, "unhandled, unknown C++ exception" );
      result = nullptr;
   }
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// Helper code to prevent some code duplication; this code embeds a ROOT "try/catch"
/// block that saves the stack for restoration in case of an otherwise fatal signal.

inline PyObject* PyROOT::TMethodHolder::CallSafe( void* self, ptrdiff_t offset, TCallContext* ctxt )
{
   PyObject* result = 0;

   TRY {       // ROOT "try block"
      result = CallFast( self, offset, ctxt );
   } CATCH( excode ) {
      PyErr_SetString( PyExc_SystemError, "problem in C++; program state has been reset" );
      result = 0;
      Throw( excode );
   } ENDTRY;

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// build buffers for argument dispatching

Bool_t PyROOT::TMethodHolder::InitConverters_()
{
   const size_t nArgs = Cppyy::GetMethodNumArgs( fMethod );
   fConverters.resize( nArgs );

// setup the dispatch cache
   for ( size_t iarg = 0; iarg < nArgs; ++iarg ) {
      const std::string& fullType = Cppyy::GetMethodArgType( fMethod, iarg );

   // CLING WORKAROUND -- std::string can not use kExactMatch as that will
   //                     fail, but if no exact match is used, the const-ref
   //                     std::string arguments will mask the const char* ones,
   //                     even though the extra default arguments differ
      if ( Cppyy::GetFinalName( fScope ) == "string" && Cppyy::GetMethodName( fMethod ) == "string" &&
           // Note with the improve naming normalization we should see only
           // the spelling "const string&" (but soon it will be "const std::string&")
           ( fullType == "const std::string&" || fullType == "const std::string &"
             || fullType == "const string&" || fullType == "const string &" ) ) {
         fConverters[ iarg ] = new TStrictCppObjectConverter( Cppyy::GetScope( "string" ), kFALSE ); // TODO: this is sooo wrong
   // -- CLING WORKAROUND
      } else
         fConverters[ iarg ] = CreateConverter( fullType );

      if ( ! fConverters[ iarg ] ) {
         PyErr_Format( PyExc_TypeError, "argument type %s not handled", fullType.c_str() );
         return kFALSE;
      }

   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// install executor conform to the return type

Bool_t PyROOT::TMethodHolder::InitExecutor_( TExecutor*& executor, TCallContext* ctxt )
{
   executor = CreateExecutor( (Bool_t)fMethod == true ?
      Cppyy::ResolveName( Cppyy::GetMethodResultType( fMethod ) ) : Cppyy::GetScopedFinalName( fScope ),
      ctxt ? ManagesSmartPtr( ctxt ) : kFALSE );

   if ( ! executor )
      return kFALSE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// built a signature representation (used for doc strings)

std::string PyROOT::TMethodHolder::GetSignatureString()
{
   std::stringstream sig; sig << "(";
   Int_t ifirst = 0;
   const size_t nArgs = Cppyy::GetMethodNumArgs( fMethod );
   for ( size_t iarg = 0; iarg < nArgs; ++iarg ) {
      if ( ifirst ) sig << ", ";

      sig << Cppyy::GetMethodArgType( fMethod, iarg );

      const std::string& parname = Cppyy::GetMethodArgName( fMethod, iarg );
      if ( ! parname.empty() )
         sig << " " << parname;

      const std::string& defvalue = Cppyy::GetMethodArgDefault( fMethod, iarg );
      if ( ! defvalue.empty() )
         sig << " = " << defvalue;
      ifirst++;
   }
   sig << ")";
   return sig.str();
}

////////////////////////////////////////////////////////////////////////////////
/// helper to report errors in a consistent format (derefs msg)

void PyROOT::TMethodHolder::SetPyError_( PyObject* msg )
{
   PyObject *etype, *evalue, *etrace;
   PyErr_Fetch( &etype, &evalue, &etrace );

   std::string details = "";
   if ( evalue ) {
      PyObject* descr = PyObject_Str( evalue );
      if ( descr ) {
         details = PyROOT_PyUnicode_AsString( descr );
         Py_DECREF( descr );
      }
   }

   Py_XDECREF( evalue ); Py_XDECREF( etrace );

   PyObject* doc = GetDocString();
   PyObject* errtype = etype ? etype : PyExc_TypeError;
   if ( details.empty() ) {
      PyErr_Format( errtype, "%s =>\n    %s", PyROOT_PyUnicode_AsString( doc ),
         msg ? PyROOT_PyUnicode_AsString( msg ) : ""  );
   } else if ( msg ) {
      PyErr_Format( errtype, "%s =>\n    %s (%s)",
         PyROOT_PyUnicode_AsString( doc ), PyROOT_PyUnicode_AsString( msg ), details.c_str() );
   } else {
      PyErr_Format( errtype, "%s =>\n    %s",
         PyROOT_PyUnicode_AsString( doc ), details.c_str() );
   }

   Py_XDECREF( etype );
   Py_DECREF( doc );
   Py_XDECREF( msg );
}

//- constructors and destructor ----------------------------------------------
PyROOT::TMethodHolder::TMethodHolder(
      Cppyy::TCppScope_t scope, Cppyy::TCppMethod_t method ) :
   fMethod( method ), fScope( scope ), fExecutor( nullptr ), fArgsRequired( -1 ),
   fIsInitialized( kFALSE )
{
   // empty
}

PyROOT::TMethodHolder::TMethodHolder( const TMethodHolder& other ) :
      PyCallable( other ), fMethod( other.fMethod ), fScope( other.fScope )
{
// copy constructor
   Copy_( other );
}

////////////////////////////////////////////////////////////////////////////////
/// assignment operator

PyROOT::TMethodHolder& PyROOT::TMethodHolder::operator=( const TMethodHolder& other )
{
   if ( this != &other ) {
      Destroy_();
      Copy_( other );
      fScope  = other.fScope;
      fMethod = other.fMethod;
   }

   return *this;
}

////////////////////////////////////////////////////////////////////////////////
/// destructor

PyROOT::TMethodHolder::~TMethodHolder()
{
   Destroy_();
}


//- public members -----------------------------------------------------------
PyObject* PyROOT::TMethodHolder::GetPrototype()
{
// construct python string from the method's prototype
   return PyROOT_PyUnicode_FromFormat( "%s%s %s::%s%s",
      ( Cppyy::IsStaticMethod( fMethod ) ? "static " : "" ),
      Cppyy::GetMethodResultType( fMethod ).c_str(),
      Cppyy::GetFinalName( fScope ).c_str(), Cppyy::GetMethodName( fMethod ).c_str(),
      GetSignatureString().c_str() );
}

////////////////////////////////////////////////////////////////////////////////
/// Method priorities exist (in lieu of true overloading) there to prevent
/// void* or <unknown>* from usurping otherwise valid calls. TODO: extend this
/// to favour classes that are not bases.

Int_t PyROOT::TMethodHolder::GetPriority()
{
   Int_t priority = 0;

   const size_t nArgs = Cppyy::GetMethodNumArgs( fMethod );
   for ( size_t iarg = 0; iarg < nArgs; ++iarg ) {
      const std::string aname = Cppyy::GetMethodArgType( fMethod, iarg );

   // the following numbers are made up and may cause problems in specific
   // situations: use <obj>.<meth>.disp() for choice of exact dispatch
      if ( Cppyy::IsBuiltin( aname ) ) {
      // happens for builtin types (and namespaces, but those can never be an
      // argument), NOT for unknown classes as that concept no longer exists
         if ( strstr( aname.c_str(), "void*" ) )
         // TODO: figure out in general all void* converters
            priority -= 10000;     // void*/void** shouldn't be too greedy
         else if ( strstr( aname.c_str(), "float" ) )
            priority -= 1000;      // double preferred (no float in python)
         else if ( strstr( aname.c_str(), "long double" ) )
            priority -= 100;       // id, but better than float
         else if ( strstr( aname.c_str(), "double" ) )
            priority -= 10;        // char, int, long can't convert float,
                                   // but vv. works, so prefer the int types
         else if ( strstr( aname.c_str(), "bool" ) )
            priority += 1;         // bool over int (does accept 1 and 0)

      } else if ( !aname.empty() && !Cppyy::IsComplete( aname ) ) {
      // class is known, but no dictionary available, 2 more cases: * and &
         if ( aname[ aname.size() - 1 ] == '&' )
            priority -= 1000000;
         else
            priority -= 100000; // prefer pointer passing over reference

      } else {
      // resolve a few special cases (these are valid & known, but are lined up
      // with derived classes in there interface that should have preference
         if ( aname == "IBaseFunctionMultiDim")
            priority -= 1;
         else if ( aname == "RooAbsReal" )
            priority -= 1;
      }

   }

// add a small penalty to prefer non-const methods over const ones for
// getitem/setitem
   if ( Cppyy::IsConstMethod( fMethod ) && Cppyy::GetMethodName( fMethod ) == "operator[]" )
       priority -= 1;

// another special case for RooFit, as it is inconsistent on base <-> derived
   if ( Cppyy::GetMethodName( fMethod ) == "import" &&
        nArgs != 0 && Cppyy::GetMethodArgType( fMethod, 0 ) == "TObject&" )
      priority -= 1000;

   return priority;
}

////////////////////////////////////////////////////////////////////////////////

Int_t PyROOT::TMethodHolder::GetMaxArgs()
{
   return Cppyy::GetMethodNumArgs( fMethod );
}

////////////////////////////////////////////////////////////////////////////////
/// Build a tuple of the argument types/names.

PyObject* PyROOT::TMethodHolder::GetCoVarNames()
{
   int co_argcount = (int)GetMaxArgs() /* +1 for self */;

// TODO: static methods need no 'self' (but is harmless otherwise)

   PyObject* co_varnames = PyTuple_New( co_argcount + 1 /* self */ );
   PyTuple_SET_ITEM( co_varnames, 0, PyROOT_PyUnicode_FromString( "self" ) );
   for ( int iarg = 0; iarg < co_argcount; ++iarg ) {
      std::string argrep = Cppyy::GetMethodArgType( fMethod, iarg );
      const std::string& parname = Cppyy::GetMethodArgName( fMethod, iarg );
      if ( ! parname.empty() ) {
         argrep += " ";
         argrep += parname;
      }

      PyObject* pyspec = PyROOT_PyUnicode_FromString( argrep.c_str() );

      PyTuple_SET_ITEM( co_varnames, iarg + 1, pyspec );
   }

   return co_varnames;
}

////////////////////////////////////////////////////////////////////////////////
/// get the default value (if any) of argument iarg of this method

PyObject* PyROOT::TMethodHolder::GetArgDefault( Int_t iarg )
{
   if ( iarg >= (int)GetMaxArgs() )
      return 0;

   const std::string& defvalue = Cppyy::GetMethodArgDefault( fMethod, iarg );
   if ( ! defvalue.empty() ) {

   // attempt to evaluate the string representation (will work for all builtin types)
      PyObject* pyval = (PyObject*)PyRun_String(
          (char*)defvalue.c_str(), Py_eval_input, gRootModule, gRootModule );
      if ( ! pyval && PyErr_Occurred() ) {
         PyErr_Clear();
         return PyROOT_PyUnicode_FromString( defvalue.c_str() );
      }

      return pyval;
   }

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// Get or build the scope of this method.

PyObject* PyROOT::TMethodHolder::GetScopeProxy()
{
   return CreateScopeProxy( fScope );
}

////////////////////////////////////////////////////////////////////////////////
/// done if cache is already setup

Bool_t PyROOT::TMethodHolder::Initialize( TCallContext* ctxt )
{
   if ( fIsInitialized == kTRUE )
      return kTRUE;

   if ( ! InitConverters_() )
      return kFALSE;

   if ( ! InitExecutor_( fExecutor, ctxt ) )
      return kFALSE;

// minimum number of arguments when calling
   fArgsRequired = (Bool_t)fMethod == true ? Cppyy::GetMethodReqArgs( fMethod ) : 0;

// init done
   fIsInitialized = kTRUE;

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// verify existence of self, return if ok

PyObject* PyROOT::TMethodHolder::PreProcessArgs( ObjectProxy*& self, PyObject* args, PyObject* )
{
   if ( self != 0 ) {
      Py_INCREF( args );
      return args;
   }

// otherwise, check for a suitable 'self' in args and update accordingly
   if ( PyTuple_GET_SIZE( args ) != 0 ) {
      ObjectProxy* pyobj = (ObjectProxy*)PyTuple_GET_ITEM( args, 0 );

   // demand PyROOT object, and an argument that may match down the road
      if ( ObjectProxy_Check( pyobj ) &&
           ( fScope == Cppyy::gGlobalScope ||               // free global
           ( pyobj->ObjectIsA() == 0 )     ||               // null pointer or ctor call
           ( Cppyy::IsSubtype( pyobj->ObjectIsA(), fScope ) ) ) // matching types
         ) {
      // reset self
         self = pyobj;
         Py_INCREF( self );        // corresponding Py_DECREF is in MethodProxy

      // offset args by 1 (new ref)
         return PyTuple_GetSlice( args, 1, PyTuple_GET_SIZE( args ) );
      }
   }

// no self, set error and lament
   SetPyError_( PyROOT_PyUnicode_FromFormat(
      "unbound method %s::%s must be called with a %s instance as first argument",
      Cppyy::GetFinalName( fScope ).c_str(), Cppyy::GetMethodName( fMethod ).c_str(),
      Cppyy::GetFinalName( fScope ).c_str() ) );
   return 0;
}

////////////////////////////////////////////////////////////////////////////////

Bool_t PyROOT::TMethodHolder::ConvertAndSetArgs( PyObject* args, TCallContext* ctxt )
{
   int argc = PyTuple_GET_SIZE( args );
   int argMax = fConverters.size();

// argc must be between min and max number of arguments
   if ( argc < fArgsRequired ) {
      SetPyError_( PyROOT_PyUnicode_FromFormat(
         "takes at least %d arguments (%d given)", fArgsRequired, argc ) );
      return kFALSE;
   } else if ( argMax < argc ) {
      SetPyError_( PyROOT_PyUnicode_FromFormat(
         "takes at most %d arguments (%d given)", argMax, argc ) );
      return kFALSE;
   }

// convert the arguments to the method call array
   ctxt->fArgs.resize( argc );
   for ( int i = 0; i < argc; ++i ) {
      if ( ! fConverters[ i ]->SetArg(
              PyTuple_GET_ITEM( args, i ), ctxt->fArgs[i], ctxt ) ) {
         SetPyError_( PyROOT_PyUnicode_FromFormat( "could not convert argument %d", i+1 ) );
         return kFALSE;
      }
   }

   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// call the interface method

PyObject* PyROOT::TMethodHolder::Execute( void* self, ptrdiff_t offset, TCallContext* ctxt )
{
   PyObject* result = 0;

   if ( TCallContext::sSignalPolicy == TCallContext::kFast ) {
   // bypasses ROOT try block (i.e. segfaults will abort)
      result = CallFast( self, offset, ctxt );
   } else {
   // at the cost of ~10% performance, don't abort the interpreter on any signal
      result = CallSafe( self, offset, ctxt );
   }

   if ( result && Utility::PyErr_Occurred_WithGIL() ) {
   // can happen in the case of a CINT error: trigger exception processing
      Py_DECREF( result );
      result = 0;
   } else if ( ! result && PyErr_Occurred() )
      SetPyError_( 0 );

   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// preliminary check in case keywords are accidently used (they are ignored otherwise)

PyObject* PyROOT::TMethodHolder::Call(
      ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt )
{
   if ( kwds != 0 && PyDict_Size( kwds ) ) {
      PyErr_SetString( PyExc_TypeError, "keyword arguments are not yet supported" );
      return 0;
   }

// setup as necessary
   if ( ! Initialize( ctxt ) )
      return 0;                              // important: 0, not Py_None

// fetch self, verify, and put the arguments in usable order
   if ( ! ( args = PreProcessArgs( self, args, kwds ) ) )
      return 0;

// translate the arguments
   Bool_t bConvertOk = ConvertAndSetArgs( args, ctxt );
   Py_DECREF( args );

   if ( bConvertOk == kFALSE )
      return 0;                              // important: 0, not Py_None

// get the ROOT object that this object proxy is a handle for
   void* object = self->GetObject();

// validity check that should not fail
   if ( ! object ) {
      PyErr_SetString( PyExc_ReferenceError, "attempt to access a null-pointer" );
      return 0;
   }

// get its class
   Cppyy::TCppType_t derived = self->ObjectIsA();

// calculate offset (the method expects 'this' to be an object of fScope)
   ptrdiff_t offset = 0;
   if ( derived && derived != fScope )
      offset = Cppyy::GetBaseOffset( derived, fScope, object, 1 /* up-cast */ );

// actual call; recycle self instead of returning new object for same address objects
   ObjectProxy* pyobj = (ObjectProxy*)Execute( object, offset, ctxt );

   if ( ObjectProxy_Check( pyobj ) &&
        derived && pyobj->ObjectIsA() == derived &&
        pyobj->GetObject() == object ) {
      Py_INCREF( (PyObject*)self );
      Py_DECREF( pyobj );
      return (PyObject*)self;
   }

   return (PyObject*)pyobj;
}

//- protected members --------------------------------------------------------
PyObject* PyROOT::TMethodHolder::GetSignature()
{
// construct python string from the method's signature
   return PyROOT_PyUnicode_FromString( GetSignatureString().c_str() );
}

////////////////////////////////////////////////////////////////////////////////

std::string PyROOT::TMethodHolder::GetReturnTypeName()
{
   return Cppyy::GetMethodResultType( fMethod );
}
