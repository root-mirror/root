// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

// Bindings
#include "PyROOT.h"
#include "structmember.h"    // from Python
#if PY_VERSION_HEX >= 0x02050000
#include "code.h"            // from Python
#else
#include "compile.h"         // from Python
#endif
#ifndef CO_NOFREE
// python2.2 does not have CO_NOFREE defined
#define CO_NOFREE       0x0040
#endif
#include "MethodProxy.h"
#include "ObjectProxy.h"
#include "TCallContext.h"
#include "TPyException.h"
#include "PyStrings.h"

// Standard
#include <algorithm>
#include <vector>


namespace PyROOT {

// TODO: only used here, but may be better off integrated with Pythonize.cxx callbacks
   class TPythonCallback : public PyCallable {
   public:
      PyObject* fCallable;

      TPythonCallback( PyObject* callable ):
         fCallable(nullptr)
      {
         if ( !PyCallable_Check( callable ) ) {
            PyErr_SetString(PyExc_TypeError, "parameter must be callable");
            return;
         }
         fCallable = callable;
         Py_INCREF( fCallable );
      }

      virtual ~TPythonCallback() {
         Py_DECREF( fCallable );
         fCallable = 0;
      }

      virtual PyObject* GetSignature() { return PyROOT_PyUnicode_FromString( "*args, **kwargs" ); } ;
      virtual PyObject* GetPrototype() { return PyROOT_PyUnicode_FromString( "<callback>" ); } ;
      virtual PyObject* GetDocString() {
         if ( PyObject_HasAttrString( fCallable, "__doc__" )) {
            return PyObject_GetAttrString( fCallable, "__doc__" );
         } else {
            return GetPrototype();
         }
      }

      virtual Int_t GetPriority() { return 100; };

      virtual Int_t GetMaxArgs() { return 100; };
      virtual PyObject* GetCoVarNames() { // TODO: pick these up from the callable
         Py_INCREF( Py_None );
         return Py_None;
      }
      virtual PyObject* GetArgDefault( Int_t /* iarg */ ) { // TODO: pick these up from the callable
         Py_INCREF( Py_None );
         return Py_None;
      }

      virtual PyObject* GetScopeProxy() { // should this be the module ??
         Py_INCREF( Py_None );
         return Py_None;
      }

      virtual PyCallable* Clone() { return new TPythonCallback( *this ); }

      virtual PyObject* Call(
            ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* /* ctxt = 0 */ ) {

         PyObject* newArgs = nullptr;
         if ( self ) {
            Py_ssize_t nargs = PyTuple_Size( args );
            newArgs = PyTuple_New( nargs+1 );
            Py_INCREF( self );
            PyTuple_SET_ITEM( newArgs, 0, (PyObject*)self );
            for ( Py_ssize_t iarg = 0; iarg < nargs; ++iarg ) {
               PyObject* pyarg = PyTuple_GET_ITEM( args, iarg );
               Py_INCREF( pyarg );
               PyTuple_SET_ITEM( newArgs, iarg+1, pyarg );
            }
         } else {
            Py_INCREF( args );
            newArgs = args;
         }
         return PyObject_Call( fCallable, newArgs, kwds );
      }
  };

namespace {

// helper to test whether a method is used in a pseudo-function modus
   Bool_t inline IsPseudoFunc( MethodProxy* pymeth )
   {
      return (void*)pymeth == (void*)pymeth->fSelf;
   }

// helper for collecting/maintaining exception data in overload dispatch
   struct PyError_t {
      PyError_t() { fType = fValue = fTrace = 0; }

      static void Clear( PyError_t& e )
      {
      // Remove exception information.
         Py_XDECREF( e.fType ); Py_XDECREF( e.fValue ); Py_XDECREF( e.fTrace );
         e.fType = e.fValue = e.fTrace = 0;
      }

      PyObject *fType, *fValue, *fTrace;
   };

// helper to hash tuple (using tuple hash would cause self-tailing loops)
   inline Long_t HashSignature( PyObject* args )
   {
   // Build a hash from the types of the given python function arguments.
      ULong_t hash = 0;

      Int_t nargs = PyTuple_GET_SIZE( args );
      for ( Int_t i = 0; i < nargs; ++i ) {
         hash += (ULong_t) Py_TYPE( PyTuple_GET_ITEM( args, i ) );
         hash += (hash << 10); hash ^= (hash >> 6);
      }

      hash += (hash << 3); hash ^= (hash >> 11); hash += (hash << 15);

      return hash;
   }

// helper to sort on method priority
   int PriorityCmp( PyCallable* left, PyCallable* right )
   {
      return left->GetPriority() > right->GetPriority();
   }

// return helper
   inline void ResetCallState( ObjectProxy*& selfnew, ObjectProxy* selfold, Bool_t clear ) {
      if ( selfnew != selfold ) {
         Py_XDECREF( selfnew );
         selfnew = selfold;
      }

      if ( clear )
         PyErr_Clear();
   }

// helper to factor out return logic of mp_call
   inline PyObject* HandleReturn( MethodProxy* pymeth, ObjectProxy* oldSelf, PyObject* result ) {

   // special case for python exceptions, propagated through C++ layer
      if ( result ) {

      // if this method creates new objects, always take ownership
         if ( IsCreator( pymeth->fMethodInfo->fFlags ) ) {

         // either be a constructor with a fresh object proxy self ...
            if ( IsConstructor( pymeth->fMethodInfo->fFlags ) ) {
               if ( pymeth->fSelf )
                  pymeth->fSelf->HoldOn();
            }

         // ... or be a method with an object proxy return value
            else if ( ObjectProxy_Check( result ) )
               ((ObjectProxy*)result)->HoldOn();
         }

      // if this new object falls inside self, make sure its lifetime is proper
         if ( ObjectProxy_Check( pymeth->fSelf ) && ObjectProxy_Check( result ) ) {
            Long_t ptrdiff = (Long_t)((ObjectProxy*)result)->GetObject() - (Long_t)pymeth->fSelf->GetObject();
            if ( 0 <= ptrdiff && ptrdiff < (Long_t)Cppyy::SizeOf( pymeth->fSelf->ObjectIsA() ) ) {
               if ( PyObject_SetAttr( result, PyStrings::gLifeLine, (PyObject*)pymeth->fSelf ) == -1 )
                  PyErr_Clear();     // ignored
            }
         }
      }

   // reset self as necessary to allow re-use of the MethodProxy
      ResetCallState( pymeth->fSelf, oldSelf, kFALSE );

      return result;
   }


//= PyROOT method proxy object behaviour =====================================
   PyObject* mp_name( MethodProxy* pymeth, void* )
   {
      return PyROOT_PyUnicode_FromString( pymeth->GetName().c_str() );
   }

////////////////////////////////////////////////////////////////////////////////

   PyObject* mp_module( MethodProxy* /* pymeth */, void* )
   {
      Py_INCREF( PyStrings::gROOTns );
      return PyStrings::gROOTns;
   }

////////////////////////////////////////////////////////////////////////////////
/// Build python document string ('__doc__') from all C++-side overloads.

   PyObject* mp_doc( MethodProxy* pymeth, void* )
   {
      MethodProxy::Methods_t& methods = pymeth->fMethodInfo->fMethods;

   // collect doc strings
      Int_t nMethods = methods.size();

   // from template proxy with no instantiations 
      if ( nMethods == 0 )
         return NULL;

      PyObject* doc = methods[0]->GetDocString();

   // simple case
      if ( nMethods == 1 )
         return doc;

   // overloaded method
      PyObject* separator = PyROOT_PyUnicode_FromString( "\n" );
      for ( Int_t i = 1; i < nMethods; ++i ) {
         PyROOT_PyUnicode_Append( &doc, separator );
         PyROOT_PyUnicode_AppendAndDel( &doc, methods[i]->GetDocString() );
      }
      Py_DECREF( separator );

      return doc;
   }

////////////////////////////////////////////////////////////////////////////////
/// Create a new method proxy to be returned.

   PyObject* mp_meth_func( MethodProxy* pymeth, void* )
   {
      MethodProxy* newPyMeth = (MethodProxy*)MethodProxy_Type.tp_alloc( &MethodProxy_Type, 0 );

   // method info is shared, as it contains the collected overload knowledge
      *pymeth->fMethodInfo->fRefCount += 1;
      newPyMeth->fMethodInfo = pymeth->fMethodInfo;

   // new method is unbound, use of 'meth' is for keeping track whether this
   // proxy is used in the capacity of a method or a function
      newPyMeth->fSelf = (ObjectProxy*)newPyMeth;

      return (PyObject*)newPyMeth;
   }

////////////////////////////////////////////////////////////////////////////////
/// Return the bound self, if any; in case of pseudo-function role, pretend
/// that the data member im_self does not exist.

   PyObject* mp_meth_self( MethodProxy* pymeth, void* )
   {
      if ( IsPseudoFunc( pymeth ) ) {
         PyErr_Format( PyExc_AttributeError,
            "function %s has no attribute \'im_self\'", pymeth->fMethodInfo->fName.c_str() );
         return 0;
      } else if ( pymeth->fSelf != 0 ) {
         Py_INCREF( (PyObject*)pymeth->fSelf );
         return (PyObject*)pymeth->fSelf;
      }

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Return scoping class; in case of pseudo-function role, pretend that there
/// is no encompassing class (i.e. global scope).

   PyObject* mp_meth_class( MethodProxy* pymeth, void* )
   {
      if ( ! IsPseudoFunc( pymeth ) ) {
         PyObject* pyclass = pymeth->fMethodInfo->fMethods[0]->GetScopeProxy();
         if ( ! pyclass )
            PyErr_Format( PyExc_AttributeError,
               "function %s has no attribute \'im_class\'", pymeth->fMethodInfo->fName.c_str() );
         return pyclass;
      }

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Stub only, to fill out the python function interface.

   PyObject* mp_func_closure( MethodProxy* /* pymeth */, void* )
   {
      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Code details are used in module inspect to fill out interactive help()

   PyObject* mp_func_code( MethodProxy* pymeth, void* )
   {
#if PY_VERSION_HEX < 0x03000000
      MethodProxy::Methods_t& methods = pymeth->fMethodInfo->fMethods;

   // collect arguments only if there is just 1 overload, otherwise put in a
   // fake *args (see below for co_varnames)
      PyObject* co_varnames = methods.size() == 1 ? methods[0]->GetCoVarNames() : NULL;
      if ( !co_varnames ) {
      // TODO: static methods need no 'self' (but is harmless otherwise)
         co_varnames = PyTuple_New( 1 /* self */ + 1 /* fake */ );
         PyTuple_SET_ITEM( co_varnames, 0, PyROOT_PyUnicode_FromString( "self" ) );
         PyTuple_SET_ITEM( co_varnames, 1, PyROOT_PyUnicode_FromString( "*args" ) );
      }

      int co_argcount = PyTuple_Size( co_varnames );

   // for now, code object representing the statement 'pass'
      PyObject* co_code = PyString_FromStringAndSize( "d\x00\x00S", 4 );

   // tuples with all the const literals used in the function
      PyObject* co_consts = PyTuple_New( 0 );
      PyObject* co_names = PyTuple_New( 0 );

   // names, freevars, and cellvars go unused
      PyObject* co_unused = PyTuple_New( 0 );

   // filename is made-up
      PyObject* co_filename = PyString_FromString( "ROOT.py" );

   // name is the function name, also through __name__ on the function itself
      PyObject* co_name = PyString_FromString( pymeth->GetName().c_str() );

   // firstlineno is the line number of first function code in the containing scope

   // lnotab is a packed table that maps instruction count and line number
      PyObject* co_lnotab = PyString_FromString( "\x00\x01\x0c\x01" );

      PyObject* code = (PyObject*)PyCode_New(
         co_argcount,                             // argcount
         co_argcount + 1,                         // nlocals
         2,                                       // stacksize
         CO_OPTIMIZED | CO_NEWLOCALS | CO_NOFREE, // flags
         co_code,                                 // code
         co_consts,                               // consts
         co_names,                                // names
         co_varnames,                             // varnames
         co_unused,                               // freevars
         co_unused,                               // cellvars
         co_filename,                             // filename
         co_name,                                 // name
         1,                                       // firstlineno
         co_lnotab );                             // lnotab

      Py_DECREF( co_lnotab );
      Py_DECREF( co_name );
      Py_DECREF( co_unused );
      Py_DECREF( co_filename );
      Py_DECREF( co_varnames );
      Py_DECREF( co_names );
      Py_DECREF( co_consts );
      Py_DECREF( co_code );

      return code;
#else
// not important for functioning of most code, so not implemented for p3 for now (TODO)
      pymeth = 0;
      if ( pymeth || !pymeth) Py_INCREF( Py_None );
      return Py_None;
#endif
   }

////////////////////////////////////////////////////////////////////////////////
/// Create a tuple of default values, if there is only one method (otherwise
/// leave undefined: this is only used by inspect for interactive help())

   PyObject* mp_func_defaults( MethodProxy* pymeth, void* )
   {
      MethodProxy::Methods_t& methods = pymeth->fMethodInfo->fMethods;

      if ( methods.size() != 1 )
         return PyTuple_New( 0 );

      int maxarg = methods[0]->GetMaxArgs();

      PyObject* defaults = PyTuple_New( maxarg );

      int itup = 0;
      for ( int iarg = 0; iarg < maxarg; ++iarg ) {
         PyObject* defvalue = methods[0]->GetArgDefault( iarg );
         if ( defvalue )
            PyTuple_SET_ITEM( defaults, itup++, defvalue );
      }
      _PyTuple_Resize( &defaults, itup );

      return defaults;
   }

////////////////////////////////////////////////////////////////////////////////
/// Return this function's global dict (hard-wired to be the ROOT module); used
/// for lookup of names from co_code indexing into co_names.

   PyObject* mp_func_globals( MethodProxy* /* pymeth */, void* )
   {
      PyObject* pyglobal = PyModule_GetDict( PyImport_AddModule( (char*)"ROOT" ) );
      Py_XINCREF( pyglobal );
      return pyglobal;
   }

////////////////////////////////////////////////////////////////////////////////
/// Get '_creates' boolean, which determines ownership of return values.

   PyObject* mp_getcreates( MethodProxy* pymeth, void* )
   {
      return PyInt_FromLong( (Bool_t)IsCreator( pymeth->fMethodInfo->fFlags ) );
   }

////////////////////////////////////////////////////////////////////////////////
/// Set '_creates' boolean, which determines ownership of return values.

   int mp_setcreates( MethodProxy* pymeth, PyObject* value, void* )
   {
      if ( ! value ) {        // means that _creates is being deleted
         pymeth->fMethodInfo->fFlags &= ~TCallContext::kIsCreator;
         return 0;
      }

      Long_t iscreator = PyLong_AsLong( value );
      if ( iscreator == -1 && PyErr_Occurred() ) {
         PyErr_SetString( PyExc_ValueError, "a boolean 1 or 0 is required for _creates" );
         return -1;
      }

      if ( iscreator )
         pymeth->fMethodInfo->fFlags |= TCallContext::kIsCreator;
      else
         pymeth->fMethodInfo->fFlags &= ~TCallContext::kIsCreator;

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Get '_mempolicy' enum, which determines ownership of call arguments.

   PyObject* mp_getmempolicy( MethodProxy* pymeth, void* )
   {
      if ( (Bool_t)(pymeth->fMethodInfo->fFlags & TCallContext::kUseHeuristics ) )
         return PyInt_FromLong( TCallContext::kUseHeuristics );

      if ( (Bool_t)(pymeth->fMethodInfo->fFlags & TCallContext::kUseStrict ) )
         return PyInt_FromLong( TCallContext::kUseStrict );

      return PyInt_FromLong( -1 );
   }

////////////////////////////////////////////////////////////////////////////////
/// Set '_mempolicy' enum, which determines ownership of call arguments.

   int mp_setmempolicy( MethodProxy* pymeth, PyObject* value, void* )
   {
      Long_t mempolicy = PyLong_AsLong( value );
      if ( mempolicy == TCallContext::kUseHeuristics ) {
         pymeth->fMethodInfo->fFlags |= TCallContext::kUseHeuristics;
         pymeth->fMethodInfo->fFlags &= ~TCallContext::kUseStrict;
      } else if ( mempolicy == TCallContext::kUseStrict ) {
         pymeth->fMethodInfo->fFlags |= TCallContext::kUseStrict;
         pymeth->fMethodInfo->fFlags &= ~TCallContext::kUseHeuristics;
      } else {
         PyErr_SetString( PyExc_ValueError,
            "expected kMemoryStrict or kMemoryHeuristics as value for _mempolicy" );
         return -1;
      }

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Get '_manage_smart_ptr' boolean, which determines whether or not to
/// manage returned smart pointers intelligently.

   PyObject* mp_get_manage_smart_ptr( MethodProxy* pymeth, void* )
   {
      return PyInt_FromLong(
         (Bool_t)(pymeth->fMethodInfo->fFlags & TCallContext::kManageSmartPtr) );
   }

////////////////////////////////////////////////////////////////////////////////
/// Set '_manage_smart_ptr' boolean, which determines whether or not to
/// manage returned smart pointers intelligently.

   int mp_set_manage_smart_ptr( MethodProxy* pymeth, PyObject* value, void* )
   {
      Long_t policy = PyLong_AsLong( value );
      if ( policy == -1 && PyErr_Occurred() ) {
         PyErr_SetString( PyExc_ValueError, "a boolean 1 or 0 is required for _manage_smart_ptr" );
         return -1;
      }

      pymeth->fMethodInfo->fFlags |= TCallContext::kManageSmartPtr;

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Get '_threaded' boolean, which determines whether the GIL will be released.

   PyObject* mp_getthreaded( MethodProxy* pymeth, void* )
   {
      return PyInt_FromLong(
         (Bool_t)(pymeth->fMethodInfo->fFlags & TCallContext::kReleaseGIL) );
   }

////////////////////////////////////////////////////////////////////////////////
/// Set '_threaded' boolean, which determines whether the GIL will be released.

   int mp_setthreaded( MethodProxy* pymeth, PyObject* value, void* )
   {
      Long_t isthreaded = PyLong_AsLong( value );
      if ( isthreaded == -1 && PyErr_Occurred() ) {
         PyErr_SetString( PyExc_ValueError, "a boolean 1 or 0 is required for _creates" );
         return -1;
      }

      if ( isthreaded )
         pymeth->fMethodInfo->fFlags |= TCallContext::kReleaseGIL;
      else
         pymeth->fMethodInfo->fFlags &= ~TCallContext::kReleaseGIL;

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////

   PyGetSetDef mp_getset[] = {
      { (char*)"__name__",   (getter)mp_name,   NULL, NULL, NULL },
      { (char*)"__module__", (getter)mp_module, NULL, NULL, NULL },
      { (char*)"__doc__",    (getter)mp_doc,    NULL, NULL, NULL },

   // to be more python-like, where these are duplicated as well; to actually
   // derive from the python method or function type is too memory-expensive,
   // given that most of the members of those types would not be used
      { (char*)"im_func",  (getter)mp_meth_func,  NULL, NULL, NULL },
      { (char*)"im_self",  (getter)mp_meth_self,  NULL, NULL, NULL },
      { (char*)"im_class", (getter)mp_meth_class, NULL, NULL, NULL },

      { (char*)"func_closure",  (getter)mp_func_closure,  NULL, NULL, NULL },
      { (char*)"func_code",     (getter)mp_func_code,     NULL, NULL, NULL },
      { (char*)"func_defaults", (getter)mp_func_defaults, NULL, NULL, NULL },
      { (char*)"func_globals",  (getter)mp_func_globals,  NULL, NULL, NULL },
      { (char*)"func_doc",      (getter)mp_doc,           NULL, NULL, NULL },
      { (char*)"func_name",     (getter)mp_name,          NULL, NULL, NULL },

      { (char*)"_creates", (getter)mp_getcreates, (setter)mp_setcreates,
            (char*)"For ownership rules of result: if true, objects are python-owned", NULL },
      { (char*)"__creates__", (getter)mp_getcreates, (setter)mp_setcreates,
            (char*)"For ownership rules of result: if true, objects are python-owned", NULL },
      { (char*)"_mempolicy", (getter)mp_getmempolicy, (setter)mp_setmempolicy,
            (char*)"For argument ownership rules: like global, either heuristic or strict", NULL },
      { (char*)"__mempolicy__", (getter)mp_getmempolicy, (setter)mp_setmempolicy,
            (char*)"For argument ownership rules: like global, either heuristic or strict", NULL },
      { (char*)"_manage_smart_ptr", (getter)mp_get_manage_smart_ptr, (setter)mp_set_manage_smart_ptr,
        (char*)"If a smart pointer is returned, determines management policy.", NULL },
      { (char*)"_threaded", (getter)mp_getthreaded, (setter)mp_setthreaded,
            (char*)"If true, releases GIL on call into C++", NULL },
      { (char*)"__release_gil__", (getter)mp_getthreaded, (setter)mp_setthreaded,
            (char*)"If true, releases GIL on call into C++", NULL },
      { (char*)NULL, NULL, NULL, NULL, NULL }
   };

//= PyROOT method proxy function behavior ====================================
   PyObject* mp_call( MethodProxy* pymeth, PyObject* args, PyObject* kwds )
   {
   // Call the appropriate overload of this method.

   // if called through im_func pseudo-representation (this can be gamed if the
   // user really wants to ... )
      if ( IsPseudoFunc( pymeth ) )
         pymeth->fSelf = NULL;

      ObjectProxy* oldSelf = pymeth->fSelf;

   // get local handles to proxy internals
      auto& methods     = pymeth->fMethodInfo->fMethods;
      auto& dispatchMap = pymeth->fMethodInfo->fDispatchMap;
      auto& mflags      = pymeth->fMethodInfo->fFlags;

      Int_t nMethods = methods.size();

      TCallContext ctxt = { 0 };
      ctxt.fFlags |= (mflags & TCallContext::kUseHeuristics);
      ctxt.fFlags |= (mflags & TCallContext::kUseStrict);
      ctxt.fFlags |= (mflags & TCallContext::kManageSmartPtr);
      if ( ! ctxt.fFlags ) ctxt.fFlags |= TCallContext::sMemoryPolicy;
      ctxt.fFlags |= (mflags & TCallContext::kReleaseGIL);

   // simple case
      if ( nMethods == 1 ) {
         PyObject* result = methods[0]->Call( pymeth->fSelf, args, kwds, &ctxt );
         return HandleReturn( pymeth, oldSelf, result );
      }

   // otherwise, handle overloading
      Long_t sighash = HashSignature( args );

   // look for known signatures ...
      MethodProxy::DispatchMap_t::iterator m = dispatchMap.find( sighash );
      if ( m != dispatchMap.end() ) {
         Int_t index = m->second;
         PyObject* result = methods[ index ]->Call( pymeth->fSelf, args, kwds, &ctxt );
         result = HandleReturn( pymeth, oldSelf, result );

         if ( result != 0 )
            return result;

      // fall through: python is dynamic, and so, the hashing isn't infallible
         ResetCallState( pymeth->fSelf, oldSelf, kTRUE );
      }

   // ... otherwise loop over all methods and find the one that does not fail
      if ( ! IsSorted( mflags ) ) {
         std::stable_sort( methods.begin(), methods.end(), PriorityCmp );
         mflags |= TCallContext::kIsSorted;
      }

      std::vector< PyError_t > errors;
      for ( Int_t i = 0; i < nMethods; ++i ) {
         PyObject* result = methods[i]->Call( pymeth->fSelf, args, kwds, &ctxt );

         if ( result != 0 ) {
         // success: update the dispatch map for subsequent calls
            dispatchMap[ sighash ] = i;
            std::for_each( errors.begin(), errors.end(), PyError_t::Clear );
            return HandleReturn( pymeth, oldSelf, result );
         }

      // failure: collect error message/trace (automatically clears exception, too)
         if ( ! PyErr_Occurred() ) {
         // this should not happen; set an error to prevent core dump and report
            PyObject* sig = methods[i]->GetPrototype();
            PyErr_Format( PyExc_SystemError, "%s =>\n    %s",
               PyROOT_PyUnicode_AsString( sig ), (char*)"NULL result without error in mp_call" );
            Py_DECREF( sig );
         }
         PyError_t e;
         PyErr_Fetch( &e.fType, &e.fValue, &e.fTrace );
         errors.push_back( e );
         ResetCallState( pymeth->fSelf, oldSelf, kFALSE );
      }

   // first summarize, then add details
      PyObject* value = PyROOT_PyUnicode_FromFormat(
         "none of the %d overloaded methods succeeded. Full details:", nMethods );
      PyObject* separator = PyROOT_PyUnicode_FromString( "\n  " );

   // if this point is reached, none of the overloads succeeded: notify user
      PyObject* exc_type = NULL;
      for ( std::vector< PyError_t >::iterator e = errors.begin(); e != errors.end(); ++e ) {
         if ( e->fType != PyExc_NotImplementedError ) {
            if ( ! exc_type ) exc_type = e->fType;
            else if ( exc_type != e->fType ) exc_type = PyExc_TypeError;
         }
         PyROOT_PyUnicode_Append( &value, separator );
         PyROOT_PyUnicode_Append( &value, e->fValue );
      }

      Py_DECREF( separator );
      std::for_each( errors.begin(), errors.end(), PyError_t::Clear );

   // report failure
      PyErr_SetObject( exc_type ? exc_type : PyExc_TypeError, value );
      Py_DECREF( value );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Descriptor; create and return a new bound method proxy (language requirement).

   MethodProxy* mp_descrget( MethodProxy* pymeth, ObjectProxy* pyobj, PyObject* )
   {
      MethodProxy* newPyMeth = (MethodProxy*)MethodProxy_Type.tp_alloc( &MethodProxy_Type, 0 );

   // method info is shared, as it contains the collected overload knowledge
      *pymeth->fMethodInfo->fRefCount += 1;
      newPyMeth->fMethodInfo = pymeth->fMethodInfo;

   // new method is to be bound to current object (may be NULL)
      Py_XINCREF( (PyObject*)pyobj );
      newPyMeth->fSelf = pyobj;

      return newPyMeth;
   }


//= PyROOT method proxy construction/destruction =================================
   MethodProxy* mp_new( PyTypeObject*, PyObject*, PyObject* )
   {
   // Create a new method proxy object.
      MethodProxy* pymeth = PyObject_GC_New( MethodProxy, &MethodProxy_Type );
      pymeth->fSelf = NULL;
      pymeth->fMethodInfo = new MethodProxy::MethodInfo_t;

      PyObject_GC_Track( pymeth );
      return pymeth;
   }

////////////////////////////////////////////////////////////////////////////////
/// Deallocate memory held by method proxy object.

   void mp_dealloc( MethodProxy* pymeth )
   {
      PyObject_GC_UnTrack( pymeth );

      if ( ! IsPseudoFunc( pymeth ) )
         Py_CLEAR( pymeth->fSelf );
      pymeth->fSelf = NULL;

      if ( --(*pymeth->fMethodInfo->fRefCount) <= 0 ) {
         delete pymeth->fMethodInfo;
      }

      PyObject_GC_Del( pymeth );
   }


////////////////////////////////////////////////////////////////////////////////
/// Hash of method proxy object for insertion into dictionaries; with actual
/// method (fMethodInfo) shared, its address is best suited.

   Long_t mp_hash( MethodProxy* pymeth )
   {
      return _Py_HashPointer( pymeth->fMethodInfo );
   }

////////////////////////////////////////////////////////////////////////////////
/// Garbage collector traverse of held python member objects.

   int mp_traverse( MethodProxy* pymeth, visitproc visit, void* args )
   {
      if ( pymeth->fSelf && ! IsPseudoFunc( pymeth ) )
         return visit( (PyObject*)pymeth->fSelf, args );

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Garbage collector clear of held python member objects.

   int mp_clear( MethodProxy* pymeth )
   {
      if ( ! IsPseudoFunc( pymeth ) )
         Py_CLEAR( pymeth->fSelf );
      pymeth->fSelf = NULL;

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Rich set of comparison objects; only equals is defined.

   PyObject* mp_richcompare( MethodProxy* self, MethodProxy* other, int op )
   {
      if ( op != Py_EQ )
         return PyType_Type.tp_richcompare( (PyObject*)self, (PyObject*)other, op );

      // defined by type + (shared) MethodInfo + bound self, with special case for fSelf (i.e. pseudo-function)
      if ( ( Py_TYPE(self) == Py_TYPE(other) && self->fMethodInfo == other->fMethodInfo ) && \
           ( ( IsPseudoFunc( self ) && IsPseudoFunc( other ) ) || self->fSelf == other->fSelf ) ) {
         Py_INCREF( Py_True );
         return Py_True;
      }
      Py_INCREF( Py_False );
      return Py_False;
   }


//= PyROOT method proxy access to internals =================================
   PyObject* mp_disp( MethodProxy* pymeth, PyObject* sigarg )
   {
   // Select and call a specific C++ overload, based on its signature.
      if ( ! PyROOT_PyUnicode_Check( sigarg ) ) {
         PyErr_Format( PyExc_TypeError, "disp() argument 1 must be string, not %.50s",
                       sigarg == Py_None ? "None" : Py_TYPE(sigarg)->tp_name );
         return 0;
      }

      PyObject* sig1 = PyROOT_PyUnicode_FromFormat( "(%s)", PyROOT_PyUnicode_AsString( sigarg ) );

      MethodProxy::Methods_t& methods = pymeth->fMethodInfo->fMethods;
      for ( Int_t i = 0; i < (Int_t)methods.size(); ++i ) {

         PyObject* sig2 = methods[ i ]->GetSignature();
         if ( PyObject_RichCompareBool( sig1, sig2, Py_EQ ) ) {
            Py_DECREF( sig2 );

            MethodProxy* newmeth = mp_new( NULL, NULL, NULL );
            MethodProxy::Methods_t vec; vec.push_back( methods[ i ]->Clone() );
            newmeth->Set( pymeth->fMethodInfo->fName, vec );

            if ( pymeth->fSelf && ! IsPseudoFunc( pymeth ) ) {
               Py_INCREF( pymeth->fSelf );
               newmeth->fSelf = pymeth->fSelf;
            }

            Py_DECREF( sig1 );
            return (PyObject*)newmeth;
         }

         Py_DECREF( sig2 );
      }

      Py_DECREF( sig1 );
      PyErr_Format( PyExc_LookupError, "signature \"%s\" not found", PyROOT_PyUnicode_AsString( sigarg ) );
      return 0;
   }

//= PyROOT method proxy access to internals =================================
   PyObject* mp_add_overload( MethodProxy* pymeth, PyObject* new_overload )
   {
      TPythonCallback* cb = new TPythonCallback(new_overload);
      pymeth->AddMethod( cb );
      Py_INCREF( Py_None );
      return Py_None;
   }

   PyMethodDef mp_methods[] = {
      { (char*)"disp",             (PyCFunction)mp_disp, METH_O, (char*)"select overload for dispatch" },
      { (char*)"__add_overload__", (PyCFunction)mp_add_overload, METH_O, (char*)"add a new overload" },
      { (char*)NULL, NULL, 0, NULL }
   };

} // unnamed namespace

////////////////////////////////////////////////////////////////////////////////


//= PyROOT method proxy type =================================================
PyTypeObject MethodProxy_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.MethodProxy", // tp_name
   sizeof(MethodProxy),       // tp_basicsize
   0,                         // tp_itemsize
   (destructor)mp_dealloc,    // tp_dealloc
   0,                         // tp_print (python < 3.8)
                              // tp_vectorcall_offset (python >= 3.8)
   0,                         // tp_getattr
   0,                         // tp_setattr
   0,                         // tp_compare
   0,                         // tp_repr
   0,                         // tp_as_number
   0,                         // tp_as_sequence
   0,                         // tp_as_mapping
   (hashfunc)mp_hash,         // tp_hash
   (ternaryfunc)mp_call,      // tp_call
   0,                         // tp_str
   0,                         // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,      // tp_flags
   (char*)"PyROOT method proxy (internal)",      // tp_doc
   (traverseproc)mp_traverse, // tp_traverse
   (inquiry)mp_clear,         // tp_clear
   (richcmpfunc)mp_richcompare,                  // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   mp_methods,                // tp_methods
   0,                         // tp_members
   mp_getset,                 // tp_getset
   0,                         // tp_base
   0,                         // tp_dict
   (descrgetfunc)mp_descrget, // tp_descr_get
   0,                         // tp_descr_set
   0,                         // tp_dictoffset
   0,                         // tp_init
   0,                         // tp_alloc
   (newfunc)mp_new,           // tp_new
   0,                         // tp_free
   0,                         // tp_is_gc
   0,                         // tp_bases
   0,                         // tp_mro
   0,                         // tp_cache
   0,                         // tp_subclasses
   0                          // tp_weaklist
#if PY_VERSION_HEX >= 0x02030000
   , 0                        // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
   , 0                        // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x03040000
   , 0                        // tp_finalize
#endif
#if PY_VERSION_HEX >= 0x03080000
   , 0                        // tp_vectorcall
#if PY_VERSION_HEX < 0x03090000
   , 0                        // tp_print (python 3.8 only)
#endif
#endif
};

} // namespace PyROOT


//- public members -----------------------------------------------------------
void PyROOT::MethodProxy::Set( const std::string& name, std::vector< PyCallable* >& methods )
{
// Fill in the data of a freshly created method proxy.
   fMethodInfo->fName = name;
   fMethodInfo->fMethods.swap( methods );
   fMethodInfo->fFlags &= ~TCallContext::kIsSorted;
   fMethodInfo->fFlags |= TCallContext::kManageSmartPtr;

// special case: all constructors are considered creators by default
   if ( name == "__init__" )
      fMethodInfo->fFlags |= (TCallContext::kIsCreator | TCallContext::kIsConstructor);

// special case, in heuristics mode also tag *Clone* methods as creators
   if ( TCallContext::sMemoryPolicy == TCallContext::kUseHeuristics && \
        name.find( "Clone" ) != std::string::npos )
      fMethodInfo->fFlags |= TCallContext::kIsCreator;
}

////////////////////////////////////////////////////////////////////////////////
/// Fill in the data of a freshly created method proxy.

void PyROOT::MethodProxy::AddMethod( PyCallable* pc )
{
   fMethodInfo->fMethods.push_back( pc );
   fMethodInfo->fFlags &= ~TCallContext::kIsSorted;
}

////////////////////////////////////////////////////////////////////////////////

void PyROOT::MethodProxy::AddMethod( MethodProxy* meth )
{
   fMethodInfo->fMethods.insert( fMethodInfo->fMethods.end(),
      meth->fMethodInfo->fMethods.begin(), meth->fMethodInfo->fMethods.end() );
   fMethodInfo->fFlags &= ~TCallContext::kIsSorted;
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor (this object is reference counted).

PyROOT::MethodProxy::MethodInfo_t::~MethodInfo_t()
{
   for ( Methods_t::iterator it = fMethods.begin(); it != fMethods.end(); ++it ) {
      delete *it;
   }
   fMethods.clear();
   delete fRefCount;
}
