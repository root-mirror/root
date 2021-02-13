// Author: Wim Lavrijsen, Jan 2005

// Bindings
#include "PyROOT.h"
#include "TemplateProxy.h"
#include "MethodProxy.h"
#include "TFunctionHolder.h"
#include "TMethodHolder.h"
#include "PyCallable.h"
#include "PyStrings.h"
#include "Utility.h"
#include "PyRootType.h"

// ROOT
#include "TClass.h"
#include "TMethod.h"


namespace PyROOT {

////////////////////////////////////////////////////////////////////////////////
/// Initialize the proxy for the given 'pyclass.'

void TemplateProxy::Set( const std::string& name, PyObject* pyclass )
{
   fPyName       = PyROOT_PyUnicode_FromString( const_cast< char* >( name.c_str() ) );
   Py_XINCREF( pyclass );
   fPyClass      = pyclass;
   fSelf         = NULL;
   std::vector< PyCallable* > dummy;
   fNonTemplated = MethodProxy_New( name, dummy );
   fTemplated    = MethodProxy_New( name, dummy );
}

////////////////////////////////////////////////////////////////////////////////
/// Store overloads of this templated method.

void TemplateProxy::AddOverload( MethodProxy* mp ) {
   fNonTemplated->AddMethod( mp );
}

void TemplateProxy::AddOverload( PyCallable* pc ) {
// Store overload of this templated method.
   fNonTemplated->AddMethod( pc );
}

void TemplateProxy::AddTemplate( PyCallable* pc )
{
// Store know template methods.
   fTemplated->AddMethod( pc );
}


namespace {

//= PyROOT template proxy construction/destruction ===========================
   TemplateProxy* tpp_new( PyTypeObject*, PyObject*, PyObject* )
   {
   // Create a new empty template method proxy.
      TemplateProxy* pytmpl = PyObject_GC_New( TemplateProxy, &TemplateProxy_Type );
      pytmpl->fPyName       = NULL;
      pytmpl->fPyClass      = NULL;
      pytmpl->fSelf         = NULL;
      pytmpl->fNonTemplated = NULL;
      pytmpl->fTemplated    = NULL;

      PyObject_GC_Track( pytmpl );
      return pytmpl;
   }

////////////////////////////////////////////////////////////////////////////////
/// Garbage collector clear of held python member objects.

   int tpp_clear( TemplateProxy* pytmpl )
   {
      Py_CLEAR( pytmpl->fPyName );
      Py_CLEAR( pytmpl->fPyClass );
      Py_CLEAR( pytmpl->fSelf );
      Py_CLEAR( pytmpl->fNonTemplated );
      Py_CLEAR( pytmpl->fTemplated );

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Destroy the given template method proxy.

   void tpp_dealloc( TemplateProxy* pytmpl )
   {
      PyObject_GC_UnTrack( pytmpl );
      tpp_clear( pytmpl );
      PyObject_GC_Del( pytmpl );
   }

////////////////////////////////////////////////////////////////////////////////
/// Forward to method proxies to doc all overloads

   PyObject* tpp_doc( TemplateProxy* pytmpl, void* )
   {
      PyObject* doc = nullptr;
      if ( pytmpl->fNonTemplated )
         doc = PyObject_GetAttrString( (PyObject*)pytmpl->fNonTemplated, "__doc__" );
      if ( pytmpl->fTemplated ) {
         PyObject* doc2 = PyObject_GetAttrString( (PyObject*)pytmpl->fTemplated, "__doc__" );
         if ( doc && doc2 ) {
            PyROOT_PyUnicode_AppendAndDel( &doc, PyROOT_PyUnicode_FromString( "\n" ));
            PyROOT_PyUnicode_AppendAndDel( &doc, doc2 );
         } else if ( !doc && doc2 ) {
            doc = doc2;
         }
      }

      if ( doc )
         return doc;

      return PyROOT_PyUnicode_FromString( TemplateProxy_Type.tp_doc );
   }

////////////////////////////////////////////////////////////////////////////////
/// Garbage collector traverse of held python member objects.

   int tpp_traverse( TemplateProxy* pytmpl, visitproc visit, void* arg )
   {
      Py_VISIT( pytmpl->fPyName );
      Py_VISIT( pytmpl->fPyClass );
      Py_VISIT( pytmpl->fSelf );
      Py_VISIT( pytmpl->fNonTemplated );
      Py_VISIT( pytmpl->fTemplated );

      return 0;
   }

//= PyROOT template proxy callable behavior ==================================
   PyObject* tpp_call( TemplateProxy* pytmpl, PyObject* args, PyObject* kwds )
   {
   // Dispatcher to the actual member method, several uses possible; in order:
   //
   // case 1:
   //
   //    obj.method( a0, a1, ... )              # non-template
   //       => obj->method( a0, a1, ... )      // non-template
   //
   // case 2:
   //
   //    obj.method( t0, t1, ... )( a0, a1, ... )
   //       -> getattr( obj, 'method< t0, t1, ... >' )( a0, a1, ... )
   //       => obj->method< t0, t1, ... >( a0, a1, ... )
   //
   // case 3:
   //
   //    obj.method( a0, a1, ... )              # all known templates
   //       => obj->method( a0, a1, ... )      // all known templates
   //
   // case 4:
   //
   //    collect types of arguments unless they are a type themselves (the idea
   //    here is it is more likely for e.g. (1, 3.14) to be real arguments and
   //    e.g. (int, 'double') to be template arguments:
   //       a) if ! is_type(arg)
   //             template<> method< T0, T1, ... >( type(a0), type(a1), ... )
   //       b) else
   //             drop to case 5
   //    TODO: there is ambiguity in the above if the arguments are strings
   //
   // case 5:
   //
   //    instantiate template<> method< t0, t1, ... >

   // case 1: obj->method( a0, a1, ... )

   // simply forward the call: all non-templated methods are defined on class definition
   // and thus already available
      PyObject* pymeth = MethodProxy_Type.tp_descr_get(
         (PyObject*)pytmpl->fNonTemplated, pytmpl->fSelf, (PyObject*)&MethodProxy_Type );
      if ( MethodProxy_Check( pymeth ) ) {
      // now call the method with the arguments
         PyObject* result = MethodProxy_Type.tp_call( pymeth, args, kwds );
         Py_DECREF( pymeth ); pymeth = 0;
         if ( result )
            return result;
      // TODO: collect error here, as the failure may be either an overload
      // failure after which we should continue; or a real failure, which should
      // be reported.
      }
      Py_XDECREF( pymeth ); pymeth = 0;
      PyErr_Clear();

   // error check on method() which can not be derived if non-templated case fails
      Py_ssize_t nArgs = PyTuple_GET_SIZE( args );
      if ( nArgs == 0 ) {
         PyErr_Format( PyExc_TypeError, "template method \'%s\' with no arguments must be explicit",
            PyROOT_PyUnicode_AsString( pytmpl->fPyName ) );
         return 0;
      }

   // case 2: non-instantiating obj->method< t0, t1, ... >( a0, a1, ... )

      Bool_t isType = kFALSE;
      Int_t nStrings = 0;
      PyObject* tpArgs = PyTuple_New( nArgs );
      for ( Int_t i = 0; i < nArgs; ++i ) {
         PyObject* itemi = PyTuple_GET_ITEM( args, i );
         if ( PyType_Check( itemi ) ) isType = kTRUE;
#if PY_VERSION_HEX >= 0x03000000
         else if ( ! isType && PyUnicode_Check( itemi ) ) nStrings += 1;
#else
         else if ( ! isType && PyBytes_Check( itemi ) ) nStrings += 1;
#endif
      // special case for arrays
         PyObject* pytc = PyObject_GetAttr( itemi, PyStrings::gTypeCode );
         if ( ! ( pytc && PyROOT_PyUnicode_Check( pytc ) ) ) {
         // normal case (not an array)
            PyErr_Clear();
            PyObject* tp = (PyObject*)Py_TYPE( itemi );
            Py_INCREF( tp );
            PyTuple_SET_ITEM( tpArgs, i, tp );
         } else {
         // array, build up a pointer type
            char tc = ((char*)PyROOT_PyUnicode_AsString( pytc ))[0];
            const char* ptrname = 0;
            switch ( tc ) {
               case 'b': ptrname = "char*";           break;
               case 'h': ptrname = "short*";          break;
               case 'H': ptrname = "unsigned short*"; break;
               case 'i': ptrname = "int*";            break;
               case 'I': ptrname = "unsigned int*";   break;
               case 'l': ptrname = "long*";           break;
               case 'L': ptrname = "unsigned long*";  break;
               case 'f': ptrname = "float*";          break;
               case 'd': ptrname = "double*";         break;
               default:  ptrname = "void*";  // TODO: verify if this is right
            }
            if ( ptrname ) {
               //PyObject* pyptrname = PyBytes_FromString( ptrname );
               PyObject *pyptrname = PyROOT_PyUnicode_FromString(ptrname);
               PyTuple_SET_ITEM( tpArgs, i, pyptrname );
            // string added, but not counted towards nStrings
            } else {
            // this will cleanly fail instantiation
               Py_INCREF( pytc );
               PyTuple_SET_ITEM( tpArgs, i, pytc );
            }
         }
         Py_XDECREF( pytc );
      }

   // build "< type, type, ... >" part of method name
      PyObject* pyname_v1 = Utility::BuildTemplateName( pytmpl->fPyName, args, 0 );
      if ((isType || nStrings == nArgs) && pyname_v1) {  // types in args or all strings
      // lookup method on self (to make sure it propagates), which is readily callable
         pymeth = PyObject_GetAttr( pytmpl->fSelf ? pytmpl->fSelf : pytmpl->fPyClass, pyname_v1 );
         if ( pymeth ) { // overloads stop here, as this is an explicit match
            Py_DECREF( pyname_v1 );
            if (PyErr_WarnEx(PyExc_FutureWarning,
                       "Instantiating a function template with parentheses ( f(type1, ..., typeN) ) "
                       "is deprecated and will not be supported in a future version of ROOT. "
                       "Instead, use square brackets: f[type1, ..., typeN]", 1) < 0) {
               return nullptr;
            }
            return pymeth;         // callable method, next step is by user
         }
      }
      PyErr_Clear();

   // case 3: loop over all previously instantiated templates
      pymeth = MethodProxy_Type.tp_descr_get(
         (PyObject*)pytmpl->fTemplated, pytmpl->fSelf, (PyObject*)&MethodProxy_Type );
      if ( MethodProxy_Check( pymeth ) ) {
      // now call the method with the arguments
         PyObject* result = MethodProxy_Type.tp_call( pymeth, args, kwds );
         Py_DECREF( pymeth ); pymeth = 0;
         if ( result ) {
            Py_XDECREF( pyname_v1 );
            return result;
         }
      // TODO: collect error here, as the failure may be either an overload
      // failure after which we should continue; or a real failure, which should
      // be reported.
      }
      Py_XDECREF( pymeth ); pymeth = 0;
      PyErr_Clear();

   // still here? try instantiating methods

      PyObject* clName = PyObject_GetAttr( pytmpl->fPyClass, PyStrings::gCppName );
      if (!clName) {
         PyErr_Clear();
         clName = PyObject_GetAttr(pytmpl->fPyClass, PyStrings::gName);
      }
      auto clNameStr = std::string(PyROOT_PyUnicode_AsString(clName));
      if (clNameStr == "_global_cpp")
         clNameStr = ""; // global namespace
      TClass* klass = TClass::GetClass(clNameStr.c_str());
      Py_DECREF( clName );
      const std::string& tmplname = pytmpl->fNonTemplated->fMethodInfo->fName;

    // case 4a: instantiating obj->method< T0, T1, ... >( type(a0), type(a1), ... )( a0, a1, ... )
      if ( ! isType && ! ( nStrings == nArgs ) ) {    // no types among args and not all strings
         for (auto pref : {Utility::kReference, Utility::kPointer, Utility::kValue}) {
            int pcnt = 0;
            PyObject* pyname_v2 = Utility::BuildTemplateName( NULL, tpArgs, 0, args, pref, &pcnt, true );
            if ( pyname_v2 ) {
               std::string mname = PyROOT_PyUnicode_AsString( pyname_v2 );
               Py_DECREF( pyname_v2 );
               std::string proto = mname.substr( 1, mname.size() - 2 );
            // the following causes instantiation as necessary
               auto scope = Cppyy::GetScope(clNameStr);
               auto cppmeth = Cppyy::GetMethodTemplate(scope, tmplname, proto);
               if ( cppmeth ) {    // overload stops here
                  Py_XDECREF( pyname_v1 );
                  if (Cppyy::IsNamespace(scope) || Cppyy::IsStaticMethod(cppmeth)) {
                     pytmpl->fTemplated->AddMethod( new TFunctionHolder( scope, cppmeth ) );
                     pymeth = (PyObject*)MethodProxy_New(
                        Cppyy::GetMethodName(cppmeth).c_str(), new TFunctionHolder( scope, cppmeth ) );
                  } else {
                     pytmpl->fTemplated->AddMethod( new TMethodHolder( scope, cppmeth ) );
                     pymeth = (PyObject*)MethodProxy_New(
                        Cppyy::GetMethodName(cppmeth).c_str(), new TMethodHolder( scope, cppmeth ) );
                  }
                  PyObject_SetAttrString( pytmpl->fPyClass, (char*)Cppyy::GetMethodName(cppmeth).c_str(), (PyObject*)pymeth );
                  Py_DECREF( pymeth );
                  pymeth = PyObject_GetAttrString(
                     pytmpl->fSelf ? pytmpl->fSelf : pytmpl->fPyClass, (char*)Cppyy::GetMethodName(cppmeth).c_str() );
                  PyObject* result = MethodProxy_Type.tp_call( pymeth, args, kwds );
                  Py_DECREF( pymeth );
                  return result;
               }
            }
            if (!pcnt) break; // preference never used; no point trying others
         }
      }

   // case 4b/5: instantiating obj->method< t0, t1, ... >( a0, a1, ... )
      if ( pyname_v1 ) {
         std::string mname = PyROOT_PyUnicode_AsString( pyname_v1 );
       // the following causes instantiation as necessary
         TMethod* cppmeth = klass ? klass->GetMethodAny( mname.c_str() ) : 0;
         if ( cppmeth ) {    // overload stops here
            pymeth = (PyObject*)MethodProxy_New(
               mname, new TMethodHolder( Cppyy::GetScope( klass->GetName() ), (Cppyy::TCppMethod_t)cppmeth ) );
            PyObject_SetAttr( pytmpl->fPyClass, pyname_v1, (PyObject*)pymeth );
            if ( mname != cppmeth->GetName() ) // happens with typedefs and template default arguments
               PyObject_SetAttrString( pytmpl->fPyClass, (char*)mname.c_str(), (PyObject*)pymeth );
            Py_DECREF( pymeth );
            pymeth = PyObject_GetAttr( pytmpl->fSelf ? pytmpl->fSelf : pytmpl->fPyClass, pyname_v1 );
            Py_DECREF( pyname_v1 );
            if (PyErr_WarnEx(PyExc_FutureWarning,
                       "Instantiating a function template with parentheses ( f(type1, ..., typeN) ) "
                       "is deprecated and will not be supported in a future version of ROOT. "
                       "Instead, use square brackets: f[type1, ..., typeN]", 1) < 0) {
               return nullptr;
            }
            return pymeth;         // callable method, next step is by user
         }
         Py_DECREF( pyname_v1 );
      }

   // moderately generic error message, but should be clear enough
      PyErr_Format( PyExc_TypeError, "can not resolve method template call for \'%s\'",
         PyROOT_PyUnicode_AsString( pytmpl->fPyName ) );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// create and use a new template proxy (language requirement)

   TemplateProxy* tpp_descrget( TemplateProxy* pytmpl, PyObject* pyobj, PyObject* )
   {
      TemplateProxy* newPyTmpl = (TemplateProxy*)TemplateProxy_Type.tp_alloc( &TemplateProxy_Type, 0 );

   // copy name and class pointers
      Py_INCREF( pytmpl->fPyName );
      newPyTmpl->fPyName = pytmpl->fPyName;

      Py_XINCREF( pytmpl->fPyClass );
      newPyTmpl->fPyClass = pytmpl->fPyClass;

   // copy non-templated method proxy pointer
      Py_INCREF( pytmpl->fNonTemplated );
      newPyTmpl->fNonTemplated = pytmpl->fNonTemplated;

   // copy templated method proxy pointer
      Py_INCREF( pytmpl->fTemplated );
      newPyTmpl->fTemplated = pytmpl->fTemplated;

   // new method is to be bound to current object (may be NULL)
      Py_XINCREF( pyobj );
      newPyTmpl->fSelf = pyobj;

      return newPyTmpl;
   }

////////////////////////////////////////////////////////////////////////////////
/// Explicit instantiation with square bracket syntax.
/// Implementation is equivalent to case 2&4b in tpp_call (parenthesis syntax)

   PyObject *tpp_subscript(TemplateProxy *pytmpl, PyObject *args)
   {
      bool justOne = !PyTuple_CheckExact(args);
      Py_ssize_t nArgs;
      if (justOne) {
         nArgs = 1;
         auto item = args;
         args = PyTuple_New(nArgs);
         Py_INCREF(item);
         PyTuple_SET_ITEM(args, 0, item);
      } else {
         nArgs = PyTuple_GET_SIZE(args);
      }

      Bool_t isType = false;
      Int_t nStrings = 0;
      for (int i = 0; i < nArgs; ++i) {
         PyObject* itemi = PyTuple_GET_ITEM(args, i);
         if (PyType_Check(itemi)) isType = kTRUE;
#if PY_VERSION_HEX >= 0x03000000
         else if (! isType && PyUnicode_Check(itemi)) nStrings += 1;
#else
         else if (! isType && PyBytes_Check(itemi)) nStrings += 1;
#endif
      }

      // Build "< type, type, ... >" part of method name
      PyObject* pyname = Utility::BuildTemplateName(pytmpl->fPyName, args, 0);
      if (justOne) Py_DECREF(args);

      // Non-instantiating obj->method< t0, t1, ... >( a0, a1, ... )
      if ((isType || nStrings == nArgs) && pyname) {  // types in args or all strings
         // Lookup method on self (to make sure it propagates), which is readily callable
         PyObject* pymeth = PyObject_GetAttr(pytmpl->fSelf ? pytmpl->fSelf : pytmpl->fPyClass, pyname);
         if (pymeth) { // Overloads stop here, as this is an explicit match
            Py_DECREF(pyname);
            return pymeth; // Callable method, next step is by user
         }
         PyErr_Clear();
      }

      // Still here? try instantiating methods
      PyObject* clName = PyObject_GetAttr(pytmpl->fPyClass, PyStrings::gCppName);
      if (!clName) {
         PyErr_Clear();
         clName = PyObject_GetAttr(pytmpl->fPyClass, PyStrings::gName);
      }
      auto clNameStr = std::string(PyROOT_PyUnicode_AsString(clName));
      if (clNameStr == "_global_cpp")
         clNameStr = ""; // global namespace
      auto klass = TClass::GetClass(clNameStr.c_str());
      Py_DECREF(clName);

      // Instantiating obj->method< t0, t1, ... >( a0, a1, ... )
      if (pyname) {
         std::string mname = PyROOT_PyUnicode_AsString(pyname);
         // The following causes instantiation as necessary
         TMethod *cppmeth = klass ? klass->GetMethodAny(mname.c_str()) : nullptr;
         if (cppmeth) {    // overload stops here
            PyObject *pymeth = (PyObject*)MethodProxy_New(
               mname, new TMethodHolder(Cppyy::GetScope(klass->GetName()),(Cppyy::TCppMethod_t)cppmeth));
            PyObject_SetAttr(pytmpl->fPyClass, pyname, (PyObject*)pymeth);
            if (mname != cppmeth->GetName()) // happens with typedefs and template default arguments
               PyObject_SetAttrString(pytmpl->fPyClass, (char*)mname.c_str(), (PyObject*)pymeth);
            Py_DECREF(pymeth);
            pymeth = PyObject_GetAttr(pytmpl->fSelf ? pytmpl->fSelf : pytmpl->fPyClass, pyname);
            Py_DECREF(pyname);
            return pymeth;         // callable method, next step is by user
         }
         Py_DECREF(pyname);
      }

      PyErr_Format(PyExc_TypeError, "cannot resolve method template instantiation for \'%s\'",
         PyROOT_PyUnicode_AsString(pytmpl->fPyName));
      return nullptr;
   }

////////////////////////////////////////////////////////////////////////////////

   static PyMappingMethods tpp_as_mapping = {
      nullptr, (binaryfunc)tpp_subscript, nullptr
   };

   PyGetSetDef tpp_getset[] = {
      { (char*)"__doc__",    (getter)tpp_doc,    NULL, NULL, NULL },
      { (char*)NULL, NULL, NULL, NULL, NULL }
   };

} // unnamed namespace


//= PyROOT template proxy type ===============================================
PyTypeObject TemplateProxy_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.TemplateProxy", // tp_name
   sizeof(TemplateProxy),     // tp_basicsize
   0,                         // tp_itemsize
   (destructor)tpp_dealloc,   // tp_dealloc
   0,                         // tp_print (python < 3.8)
                              // tp_vectorcall_offset (python >= 3.8)
   0,                         // tp_getattr
   0,                         // tp_setattr
   0,                         // tp_compare
   0,                         // tp_repr
   0,                         // tp_as_number
   0,                         // tp_as_sequence
   &tpp_as_mapping,           // tp_as_mapping
   0,                         // tp_hash
   (ternaryfunc)tpp_call,     // tp_call
   0,                         // tp_str
   0,                         // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,      // tp_flags
   (char*)"PyROOT template proxy (internal)",    // tp_doc
   (traverseproc)tpp_traverse,// tp_traverse
   (inquiry)tpp_clear,        // tp_clear
   0,                         // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   tpp_getset,                // tp_getset
   0,                         // tp_base
   0,                         // tp_dict
   (descrgetfunc)tpp_descrget,// tp_descr_get
   0,                         // tp_descr_set
   0,                         // tp_dictoffset
   0,                         // tp_init
   0,                         // tp_alloc
   (newfunc)tpp_new,          // tp_new
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
