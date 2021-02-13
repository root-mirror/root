// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

// Bindings
#include "PyROOT.h"
#include "PyRootType.h"
#include "MethodProxy.h"
#include "PropertyProxy.h"
#include "RootWrapper.h"
#include "TFunctionHolder.h"
#include "TemplateProxy.h"
#include "PyStrings.h"

// ROOT
#include "TClass.h"     // for method and enum finding
#include "TList.h"      // id.

// Standard
#include <string.h>
#include <string>
#include <vector>

namespace PyROOT {

namespace {

//= PyROOT type proxy construction/destruction ===============================
   PyObject* meta_alloc( PyTypeObject* metatype, Py_ssize_t nitems )
   {
   // specialized allocator, fitting in a few extra bytes for a TClassRef
      PyObject* pyclass = PyType_Type.tp_alloc( metatype, nitems );

      return pyclass;
   }

////////////////////////////////////////////////////////////////////////////////

   void meta_dealloc( PyRootClass* pytype )
   {
      return PyType_Type.tp_dealloc( (PyObject*)pytype );
   }

////////////////////////////////////////////////////////////////////////////////
/// Called when PyRootType acts as a metaclass; since type_new always resets
/// tp_alloc, and since it does not call tp_init on types, the metaclass is
/// being fixed up here, and the class is initialized here as well.

   PyObject* pt_new( PyTypeObject* subtype, PyObject* args, PyObject* kwds )
   {
   // fixup of metaclass (left permanent, and in principle only called once b/c
   // PyROOT caches python classes)
      subtype->tp_alloc   = (allocfunc)meta_alloc;
      subtype->tp_dealloc = (destructor)meta_dealloc;

   // creation of the python-side class
      PyRootClass* result = (PyRootClass*)PyType_Type.tp_new( subtype, args, kwds );

   // initialization of class (based on name only, initially, which is lazy)

   // there's a snag here: if a python class is derived from the bound class,
   // the name will not be known by TClassRef, hence we'll use the meta class
   // name from the subtype, rather than given class name

      const char* mp = strstr( subtype->tp_name, "_meta" );
      if ( ! mp ) {
      // there has been a user meta class override in a derived class, so do
      // the consistent thing, thus allowing user control over naming
         result->fCppType = Cppyy::GetScope(
            PyROOT_PyUnicode_AsString( PyTuple_GET_ITEM( args, 0 ) ) );
      } else {
      // coming here from PyROOT, use meta class name instead of given name,
      // so that it is safe to inherit python classes from the bound class
         result->fCppType = Cppyy::GetScope(
            std::string( subtype->tp_name ).substr( 0, mp-subtype->tp_name ).c_str() );
      }

      return (PyObject*)result;
   }


//= PyROOT type metaclass behavior ===========================================
   PyObject* pt_getattro( PyObject* pyclass, PyObject* pyname )
   {
   // normal type lookup
      PyObject* attr = PyType_Type.tp_getattro( pyclass, pyname );

   // extra ROOT lookup in case of failure (e.g. for inner classes on demand)
      if ( ! attr && PyROOT_PyUnicode_CheckExact( pyname ) ) {
         PyObject *etype, *value, *trace;
         PyErr_Fetch( &etype, &value, &trace );         // clears current exception

      // filter for python specials and lookup qualified class or function
         std::string name = PyROOT_PyUnicode_AsString( pyname );

         if ( name.size() <= 2 || name.substr( 0, 2 ) != "__" ) {
            attr = CreateScopeProxy( name, pyclass );

         // namespaces may have seen updates in their list of global functions, which
         // are available as "methods" even though they're not really that
            if ( ! attr && ! PyRootType_CheckExact( pyclass ) && PyType_Check( pyclass ) ) {
               PyErr_Clear();
               PyObject* pycppname = PyObject_GetAttr( pyclass, PyStrings::gCppName );
               const char* cppname = PyROOT_PyUnicode_AsString(pycppname);
               Py_DECREF(pycppname);
               Cppyy::TCppScope_t scope = Cppyy::GetScope( cppname );
               TClass* klass = TClass::GetClass( cppname );
               if ( Cppyy::IsNamespace( scope ) ) {
               // tickle lazy lookup of functions
                  if ( ! attr ) {
                     TObject *methObj = nullptr;
                     if ( (methObj = klass->GetListOfMethods()->FindObject(name.c_str())) ) {
                        // function exists, now collect overloads
                        auto completeName = methObj->GetName();
                        std::vector< PyCallable* > overloads;
                        const size_t nmeth = Cppyy::GetNumMethods( scope );
                        for ( size_t imeth = 0; imeth < nmeth; ++imeth ) {
                           Cppyy::TCppMethod_t method = Cppyy::GetMethod( scope, imeth );
                           auto currentName = Cppyy::GetMethodName(method);
                           // We need to compare with a final complete name, where:
                           // - Typedefs are resolved (e.g. Float_t -> float)
                           // - Namespaces are added (e.g. vector<float> -> std::vector<float>
                           // - All template types are included (e.g. std::vector<float> ->
                           //   std::vector<float, std::allocator<float>>)
                           if (currentName == completeName)
                              overloads.push_back( new TFunctionHolder( scope, method ) );
                        }

                     // Note: can't re-use Utility::AddClass here, as there's the risk of
                     // a recursive call. Simply add method directly, as we're guaranteed
                     // that it doesn't exist yet.
                        attr = (PyObject*)MethodProxy_New( name.c_str(), overloads );
                     }
                  }

               // tickle lazy lookup of data members
                  if ( ! attr ) {
                      Cppyy::TCppIndex_t dmi = Cppyy::GetDatamemberIndex( scope, name );
                      if ( 0 <= dmi ) attr = (PyObject*)PropertyProxy_New( scope, dmi );
                  }
               }

            // function templates that have not been instantiated
               if ( ! attr && klass ) {
                  TFunctionTemplate* tmpl = klass->GetFunctionTemplate( name.c_str() );
                  if ( tmpl )
                     attr = (PyObject*)TemplateProxy_New( name, pyclass );
               }

            // enums types requested as type (rather than the constants)
               if (!attr) {
                  if (Cppyy::IsEnum(Cppyy::GetScopedFinalName(scope)+"::"+name)) {
                     // enum types (incl. named and class enums)
                     Cppyy::TCppEnum_t enumtype = Cppyy::GetEnum(scope, name);
                     if (enumtype) {
                        // collect the enum values
                        Cppyy::TCppIndex_t ndata = Cppyy::GetNumEnumData(enumtype);
                        PyObject* dct = PyDict_New();
                        for (Cppyy::TCppIndex_t idata = 0; idata < ndata; ++idata) {
                           PyObject* val = PyLong_FromLongLong(Cppyy::GetEnumDataValue(enumtype, idata));
                           PyDict_SetItemString(dct, Cppyy::GetEnumDataName(enumtype, idata).c_str(), val);
                           Py_DECREF(val);
                        }

                        // add the __cppname__ for templates
                        PyObject* cppnamepy = PyROOT_PyUnicode_FromString((Cppyy::GetScopedFinalName(scope)+"::"+name).c_str());
                        PyDict_SetItem(dct, PyStrings::gCppName, cppnamepy);
                         // add also __cpp_name__ for forward compatibility
                        PyDict_SetItem(dct, PyStrings::gCppNameNew, cppnamepy);
                        Py_DECREF(cppnamepy);

                        // create new type with labeled values in place
                        PyObject* pybases = PyTuple_New(1);
                        Py_INCREF(&PyInt_Type);
                        PyTuple_SET_ITEM(pybases, 0, (PyObject*)&PyInt_Type);
                        PyObject* args = Py_BuildValue((char*)"sOO", name.c_str(), pybases, dct);
                        attr = Py_TYPE(&PyInt_Type)->tp_new(Py_TYPE(&PyInt_Type), args, nullptr);
                        Py_DECREF(args);
                        Py_DECREF(pybases);
                        Py_DECREF(dct);
                     } else {
                        // presumably not a class enum; simply pretend int
                        Py_INCREF(&PyInt_Type);
                        attr = (PyObject*)&PyInt_Type;
                     }
                  }
               }

               if ( attr ) {
                  PyObject_SetAttr( pyclass, pyname, attr );
                  Py_DECREF( attr );
                  attr = PyType_Type.tp_getattro( pyclass, pyname );
               }
            }

            if ( ! attr && ! PyRootType_Check( pyclass ) /* at global or module-level only */ ) {
               PyErr_Clear();
            // get class name to look up CINT tag info ...
               attr = GetCppGlobal( name /*, tag */ );
               if ( PropertyProxy_Check( attr ) ) {
                  PyObject_SetAttr( (PyObject*)Py_TYPE(pyclass), pyname, attr );
                  Py_DECREF( attr );
                  attr = PyType_Type.tp_getattro( pyclass, pyname );
               } else if ( attr )
                  PyObject_SetAttr( pyclass, pyname, attr );
            }

         }

      // if failed, then the original error is likely to be more instructive
         if ( ! attr && etype )
            PyErr_Restore( etype, value, trace );
         else if ( ! attr ) {
            PyObject* sklass = PyObject_Str( pyclass );
            PyErr_Format( PyExc_AttributeError, "%s has no attribute \'%s\'",
               PyROOT_PyUnicode_AsString( sklass ), PyROOT_PyUnicode_AsString( pyname ) );
            Py_DECREF( sklass );
         }

      // attribute is cached, if found
      }

      return attr;
   }

} // unnamed namespace


//= PyROOT object proxy type type ============================================
PyTypeObject PyRootType_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.PyRootType",  // tp_name
   sizeof(PyROOT::PyRootClass),// tp_basicsize
   0,                         // tp_itemsize
   0,                         // tp_dealloc
   0,                         // tp_print (python < 3.8)
                              // tp_vectorcall_offset (python >= 3.8)
   0,                         // tp_getattr
   0,                         // tp_setattr
   0,                         // tp_compare
   0,                         // tp_repr
   0,                         // tp_as_number
   0,                         // tp_as_sequence
   0,                         // tp_as_mapping
   0,                         // tp_hash
   0,                         // tp_call
   0,                         // tp_str
   (getattrofunc)pt_getattro, // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,     // tp_flags
   (char*)"PyROOT metatype (internal)",          // tp_doc
   0,                         // tp_traverse
   0,                         // tp_clear
   0,                         // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   0,                         // tp_getset
   &PyType_Type,              // tp_base
   0,                         // tp_dict
   0,                         // tp_descr_get
   0,                         // tp_descr_set
   0,                         // tp_dictoffset
   0,                         // tp_init
   0,                         // tp_alloc
   (newfunc)pt_new,           // tp_new
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
