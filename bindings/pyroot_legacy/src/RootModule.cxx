// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

// Bindings
#include "PyROOT.h"
#include "PyStrings.h"
#include "PyRootType.h"
#include "ObjectProxy.h"
#include "MethodProxy.h"
#include "TemplateProxy.h"
#include "PropertyProxy.h"
#include "TPyBufferFactory.h"
#include "TCustomPyTypes.h"
#include "TTupleOfInstances.h"
#include "RootWrapper.h"
#include "TCallContext.h"
#include "Utility.h"

// ROOT
#include "TObject.h"         // for FindObject
#include "TROOT.h"           // for ProcessLine and FindObject


#include "TBufferFile.h"     // for pickling

// Standard
#include <string>
#include <sstream>
#include <utility>
#include <vector>


//- from Python's dictobject.c -------------------------------------------------
#if PY_VERSION_HEX >= 0x03030000
   typedef struct PyDictKeyEntry {
      /* Cached hash code of me_key. */
      Py_hash_t me_hash;
      PyObject *me_key;
      PyObject *me_value; /* This field is only meaningful for combined tables */
   } PyDictEntry;

   typedef struct _dictkeysobject {
      Py_ssize_t dk_refcnt;
      Py_ssize_t dk_size;
      dict_lookup_func dk_lookup;
      Py_ssize_t dk_usable;
      PyDictKeyEntry dk_entries[1];
   } PyDictKeysObject;

#define PYROOT_GET_DICT_LOOKUP( mp )\
   ((dict_lookup_func&)mp->ma_keys->dk_lookup)

#else

#define PYROOT_GET_DICT_LOOKUP( mp )\
   ((dict_lookup_func&)mp->ma_lookup)

#endif

//- data -----------------------------------------------------------------------
static PyObject* nullptr_repr( PyObject* )
{
   return PyBytes_FromString( "nullptr" );
}

static void nullptr_dealloc( PyObject* )
{
   Py_FatalError( "deallocating nullptr" );
}

static int nullptr_nonzero( PyObject* )
{
   return 0;
}

static PyNumberMethods nullptr_as_number = {
   0, 0, 0,
#if PY_VERSION_HEX < 0x03000000
   0,
#endif
   0, 0, 0, 0, 0, 0,
   (inquiry)nullptr_nonzero,          // tp_nonzero (nb_bool in p3)
   0, 0, 0, 0, 0, 0,
#if PY_VERSION_HEX < 0x03000000
   0,                                 // nb_coerce
#endif
   0, 0, 0,
#if PY_VERSION_HEX < 0x03000000
   0, 0,
#endif
   0, 0, 0,
#if PY_VERSION_HEX < 0x03000000
   0,                                 // nb_inplace_divide
#endif
   0, 0, 0, 0, 0, 0, 0
#if PY_VERSION_HEX >= 0x02020000
   , 0                                // nb_floor_divide
#if PY_VERSION_HEX < 0x03000000
   , 0                                // nb_true_divide
#else
   , 0                                // nb_true_divide
#endif
   , 0, 0
#endif
#if PY_VERSION_HEX >= 0x02050000
   , 0                                // nb_index
#endif
#if PY_VERSION_HEX >= 0x03050000
   , 0                                // nb_matrix_multiply
   , 0                                // nb_inplace_matrix_multiply
#endif

   };

static PyTypeObject PyNullPtr_t_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   "nullptr_t",        // tp_name
   sizeof(PyObject),   // tp_basicsize
   0,                  // tp_itemsize
   nullptr_dealloc,    // tp_dealloc (never called)
   0, 0, 0, 0,
   nullptr_repr,       // tp_repr
   &nullptr_as_number, // tp_as_number
   0, 0,
   (hashfunc)_Py_HashPointer, // tp_hash
   0, 0, 0, 0, 0, Py_TPFLAGS_DEFAULT, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
#if PY_VERSION_HEX >= 0x02030000
   , 0                 // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
   , 0                 // tp_version_tag
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

PyObject _PyROOT_NullPtrStruct = {
  _PyObject_EXTRA_INIT
  1, &PyNullPtr_t_Type
};

namespace PyROOT {
   PyObject* gRootModule = 0;
   PyObject* gNullPtrObject = 0;
   std::vector<std::pair<Cppyy::TCppType_t, Cppyy::TCppType_t> > gPinnedTypes;
   std::vector<Cppyy::TCppType_t> gIgnorePinnings;
}


//- private helpers ------------------------------------------------------------
namespace {

   using namespace PyROOT;

////////////////////////////////////////////////////////////////////////////////

   PyObject* RootModuleResetCallback( PyObject*, PyObject* )
   {
      gRootModule = 0;   // reference was borrowed
      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Find a match within the ROOT module for something with name 'pyname'.

   PyObject* LookupCppEntity( PyObject* pyname, PyObject* args )
   {
      const char* cname = 0; long macro_ok = 0;
      if ( pyname && PyROOT_PyUnicode_CheckExact( pyname ) )
         cname = PyROOT_PyUnicode_AsString( pyname );
      else if ( ! ( args && PyArg_ParseTuple( args, const_cast< char* >( "s|l" ), &cname, &macro_ok ) ) )
         return 0;

   // we may have been destroyed if this code is called during shutdown
      if ( !gRootModule ) {
         PyErr_Format( PyExc_AttributeError, "%s", cname );
         return 0;
      }

      std::string name = cname;

   // block search for privates
      if ( name.size() <= 2 || name.substr( 0, 2 ) != "__" ) {
      // 1st attempt: look in myself
         PyObject* attr = PyObject_GetAttrString( gRootModule, const_cast< char* >( cname ) );
         if ( attr != 0 )
            return attr;

      // 2nd attempt: construct name as a class
         PyErr_Clear();
         attr = CreateScopeProxy( name, 0 /* parent */);
         if ( attr != 0 )
            return attr;

      // 3rd attempt: lookup name as global variable
         PyErr_Clear();
         attr = GetCppGlobal( name );
         if ( attr != 0 )
            return attr;

      // 4th attempt: find existing object (e.g. from file)
         PyErr_Clear();
         TObject* object = gROOT->FindObject( name.c_str() );
         if ( object != 0 )
            return BindCppObject( object, object->IsA()->GetName() );

      // 5th attempt: global enum
         if (Cppyy::IsEnum(name)) {
            // enum types (incl. named and class enums)
            Cppyy::TCppEnum_t enumtype = Cppyy::GetEnum(Cppyy::gGlobalScope, name);
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
               PyObject* cppnamepy = PyROOT_PyUnicode_FromString(cname);
               PyDict_SetItem(dct, PyStrings::gCppName, cppnamepy);
               // add also __cpp_name__ for forward compatibility
               PyDict_SetItem(dct, PyStrings::gCppNameNew, cppnamepy);
               Py_DECREF(cppnamepy);

               // create new type with labeled values in place
               PyObject* pybases = PyTuple_New(1);
               Py_INCREF(&PyInt_Type);
               PyTuple_SET_ITEM(pybases, 0, (PyObject*)&PyInt_Type);
               PyObject* argsnt = Py_BuildValue((char*)"sOO", name.c_str(), pybases, dct);
               attr = Py_TYPE(&PyInt_Type)->tp_new(Py_TYPE(&PyInt_Type), argsnt, nullptr);
               Py_DECREF(argsnt);
               Py_DECREF(pybases);
               Py_DECREF(dct);
            } else {
               // presumably not a class enum; simply pretend int
               Py_INCREF(&PyInt_Type);
               attr = (PyObject*)&PyInt_Type;
            }
            return attr;
         }

      // 6th attempt: check macro's (debatable, but this worked in CINT)
         if ( macro_ok ) {
            PyErr_Clear();
            std::ostringstream ismacro;
            ismacro << "#ifdef " << name << "\n_pyroot_" << name << "=" << name
                    << ";true;\n#else\nfalse;\n#endif";
            if ( gROOT->ProcessLine( ismacro.str().c_str() ) ) {
            // can now retrieve this as a global
               attr = GetCppGlobal( "_pyroot_"+name );
               if ( attr != 0 )
                  return attr;
            }
         }
      }

   // still here? raise attribute error
      PyErr_Format( PyExc_AttributeError, "%s", name.c_str() );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////

#if PY_VERSION_HEX >= 0x03030000
   inline PyDictKeyEntry* OrgDictLookup(
         PyDictObject* mp, PyObject* key, Py_hash_t hash, PyObject*** value_addr )
   {
      return (*gDictLookupOrg)( mp, key, hash, value_addr );
   }

#define PYROOT_ORGDICT_LOOKUP( mp, key, hash, value_addr )\
   OrgDictLookup( mp, key, hash, value_addr )

   PyDictKeyEntry* RootLookDictString(
         PyDictObject* mp, PyObject* key, Py_hash_t hash, PyObject*** value_addr )
#else
   inline PyDictEntry* OrgDictLookup( PyDictObject* mp, PyObject* key, Long_t hash )
   {
       return (*gDictLookupOrg)( mp, key, hash );
   }

#define PYROOT_ORGDICT_LOOKUP( mp, key, hash, value_addr )\
   OrgDictLookup( mp, key, hash )

   PyDictEntry* RootLookDictString( PyDictObject* mp, PyObject* key, Long_t hash )
#endif
   {
   // first search dictionary itself
      PyDictEntry* ep = PYROOT_ORGDICT_LOOKUP( mp, key, hash, value_addr );
      if ( ! ep || (ep->me_key && ep->me_value) || gDictLookupActive )
         return ep;

   // filter for builtins
      if ( PyDict_GetItem( PyEval_GetBuiltins(), key ) != 0 ) {
         return ep;
      }

   // all failed, start calling into ROOT
      gDictLookupActive = kTRUE;

   // ROOT globals (the round-about lookup is to prevent recursion)
      PyObject* gval = PyDict_GetItem( PyModule_GetDict( gRootModule ), key );
      if ( gval ) {
         Py_INCREF( gval );
         ep->me_value = gval;
         ep->me_key   = key;
         ep->me_hash  = hash;
#if PY_VERSION_HEX >= 0x03030000
         *value_addr  = &gval;
#endif
         gDictLookupActive = kFALSE;
         return ep;
      }

   // attempt to get ROOT enum/global/class
      PyObject* val = LookupCppEntity( key, 0 );

      if ( val != 0 ) {
      // success ...

         if ( PropertyProxy_CheckExact( val ) ) {
         // don't want to add to dictionary (the proper place would be the
         // dictionary of the (meta)class), but modifying ep will be noticed no
         // matter what; just return the actual value and live with the copy in
         // the dictionary (mostly, this is correct)
            PyObject* actual_val = Py_TYPE(val)->tp_descr_get( val, NULL, NULL );
            Py_DECREF( val );
            val = actual_val;
         }

      // add reference to ROOT entity in the given dictionary
         PYROOT_GET_DICT_LOOKUP( mp ) = gDictLookupOrg;     // prevent recursion
         if ( PyDict_SetItem( (PyObject*)mp, key, val ) == 0 ) {
            ep = PYROOT_ORGDICT_LOOKUP( mp, key, hash, value_addr );
         } else {
            ep->me_key   = 0;
            ep->me_value = 0;
         }
         PYROOT_GET_DICT_LOOKUP( mp ) = RootLookDictString; // restore

      // done with val
         Py_DECREF( val );
      } else
         PyErr_Clear();

#if PY_VERSION_HEX >= 0x03030000
      if ( mp->ma_keys->dk_usable <= 0 ) {
      // big risk that this lookup will result in a resize, so force it here
      // to be able to reset the lookup function; of course, this is nowhere
      // near fool-proof, but should cover interactive usage ...
         PYROOT_GET_DICT_LOOKUP( mp ) = gDictLookupOrg;
         const int maxinsert = 5;
         PyObject* buf[maxinsert];
         for ( int varmax = 1; varmax <= maxinsert; ++varmax ) {
            for ( int ivar = 0; ivar < varmax; ++ivar ) {
               buf[ivar] = PyROOT_PyUnicode_FromFormat( "__ROOT_FORCE_RESIZE_%d", ivar );
               PyDict_SetItem( (PyObject*)mp, buf[ivar], Py_None);
            }
            for ( int ivar = 0; ivar < varmax; ++ivar ) {
               PyDict_DelItem( (PyObject*)mp, buf[ivar] );
               Py_DECREF( buf[ivar] );
            }
            if ( 0 < mp->ma_keys->dk_usable )
               break;
         }

      // make sure the entry pointer is still valid by re-doing the lookup
         ep = PYROOT_ORGDICT_LOOKUP( mp, key, hash, value_addr );

      // full reset of all lookup functions
         gDictLookupOrg = PYROOT_GET_DICT_LOOKUP( mp );
         PYROOT_GET_DICT_LOOKUP( mp ) = RootLookDictString; // restore
      }
#endif

   // stopped calling into ROOT
      gDictLookupActive = kFALSE;

      return ep;
   }

////////////////////////////////////////////////////////////////////////////////
/// Modify the given dictionary to install the lookup function that also
/// tries the ROOT namespace before failing. Called on a module's dictionary,
/// this allows for lazy lookups.

   PyObject* SetRootLazyLookup( PyObject*, PyObject* args )
   {
      PyDictObject* dict = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!" ), &PyDict_Type, &dict ) )
         return 0;

   // Notwithstanding the code changes, the following does not work for p3.3 and
   // later: once the dictionary is resized for anything other than an insert (see
   // hack in RootLookDictString), its lookup function on its keys will revert to
   // the default (lookdict_unicode_nodummy) and only if the resizing dictionary
   // has the generic lookdict function as dk_lookup for its keys, will this be
   // set on the new keys.
      PYROOT_GET_DICT_LOOKUP( dict ) = RootLookDictString;

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Create a binding for a templated class instantiation.

   PyObject* MakeRootTemplateClass( PyObject*, PyObject* args )
   {
   // args is class name + template arguments; build full instantiation
      Py_ssize_t nArgs = PyTuple_GET_SIZE( args );
      if ( nArgs < 2 ) {
         PyErr_Format( PyExc_TypeError, "too few arguments for template instantiation" );
         return 0;
      }

   // build "< type, type, ... >" part of class name (modifies pyname)
      PyObject* pyname = Utility::BuildTemplateName( PyTuple_GET_ITEM( args, 0 ), args, 1 );
      if ( ! pyname )
         return 0;

      std::string name = PyROOT_PyUnicode_AsString( pyname );
      Py_DECREF( pyname );

      return CreateScopeProxy( name );
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper to get the address (address-of-address) of various object proxy types.

   void* GetObjectProxyAddress( PyObject*, PyObject* args )
   {
      ObjectProxy* pyobj = 0;
      PyObject* pyname = 0;
      if ( PyArg_ParseTuple( args, const_cast< char* >( "O|O!" ), &pyobj,
             &PyROOT_PyUnicode_Type, &pyname ) &&
           ObjectProxy_Check( pyobj ) && pyobj->fObject ) {

         if ( pyname != 0 ) {
         // locate property proxy for offset info
            PropertyProxy* pyprop = 0;

            PyObject* pyclass = PyObject_GetAttr( (PyObject*)pyobj, PyStrings::gClass );

            if ( pyclass ) {
               PyObject* dict = PyObject_GetAttr( pyclass, PyStrings::gDict );
               pyprop = (PropertyProxy*)PyObject_GetItem( dict, pyname );
               Py_DECREF( dict );
            }
            Py_XDECREF( pyclass );

            if ( PropertyProxy_Check( pyprop ) ) {
            // this is an address of a value (i.e. &myobj->prop)
               void* addr = (void*)pyprop->GetAddress( pyobj );
               Py_DECREF( pyprop );
               return addr;
            }

            Py_XDECREF( pyprop );

            PyErr_Format( PyExc_TypeError,
               "%s is not a valid data member", PyROOT_PyUnicode_AsString( pyname ) );
            return 0;
         }

      // this is an address of an address (i.e. &myobj, with myobj of type MyObj*)
         return (void*)&pyobj->fObject;
      }

      PyErr_SetString( PyExc_ValueError, "invalid argument for AddressOf()" );
      return 0;
   }

   PyObject* _addressof_common( PyObject* dummy ) {
      if ( dummy == Py_None || dummy == gNullPtrObject ) {
         Py_INCREF( gNullPtrObject );
         return gNullPtrObject;
      }
      if ( !PyErr_Occurred() ) {
         PyObject* str = PyObject_Str( dummy );
         if ( str && PyROOT_PyUnicode_Check( str ) )
            PyErr_Format( PyExc_ValueError, "unknown object %s", PyBytes_AS_STRING( str ) );
         else
            PyErr_Format( PyExc_ValueError, "unknown object at %p", (void*)dummy );
         Py_XDECREF( str );
      }
      return 0;
   }

   PyObject* AddressOf( PyObject* dummy, PyObject* args )
   {
   // Return object proxy address as an indexable buffer.
      void* addr = GetObjectProxyAddress( dummy, args );
      if ( addr )
         return BufFac_t::Instance()->PyBuffer_FromMemory( (Long_t*)addr, sizeof(Long_t) );
      if ( ! addr && PyTuple_Size( args ) ) {
         Utility::GetBuffer( PyTuple_GetItem( args, 0 ), '*', 1, addr, kFALSE );
         if ( addr )
            return BufFac_t::Instance()->PyBuffer_FromMemory( (Long_t*)&addr, sizeof(Long_t) );
      }
      return 0;//_addressof_common( dummy );
   }

   PyObject* addressof( PyObject* dummy, PyObject* args )
   {
   // Return object proxy address as a value (cppyy-style), or the same for an array.
      void* addr = GetObjectProxyAddress( dummy, args );
      if ( addr )
         return PyLong_FromLong( *(Long_t*)addr );
      else if ( PyTuple_Size( args ) ) {
         PyErr_Clear();
         Utility::GetBuffer( PyTuple_GetItem( args, 0 ), '*', 1, addr, kFALSE );
         if ( addr ) return PyLong_FromLong( (Long_t)addr );
      }
      return _addressof_common( dummy );
   }

   PyObject* AsCObject( PyObject* dummy, PyObject* args )
   {
   // Return object proxy as an opaque CObject.
      void* addr = GetObjectProxyAddress( dummy, args );
      if ( addr )
         return PyROOT_PyCapsule_New( (void*)(*(Long_t*)addr), NULL, NULL );

      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper to factorize the common code between MakeNullPointer and BindObject.

   PyObject* BindObject_( void* addr, PyObject* pyname )
   {
      if ( ! PyROOT_PyUnicode_Check( pyname ) ) {     // name given as string
         PyObject* nattr = PyObject_GetAttr( pyname, PyStrings::gCppName );
         if ( ! nattr ) nattr = PyObject_GetAttr( pyname, PyStrings::gName );
         if ( nattr )                        // object is actually a class
            pyname = nattr;
         pyname = PyObject_Str( pyname );
         Py_XDECREF( nattr );
      } else {
         Py_INCREF( pyname );
      }

      Cppyy::TCppType_t klass = (Cppyy::TCppType_t)Cppyy::GetScope( PyROOT_PyUnicode_AsString( pyname ) );
      Py_DECREF( pyname );

      if ( ! klass ) {
         PyErr_SetString( PyExc_TypeError,
            "BindObject expects a valid class or class name as an argument" );
         return 0;
      }

      return BindCppObjectNoCast( addr, klass, kFALSE );
   }

////////////////////////////////////////////////////////////////////////////////
/// From a long representing an address or a PyCapsule/CObject, bind to a class.

   PyObject* BindObject( PyObject*, PyObject* args )
   {
      Py_ssize_t argc = PyTuple_GET_SIZE( args );
      if ( argc != 2 ) {
         PyErr_Format( PyExc_TypeError,
           "BindObject takes exactly 2 argumenst (" PY_SSIZE_T_FORMAT " given)", argc );
         return 0;
      }

   // try to convert first argument: either PyCapsule/CObject or long integer
      PyObject* pyaddr = PyTuple_GET_ITEM( args, 0 );
      void* addr = PyROOT_PyCapsule_GetPointer( pyaddr, NULL );
      if ( PyErr_Occurred() ) {
         PyErr_Clear();

         addr = PyLong_AsVoidPtr( pyaddr );
         if ( PyErr_Occurred() ) {
            PyErr_Clear();

         // last chance, perhaps it's a buffer/array (return from void*)
            int buflen = Utility::GetBuffer( PyTuple_GetItem( args, 0 ), '*', 1, addr, kFALSE );
            if ( ! addr || ! buflen ) {
               PyErr_SetString( PyExc_TypeError,
                  "BindObject requires a CObject or long integer as first argument" );
               return 0;
            }
         }
      }

      return BindObject_( addr, PyTuple_GET_ITEM( args, 1 ) );
   }

////////////////////////////////////////////////////////////////////////////////
/// Create an object of the given type point to NULL (historic note: this
/// function is older than BindObject(), which can be used instead).

   PyObject* MakeNullPointer( PyObject*, PyObject* args )
   {
      Py_ssize_t argc = PyTuple_GET_SIZE( args );
      if ( argc != 0 && argc != 1 ) {
         PyErr_Format( PyExc_TypeError,
            "MakeNullPointer takes at most 1 argument (" PY_SSIZE_T_FORMAT " given)", argc );
         return 0;
      }

   // no class given, use None as generic
      if ( argc == 0 ) {
         Py_INCREF( Py_None );
         return Py_None;
      }

      return BindObject_( 0, PyTuple_GET_ITEM( args, 0 ) );
   }

////////////////////////////////////////////////////////////////////////////////
/// This method is a helper for (un)pickling of ObjectProxy instances.

   PyObject* ObjectProxyExpand( PyObject*, PyObject* args )
   {
      PyObject* pybuf = 0, *pyname = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!:__expand__" ),
               &PyBytes_Type, &pybuf, &PyBytes_Type, &pyname ) )
         return 0;

      const char* clname = PyBytes_AS_STRING(pyname);

   // make sure that ROOT.py is loaded and fully initialized by accessing on it
      PyObject* mod = PyImport_ImportModule( (char*)"ROOT" );
      if ( mod ) {
         PyObject* dummy = PyObject_GetAttrString( mod, (char*)"kRed" );
         Py_XDECREF( dummy );
         Py_DECREF( mod );
      }

   // TBuffer and its derived classes can't write themselves, but can be created
   // directly from the buffer, so handle them in a special case
      void* newObj = 0;
      if ( strcmp( clname, "TBufferFile" ) == 0 ) {
         TBufferFile* buf = new TBufferFile( TBuffer::kWrite );
         buf->WriteFastArray( PyBytes_AS_STRING(pybuf), PyBytes_GET_SIZE( pybuf ) );
         newObj = buf;
      } else {
      // use the PyString macro's to by-pass error checking; do not adopt the buffer,
      // as the local TBufferFile can go out of scope (there is no copying)
         TBufferFile buf( TBuffer::kRead,
            PyBytes_GET_SIZE( pybuf ), PyBytes_AS_STRING( pybuf ), kFALSE );
         newObj = buf.ReadObjectAny( 0 );
      }

      PyObject* result = BindCppObject( newObj, clname );
      if ( result ) {
      // this object is to be owned by the interpreter, assuming that the call
      // originated from there
         ((ObjectProxy*)result)->HoldOn();
      }

      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Set the global memory policy, which affects object ownership when objects
/// are passed as function arguments.

   PyObject* SetMemoryPolicy( PyObject*, PyObject* args )
   {
      PyObject* policy = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!" ), &PyInt_Type, &policy ) )
         return 0;

      Long_t l = PyInt_AS_LONG( policy );
      if ( TCallContext::SetMemoryPolicy( (TCallContext::ECallFlags)l ) ) {
         Py_INCREF( Py_None );
         return Py_None;
      }

      PyErr_Format( PyExc_ValueError, "Unknown policy %ld", l );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Set the global signal policy, which determines whether a jmp address
/// should be saved to return to after a C++ segfault.

   PyObject* SetSignalPolicy( PyObject*, PyObject* args )
   {
      PyObject* policy = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!" ), &PyInt_Type, &policy ) )
         return 0;

      Long_t l = PyInt_AS_LONG( policy );
      if ( TCallContext::SetSignalPolicy( (TCallContext::ECallFlags)l ) ) {
         Py_INCREF( Py_None );
         return Py_None;
      }

      PyErr_Format( PyExc_ValueError, "Unknown policy %ld", l );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Set the ownership (True is python-owns) for the given object.

   PyObject* SetOwnership( PyObject*, PyObject* args )
   {
      ObjectProxy* pyobj = 0; PyObject* pykeep = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!" ),
                &ObjectProxy_Type, (void*)&pyobj, &PyInt_Type, &pykeep ) )
         return 0;

      (Bool_t)PyLong_AsLong( pykeep ) ? pyobj->HoldOn() : pyobj->Release();

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Add a smart pointer to the list of known smart pointer types.

   PyObject* AddSmartPtrType( PyObject*, PyObject* args )
   {
      const char* type_name;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "s" ), &type_name ) )
         return nullptr;

      Cppyy::AddSmartPtrType( type_name );

      Py_RETURN_NONE;
   }


////////////////////////////////////////////////////////////////////////////////
/// Add a pinning so that objects of type `derived' are interpreted as
/// objects of type `base'.

   PyObject* SetTypePinning( PyObject*, PyObject* args )
   {
      PyRootClass* derived = nullptr, *base = nullptr;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!" ),
                               &PyRootType_Type, &derived,
                               &PyRootType_Type, &base ) )
         return nullptr;
      gPinnedTypes.push_back( std::make_pair( derived->fCppType, base->fCppType ) );

      Py_RETURN_NONE;
   }

////////////////////////////////////////////////////////////////////////////////
/// Add an exception to the type pinning for objects of type `derived'.

   PyObject* IgnoreTypePinning( PyObject*, PyObject* args )
   {
      PyRootClass* derived = nullptr;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!" ),
                               &PyRootType_Type, &derived ) )
         return nullptr;
      gIgnorePinnings.push_back( derived->fCppType );

      Py_RETURN_NONE;
   }

////////////////////////////////////////////////////////////////////////////////
/// Cast `obj' to type `type'.

   PyObject* Cast( PyObject*, PyObject* args )
   {
      ObjectProxy* obj = nullptr;
      PyRootClass* type = nullptr;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!" ),
                               &ObjectProxy_Type, &obj,
                               &PyRootType_Type, &type ) )
         return nullptr;
      // TODO: this misses an offset calculation, and reference type must not
      // be cast ...
      return BindCppObjectNoCast( obj->GetObject(), type->fCppType,
                                  obj->fFlags & ObjectProxy::kIsReference );
   }

} // unnamed namespace


//- data -----------------------------------------------------------------------
static PyMethodDef gPyROOTMethods[] = {
   { (char*) "CreateScopeProxy", (PyCFunction)PyROOT::CreateScopeProxy,
     METH_VARARGS, (char*) "PyROOT internal function" },
   { (char*) "GetCppGlobal", (PyCFunction)PyROOT::GetCppGlobal,
     METH_VARARGS, (char*) "PyROOT internal function" },
   { (char*) "LookupCppEntity", (PyCFunction)LookupCppEntity,
     METH_VARARGS, (char*) "PyROOT internal function" },
   { (char*) "SetRootLazyLookup", (PyCFunction)SetRootLazyLookup,
     METH_VARARGS, (char*) "PyROOT internal function" },
   { (char*) "MakeRootTemplateClass", (PyCFunction)MakeRootTemplateClass,
     METH_VARARGS, (char*) "PyROOT internal function" },
   { (char*) "_DestroyPyStrings", (PyCFunction)PyROOT::DestroyPyStrings,
     METH_NOARGS, (char*) "PyROOT internal function" },
   { (char*) "_ResetRootModule", (PyCFunction)RootModuleResetCallback,
     METH_NOARGS, (char*) "PyROOT internal function" },
   { (char*) "ClearProxiedObjects", (PyCFunction)ClearProxiedObjects,
     METH_NOARGS, (char*) "PyROOT internal function" },
   { (char*) "AddressOf", (PyCFunction)AddressOf,
     METH_VARARGS, (char*) "Retrieve address of held object in a buffer" },
   { (char*) "addressof", (PyCFunction)addressof,
     METH_VARARGS, (char*) "Retrieve address of held object as a value" },
   { (char*) "AsCObject", (PyCFunction)AsCObject,
     METH_VARARGS, (char*) "Retrieve held object in a CObject" },
   { (char*) "as_cobject", (PyCFunction)AsCObject,
     METH_VARARGS, (char*) "Retrieve held object in a CObject" },
   { (char*) "BindObject", (PyCFunction)BindObject,
     METH_VARARGS, (char*) "Create an object of given type, from given address" },
   { (char*) "bind_object", (PyCFunction)BindObject,
     METH_VARARGS, (char*) "Create an object of given type, from given address" },
   { (char*) "MakeNullPointer", (PyCFunction)MakeNullPointer,
     METH_VARARGS, (char*) "Create a NULL pointer of the given type" },
   { (char*) "_ObjectProxy__expand__", (PyCFunction)ObjectProxyExpand,
     METH_VARARGS, (char*) "Helper method for pickling" },
   { (char*) "SetMemoryPolicy", (PyCFunction)SetMemoryPolicy,
     METH_VARARGS, (char*) "Determines object ownership model" },
   { (char*) "SetSignalPolicy", (PyCFunction)SetSignalPolicy,
     METH_VARARGS, (char*) "Trap signals in safe mode to prevent interpreter abort" },
   { (char*) "SetOwnership", (PyCFunction)SetOwnership,
     METH_VARARGS, (char*) "Modify held C++ object ownership" },
   { (char*) "AddSmartPtrType", (PyCFunction)AddSmartPtrType,
     METH_VARARGS, (char*) "Add a smart pointer to the list of known smart pointer types" },
   { (char*) "InstallGUIEventInputHook", (PyCFunction)PyROOT::Utility::InstallGUIEventInputHook,
     METH_NOARGS, (char*) "Install input hook to sent GUI events" },
   { (char*) "RemoveGUIEventInputHook", (PyCFunction)PyROOT::Utility::RemoveGUIEventInputHook,
     METH_NOARGS, (char*) "Remove input hook to sent GUI events" },
   { (char*) "SetTypePinning", (PyCFunction)SetTypePinning,
     METH_VARARGS, (char*) "Install a type pinning" },
   { (char*) "IgnoreTypePinning", (PyCFunction)IgnoreTypePinning,
     METH_VARARGS, (char*) "Don't pin the given type" },
   { (char*) "Cast", (PyCFunction)Cast,
     METH_VARARGS, (char*) "Cast the given object to the given type" },
   { NULL, NULL, 0, NULL }
};


#if PY_VERSION_HEX >= 0x03000000
struct module_state {
    PyObject *error;
};

#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))

static int rootmodule_traverse( PyObject* m, visitproc visit, void* arg )
{
    Py_VISIT( GETSTATE( m )->error );
    return 0;
}

static int rootmodule_clear( PyObject* m )
{
    Py_CLEAR( GETSTATE( m )->error );
    return 0;
}


static struct PyModuleDef moduledef = {
   PyModuleDef_HEAD_INIT,
   "libPyROOT",
   NULL,
   sizeof(struct module_state),
   gPyROOTMethods,
   NULL,
   rootmodule_traverse,
   rootmodule_clear,
   NULL
};

////////////////////////////////////////////////////////////////////////////////
/// Initialization of extension module libPyROOT.

#define PYROOT_INIT_ERROR return NULL
extern "C" PyObject* PyInit_libPyROOT()
#else
#define PYROOT_INIT_ERROR return
extern "C" void initlibPyROOT()
#endif
{
   using namespace PyROOT;

// load commonly used python strings
   if ( ! PyROOT::CreatePyStrings() )
      PYROOT_INIT_ERROR;

// prepare for lazyness
   PyObject* dict = PyDict_New();
#if PY_VERSION_HEX >= 0x03030000
   gDictLookupOrg = (dict_lookup_func)((PyDictObject*)dict)->ma_keys->dk_lookup;
#else
   gDictLookupOrg = (dict_lookup_func)((PyDictObject*)dict)->ma_lookup;
#endif
   Py_DECREF( dict );

// setup PyROOT
#if PY_VERSION_HEX >= 0x03000000
   gRootModule = PyModule_Create( &moduledef );
#else
   gRootModule = Py_InitModule( const_cast< char* >( "libPyROOT" ), gPyROOTMethods );
#endif
   if ( ! gRootModule )
      PYROOT_INIT_ERROR;

// keep gRootModule, but do not increase its reference count even as it is borrowed,
// or a self-referencing cycle would be created

// Pythonizations ...
   PyObject* userPythonizations = PyDict_New();
   PyObject* gblList = PyList_New( 0 );
   PyDict_SetItemString( userPythonizations, "__global__", gblList );
   Py_DECREF( gblList );
   PyModule_AddObject( gRootModule, "UserPythonizations", userPythonizations );
   PyModule_AddObject( gRootModule, "UserExceptions",     PyDict_New() );
   PyModule_AddObject( gRootModule, "PythonizationScope", PyROOT_PyUnicode_FromString( "__global__" ) );

// inject meta type
   if ( ! Utility::InitProxy( gRootModule, &PyRootType_Type, "PyRootType" ) )
      PYROOT_INIT_ERROR;

// inject object proxy type
   if ( ! Utility::InitProxy( gRootModule, &ObjectProxy_Type, "ObjectProxy" ) )
      PYROOT_INIT_ERROR;

// inject method proxy type
   if ( ! Utility::InitProxy( gRootModule, &MethodProxy_Type, "MethodProxy" ) )
      PYROOT_INIT_ERROR;

// inject template proxy type
   if ( ! Utility::InitProxy( gRootModule, &TemplateProxy_Type, "TemplateProxy" ) )
      PYROOT_INIT_ERROR;

// inject property proxy type
   if ( ! Utility::InitProxy( gRootModule, &PropertyProxy_Type, "PropertyProxy" ) )
      PYROOT_INIT_ERROR;

// inject custom data types
   if ( ! Utility::InitProxy( gRootModule, &TCustomFloat_Type, "Double" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &TCustomInt_Type, "Long" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &TCustomFloat_Type, "double" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &TCustomInt_Type, "long" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &TCustomInstanceMethod_Type, "InstanceMethod" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &TTupleOfInstances_Type, "InstancesArray" ) )
      PYROOT_INIT_ERROR;

   if ( ! Utility::InitProxy( gRootModule, &PyNullPtr_t_Type, "nullptr_t" ) )
      PYROOT_INIT_ERROR;

// inject identifiable nullptr
   gNullPtrObject = (PyObject*)&_PyROOT_NullPtrStruct;
   Py_INCREF( gNullPtrObject );
   PyModule_AddObject( gRootModule, (char*)"nullptr", gNullPtrObject );

// policy labels
   PyModule_AddObject( gRootModule, (char*)"kMemoryHeuristics",
      PyInt_FromLong( (int)TCallContext::kUseHeuristics ) );
   PyModule_AddObject( gRootModule, (char*)"kMemoryStrict",
      PyInt_FromLong( (int)TCallContext::kUseStrict ) );
   PyModule_AddObject( gRootModule, (char*)"kSignalFast",
      PyInt_FromLong( (int)TCallContext::kFast ) );
   PyModule_AddObject( gRootModule, (char*)"kSignalSafe",
      PyInt_FromLong( (int)TCallContext::kSafe ) );

// setup ROOT
   PyROOT::InitRoot();

// signal policy: don't abort interpreter in interactive mode
   TCallContext::SetSignalPolicy( gROOT->IsBatch() ? TCallContext::kFast : TCallContext::kSafe );

// inject ROOT namespace for convenience
   PyModule_AddObject( gRootModule, (char*)"ROOT", CreateScopeProxy( "ROOT" ) );

#if PY_VERSION_HEX >= 0x03000000
   Py_INCREF( gRootModule );
   return gRootModule;
#endif
}
