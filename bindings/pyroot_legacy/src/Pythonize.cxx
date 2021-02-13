// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jul 2004

// Bindings
#include "PyROOT.h"
#include "PyStrings.h"
#include "Pythonize.h"
#include "ObjectProxy.h"
#include "MethodProxy.h"
#include "RootWrapper.h"
#include "Utility.h"
#include "PyCallable.h"
#include "TPyBufferFactory.h"
#include "TFunctionHolder.h"
#include "Converters.h"
#include "TMemoryRegulator.h"
#include "Utility.h"

// ROOT
#include "TClass.h"
#include "TFunction.h"
#include "TInterpreter.h"
#include "TMethod.h"

#include "TClonesArray.h"
#include "TCollection.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSeqCollection.h"

#include "TTree.h"
#include "TBranch.h"
#include "TBranchElement.h"
#include "TBranchObject.h"
#include "TLeaf.h"
#include "TLeafElement.h"
#include "TLeafObject.h"
#include "TStreamerElement.h"
#include "TStreamerInfo.h"
#include "TInterpreterValue.h"

#include "ROOT/RVec.hxx"

// Standard
#include <stdexcept>
#include <string>
#include <utility>
#include <sstream>

#include <stdio.h>
#include <string.h>     // only needed for Cling TMinuit workaround


// temp (?)
static inline TClass* OP2TCLASS( PyROOT::ObjectProxy* pyobj ) {
   return TClass::GetClass( Cppyy::GetFinalName( pyobj->ObjectIsA() ).c_str());
}
//-- temp

//- data and local helpers ---------------------------------------------------
namespace PyROOT {
   R__EXTERN PyObject* gRootModule;
}

namespace {

// for convenience
   using namespace PyROOT;

////////////////////////////////////////////////////////////////////////////////
/// prevents calls to Py_TYPE(pyclass)->tp_getattr, which is unnecessary for our
/// purposes here and could tickle problems w/ spurious lookups into ROOT meta

   Bool_t HasAttrDirect( PyObject* pyclass, PyObject* pyname, Bool_t mustBePyROOT = kFALSE ) {
      PyObject* attr = PyType_Type.tp_getattro( pyclass, pyname );
      if ( attr != 0 && ( ! mustBePyROOT || MethodProxy_Check( attr ) ) ) {
         Py_DECREF( attr );
         return kTRUE;
      }

      PyErr_Clear();
      return kFALSE;
   }

////////////////////////////////////////////////////////////////////////////////
/// prevents calls to descriptors

   PyObject* PyObject_GetAttrFromDict( PyObject* pyclass, PyObject* pyname ) {
      PyObject* dict = PyObject_GetAttr( pyclass, PyStrings::gDict );
      PyObject* attr = PyObject_GetItem( dict, pyname );
      Py_DECREF( dict );
      return attr;
   }

////////////////////////////////////////////////////////////////////////////////
/// Scan the name of the class and determine whether it is a template instantiation.

   inline Bool_t IsTemplatedSTLClass( const std::string& name, const std::string& klass ) {
      const int nsize = (int)name.size();
      const int ksize = (int)klass.size();

      return ( ( ksize   < nsize && name.substr(0,ksize) == klass ) ||
               ( ksize+5 < nsize && name.substr(5,ksize) == klass ) ) &&
             name.find( "::", name.find( ">" ) ) == std::string::npos;
   }

// to prevent compiler warnings about const char* -> char*
   inline PyObject* CallPyObjMethod( PyObject* obj, const char* meth )
   {
   // Helper; call method with signature: obj->meth().
      Py_INCREF( obj );
      PyObject* result = PyObject_CallMethod( obj, const_cast< char* >( meth ), const_cast< char* >( "" ) );
      Py_DECREF( obj );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper; call method with signature: obj->meth( arg1 ).

   inline PyObject* CallPyObjMethod( PyObject* obj, const char* meth, PyObject* arg1 )
   {
      Py_INCREF( obj );
      PyObject* result = PyObject_CallMethod(
         obj, const_cast< char* >( meth ), const_cast< char* >( "O" ), arg1 );
      Py_DECREF( obj );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper; call method with signature: obj->meth( arg1, arg2 ).

   inline PyObject* CallPyObjMethod(
      PyObject* obj, const char* meth, PyObject* arg1, PyObject* arg2 )
   {
      Py_INCREF( obj );
      PyObject* result = PyObject_CallMethod(
         obj, const_cast< char* >( meth ), const_cast< char* >( "OO" ), arg1, arg2 );
      Py_DECREF( obj );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper; call method with signature: obj->meth( arg1, int ).

   inline PyObject* CallPyObjMethod( PyObject* obj, const char* meth, PyObject* arg1, int arg2 )
   {
      Py_INCREF( obj );
      PyObject* result = PyObject_CallMethod(
         obj, const_cast< char* >( meth ), const_cast< char* >( "Oi" ), arg1, arg2 );
      Py_DECREF( obj );
      return result;
   }


//- helpers --------------------------------------------------------------------
   PyObject* PyStyleIndex( PyObject* self, PyObject* index )
   {
   // Helper; converts python index into straight C index.
      Py_ssize_t idx = PyInt_AsSsize_t( index );
      if ( idx == (Py_ssize_t)-1 && PyErr_Occurred() )
         return 0;

      Py_ssize_t size = PySequence_Size( self );
      if ( idx >= size || ( idx < 0 && idx < -size ) ) {
         PyErr_SetString( PyExc_IndexError, "index out of range" );
         return 0;
      }

      PyObject* pyindex = 0;
      if ( idx >= 0 ) {
         Py_INCREF( index );
         pyindex = index;
      } else
         pyindex = PyLong_FromLong( size + idx );

      return pyindex;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper; call method with signature: meth( pyindex ).

   inline PyObject* CallSelfIndex( ObjectProxy* self, PyObject* idx, const char* meth )
   {
      Py_INCREF( (PyObject*)self );
      PyObject* pyindex = PyStyleIndex( (PyObject*)self, idx );
      if ( ! pyindex ) {
         Py_DECREF( (PyObject*)self );
         return 0;
      }

      PyObject* result = CallPyObjMethod( (PyObject*)self, meth, pyindex );
      Py_DECREF( pyindex );
      Py_DECREF( (PyObject*)self );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Helper; convert generic python object into a boolean value.

   inline PyObject* BoolNot( PyObject* value )
   {
      if ( PyObject_IsTrue( value ) == 1 ) {
         Py_INCREF( Py_False );
         Py_DECREF( value );
         return Py_False;
      } else {
         Py_INCREF( Py_True );
         Py_XDECREF( value );
         return Py_True;
      }
   }

//- "smart pointer" behavior ---------------------------------------------------
   PyObject* DeRefGetAttr( PyObject* self, PyObject* name )
   {
   // Follow operator*() if present (available in python as __deref__), so that
   // smart pointers behave as expected.
      if ( ! PyROOT_PyUnicode_Check( name ) )
         PyErr_SetString( PyExc_TypeError, "getattr(): attribute name must be string" );

      PyObject* pyptr = CallPyObjMethod( self, "__deref__" );
      if ( ! pyptr )
         return 0;

   // prevent a potential infinite loop
      if ( Py_TYPE(pyptr) == Py_TYPE(self) ) {
         PyObject* val1 = PyObject_Str( self );
         PyObject* val2 = PyObject_Str( name );
         PyErr_Format( PyExc_AttributeError, "%s has no attribute \'%s\'",
            PyROOT_PyUnicode_AsString( val1 ), PyROOT_PyUnicode_AsString( val2 ) );
         Py_DECREF( val2 );
         Py_DECREF( val1 );

         Py_DECREF( pyptr );
         return 0;
      }

      PyObject* result = PyObject_GetAttr( pyptr, name );
      Py_DECREF( pyptr );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Follow operator->() if present (available in python as __follow__), so that
/// smart pointers behave as expected.

   PyObject* FollowGetAttr( PyObject* self, PyObject* name )
   {
      if ( ! PyROOT_PyUnicode_Check( name ) )
         PyErr_SetString( PyExc_TypeError, "getattr(): attribute name must be string" );

      PyObject* pyptr = CallPyObjMethod( self, "__follow__" );
      if ( ! pyptr )
         return 0;

      PyObject* result = PyObject_GetAttr( pyptr, name );
      Py_DECREF( pyptr );
      return result;
   }

//- TObject behavior -----------------------------------------------------------
   PyObject* TObjectContains( PyObject* self, PyObject* obj )
   {
   // Implement python's __contains__ with TObject::FindObject.
      if ( ! ( ObjectProxy_Check( obj ) || PyROOT_PyUnicode_Check( obj ) ) )
         return PyInt_FromLong( 0l );

      PyObject* found = CallPyObjMethod( self, "FindObject", obj );
      PyObject* result = PyInt_FromLong( PyObject_IsTrue( found ) );
      Py_DECREF( found );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __cmp__ with TObject::Compare.

   PyObject* TObjectCompare( PyObject* self, PyObject* obj )
   {
      if ( ! ObjectProxy_Check( obj ) )
         return PyInt_FromLong( -1l );

      return CallPyObjMethod( self, "Compare", obj );
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __eq__ with TObject::IsEqual.

   PyObject* TObjectIsEqual( PyObject* self, PyObject* obj )
   {
      if ( ! ObjectProxy_Check( obj ) || ! ((ObjectProxy*)obj)->fObject )
         return ObjectProxy_Type.tp_richcompare( self, obj, Py_EQ );

      return CallPyObjMethod( self, "IsEqual", obj );
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __ne__ in terms of not TObject::IsEqual.

   PyObject* TObjectIsNotEqual( PyObject* self, PyObject* obj )
   {
      if ( ! ObjectProxy_Check( obj ) || ! ((ObjectProxy*)obj)->fObject )
         return ObjectProxy_Type.tp_richcompare( self, obj, Py_NE );

      return BoolNot( CallPyObjMethod( self, "IsEqual", obj ) );
   }

////////////////////////////////////////////////////////////////////////////////
/// Contrary to TObjectIsEqual, it can now not be relied upon that the only
/// non-ObjectProxy obj is None, as any operator==(), taking any object (e.g.
/// an enum) can be implemented. However, those cases will yield an exception
/// if presented with None.

   PyObject* GenObjectIsEqual( PyObject* self, PyObject* obj )
   {
      PyObject* result = CallPyObjMethod( self, "__cpp_eq__", obj );
      if ( ! result ) {
         PyErr_Clear();
         result = ObjectProxy_Type.tp_richcompare( self, obj, Py_EQ );
      }

      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Reverse of GenObjectIsEqual, if operator!= defined.

   PyObject* GenObjectIsNotEqual( PyObject* self, PyObject* obj )
   {
      PyObject* result = CallPyObjMethod( self, "__cpp_ne__", obj );
      if ( ! result ) {
         PyErr_Clear();
         result = ObjectProxy_Type.tp_richcompare( self, obj, Py_NE );
      }

      return result;
   }

//- TClass behavior ------------------------------------------------------------
   PyObject* TClassStaticCast( ObjectProxy* self, PyObject* args )
   {
   // Implemented somewhat different than TClass::DynamicClass, in that "up" is
   // chosen automatically based on the relationship between self and arg pyclass.
      ObjectProxy* pyclass = 0; PyObject* pyobject = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O:StaticCast" ),
              &ObjectProxy_Type, &pyclass, &pyobject ) )
         return 0;

   // check the given arguments (dcasts are necessary b/c of could be a TQClass
      TClass* from = (TClass*)OP2TCLASS(self)->DynamicCast( TClass::Class(), self->GetObject() );
      TClass* to   = (TClass*)OP2TCLASS(self)->DynamicCast( TClass::Class(), pyclass->GetObject() );

      if ( ! from ) {
         PyErr_SetString( PyExc_TypeError, "unbound method TClass::StaticCast "
            "must be called with a TClass instance as first argument" );
         return 0;
      }

      if ( ! to ) {
         PyErr_SetString( PyExc_TypeError, "could not convert argument 1 (TClass* expected)" );
         return 0;
      }

   // retrieve object address
      void* address = 0;
      if ( ObjectProxy_Check( pyobject ) ) address = ((ObjectProxy*)pyobject)->GetObject();
      else if ( PyInt_Check( pyobject ) || PyLong_Check( pyobject ) ) address = (void*)PyLong_AsLong( pyobject );
      else Utility::GetBuffer( pyobject, '*', 1, address, kFALSE );

      if ( ! address ) {
         PyErr_SetString( PyExc_TypeError, "could not convert argument 2 (void* expected)" );
         return 0;
      }

   // determine direction of cast
      int up = -1;
      if ( from->InheritsFrom( to ) ) up = 1;
      else if ( to->InheritsFrom( from ) ) {
         TClass* tmp = to; to = from; from = tmp;
         up = 0;
      }

      if ( up == -1 ) {
         PyErr_Format( PyExc_TypeError, "unable to cast %s to %s", from->GetName(), to->GetName() );
         return 0;
      }

   // perform actual cast
      void* result = from->DynamicCast( to, address, (Bool_t)up );

   // at this point, "result" can't be null (but is still safe if it is)
      return BindCppObjectNoCast( result, Cppyy::GetScope( to->GetName() ) );
   }

////////////////////////////////////////////////////////////////////////////////
/// TClass::DynamicCast returns a void* that the user still has to cast (it
/// will have the proper offset, though). Fix this by providing the requested
/// binding if the cast succeeded.

   PyObject* TClassDynamicCast( ObjectProxy* self, PyObject* args )
   {
      ObjectProxy* pyclass = 0; PyObject* pyobject = 0;
      Long_t up = 1;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O|l:DynamicCast" ),
              &ObjectProxy_Type, &pyclass, &pyobject, &up ) )
         return 0;

   // perform actual cast
      PyObject* meth = PyObject_GetAttr( (PyObject*)self, PyStrings::gTClassDynCast );
      PyObject* ptr = meth ? PyObject_Call( meth, args, 0 ) : 0;
      Py_XDECREF( meth );

   // simply forward in case of call failure
      if ( ! ptr )
         return ptr;

   // retrieve object address
      void* address = 0;
      if ( ObjectProxy_Check( pyobject ) ) address = ((ObjectProxy*)pyobject)->GetObject();
      else if ( PyInt_Check( pyobject ) || PyLong_Check( pyobject ) ) address = (void*)PyLong_AsLong( pyobject );
      else Utility::GetBuffer( pyobject, '*', 1, address, kFALSE );

      if ( PyErr_Occurred() ) {
         PyErr_Clear();
         return ptr;
      }

   // now use binding to return a usable class
      TClass* klass = 0;
      if ( up ) {                  // up-cast: result is a base
         klass = (TClass*)OP2TCLASS(pyclass)->DynamicCast( TClass::Class(), pyclass->GetObject() );
      } else {                     // down-cast: result is a derived
         klass = (TClass*)OP2TCLASS(self)->DynamicCast( TClass::Class(), self->GetObject() );
      }

      PyObject* result = BindCppObjectNoCast( (void*)address, Cppyy::GetScope( klass->GetName() ) );
      Py_DECREF( ptr );

      return result;
   }

//- TCollection behavior -------------------------------------------------------
   PyObject* TCollectionExtend( PyObject* self, PyObject* obj )
   {
   // Implement a python-style extend with TCollection::Add.
      for ( Py_ssize_t i = 0; i < PySequence_Size( obj ); ++i ) {
         PyObject* item = PySequence_GetItem( obj, i );
         PyObject* result = CallPyObjMethod( self, "Add", item );
         Py_XDECREF( result );
         Py_DECREF( item );
      }

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style remove with TCollection::Add.

   PyObject* TCollectionRemove( PyObject* self, PyObject* obj )
   {
      PyObject* result = CallPyObjMethod( self, "Remove", obj );
      if ( ! result )
         return 0;

      if ( ! PyObject_IsTrue( result ) ) {
         Py_DECREF( result );
         PyErr_SetString( PyExc_ValueError, "list.remove(x): x not in list" );
         return 0;
      }

      Py_DECREF( result );
      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __add__ with the pythonized extend for TCollections.

   PyObject* TCollectionAdd( PyObject* self, PyObject* other )
   {
      PyObject* l = CallPyObjMethod( self, "Clone" );
      if ( ! l )
         return 0;

      PyObject* result = CallPyObjMethod( l, "extend", other );
      if ( ! result ) {
         Py_DECREF( l );
         return 0;
      }

      return l;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __mul__ with the pythonized extend for TCollections.

   PyObject* TCollectionMul( ObjectProxy* self, PyObject* pymul )
   {
      Long_t imul = PyLong_AsLong( pymul );
      if ( imul == -1 && PyErr_Occurred() )
         return 0;

      if ( ! self->GetObject() ) {
         PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
         return 0;
      }

      PyObject* nseq = BindCppObject(
         Cppyy::Construct( self->ObjectIsA() ), self->ObjectIsA() );

      for ( Long_t i = 0; i < imul; ++i ) {
         PyObject* result = CallPyObjMethod( nseq, "extend", (PyObject*)self );
         Py_DECREF( result );
      }

      return nseq;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __imul__ with the pythonized extend for TCollections.

   PyObject* TCollectionIMul( PyObject* self, PyObject* pymul )
   {
      Long_t imul = PyLong_AsLong( pymul );
      if ( imul == -1 && PyErr_Occurred() )
         return 0;

      PyObject* l = PySequence_List( self );

      for ( Long_t i = 0; i < imul - 1; ++i ) {
         CallPyObjMethod( self, "extend", l );
      }

      Py_INCREF( self );
      return self;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style count for TCollections.

   PyObject* TCollectionCount( PyObject* self, PyObject* obj )
   {
      Py_ssize_t count = 0;
      for ( Py_ssize_t i = 0; i < PySequence_Size( self ); ++i ) {
         PyObject* item = PySequence_GetItem( self, i );
         PyObject* found = PyObject_RichCompare( item, obj, Py_EQ );

         Py_DECREF( item );

         if ( ! found )
            return 0;                        // internal problem

         if ( PyObject_IsTrue( found ) )
            count += 1;
         Py_DECREF( found );
      }

      return PyInt_FromSsize_t( count );
   }

////////////////////////////////////////////////////////////////////////////////
/// Python __iter__ protocol for TCollections.

   PyObject* TCollectionIter( ObjectProxy* self ) {
      if ( ! self->GetObject() ) {
         PyErr_SetString( PyExc_TypeError, "iteration over non-sequence" );
         return 0;
      }

      TCollection* col =
         (TCollection*)OP2TCLASS(self)->DynamicCast( TCollection::Class(), self->GetObject() );

      PyObject* pyobject = BindCppObject( (void*) new TIter( col ), "TIter" );
      ((ObjectProxy*)pyobject)->HoldOn();
      return pyobject;
   }


//- TSeqCollection behavior ----------------------------------------------------
   PyObject* TSeqCollectionGetItem( ObjectProxy* self, PySliceObject* index )
   {
   // Python-style indexing and size checking for getting objects from a TCollection.
      if ( PySlice_Check( index ) ) {
         if ( ! self->GetObject() ) {
            PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
            return 0;
         }

         TClass* clSeq = OP2TCLASS(self);
         TSeqCollection* oseq =
            (TSeqCollection*)clSeq->DynamicCast( TSeqCollection::Class(), self->GetObject() );
         TSeqCollection* nseq = (TSeqCollection*)clSeq->New();

         Py_ssize_t start, stop, step;
         PySlice_GetIndices( (PyROOT_PySliceCast)index, oseq->GetSize(), &start, &stop, &step );

         for ( Py_ssize_t i = start; i < stop; i += step ) {
            nseq->Add( oseq->At( (Int_t)i ) );
         }

         return BindCppObject( (void*) nseq, clSeq->GetName() );
      }

      return CallSelfIndex( self, (PyObject*)index, "At" );
   }

////////////////////////////////////////////////////////////////////////////////
/// Python-style indexing and size checking for setting objects in a TCollection.

   PyObject* TSeqCollectionSetItem( ObjectProxy* self, PyObject* args )
   {
      PyObject* index = 0, *obj = 0;
      if ( ! PyArg_ParseTuple( args,
                const_cast< char* >( "OO:__setitem__" ), &index, &obj ) )
         return 0;

      if ( PySlice_Check( index ) ) {
         if ( ! self->GetObject() ) {
            PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
            return 0;
         }

         TSeqCollection* oseq = (TSeqCollection*)OP2TCLASS(self)->DynamicCast(
            TSeqCollection::Class(), self->GetObject() );

         Py_ssize_t start, stop, step;
         PySlice_GetIndices( (PyROOT_PySliceCast)index, oseq->GetSize(), &start, &stop, &step );
         for ( Py_ssize_t i = stop - step; i >= start; i -= step ) {
            oseq->RemoveAt( (Int_t)i );
         }

         for ( Py_ssize_t i = 0; i < PySequence_Size( obj ); ++i ) {
            ObjectProxy* item = (ObjectProxy*)PySequence_GetItem( obj, i );
            item->Release();
            oseq->AddAt( (TObject*) item->GetObject(), (Int_t)(i + start) );
            Py_DECREF( item );
         }

         Py_INCREF( Py_None );
         return Py_None;
      }

      PyObject* pyindex = PyStyleIndex( (PyObject*)self, index );
      if ( ! pyindex )
         return 0;

      PyObject* result  = CallPyObjMethod( (PyObject*)self, "RemoveAt", pyindex );
      if ( ! result ) {
         Py_DECREF( pyindex );
         return 0;
      }

      Py_DECREF( result );
      result = CallPyObjMethod( (PyObject*)self, "AddAt", obj, pyindex );
      Py_DECREF( pyindex );
      return result;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement python's __del__ with TCollection::RemoveAt.

   PyObject* TSeqCollectionDelItem( ObjectProxy* self, PySliceObject* index )
   {
      if ( PySlice_Check( index ) ) {
         if ( ! self->GetObject() ) {
            PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
            return 0;
         }

         TSeqCollection* oseq = (TSeqCollection*)OP2TCLASS(self)->DynamicCast(
            TSeqCollection::Class(), self->GetObject() );

         Py_ssize_t start, stop, step;
         PySlice_GetIndices( (PyROOT_PySliceCast)index, oseq->GetSize(), &start, &stop, &step );
         for ( Py_ssize_t i = stop - step; i >= start; i -= step ) {
            oseq->RemoveAt( (Int_t)i );
         }

         Py_INCREF( Py_None );
         return Py_None;
      }

      PyObject* result = CallSelfIndex( self, (PyObject*)index, "RemoveAt" );
      if ( ! result )
         return 0;

      Py_DECREF( result );
      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Python-style insertion implemented with TCollection::AddAt.

   PyObject* TSeqCollectionInsert( PyObject* self, PyObject* args )
   {
      PyObject* obj = 0; Long_t idx = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "lO:insert" ), &idx, &obj ) )
         return 0;

      Py_ssize_t size = PySequence_Size( self );
      if ( idx < 0 )
         idx = 0;
      else if ( size < idx )
         idx = size;

      return CallPyObjMethod( self, "AddAt", obj, idx );
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style pop for TCollections.

   PyObject* TSeqCollectionPop( ObjectProxy* self, PyObject* args )
   {
      int nArgs = PyTuple_GET_SIZE( args );
      if ( nArgs == 0 ) {
      // create the default argument 'end of sequence'
         PyObject* index = PyInt_FromSsize_t( PySequence_Size( (PyObject*)self ) - 1 );
         PyObject* result = CallSelfIndex( self, index, "RemoveAt" );
         Py_DECREF( index );
         return result;
      } else if ( nArgs != 1 ) {
         PyErr_Format( PyExc_TypeError,
            "pop() takes at most 1 argument (%d given)", nArgs );
         return 0;
      }

      return CallSelfIndex( self, PyTuple_GET_ITEM( args, 0 ), "RemoveAt" );
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style reverse for TCollections.

   PyObject* TSeqCollectionReverse( PyObject* self )
   {
      PyObject* tup = PySequence_Tuple( self );
      if ( ! tup )
         return 0;

      PyObject* result = CallPyObjMethod( self, "Clear" );
      Py_XDECREF( result );

      for ( Py_ssize_t i = 0; i < PySequence_Size( tup ); ++i ) {
         PyObject* retval = CallPyObjMethod( self, "AddAt", PyTuple_GET_ITEM( tup, i ), 0 );
         Py_XDECREF( retval );
      }

      Py_INCREF( Py_None );
      return Py_None;
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style sort for TCollections.

   PyObject* TSeqCollectionSort( PyObject* self, PyObject* args, PyObject* kw )
   {
      if ( PyTuple_GET_SIZE( args ) == 0 && ! kw ) {
      // no specialized sort, use ROOT one
         return CallPyObjMethod( self, "Sort" );
      } else {
      // sort in a python list copy
         PyObject* l = PySequence_List( self );
         PyObject* result = 0;
         if ( PyTuple_GET_SIZE( args ) == 1 )
            result = CallPyObjMethod( l, "sort", PyTuple_GET_ITEM( args, 0 ) );
         else {
            PyObject* pymeth = PyObject_GetAttrString( l, const_cast< char* >( "sort" ) );
            result = PyObject_Call( pymeth, args, kw );
            Py_DECREF( pymeth );
         }

         Py_XDECREF( result );
         if ( PyErr_Occurred() ) {
            Py_DECREF( l );
            return 0;
         }

         result = CallPyObjMethod( self, "Clear" );
         Py_XDECREF( result );
         result = CallPyObjMethod( self, "extend", l );
         Py_XDECREF( result );
         Py_DECREF( l );

         Py_INCREF( Py_None );
         return Py_None;
      }
   }

////////////////////////////////////////////////////////////////////////////////
/// Implement a python-style index with TCollection::IndexOf.

   PyObject* TSeqCollectionIndex( PyObject* self, PyObject* obj )
   {
      PyObject* index = CallPyObjMethod( self, "IndexOf", obj );
      if ( ! index )
         return 0;

      if ( PyLong_AsLong( index ) < 0 ) {
         Py_DECREF( index );
         PyErr_SetString( PyExc_ValueError, "list.index(x): x not in list" );
         return 0;
      }

      return index;
   }

//- TObjArray behavior ---------------------------------------------------------
   PyObject* TObjArrayLen( PyObject* self )
   {
   // GetSize on a TObjArray returns its capacity, not size in use
      PyObject* size = CallPyObjMethod( self, "GetLast" );
      if ( ! size )
         return 0;

      long lsize = PyLong_AsLong( size );
      if ( lsize == -1 && PyErr_Occurred() )
         return 0;

      Py_DECREF( size );
      return PyInt_FromLong( lsize + 1 );
   }


//- TClonesArray behavior ------------------------------------------------------
   PyObject* TClonesArraySetItem( ObjectProxy* self, PyObject* args )
   {
   // TClonesArray sets objects by constructing them in-place; which is impossible
   // to support as the python object given as value must exist a priori. It can,
   // however, be memcpy'd and stolen, caveat emptor.
      ObjectProxy* pyobj = 0; PyObject* idx = 0;
      if ( ! PyArg_ParseTuple( args,
               const_cast< char* >( "OO!:__setitem__" ), &idx, &ObjectProxy_Type, &pyobj ) )
         return 0;

      if ( ! self->GetObject() ) {
         PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
         return 0;
      }

      PyObject* pyindex = PyStyleIndex( (PyObject*)self, idx );
      if ( ! pyindex )
         return 0;
      int index = (int)PyLong_AsLong( pyindex );
      Py_DECREF( pyindex );

   // get hold of the actual TClonesArray
      TClonesArray* cla =
         (TClonesArray*)OP2TCLASS(self)->DynamicCast( TClonesArray::Class(), self->GetObject() );

      if ( ! cla ) {
         PyErr_SetString( PyExc_TypeError, "attempt to call with null object" );
         return 0;
      }

      if ( Cppyy::GetScope( cla->GetClass()->GetName() ) != pyobj->ObjectIsA() ) {
         PyErr_Format( PyExc_TypeError, "require object of type %s, but %s given",
            cla->GetClass()->GetName(), Cppyy::GetFinalName( pyobj->ObjectIsA() ).c_str() );
      }

   // destroy old stuff, if applicable
      if ( ((const TClonesArray&)*cla)[index] ) {
         cla->RemoveAt( index );
      }

      if ( pyobj->GetObject() ) {
      // accessing an entry will result in new, unitialized memory (if properly used)
         TObject* object = (*cla)[index];
         pyobj->Release();
         TMemoryRegulator::RegisterObject( pyobj, object );
         memcpy( (void*)object, pyobj->GetObject(), cla->GetClass()->Size() );
      }

      Py_INCREF( Py_None );
      return Py_None;
   }

//- vector behavior as primitives ----------------------------------------------
   typedef struct {
      PyObject_HEAD
      PyObject*           vi_vector;
      void*               vi_data;
      PyROOT::TConverter* vi_converter;
      Py_ssize_t          vi_pos;
      Py_ssize_t          vi_len;
      Py_ssize_t          vi_stride;
   } vectoriterobject;

   static void vectoriter_dealloc( vectoriterobject* vi ) {
      Py_XDECREF( vi->vi_vector );
      delete vi->vi_converter;
      PyObject_GC_Del( vi );
   }

   static int vectoriter_traverse( vectoriterobject* vi, visitproc visit, void* arg ) {
      Py_VISIT( vi->vi_vector );
      return 0;
   }

   static PyObject* vectoriter_iternext( vectoriterobject* vi ) {
      if ( vi->vi_pos >= vi->vi_len )
         return nullptr;

      PyObject* result = nullptr;

      if ( vi->vi_data && vi->vi_converter ) {
         void* location  = (void*)((ptrdiff_t)vi->vi_data + vi->vi_stride * vi->vi_pos );
         result = vi->vi_converter->FromMemory( location );
      } else {
         PyObject* pyindex = PyLong_FromLong( vi->vi_pos );
         result = CallPyObjMethod( (PyObject*)vi->vi_vector, "_vector__at", pyindex );
         Py_DECREF( pyindex );
      }

      vi->vi_pos += 1;
      return result;
   }

   PyTypeObject VectorIter_Type = {
      PyVarObject_HEAD_INIT( &PyType_Type, 0 )
      (char*)"ROOT.vectoriter",  // tp_name
      sizeof(vectoriterobject),  // tp_basicsize
      0,
      (destructor)vectoriter_dealloc,            // tp_dealloc
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      Py_TPFLAGS_DEFAULT |
         Py_TPFLAGS_HAVE_GC,     // tp_flags
      0,
      (traverseproc)vectoriter_traverse,         // tp_traverse
      0, 0, 0,
      PyObject_SelfIter,         // tp_iter
      (iternextfunc)vectoriter_iternext,         // tp_iternext
      0,                         // tp_methods
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
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

   static PyObject* vector_iter( PyObject* v ) {
      vectoriterobject* vi = PyObject_GC_New( vectoriterobject, &VectorIter_Type );
      if ( ! vi ) return NULL;

      Py_INCREF( v );
      vi->vi_vector = v;

      PyObject* pyvalue_type = PyObject_GetAttrString( (PyObject*)Py_TYPE(v), "value_type" );
      PyObject* pyvalue_size = PyObject_GetAttrString( (PyObject*)Py_TYPE(v), "value_size" );

      if ( pyvalue_type && pyvalue_size ) {
         PyObject* pydata = CallPyObjMethod( v, "data" );
         if ( !pydata || Utility::GetBuffer( pydata, '*', 1, vi->vi_data, kFALSE ) == 0 )
            vi->vi_data = nullptr;
         Py_XDECREF( pydata );

         vi->vi_converter = PyROOT::CreateConverter( PyROOT_PyUnicode_AsString( pyvalue_type ) );
         vi->vi_stride    = PyLong_AsLong( pyvalue_size );
      } else {
         PyErr_Clear();
         vi->vi_data      = nullptr;
         vi->vi_converter = nullptr;
         vi->vi_stride    = 0;
      }

      Py_XDECREF( pyvalue_size );
      Py_XDECREF( pyvalue_type );

      vi->vi_len = vi->vi_pos = 0;
      vi->vi_len = PySequence_Size( v );

      PyObject_GC_Track( vi );
      return (PyObject*)vi;
   }


   PyObject* VectorGetItem( ObjectProxy* self, PySliceObject* index )
   {
   // Implement python's __getitem__ for std::vector<>s.
      if ( PySlice_Check( index ) ) {
         if ( ! self->GetObject() ) {
            PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
            return 0;
         }

         PyObject* pyclass = PyObject_GetAttr( (PyObject*)self, PyStrings::gClass );
         PyObject* nseq = PyObject_CallObject( pyclass, NULL );
         Py_DECREF( pyclass );

         Py_ssize_t start, stop, step;
         PySlice_GetIndices( (PyROOT_PySliceCast)index, PyObject_Length( (PyObject*)self ), &start, &stop, &step );
         for ( Py_ssize_t i = start; i < stop; i += step ) {
            PyObject* pyidx = PyInt_FromSsize_t( i );
            CallPyObjMethod( nseq, "push_back", CallPyObjMethod( (PyObject*)self, "_vector__at", pyidx ) );
            Py_DECREF( pyidx );
         }

         return nseq;
      }

      return CallSelfIndex( self, (PyObject*)index, "_vector__at" );
   }

   PyObject* VectorBoolSetItem( ObjectProxy* self, PyObject* args )
   {
   // std::vector<bool> is a special-case in C++, and its return type depends on
   // the compiler: treat it special here as well
      int bval = 0; PyObject* idx = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "Oi:__setitem__" ), &idx, &bval ) )
         return 0;

      if ( ! self->GetObject() ) {
         PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
         return 0;
      }

      PyObject* pyindex = PyStyleIndex( (PyObject*)self, idx );
      if ( ! pyindex )
         return 0;
      int index = (int)PyLong_AsLong( pyindex );
      Py_DECREF( pyindex );

      std::string clName = Cppyy::GetFinalName( self->ObjectIsA() );
      std::string::size_type pos = clName.find( "vector<bool" );
      if ( pos != 0 && pos != 5 /* following std:: */ ) {
         PyErr_Format( PyExc_TypeError,
                       "require object of type std::vector<bool>, but %s given",
                       Cppyy::GetFinalName( self->ObjectIsA() ).c_str() );
         return 0;
      }

   // get hold of the actual std::vector<bool> (no cast, as vector is never a base)
      std::vector<bool>* vb = (std::vector<bool>*)self->GetObject();

   // finally, set the value
      (*vb)[ index ] = (bool)bval;

      Py_INCREF( Py_None );
      return Py_None;
   }

//- map behavior as primitives ------------------------------------------------
   PyObject* MapContains( PyObject* self, PyObject* obj )
   {
   // Implement python's __contains__ for std::map<>s.
      PyObject* result = 0;

      PyObject* iter = CallPyObjMethod( self, "find", obj );
      if ( ObjectProxy_Check( iter ) ) {
         PyObject* end = CallPyObjMethod( self, "end" );
         if ( ObjectProxy_Check( end ) ) {
            if ( ! PyObject_RichCompareBool( iter, end, Py_EQ ) ) {
               Py_INCREF( Py_True );
               result = Py_True;
            }
         }
         Py_XDECREF( end );
      }
      Py_XDECREF( iter );

      if ( ! result ) {
         PyErr_Clear();            // e.g. wrong argument type, which should always lead to False
         Py_INCREF( Py_False );
         result = Py_False;
      }

      return result;
   }

//- STL container iterator support --------------------------------------------
   PyObject* StlSequenceIter( PyObject* self )
   {
   // Implement python's __iter__ for std::iterator<>s.
      PyObject* iter = CallPyObjMethod( self, "begin" );
      if ( iter ) {
         PyObject* end = CallPyObjMethod( self, "end" );
         if ( end )
            PyObject_SetAttr( iter, PyStrings::gEnd, end );
         Py_XDECREF( end );

         // add iterated collection as attribute so its refcount stays >= 1 while it's being iterated over
         PyObject_SetAttr( iter, PyUnicode_FromString("_collection"), self );
      }
      return iter;
   }

//- safe indexing for STL-like vector w/o iterator dictionaries ---------------
   PyObject* CheckedGetItem( PyObject* self, PyObject* obj )
   {
   // Implement a generic python __getitem__ for std::vector<>s that are missing
   // their std::vector<>::iterator dictionary. This is then used for iteration
   // by means of consecutive index.
      Bool_t inbounds = kFALSE;
      Py_ssize_t size = PySequence_Size( self );
      Py_ssize_t idx  = PyInt_AsSsize_t( obj );
      if ( 0 <= idx && 0 <= size && idx < size )
         inbounds = kTRUE;

      if ( inbounds ) {
         return CallPyObjMethod( self, "_getitem__unchecked", obj );
      } else if ( PyErr_Occurred() ) {
      // argument conversion problem: let method itself resolve anew and report
         PyErr_Clear();
         return CallPyObjMethod( self, "_getitem__unchecked", obj );
      } else {
         PyErr_SetString( PyExc_IndexError, "index out of range" );
      }

      return 0;
   }

//- pair as sequence to allow tuple unpacking ---------------------------------
   PyObject* PairUnpack( PyObject* self, PyObject* pyindex )
   {
   // For std::map<> iteration, unpack std::pair<>s into tuples for the loop.
      Long_t idx = PyLong_AsLong( pyindex );
      if ( idx == -1 && PyErr_Occurred() )
         return 0;

      if ( ! ObjectProxy_Check( self ) || ! ((ObjectProxy*)self)->GetObject() ) {
         PyErr_SetString( PyExc_TypeError, "unsubscriptable object" );
         return 0;
      }

      if ( (int)idx == 0 )
         return PyObject_GetAttr( self, PyStrings::gFirst );
      else if ( (int)idx == 1 )
         return PyObject_GetAttr( self, PyStrings::gSecond );

   // still here? Trigger stop iteration
      PyErr_SetString( PyExc_IndexError, "out of bounds" );
      return 0;
   }

//- string behavior as primitives ----------------------------------------------
#if PY_VERSION_HEX >= 0x03000000
// TODO: this is wrong, b/c it doesn't order
static int PyObject_Compare( PyObject* one, PyObject* other ) {
   return ! PyObject_RichCompareBool( one, other, Py_EQ );
}
#endif
   static inline PyObject* PyROOT_PyString_FromCppString( std::string* s ) {
      return PyROOT_PyUnicode_FromStringAndSize( s->c_str(), s->size() );
   }

   static inline PyObject* PyROOT_PyString_FromCppString( TString* s ) {
      return PyROOT_PyUnicode_FromStringAndSize( s->Data(), s->Length() );
   }

   static inline PyObject* PyROOT_PyString_FromCppString( TObjString* s ) {
      return PyROOT_PyUnicode_FromStringAndSize( s->GetString().Data(), s->GetString().Length() );
   }

#define PYROOT_IMPLEMENT_STRING_PYTHONIZATION( type, name )                   \
   inline PyObject* name##GetData( PyObject* self ) {                         \
      if ( PyROOT::ObjectProxy_Check( self ) ) {                              \
         type* obj = ((type*)((ObjectProxy*)self)->GetObject());              \
         if ( obj ) {                                                         \
            return PyROOT_PyString_FromCppString( obj );                      \
         } else {                                                             \
            return ObjectProxy_Type.tp_str( self );                           \
         }                                                                    \
      }                                                                       \
      PyErr_Format( PyExc_TypeError, "object mismatch (%s expected)", #type );\
      return 0;                                                               \
   }                                                                          \
                                                                              \
   PyObject* name##StringRepr( PyObject* self )                               \
   {                                                                          \
      PyObject* data = name##GetData( self );                                 \
      if ( data ) {                                                           \
         PyObject* repr = PyROOT_PyUnicode_FromFormat( "\'%s\'", PyROOT_PyUnicode_AsString( data ) ); \
         Py_DECREF( data );                                                   \
         return repr;                                                         \
      }                                                                       \
      return 0;                                                               \
   }                                                                          \
                                                                              \
   PyObject* name##StringIsEqual( PyObject* self, PyObject* obj )             \
   {                                                                          \
      PyObject* data = name##GetData( self );                                 \
      if ( data ) {                                                           \
         PyObject* result = PyObject_RichCompare( data, obj, Py_EQ );         \
         Py_DECREF( data );                                                   \
         return result;                                                       \
      }                                                                       \
      return 0;                                                               \
   }                                                                          \
                                                                              \
   PyObject* name##StringIsNotEqual( PyObject* self, PyObject* obj )          \
   {                                                                          \
      PyObject* data = name##GetData( self );                                 \
      if ( data ) {                                                           \
         PyObject* result = PyObject_RichCompare( data, obj, Py_NE );         \
         Py_DECREF( data );                                                   \
         return result;                                                       \
      }                                                                       \
      return 0;                                                               \
   }

   // Only define StlStringCompare:
   // TStringCompare is unused and generates a warning;
#define PYROOT_IMPLEMENT_STRING_PYTHONIZATION_CMP( type, name )               \
   PYROOT_IMPLEMENT_STRING_PYTHONIZATION( type, name )                        \
   PyObject* name##StringCompare( PyObject* self, PyObject* obj )             \
   {                                                                          \
      PyObject* data = name##GetData( self );                                 \
      int result = 0;                                                         \
      if ( data ) {                                                           \
         result = PyObject_Compare( data, obj );                              \
         Py_DECREF( data );                                                   \
      }                                                                       \
      if ( PyErr_Occurred() )                                                 \
         return 0;                                                            \
      return PyInt_FromLong( result );                                        \
   }

   PYROOT_IMPLEMENT_STRING_PYTHONIZATION_CMP( std::string, Stl )
   PYROOT_IMPLEMENT_STRING_PYTHONIZATION( TString, T )


//- TObjString behavior --------------------------------------------------------
   PYROOT_IMPLEMENT_STRING_PYTHONIZATION_CMP( TObjString, TObj )

////////////////////////////////////////////////////////////////////////////////
/// Implementation of python __len__ for TObjString.

   PyObject* TObjStringLength( PyObject* self )
   {
      PyObject* data = CallPyObjMethod( self, "GetName" );
      Py_ssize_t size = PySequence_Size( data );
      Py_DECREF( data );
      return PyInt_FromSsize_t( size );
   }


//- TIter behavior -------------------------------------------------------------
   PyObject* TIterNext( PyObject* self )
   {
   // Implementation of python __next__ (iterator protocol) for TIter.
      PyObject* next = CallPyObjMethod( self, "Next" );

      if ( ! next )
         return 0;

      if ( ! PyObject_IsTrue( next ) ) {
         Py_DECREF( next );
         PyErr_SetString( PyExc_StopIteration, "" );
         return 0;
      }

      return next;
   }


//- STL iterator behavior ------------------------------------------------------
   PyObject* StlIterNext( PyObject* self )
   {
   // Python iterator protocol __next__ for STL forward iterators.
      PyObject* next = 0;
      PyObject* last = PyObject_GetAttr( self, PyStrings::gEnd );

      if ( last != 0 ) {
      // handle special case of empty container (i.e. self is end)
         if ( PyObject_RichCompareBool( last, self, Py_EQ ) ) {
            PyErr_SetString( PyExc_StopIteration, "" );
         } else {
            PyObject* dummy = PyInt_FromLong( 1l );
            PyObject* iter = CallPyObjMethod( self, "__postinc__", dummy );
            Py_DECREF( dummy );
            if ( iter != 0 ) {
               if ( PyObject_RichCompareBool( last, iter, Py_EQ ) )
                  PyErr_SetString( PyExc_StopIteration, "" );
               else
                  next = CallPyObjMethod( iter, "__deref__" );
            } else {
               PyErr_SetString( PyExc_StopIteration, "" );
            }
            Py_XDECREF( iter );
         }
      } else {
         PyErr_SetString( PyExc_StopIteration, "" );
      }

      Py_XDECREF( last );
      return next;
   }

////////////////////////////////////////////////////////////////////////////////
/// Called if operator== not available (e.g. if a global overload as under gcc).
/// An exception is raised as the user should fix the dictionary.

   PyObject* StlIterIsEqual( PyObject* self, PyObject* other )
   {
      if (other != Py_None) {
         if (Utility::AddBinaryOperator(self, other, "==", "__eq__", nullptr, true))
            return PyObject_CallMethodObjArgs(self, PyStrings::gEq, other, nullptr);
      }

      return PyErr_Format( PyExc_LookupError,
         "No operator==(const %s&, const %s&) available in the dictionary!",
         Utility::ClassName( self ).c_str(), Utility::ClassName( other ).c_str()  );
   }

////////////////////////////////////////////////////////////////////////////////
/// Called if operator!= not available (e.g. if a global overload as under gcc).
/// An exception is raised as the user should fix the dictionary.

   PyObject* StlIterIsNotEqual( PyObject* self, PyObject* other )
   {
      if (other != Py_None) {
         if (Utility::AddBinaryOperator(self, other, "!=", "__ne__", nullptr, true))
            return PyObject_CallMethodObjArgs(self, PyStrings::gNe, other, nullptr);
      }

      return PyErr_Format( PyExc_LookupError,
         "No operator!=(const %s&, const %s&) available in the dictionary!",
         Utility::ClassName( self ).c_str(), Utility::ClassName( other ).c_str()  );
   }


//- TDirectory member templates ----------------------------------------------
   PyObject* TDirectoryGetObject( ObjectProxy* self, PyObject* args )
   {
   // Pythonization of TDirector::GetObject().
      PyObject* name = 0; ObjectProxy* ptr = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!:TDirectory::GetObject" ),
               &PyROOT_PyUnicode_Type, &name, &ObjectProxy_Type, &ptr ) )
         return 0;

      TDirectory* dir =
         (TDirectory*)OP2TCLASS(self)->DynamicCast( TDirectory::Class(), self->GetObject() );

      if ( ! dir ) {
         PyErr_SetString( PyExc_TypeError,
           "TDirectory::GetObject must be called with a TDirectory instance as first argument" );
         return 0;
      }

      void* address = dir->GetObjectChecked( PyROOT_PyUnicode_AsString( name ), OP2TCLASS(ptr) );
      if ( address ) {
         ptr->Set( address );

         Py_INCREF( Py_None );
         return Py_None;
      }

      PyErr_Format( PyExc_LookupError, "no such object, \"%s\"", PyROOT_PyUnicode_AsString( name ) );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////
/// Type-safe version of TDirectory::WriteObjectAny, which is a template for
/// the same reason on the C++ side.

   PyObject* TDirectoryWriteObject( ObjectProxy* self, PyObject* args )
   {
      ObjectProxy *wrt = 0; PyObject *name = 0, *option = 0;
      Int_t bufsize = 0;
      if ( ! PyArg_ParseTuple( args, const_cast< char* >( "O!O!|O!i:TDirectory::WriteObject" ),
               &ObjectProxy_Type, &wrt, &PyROOT_PyUnicode_Type, &name,
               &PyROOT_PyUnicode_Type, &option, &bufsize ) )
         return 0;

      TDirectory* dir =
         (TDirectory*)OP2TCLASS(self)->DynamicCast( TDirectory::Class(), self->GetObject() );

      if ( ! dir ) {
         PyErr_SetString( PyExc_TypeError,
           "TDirectory::WriteObject must be called with a TDirectory instance as first argument" );
         return 0;
      }

      Int_t result = 0;
      if ( option != 0 ) {
         result = dir->WriteObjectAny( wrt->GetObject(), OP2TCLASS(wrt),
            PyROOT_PyUnicode_AsString( name ), PyROOT_PyUnicode_AsString( option ), bufsize );
      } else {
         result = dir->WriteObjectAny(
            wrt->GetObject(), OP2TCLASS(wrt), PyROOT_PyUnicode_AsString( name ) );
      }

      return PyInt_FromLong( (Long_t)result );
   }

}


namespace PyROOT {      // workaround for Intel icc on Linux

//- TTree behavior ------------------------------------------------------------
   PyObject* TTreeGetAttr( ObjectProxy* self, PyObject* pyname )
   {
   // allow access to branches/leaves as if they are data members
      const char* name1 = PyROOT_PyUnicode_AsString( pyname );
      if ( ! name1 )
         return 0;

   // get hold of actual tree
      TTree* tree =
         (TTree*)OP2TCLASS(self)->DynamicCast( TTree::Class(), self->GetObject() );

      if ( ! tree ) {
         PyErr_SetString( PyExc_ReferenceError, "attempt to access a null-pointer" );
         return 0;
      }

   // deal with possible aliasing
      const char* name = tree->GetAlias( name1 );
      if ( ! name ) name = name1;

   // search for branch first (typical for objects)
      TBranch* branch = tree->GetBranch( name );
      if ( ! branch ) {
      // for benefit of naming of sub-branches, the actual name may have a trailing '.'
         branch = tree->GetBranch( (std::string( name ) + '.' ).c_str() );
      }

      if ( branch ) {
      // found a branched object, wrap its address for the object it represents

      // for partial return of a split object
         if ( branch->InheritsFrom(TBranchElement::Class()) ) {
            TBranchElement* be = (TBranchElement*)branch;
            if ( be->GetCurrentClass() && (be->GetCurrentClass() != be->GetTargetClass()) && (0 <= be->GetID()) ) {
               Long_t offset = ((TStreamerElement*)be->GetInfo()->GetElements()->At(be->GetID()))->GetOffset();
               return BindCppObjectNoCast( be->GetObject() + offset, Cppyy::GetScope( be->GetCurrentClass()->GetName() ) );
            }
         }

      // for return of a full object
         if ( branch->IsA() == TBranchElement::Class() || branch->IsA() == TBranchObject::Class() ) {
            TClass* klass = TClass::GetClass( branch->GetClassName() );
            if ( klass && branch->GetAddress() )
               return BindCppObjectNoCast( *(void**)branch->GetAddress(), Cppyy::GetScope( branch->GetClassName() ) );

         // try leaf, otherwise indicate failure by returning a typed null-object
            TObjArray* leaves = branch->GetListOfLeaves();
            if ( klass && ! tree->GetLeaf( name ) &&
                 ! (leaves->GetSize() && ( leaves->First() == leaves->Last() ) ) )
               return BindCppObjectNoCast( NULL, Cppyy::GetScope( branch->GetClassName() ) );
         }
      }

   // if not, try leaf
      TLeaf* leaf = tree->GetLeaf( name );
      if ( branch && ! leaf ) {
         leaf = branch->GetLeaf( name );
         if ( ! leaf ) {
            TObjArray* leaves = branch->GetListOfLeaves();
            if ( leaves->GetSize() && ( leaves->First() == leaves->Last() ) ) {
            // i.e., if unambiguously only this one
               leaf = (TLeaf*)leaves->At( 0 );
            }
         }
      }

      if ( leaf ) {
      // found a leaf, extract value and wrap
         if ( 1 < leaf->GetLenStatic() || leaf->GetLeafCount() ) {
         // array types
            std::string typeName = leaf->GetTypeName();
            TConverter* pcnv = CreateConverter( typeName + '*', leaf->GetNdata() );

            void* address = 0;
            if ( leaf->GetBranch() ) address = (void*)leaf->GetBranch()->GetAddress();
            if ( ! address ) address = (void*)leaf->GetValuePointer();

            PyObject* value = pcnv->FromMemory( &address );
            delete pcnv;

            return value;
         } else if ( leaf->GetValuePointer() ) {
         // value types
            TConverter* pcnv = CreateConverter( leaf->GetTypeName() );
            PyObject* value = 0;
            if ( leaf->IsA() == TLeafElement::Class() || leaf->IsA() == TLeafObject::Class() )
               value = pcnv->FromMemory( (void*)*(void**)leaf->GetValuePointer() );
            else
               value = pcnv->FromMemory( (void*)leaf->GetValuePointer() );
            delete pcnv;

            return value;
         }
      }

   // confused
      PyErr_Format( PyExc_AttributeError,
          "\'%s\' object has no attribute \'%s\'", tree->IsA()->GetName(), name );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////

   class TTreeMemberFunction : public PyCallable {
   protected:
      TTreeMemberFunction( MethodProxy* org ) { Py_INCREF( org ); fOrg = org; }
      TTreeMemberFunction( const TTreeMemberFunction& t ) : PyCallable( t )
      {
      // Copy constructor; conform to python reference counting.
         Py_INCREF( t.fOrg );
         fOrg = t.fOrg;
      }
      TTreeMemberFunction& operator=( const TTreeMemberFunction& t )
      {
      // Assignment operator; conform to python reference counting.
         if ( &t != this ) {
            Py_INCREF( t.fOrg );
            Py_XDECREF( fOrg );
            fOrg = t.fOrg;
         }
         return *this;
      }
      ~TTreeMemberFunction() { Py_DECREF( fOrg ); fOrg = 0; }

   public:
      virtual PyObject* GetSignature() { return PyROOT_PyUnicode_FromString( "(...)" ); }
      virtual PyObject* GetPrototype() { return PyObject_GetAttrString( (PyObject*)fOrg, (char*)"__doc__" ); }
      virtual Int_t GetPriority() { return 100; }
      virtual PyObject* GetCoVarNames() {
         PyObject* co_varnames = PyTuple_New( 1 /* self */ + 1 /* fake */ );
         PyTuple_SET_ITEM( co_varnames, 0, PyROOT_PyUnicode_FromString( "self" ) );
         PyTuple_SET_ITEM( co_varnames, 1, PyROOT_PyUnicode_FromString( "*args" ) );
         return co_varnames;
      }
      virtual PyObject* GetArgDefault( Int_t ) { return NULL; }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TTree" ); }

   protected:
      MethodProxy* fOrg;
   };

////////////////////////////////////////////////////////////////////////////////

   class TTreeBranch : public TTreeMemberFunction {
   public:
      TTreeBranch( MethodProxy* org ) : TTreeMemberFunction( org ) {}

   public:
      virtual Int_t GetMaxArgs() { return 5; }
      virtual PyCallable* Clone() { return new TTreeBranch( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* /* ctxt */ )
      {
      // acceptable signatures:
      //   ( const char*, void*, const char*, Int_t = 32000 )
      //   ( const char*, const char*, T**, Int_t = 32000, Int_t = 99 )
      //   ( const char*, T**, Int_t = 32000, Int_t = 99 )
         int argc = PyTuple_GET_SIZE( args );

         if ( 2 <= argc ) {
            TTree* tree =
               (TTree*)OP2TCLASS(self)->DynamicCast( TTree::Class(), self->GetObject() );

            if ( ! tree ) {
               PyErr_SetString( PyExc_TypeError,
                  "TTree::Branch must be called with a TTree instance as first argument" );
               return 0;
            }

            PyObject *name = 0, *clName = 0, *leaflist = 0;
            PyObject *address = 0;
            PyObject *bufsize = 0, *splitlevel = 0;

         // try: ( const char*, void*, const char*, Int_t = 32000 )
            if ( PyArg_ParseTuple( args, const_cast< char* >( "O!OO!|O!:Branch" ),
                   &PyROOT_PyUnicode_Type, &name, &address, &PyROOT_PyUnicode_Type,
                   &leaflist, &PyInt_Type, &bufsize ) ) {

               void* buf = 0;
               if ( ObjectProxy_Check( address ) )
                  buf = (void*)((ObjectProxy*)address)->GetObject();
               else
                  Utility::GetBuffer( address, '*', 1, buf, kFALSE );

               if ( buf != 0 ) {
                  TBranch* branch = 0;
                  if ( argc == 4 ) {
                     branch = tree->Branch( PyROOT_PyUnicode_AsString( name ), buf,
                        PyROOT_PyUnicode_AsString( leaflist ), PyInt_AS_LONG( bufsize ) );
                  } else {
                     branch = tree->Branch( PyROOT_PyUnicode_AsString( name ), buf,
                        PyROOT_PyUnicode_AsString( leaflist ) );
                  }

                  return BindCppObject( branch, "TBranch" );
               }

            }
            PyErr_Clear();

         // try: ( const char*, const char*, T**, Int_t = 32000, Int_t = 99 )
         //  or: ( const char*,              T**, Int_t = 32000, Int_t = 99 )
            Bool_t bIsMatch = kFALSE;
            if ( PyArg_ParseTuple( args, const_cast< char* >( "O!O!O|O!O!:Branch" ),
                   &PyROOT_PyUnicode_Type, &name, &PyROOT_PyUnicode_Type, &clName, &address,
                   &PyInt_Type, &bufsize, &PyInt_Type, &splitlevel ) ) {
               bIsMatch = kTRUE;
            } else {
               PyErr_Clear(); clName = 0;    // clName no longer used
               if ( PyArg_ParseTuple( args, const_cast< char* >( "O!O|O!O!" ),
                      &PyROOT_PyUnicode_Type, &name, &address,
                      &PyInt_Type, &bufsize, &PyInt_Type, &splitlevel ) ) {
                  bIsMatch = kTRUE;
               } else
                  PyErr_Clear();
            }

            if ( bIsMatch == kTRUE ) {
               std::string klName = clName ? PyROOT_PyUnicode_AsString( clName ) : "";
               void* buf = 0;

               if ( ObjectProxy_Check( address ) ) {
                  if ( ((ObjectProxy*)address)->fFlags & ObjectProxy::kIsReference )
                     buf = (void*)((ObjectProxy*)address)->fObject;
                  else
                     buf = (void*)&((ObjectProxy*)address)->fObject;

                  if ( ! clName ) {
                     klName = OP2TCLASS((ObjectProxy*)address)->GetName();
                     argc += 1;
                  }
               } else
                  Utility::GetBuffer( address, '*', 1, buf, kFALSE );

               if ( buf != 0 && klName != "" ) {
                  TBranch* branch = 0;
                  if ( argc == 3 ) {
                     branch = tree->Branch( PyROOT_PyUnicode_AsString( name ), klName.c_str(), buf );
                  } else if ( argc == 4 ) {
                     branch = tree->Branch( PyROOT_PyUnicode_AsString( name ), klName.c_str(), buf,
                        PyInt_AS_LONG( bufsize ) );
                  } else if ( argc == 5 ) {
                     branch = tree->Branch( PyROOT_PyUnicode_AsString( name ), klName.c_str(), buf,
                        PyInt_AS_LONG( bufsize ), PyInt_AS_LONG( splitlevel ) );
                  }

                  return BindCppObject( branch, "TBranch" );
               }
            }
         }

      // still here? Then call original Branch() to reach the other overloads:
         Py_INCREF( (PyObject*)self );
         fOrg->fSelf = self;
         PyObject* result = PyObject_Call( (PyObject*)fOrg, args, kwds );
         fOrg->fSelf = 0;
         Py_DECREF( (PyObject*)self );

         return result;
      }
   };

////////////////////////////////////////////////////////////////////////////////

   class TTreeSetBranchAddress : public TTreeMemberFunction {
   public:
      TTreeSetBranchAddress( MethodProxy* org ) : TTreeMemberFunction( org ) {}

   public:
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString( "TBranch* TTree::SetBranchAddress( ... )" );
      }

      virtual Int_t GetMaxArgs() { return 2; }
      virtual PyCallable* Clone() { return new TTreeSetBranchAddress( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* /* ctxt */ )
      {
      // acceptable signature:
      //   ( const char*, void* )
         int argc = PyTuple_GET_SIZE( args );

         if ( 2 == argc ) {
            TTree* tree =
               (TTree*)OP2TCLASS(self)->DynamicCast( TTree::Class(), self->GetObject() );

            if ( ! tree ) {
               PyErr_SetString( PyExc_TypeError,
                  "TTree::SetBranchAddress must be called with a TTree instance as first argument" );
               return 0;
            }

            PyObject *name = 0, *address = 0;

         // try: ( const char*, void* )
            if ( PyArg_ParseTuple( args, const_cast< char* >( "SO:SetBranchAddress" ),
                    &name, &address ) ) {

               void* buf = 0;
               if ( ObjectProxy_Check( address ) ) {
                  if ( ((ObjectProxy*)address)->fFlags & ObjectProxy::kIsReference )
                     buf = (void*)((ObjectProxy*)address)->fObject;
                  else
                     buf = (void*)&((ObjectProxy*)address)->fObject;
               } else
                  Utility::GetBuffer( address, '*', 1, buf, kFALSE );

               if ( buf != 0 ) {
                  tree->SetBranchAddress( PyROOT_PyUnicode_AsString( name ), buf );

                  Py_INCREF( Py_None );
                  return Py_None;
               }
            }
         }

      // still here? Then call original Branch() to reach the other overloads:
         Py_INCREF( (PyObject*)self );
         fOrg->fSelf = self;
         PyObject* result = PyObject_Call( (PyObject*)fOrg, args, kwds );
         fOrg->fSelf = 0;
         Py_DECREF( (PyObject*)self );

         return result;
      }

   protected:
      virtual PyObject* ReportTypeError()
      {
         PyErr_SetString( PyExc_TypeError,
            "TTree::SetBranchAddress must be called with a TTree instance as first argument" );
         return 0;
      }
   };


// TChain overrides TTree's SetBranchAddress, so set it again (the python method only forwards
//   onto a TTree*, so the C++ virtual function call will make sure the right method is used)
   class TChainSetBranchAddress : public TTreeSetBranchAddress {
   public:
      TChainSetBranchAddress( MethodProxy* org ) : TTreeSetBranchAddress( org ) {}

   public:
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString( "TBranch* TChain::SetBranchAddress( ... )" );
      }

      virtual Int_t GetMaxArgs() { return 2; }
      virtual PyCallable* Clone() { return new TChainSetBranchAddress( *this ); }

   protected:
      virtual PyObject* ReportTypeError()
      {
         PyErr_SetString( PyExc_TypeError,
            "TChain::SetBranchAddress must be called with a TChain instance as first argument" );
         return 0;
      }
   };

//- TMinuit behavior ----------------------------------------------------------
   void TMinuitPyCallback( void* vpyfunc, Long_t /* npar */,
         Int_t& a0, Double_t* a1, Double_t& a2, Double_t* a3, Int_t a4 ) {
   // a void* was passed to keep the interface on builtin types only
      PyObject* pyfunc = (PyObject*)vpyfunc;

   // prepare arguments
      PyObject* pya0 = BufFac_t::Instance()->PyBuffer_FromMemory( &a0,  sizeof(Int_t) );
      PyObject* pya1 = BufFac_t::Instance()->PyBuffer_FromMemory(  a1, a0 * sizeof(Double_t) );
      PyObject* pya2 = BufFac_t::Instance()->PyBuffer_FromMemory( &a2,  sizeof(Double_t) );
      PyObject* pya3 = BufFac_t::Instance()->PyBuffer_FromMemory(  a3, -1 ); // size unknown

      if ( ! (pya0 && pya1 && pya2 && pya3) ) {
         Py_XDECREF( pya3 ); Py_XDECREF( pya2 ); Py_XDECREF( pya1 ); Py_XDECREF( pya0 );
         return;
      }

   // perform actual call
      PyObject* result = PyObject_CallFunction(
         pyfunc, (char*)"OOOOi", pya0, pya1, pya2, pya3, a4 );
      Py_DECREF( pya3 ); Py_DECREF( pya2 ); Py_DECREF( pya1 ); Py_DECREF( pya0 );

      if ( ! result ) {
         PyErr_Print();
         throw std::runtime_error( "TMinuit python fit function call failed" );
      }

      Py_XDECREF( result );
   }

//- TFN behavior --------------------------------------------------------------
   double TFNPyCallback( void* vpyfunc, Long_t npar, double* a0, double* a1 ) {
   // a void* was passed to keep the interface on builtin types only
      PyObject* pyfunc = (PyObject*)vpyfunc;

   // prepare arguments and call
      PyObject* pya0 = BufFac_t::Instance()->PyBuffer_FromMemory( a0, 4 * sizeof(double) );
      if ( ! pya0 )
         return 0.;

      PyObject* result = 0;
      if ( npar != 0 ) {
         PyObject* pya1 = BufFac_t::Instance()->PyBuffer_FromMemory( a1, npar * sizeof(double) );
         result = PyObject_CallFunction( pyfunc, (char*)"OO", pya0, pya1 );
         Py_DECREF( pya1 );
      } else
         result = PyObject_CallFunction( pyfunc, (char*)"O", pya0 );

      Py_DECREF( pya0 );

   // translate result, throw if an error has occurred
      double d = 0.;
      if ( ! result ) {
         PyErr_Print();
         throw std::runtime_error( "TFN python function call failed" );
      } else {
         d = PyFloat_AsDouble( result );
         Py_DECREF( result );
      }

      return d;
   }

} // namespace PyROOT


namespace {

// for convenience
   using namespace PyROOT;

//- THN behavior --------------------------------------------------------------
   PyObject* THNIMul( PyObject* self, PyObject* scale )
   {
   // Use THN::Scale to perform *= ... need this stub to return self.
      PyObject* result = CallPyObjMethod( self, "Scale", scale );
      if ( ! result )
         return result;

      Py_DECREF( result );
      Py_INCREF( self );
      return self;
   }


////////////////////////////////////////////////////////////////////////////////

   class TPretendInterpreted: public PyCallable {
   public:
      TPretendInterpreted( int nArgs ) : fNArgs( nArgs ) {}

   public:
      Int_t GetNArgs() { return fNArgs; }
      virtual Int_t GetPriority() { return 100; }
      virtual Int_t GetMaxArgs() { return GetNArgs()+1; }
      virtual PyObject* GetCoVarNames() {
         PyObject* co_varnames = PyTuple_New( 1 /* self */ + 1 /* fake */ );
         PyTuple_SET_ITEM( co_varnames, 0, PyROOT_PyUnicode_FromString( "self" ) );
         PyTuple_SET_ITEM( co_varnames, 1, PyROOT_PyUnicode_FromString( "*args" ) );
         return co_varnames;
      }
      virtual PyObject* GetArgDefault( Int_t ) { return NULL; }

      Bool_t IsCallable( PyObject* pyobject )
      {
      // Determine whether the given pyobject is indeed callable.
         if ( ! pyobject || ! PyCallable_Check( pyobject ) ) {
            PyObject* str = pyobject ? PyObject_Str( pyobject ) : PyROOT_PyUnicode_FromString( "null pointer" );
            PyErr_Format( PyExc_ValueError,
               "\"%s\" is not a valid python callable", PyROOT_PyUnicode_AsString( str ) );
            Py_DECREF( str );
            return kFALSE;
         }

         return kTRUE;
      }

   private:
      Int_t fNArgs;
   };

////////////////////////////////////////////////////////////////////////////////

   class TF1InitWithPyFunc : public TPretendInterpreted {
   public:
      TF1InitWithPyFunc( int ntf = 1 ) : TPretendInterpreted( 2 + 2*ntf ) {}

   public:
      virtual PyObject* GetSignature() { return PyROOT_PyUnicode_FromString( "(...)" ); }
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TF1::TF1(const char* name, PyObject* callable, "
            "Double_t xmin, Double_t xmax, Int_t npar = 0)" );
      }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TF1" ); }
      virtual PyCallable* Clone() { return new TF1InitWithPyFunc( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* /* kwds */, TCallContext* /* ctxt */ )
      {
      // expected signature: ( char* name, pyfunc, double xmin, double xmax, int npar = 0 )
         int argc = PyTuple_GET_SIZE( args );
         const int reqNArgs = GetNArgs();
         if ( ! ( argc == reqNArgs || argc == reqNArgs+1 ) ) {
            PyErr_Format( PyExc_TypeError,
               "TFN::TFN(const char*, PyObject* callable, ...) =>\n"
               "    takes at least %d and at most %d arguments (%d given)",
               reqNArgs, reqNArgs+1, argc );
            return 0;              // reported as an overload failure
         }

         PyObject* pyfunc = PyTuple_GET_ITEM( args, 1 );

      // verify/setup the callback parameters
         Long_t npar = 0;             // default value if not given
         if ( argc == reqNArgs+1 )
            npar = PyInt_AsLong( PyTuple_GET_ITEM( args, reqNArgs ) );

      // create signature
         std::vector<std::string> signature; signature.reserve( 2 );
         signature.push_back( "double*" );
         signature.push_back( "double*" );

      // registration with Cling
         void* fptr = Utility::CreateWrapperMethod(
            pyfunc, npar, "double", signature, "TFNPyCallback" );
         if ( ! fptr /* PyErr was set */ )
            return 0;

      // get constructor
         MethodProxy* method =
            (MethodProxy*)PyObject_GetAttr( (PyObject*)self, PyStrings::gInit );

      // build new argument array
         PyObject* newArgs = PyTuple_New( reqNArgs + 1 );

         for ( int iarg = 0; iarg < argc; ++iarg ) {
            PyObject* item = PyTuple_GET_ITEM( args, iarg );
            if ( iarg != 1 ) {
               Py_INCREF( item );
               PyTuple_SET_ITEM( newArgs, iarg, item );
            } else {
               PyTuple_SET_ITEM( newArgs, iarg, PyROOT_PyCapsule_New( fptr, NULL, NULL ) );
            }
         }

         if ( argc == reqNArgs )             // meaning: use default for last value
            PyTuple_SET_ITEM( newArgs, reqNArgs, PyInt_FromLong( 0l ) );

      // re-run constructor, will select the proper one with void* for callback
         PyObject* result = PyObject_CallObject( (PyObject*)method, newArgs );

      // done, may have worked, if not: 0 is returned
         Py_DECREF( newArgs );
         Py_DECREF( method );
         return result;
      }
   };

////////////////////////////////////////////////////////////////////////////////

   class TF2InitWithPyFunc : public TF1InitWithPyFunc {
   public:
      TF2InitWithPyFunc() : TF1InitWithPyFunc( 2 ) {}

   public:
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TF2::TF2(const char* name, PyObject* callable, "
            "Double_t xmin, Double_t xmax, "
            "Double_t ymin, Double_t ymax, Int_t npar = 0)" );
      }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TF2" ); }
      virtual PyCallable* Clone() { return new TF2InitWithPyFunc( *this ); }
   };

////////////////////////////////////////////////////////////////////////////////

   class TF3InitWithPyFunc : public TF1InitWithPyFunc {
   public:
      TF3InitWithPyFunc() : TF1InitWithPyFunc( 3 ) {}

   public:
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TF3::TF3(const char* name, PyObject* callable, "
            "Double_t xmin, Double_t xmax, "
            "Double_t ymin, Double_t ymax, "
            "Double_t zmin, Double_t zmax, Int_t npar = 0)" );
      }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TF3" ); }
      virtual PyCallable* Clone() { return new TF3InitWithPyFunc( *this ); }
   };

//- TFunction behavior ---------------------------------------------------------
   PyObject* TFunctionCall( ObjectProxy*& self, PyObject* args ) {
      return TFunctionHolder( Cppyy::gGlobalScope, (Cppyy::TCppMethod_t)self->GetObject() ).Call( self, args, 0 );
   }


//- TMinuit behavior -----------------------------------------------------------
   class TMinuitSetFCN : public TPretendInterpreted {
   public:
      TMinuitSetFCN( int nArgs = 1 ) : TPretendInterpreted( nArgs ) {}

   public:
      virtual PyObject* GetSignature() { return PyROOT_PyUnicode_FromString( "(PyObject* callable)" ); }
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TMinuit::SetFCN(PyObject* callable)" );
      }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TMinuit" ); }
      virtual PyCallable* Clone() { return new TMinuitSetFCN( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt )
      {
      // expected signature: ( pyfunc )
         int argc = PyTuple_GET_SIZE( args );
         if ( argc != 1 ) {
            PyErr_Format( PyExc_TypeError,
               "TMinuit::SetFCN(PyObject* callable, ...) =>\n"
               "    takes exactly 1 argument (%d given)", argc );
            return 0;              // reported as an overload failure
         }

         PyObject* pyfunc = PyTuple_GET_ITEM( args, 0 );
         if ( ! IsCallable( pyfunc ) )
            return 0;

      // create signature
         std::vector<std::string> signature; signature.reserve( 5 );
         signature.push_back( "Int_t&" );
         signature.push_back( "Double_t*" );
         signature.push_back( "Double_t&" );
         signature.push_back( "Double_t*" );
         signature.push_back( "Int_t" );

      // registration with Cling
         void* fptr = Utility::CreateWrapperMethod(
            pyfunc, 5, "void", signature, "TMinuitPyCallback" );
         if ( ! fptr /* PyErr was set */ )
            return 0;

      // get setter function
         MethodProxy* method =
            (MethodProxy*)PyObject_GetAttr( (PyObject*)self, PyStrings::gSetFCN );

      // CLING WORKAROUND: SetFCN(void* fun) is deprecated but for whatever reason
      // still available yet not functional; select the correct one based on its
      // signature of the full function pointer
         PyCallable* setFCN = 0;
         const MethodProxy::Methods_t& methods = method->fMethodInfo->fMethods;
         for ( MethodProxy::Methods_t::const_iterator im = methods.begin(); im != methods.end(); ++im ) {
             PyObject* sig = (*im)->GetSignature();
             if ( sig && strstr( PyROOT_PyUnicode_AsString( sig ), "Double_t&" ) ) {
             // the comparison was not exact, but this is just a workaround
                setFCN = *im;
                Py_DECREF( sig );
                break;
             }
             Py_DECREF( sig );
         }
         if ( ! setFCN )     // this never happens but Coverity insists; it can be
            return 0;        // removed with the workaround in due time
      // END CLING WORKAROUND

      // build new argument array
         PyObject* newArgs = PyTuple_New( 1 );
         PyTuple_SET_ITEM( newArgs, 0, PyROOT_PyCapsule_New( fptr, NULL, NULL ) );

      // re-run
      // CLING WORKAROUND: this is to be the call once TMinuit is fixed:
         // PyObject* result = PyObject_CallObject( (PyObject*)method, newArgs );
         PyObject* result = setFCN->Call( self, newArgs, kwds, ctxt );
      // END CLING WORKAROUND

      // done, may have worked, if not: 0 is returned
         Py_DECREF( newArgs );
         Py_DECREF( method );
         return result;
      }
   };

   class TMinuitFitterSetFCN : public TMinuitSetFCN {
   public:
      TMinuitFitterSetFCN() : TMinuitSetFCN( 1 ) {}

   public:
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TMinuitFitter::SetFCN(PyObject* callable)" );
      }

      virtual PyCallable* Clone() { return new TMinuitFitterSetFCN( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* kwds, TCallContext* ctxt )
      {
      // expected signature: ( pyfunc )
         int argc = PyTuple_GET_SIZE( args );
         if ( argc != 1 ) {
            PyErr_Format( PyExc_TypeError,
               "TMinuitFitter::SetFCN(PyObject* callable, ...) =>\n"
               "    takes exactly 1 argument (%d given)", argc );
            return 0;              // reported as an overload failure
         }

         return TMinuitSetFCN::Call( self, args, kwds, ctxt );
      }
   };

//- Fit::TFitter behavior ------------------------------------------------------
   PyObject* gFitterPyCallback = 0;

   void FitterPyCallback( int& npar, double* gin, double& f, double* u, int flag )
   {
   // Cling-callable callback for Fit::Fitter derived objects.
      PyObject* result = 0;

   // prepare arguments
      PyObject* arg1 = BufFac_t::Instance()->PyBuffer_FromMemory( &npar );

      PyObject* arg2 = BufFac_t::Instance()->PyBuffer_FromMemory( gin );

      PyObject* arg3 = PyList_New( 1 );
      PyList_SetItem( arg3, 0, PyFloat_FromDouble( f ) );

      PyObject* arg4 = BufFac_t::Instance()->PyBuffer_FromMemory( u, npar * sizeof(double) );

   // perform actual call
      result = PyObject_CallFunction(
         gFitterPyCallback, (char*)"OOOOi", arg1, arg2, arg3, arg4, flag );
      f = PyFloat_AsDouble( PyList_GetItem( arg3, 0 ) );

      Py_DECREF( arg4 ); Py_DECREF( arg3 ); Py_DECREF( arg2 ); Py_DECREF( arg1 );

      if ( ! result ) {
         PyErr_Print();
         throw std::runtime_error( "TMinuit python fit function call failed" );
      }

      Py_XDECREF( result );
   }

   class TFitterFitFCN : public TPretendInterpreted {
   public:
      TFitterFitFCN() : TPretendInterpreted( 2 ) {}

   public:
      virtual PyObject* GetSignature()
      {
         return PyROOT_PyUnicode_FromString(
            "(PyObject* callable, int npar = 0, const double* params = 0, unsigned int dataSize = 0, bool chi2fit = false)" );
      }
      virtual PyObject* GetPrototype()
      {
         return PyROOT_PyUnicode_FromString(
            "TFitter::FitFCN(PyObject* callable, int npar = 0, const double* params = 0, unsigned int dataSize = 0, bool chi2fit = false)" );
      }
      virtual PyObject* GetScopeProxy() { return CreateScopeProxy( "TFitter" ); }
      virtual PyCallable* Clone() { return new TFitterFitFCN( *this ); }

      virtual PyObject* Call(
         ObjectProxy*& self, PyObject* args, PyObject* /* kwds */, TCallContext* /* ctxt */ )
      {
      // expected signature: ( self, pyfunc, int npar = 0, const double* params = 0, unsigned int dataSize = 0, bool chi2fit = false )
         int argc = PyTuple_GET_SIZE( args );
         if ( argc < 1 ) {
            PyErr_Format( PyExc_TypeError,
               "TFitter::FitFCN(PyObject* callable, ...) =>\n"
               "    takes at least 1 argument (%d given)", argc );
            return 0;              // reported as an overload failure
         }

         PyObject* pyfunc = PyTuple_GET_ITEM( args, 0 );
         if ( ! IsCallable( pyfunc ) )
            return 0;

      // global registration
         Py_XDECREF( gFitterPyCallback );
         Py_INCREF( pyfunc );
         gFitterPyCallback = pyfunc;

      // get function
         MethodProxy* method =
            (MethodProxy*)PyObject_GetAttr( (PyObject*)self, PyStrings::gFitFCN );

      // build new argument array
         PyObject* newArgs = PyTuple_New( argc );
         PyTuple_SET_ITEM( newArgs, 0, PyROOT_PyCapsule_New( (void*)FitterPyCallback, NULL, NULL ) );
         for ( int iarg = 1; iarg < argc; ++iarg ) {
            PyObject* pyarg = PyTuple_GET_ITEM( args, iarg );
            Py_INCREF( pyarg );
            PyTuple_SET_ITEM( newArgs, iarg, pyarg );
         }

      // re-run
         PyObject* result = PyObject_CallObject( (PyObject*)method, newArgs );

      // done, may have worked, if not: 0 is returned
         Py_DECREF( newArgs );
         Py_DECREF( method );
         return result;
      }
   };


//- TFile::Get -----------------------------------------------------------------
   PyObject* TFileGetAttr( PyObject* self, PyObject* attr )
   {
   // Pythonization of TFile::Get that raises AttributeError on failure.
      PyObject* result = CallPyObjMethod( self, "Get", attr );
      if ( !result )
         return result;

      if ( !PyObject_IsTrue( result ) ) {
         PyObject* astr = PyObject_Str( attr );
         PyErr_Format( PyExc_AttributeError, "TFile object has no attribute \'%s\'",
                       PyROOT_PyUnicode_AsString( astr ) );
         Py_DECREF( astr );
         Py_DECREF( result );
         return nullptr;
      }

   // caching behavior seems to be more clear to the user; can always override said
   // behavior (i.e. re-read from file) with an explicit Get() call
      PyObject_SetAttr( self, attr, result );
      return result;
   }

// This is done for TFile, but Get() is really defined in TDirectoryFile and its base
// TDirectory suffers from a similar problem. Nevertheless, the TFile case is by far
// the most common, so we'll leave it at this until someone asks for one of the bases
// to be pythonized.
   PyObject* TDirectoryFileGet( ObjectProxy* self, PyObject* pynamecycle )
   {
   // Pythonization of TDirectoryFile::Get that handles non-TObject deriveds
      if ( ! ObjectProxy_Check( self ) ) {
         PyErr_SetString( PyExc_TypeError,
            "TDirectoryFile::Get must be called with a TDirectoryFile instance as first argument" );
         return nullptr;
      }

      TDirectoryFile* dirf =
         (TDirectoryFile*)OP2TCLASS(self)->DynamicCast( TDirectoryFile::Class(), self->GetObject() );
      if ( !dirf ) {
         PyErr_SetString( PyExc_ReferenceError, "attempt to access a null-pointer" );
         return nullptr;
      }

      const char* namecycle = PyROOT_PyUnicode_AsString( pynamecycle );
      if ( !namecycle )
         return nullptr;     // TypeError already set

      TKey* key = dirf->GetKey( namecycle );
      if ( key ) {
         void* addr = dirf->GetObjectChecked( namecycle, key->GetClassName() );
         return BindCppObjectNoCast( addr,
            (Cppyy::TCppType_t)Cppyy::GetScope( key->GetClassName() ), kFALSE );
      }

      // no key? for better or worse, call normal Get()
      void* addr = dirf->Get( namecycle );
      return BindCppObject( addr, (Cppyy::TCppType_t)Cppyy::GetScope( "TObject" ), kFALSE );
   }

   //- Pretty printing with cling::PrintValue
   PyObject *ClingPrintValue(ObjectProxy *self)
   {
      PyObject *cppname = PyObject_GetAttrString((PyObject *)self, "__cppname__");
      if (!PyROOT_PyUnicode_Check(cppname))
         return 0;
      std::string className = PyROOT_PyUnicode_AsString(cppname);
      Py_XDECREF(cppname);

      std::string printResult = gInterpreter->ToString(className.c_str(), self->GetObject());
      if (printResult.find("@0x") == 0) {
         // Fall back to __repr__ if we just get an address from cling
         auto method = PyObject_GetAttrString((PyObject*)self, "__repr__");
         auto res = PyObject_CallObject(method, nullptr);
         Py_DECREF(method);
         return res;
      } else {
         return PyROOT_PyUnicode_FromString(printResult.c_str());
      }
   }

   //- Adding array interface to classes ---------------
   void AddArrayInterface(PyObject *pyclass, PyCFunction func)
   {
      // Add a getter for the array interface dict to the class.
      Utility::AddToClass(pyclass, "_get__array_interface__", func, METH_NOARGS);
      // Add the dictionary as property to the class so that it updates automatically if accessed.
      // Since we are not able to add a property easily from C++, we do this in Python.

      // We return early if the module does not have the function to add the property, which
      // is the case if cppyy is invoked directly and not through PyROOT.
      if (!PyObject_HasAttrString(gRootModule, "_add__array_interface__")) return;

      auto f = PyObject_GetAttrString(gRootModule, "_add__array_interface__");
      auto r = PyObject_CallFunction(f, (char*)"O", pyclass);
      Py_DECREF(f);
      Py_DECREF(r);
   }

   template <typename T, char typestr>
   PyObject *ArrayInterface(ObjectProxy *self)
   {
      T *cobj = reinterpret_cast<T*>(self->GetObject());

      // Create array interface dict
      auto dict = PyDict_New();

      // Version
      auto pyversion = PyLong_FromLong(3);
      PyDict_SetItemString(dict, "version", pyversion);
      Py_DECREF(pyversion);

      // Type string
      #ifdef R__BYTESWAP
      const char endianess = '<';
#else
      const char endianess = '>';
#endif
      const UInt_t bytes = sizeof(typename T::value_type);
      auto pytypestr = PyROOT_PyUnicode_FromString(TString::Format("%c%c%i", endianess, typestr, bytes).Data());
      PyDict_SetItemString(dict, "typestr", pytypestr);
      Py_DECREF(pytypestr);

      // Shape
      auto pysize = PyLong_FromLong(cobj->size());
      auto pyshape = PyTuple_Pack(1, pysize);
      PyDict_SetItemString(dict, "shape", pyshape);
      Py_DECREF(pysize);
      Py_DECREF(pyshape);

      // Pointer
      auto ptr = reinterpret_cast<unsigned long long>(cobj->data());
      // Numpy breaks for data pointer of 0 even though the array is empty.
      // We set the pointer to 1 but the value itself is arbitrary and never accessed.
      if (cobj->empty()) ptr = 1;
      auto pyptr = PyLong_FromUnsignedLongLong(ptr);
      auto pydata = PyTuple_Pack(2, pyptr, Py_False);
      PyDict_SetItemString(dict, "data", pydata);
      Py_DECREF(pyptr);
      Py_DECREF(pydata);

      return dict;
   }

   //- simplistic len() functions -------------------------------------------------
   PyObject* ReturnThree( ObjectProxy*, PyObject* ) {
      return PyInt_FromLong( 3 );
   }

   PyObject* ReturnTwo( ObjectProxy*, PyObject* ) {
      return PyInt_FromLong( 2 );
   }

} // unnamed namespace


//- public functions -----------------------------------------------------------
Bool_t PyROOT::Pythonize( PyObject* pyclass, const std::string& name )
{
// Add pre-defined pythonizations (for STL and ROOT) to classes based on their
// signature and/or class name.
   if ( pyclass == 0 )
      return kFALSE;

   // add pretty printing
   Utility::AddToClass(pyclass, "__str__", (PyCFunction)ClingPrintValue);

   //- method name based pythonization --------------------------------------------

   // for smart pointer style classes (note fall-through)
   if ( HasAttrDirect( pyclass, PyStrings::gDeref ) ) {
      Utility::AddToClass( pyclass, "__getattr__", (PyCFunction) DeRefGetAttr, METH_O );
   } else if ( HasAttrDirect( pyclass, PyStrings::gFollow ) ) {
      Utility::AddToClass( pyclass, "__getattr__", (PyCFunction) FollowGetAttr, METH_O );
   }

// for STL containers, and user classes modeled after them
   if ( HasAttrDirect( pyclass, PyStrings::gSize ) )
      Utility::AddToClass( pyclass, "__len__", "size" );

// like-wise, some typical container sizings
   if ( HasAttrDirect( pyclass, PyStrings::gGetSize ) )
      Utility::AddToClass( pyclass, "__len__", "GetSize" );

   if ( HasAttrDirect( pyclass, PyStrings::ggetSize ) )
      Utility::AddToClass( pyclass, "__len__", "getSize" );

   if ( HasAttrDirect( pyclass, PyStrings::gBegin ) && HasAttrDirect( pyclass, PyStrings::gEnd ) ) {
   // some classes may not have dicts for their iterators, making begin/end useless
      PyObject* pyfullname = PyObject_GetAttr( pyclass, PyStrings::gCppName );
      if ( ! pyfullname ) pyfullname = PyObject_GetAttr( pyclass, PyStrings::gName );
      TClass* klass = TClass::GetClass( PyROOT_PyUnicode_AsString( pyfullname ) );
      Py_DECREF( pyfullname );

      if (!klass->InheritsFrom(TCollection::Class())) {
         // TCollection has a begin and end method so that they can be used in
         // the C++ range expression.  However, unlike any other use of TIter,
         // TCollection::begin must include the first iteration.  PyROOT is
         // handling TIter as a special case (as it should) and also does this
         // first iteration (via the first call to Next to get the first element)
         // and thus using begin in this case lead to the first element being
         // forgotten by PyROOT.
         // i.e. Don't search for begin in TCollection since we can not use.'

         TMethod* meth = klass->GetMethodAllAny( "begin" );

         TClass* iklass = 0;
         if ( meth ) {
            Int_t oldl = gErrorIgnoreLevel; gErrorIgnoreLevel = 3000;
            iklass = TClass::GetClass( meth->GetReturnTypeNormalizedName().c_str() );
            gErrorIgnoreLevel = oldl;
         }

         if ( iklass && iklass->GetClassInfo() ) {
            ((PyTypeObject*)pyclass)->tp_iter     = (getiterfunc)StlSequenceIter;
            Utility::AddToClass( pyclass, "__iter__", (PyCFunction) StlSequenceIter, METH_NOARGS );
         } else if ( HasAttrDirect( pyclass, PyStrings::gGetItem ) && HasAttrDirect( pyclass, PyStrings::gLen ) ) {
            Utility::AddToClass( pyclass, "_getitem__unchecked", "__getitem__" );
            Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) CheckedGetItem, METH_O );
         }
      }
   }

// search for global comparator overloads (may fail; not sure whether it isn't better to
// do this lazily just as is done for math operators, but this interplays nicely with the
// generic versions)
   Utility::AddBinaryOperator( pyclass, "==", "__eq__" );
   Utility::AddBinaryOperator( pyclass, "!=", "__ne__" );

// map operator==() through GenObjectIsEqual to allow comparison to None (kTRUE is to
// require that the located method is a MethodProxy; this prevents circular calls as
// GenObjectIsEqual is no MethodProxy)
   if ( HasAttrDirect( pyclass, PyStrings::gEq, kTRUE ) ) {
      Utility::AddToClass( pyclass, "__cpp_eq__",  "__eq__" );
      Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) GenObjectIsEqual, METH_O );
   }

// map operator!=() through GenObjectIsNotEqual to allow comparison to None (see note
// on kTRUE above for __eq__)
   if ( HasAttrDirect( pyclass, PyStrings::gNe, kTRUE ) ) {
      Utility::AddToClass( pyclass, "__cpp_ne__",  "__ne__" );
      Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) GenObjectIsNotEqual, METH_O );
   }


//- class name based pythonization ---------------------------------------------

   if ( name == "TObject" ) {
   // support for the 'in' operator
      Utility::AddToClass( pyclass, "__contains__", (PyCFunction) TObjectContains, METH_O );

   // comparing for lists
      Utility::AddToClass( pyclass, "__cmp__", (PyCFunction) TObjectCompare, METH_O );
      Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) TObjectIsEqual, METH_O );
      Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) TObjectIsNotEqual, METH_O );

   }

   else if ( name == "TClass" ) {
   // make DynamicCast return a usable python object, rather than void*
      Utility::AddToClass( pyclass, "_TClass__DynamicCast", "DynamicCast" );
      Utility::AddToClass( pyclass, "DynamicCast", (PyCFunction) TClassDynamicCast );

   // the following cast is easier to use (reads both ways)
      Utility::AddToClass( pyclass, "StaticCast", (PyCFunction) TClassStaticCast );

   }

   else if ( name == "TCollection" ) {
      Utility::AddToClass( pyclass, "append",   "Add" );
      Utility::AddToClass( pyclass, "extend",   (PyCFunction) TCollectionExtend, METH_O );
      Utility::AddToClass( pyclass, "remove",   (PyCFunction) TCollectionRemove, METH_O );
      Utility::AddToClass( pyclass, "__add__",  (PyCFunction) TCollectionAdd, METH_O );
      Utility::AddToClass( pyclass, "__imul__", (PyCFunction) TCollectionIMul, METH_O );
      Utility::AddToClass( pyclass, "__mul__",  (PyCFunction) TCollectionMul, METH_O );
      Utility::AddToClass( pyclass, "__rmul__", (PyCFunction) TCollectionMul, METH_O );

      Utility::AddToClass( pyclass, "count", (PyCFunction) TCollectionCount, METH_O );

      ((PyTypeObject*)pyclass)->tp_iter = (getiterfunc)TCollectionIter;
      Utility::AddToClass( pyclass, "__iter__",  (PyCFunction)TCollectionIter, METH_NOARGS );

   }

   else if ( name == "TSeqCollection" ) {
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) TSeqCollectionGetItem, METH_O );
      Utility::AddToClass( pyclass, "__setitem__", (PyCFunction) TSeqCollectionSetItem );
      Utility::AddToClass( pyclass, "__delitem__", (PyCFunction) TSeqCollectionDelItem, METH_O );

      Utility::AddToClass( pyclass, "insert",  (PyCFunction) TSeqCollectionInsert );
      Utility::AddToClass( pyclass, "pop",     (PyCFunction) TSeqCollectionPop );
      Utility::AddToClass( pyclass, "reverse", (PyCFunction) TSeqCollectionReverse, METH_NOARGS );
      Utility::AddToClass( pyclass, "sort",    (PyCFunction) TSeqCollectionSort,
                           METH_VARARGS | METH_KEYWORDS );

      Utility::AddToClass( pyclass, "index", (PyCFunction) TSeqCollectionIndex, METH_O );

   }

   else if ( name == "TObjArray" ) {
      Utility::AddToClass( pyclass, "__len__", (PyCFunction) TObjArrayLen, METH_NOARGS );
   }

   else if ( name == "TClonesArray" ) {
   // restore base TSeqCollection operator[] to prevent random object creation (it's
   // functionality is equivalent to the operator[](int) const of TClonesArray, but
   // there's no guarantee it'll be selected over the non-const version)
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) TSeqCollectionGetItem, METH_O );

   // this setitem should be used with as much care as the C++ one
      Utility::AddToClass( pyclass, "__setitem__", (PyCFunction) TClonesArraySetItem );

   }

   else if ( IsTemplatedSTLClass( name, "vector" ) || (name.find("ROOT::VecOps::RVec<") == 0) ) {

      if ( HasAttrDirect( pyclass, PyStrings::gLen ) && HasAttrDirect( pyclass, PyStrings::gAt ) ) {
         Utility::AddToClass( pyclass, "_vector__at", "at" );
      // remove iterator that was set earlier (checked __getitem__ will do the trick)
         if ( HasAttrDirect( pyclass, PyStrings::gIter ) )
            PyObject_DelAttr( pyclass, PyStrings::gIter );
      } else if ( HasAttrDirect( pyclass, PyStrings::gGetItem ) ) {
         Utility::AddToClass( pyclass, "_vector__at", "__getitem__" );   // unchecked!
      }

   // vector-optimized iterator protocol
   // Breaks for vector of bool since the "data" member function is not required in the STL.
      if ( name.find("vector<bool>") == std::string::npos && name.find("RVec<bool>") == std::string::npos ) {
         ((PyTypeObject*)pyclass)->tp_iter     = (getiterfunc)vector_iter;
      }

   // helpers for iteration
      TypedefInfo_t* ti = gInterpreter->TypedefInfo_Factory( (name+"::value_type").c_str() );
      if ( gInterpreter->TypedefInfo_IsValid( ti ) ) {
         PyObject* pyvalue_size = PyLong_FromLong( gInterpreter->TypedefInfo_Size( ti ) );
         PyObject_SetAttrString( pyclass, "value_size", pyvalue_size );
         Py_DECREF( pyvalue_size );

         PyObject* pyvalue_type = PyROOT_PyUnicode_FromString( gInterpreter->TypedefInfo_TrueName( ti ) );
         PyObject_SetAttrString( pyclass, "value_type", pyvalue_type );
         Py_DECREF( pyvalue_type );
      }
      gInterpreter->TypedefInfo_Delete( ti );

   // provide a slice-able __getitem__, if possible
      if ( HasAttrDirect( pyclass, PyStrings::gVectorAt ) )
         Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) VectorGetItem, METH_O );

   // std::vector<bool> is a special case in C++
      std::string::size_type pos = name.find( "vector<bool" ); // to cover all variations
      if ( pos == 0 /* at beginning */ || pos == 5 /* after std:: */ ) {
         Utility::AddToClass( pyclass, "__setitem__", (PyCFunction) VectorBoolSetItem );
      }

      // add array interface for STL vectors
      if (name.find("ROOT::VecOps::RVec<") == 0) {
      } else if (name == "vector<float>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<float>, 'f'>);
      } else if (name == "vector<double>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<double>, 'f'>);
      } else if (name == "vector<int>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<int>, 'i'>);
      } else if (name == "vector<unsigned int>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<unsigned int>, 'u'>);
      } else if (name == "vector<long>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<long>, 'i'>);
      } else if (name == "vector<unsigned long>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<std::vector<unsigned long>, 'u'>);
      }

      // add array interface for RVecs
      if (name.find("ROOT::VecOps::RVec<") != 0) {
      } else if (name == "ROOT::VecOps::RVec<float>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<float>, 'f'>);
      } else if (name == "ROOT::VecOps::RVec<double>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<double>, 'f'>);
      } else if (name == "ROOT::VecOps::RVec<int>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<int>, 'i'>);
      } else if (name == "ROOT::VecOps::RVec<unsigned int>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<unsigned int>, 'u'>);
      } else if (name == "ROOT::VecOps::RVec<long>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<long>, 'i'>);
      } else if (name == "ROOT::VecOps::RVec<unsigned long>") {
         AddArrayInterface(pyclass, (PyCFunction)ArrayInterface<ROOT::VecOps::RVec<unsigned long>, 'u'>);
      }
   }

   else if ( IsTemplatedSTLClass( name, "map" ) ) {
      Utility::AddToClass( pyclass, "__contains__", (PyCFunction) MapContains, METH_O );

   }

   else if ( IsTemplatedSTLClass( name, "pair" ) ) {
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) PairUnpack, METH_O );
      Utility::AddToClass( pyclass, "__len__", (PyCFunction) ReturnTwo, METH_NOARGS );

   }

   else if ( name.find( "iterator" ) != std::string::npos ) {
      ((PyTypeObject*)pyclass)->tp_iternext = (iternextfunc)StlIterNext;
      Utility::AddToClass( pyclass, PYROOT__next__, (PyCFunction) StlIterNext, METH_NOARGS );

   // special case, if operator== is a global overload and included in the dictionary
      if ( ! HasAttrDirect( pyclass, PyStrings::gCppEq, kTRUE ) )
         Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) StlIterIsEqual, METH_O );
      if ( ! HasAttrDirect( pyclass, PyStrings::gCppNe, kTRUE ) )
         Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) StlIterIsNotEqual, METH_O );

   }

   else if ( name == "string" || name == "std::string" ) {
      Utility::AddToClass( pyclass, "__repr__", (PyCFunction) StlStringRepr, METH_NOARGS );
      Utility::AddToClass( pyclass, "__str__", "c_str" );
      Utility::AddToClass( pyclass, "__cmp__", (PyCFunction) StlStringCompare, METH_O );
      Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) StlStringIsEqual, METH_O );
      Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) StlStringIsNotEqual, METH_O );

   }

   else if ( name == "TString" ) {
      Utility::AddToClass( pyclass, "__repr__", (PyCFunction) TStringRepr, METH_NOARGS );
      Utility::AddToClass( pyclass, "__str__", "Data" );
      Utility::AddToClass( pyclass, "__len__", "Length" );

      Utility::AddToClass( pyclass, "__cmp__", "CompareTo" );
      Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) TStringIsEqual, METH_O );
      Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) TStringIsNotEqual, METH_O );

   }

   else if ( name == "TObjString" ) {
      Utility::AddToClass( pyclass, "__repr__", (PyCFunction) TObjStringRepr, METH_NOARGS );
      Utility::AddToClass( pyclass, "__str__",  "GetName" );
      Utility::AddToClass( pyclass, "__len__",  (PyCFunction) TObjStringLength, METH_NOARGS );

      Utility::AddToClass( pyclass, "__cmp__", (PyCFunction) TObjStringCompare, METH_O );
      Utility::AddToClass( pyclass, "__eq__",  (PyCFunction) TObjStringIsEqual, METH_O );
      Utility::AddToClass( pyclass, "__ne__",  (PyCFunction) TObjStringIsNotEqual, METH_O );

   }

   else if ( name == "TIter" ) {
      ((PyTypeObject*)pyclass)->tp_iter     = (getiterfunc)PyObject_SelfIter;
      Utility::AddToClass( pyclass, "__iter__", (PyCFunction) PyObject_SelfIter, METH_NOARGS );

      ((PyTypeObject*)pyclass)->tp_iternext = (iternextfunc)TIterNext;
      Utility::AddToClass( pyclass, PYROOT__next__, (PyCFunction) TIterNext, METH_NOARGS );

   }

   else if ( name == "TDirectory" ) {
   // note: this replaces the already existing TDirectory::GetObject()
      Utility::AddToClass( pyclass, "GetObject", (PyCFunction) TDirectoryGetObject );

   // note: this replaces the already existing TDirectory::WriteObject()
      Utility::AddToClass( pyclass, "WriteObject", (PyCFunction) TDirectoryWriteObject );

   }

   else if ( name == "TDirectoryFile" ) {
   // add safety for non-TObject derived Get() results
      Utility::AddToClass( pyclass, "Get", (PyCFunction) TDirectoryFileGet,     METH_O );

      // Re-inject here too, since TDirectoryFile redefines GetObject
      Utility::AddToClass(pyclass, "GetObject", (PyCFunction)TDirectoryGetObject);

      return kTRUE;
   }

   else if ( name == "TTree" ) {
   // allow direct browsing of the tree
      Utility::AddToClass( pyclass, "__getattr__", (PyCFunction) TTreeGetAttr, METH_O );

   // workaround for templated member Branch()
      MethodProxy* original =
         (MethodProxy*)PyObject_GetAttrFromDict( pyclass, PyStrings::gBranch );
      MethodProxy* method = MethodProxy_New( "Branch", new TTreeBranch( original ) );
      Py_DECREF( original ); original = 0;

      PyObject_SetAttrString(
         pyclass, const_cast< char* >( method->GetName().c_str() ), (PyObject*)method );
      Py_DECREF( method ); method = 0;

   // workaround for templated member SetBranchAddress()
      original = (MethodProxy*)PyObject_GetAttrFromDict( pyclass, PyStrings::gSetBranchAddress );
      method = MethodProxy_New( "SetBranchAddress", new TTreeSetBranchAddress( original ) );
      Py_DECREF( original ); original = 0;

      PyObject_SetAttrString(
         pyclass, const_cast< char* >( method->GetName().c_str() ), (PyObject*)method );
      Py_DECREF( method ); method = 0;

   }

   else if ( name == "TChain" ) {
   // allow SetBranchAddress to take object directly, w/o needing AddressOf()
      MethodProxy* original =
         (MethodProxy*)PyObject_GetAttrFromDict( pyclass, PyStrings::gSetBranchAddress );
      MethodProxy* method = MethodProxy_New( "SetBranchAddress", new TChainSetBranchAddress( original ) );
      Py_DECREF( original ); original = 0;

      PyObject_SetAttrString(
         pyclass, const_cast< char* >( method->GetName().c_str() ), (PyObject*)method );
      Py_DECREF( method ); method = 0;

   }

   else if ( name.find("ROOT::RDataFrame") == 0 || name.find("ROOT::RDF::RInterface<") == 0 ) {
      if (PyObject_HasAttrString( gRootModule, "_RDataFrameAsNumpy" )) {
         PyObject_SetAttrString(pyclass, "AsNumpy",
                                PyObject_GetAttrString( gRootModule, "_RDataFrameAsNumpy" ));
      }
   }

   else if ( name == "TStyle" ) {
      MethodProxy* ctor = (MethodProxy*)PyObject_GetAttr( pyclass, PyStrings::gInit );
      ctor->fMethodInfo->fFlags &= ~TCallContext::kIsCreator;
      Py_DECREF( ctor );
   }

   else if ( name == "TH1" )       // allow hist *= scalar
      Utility::AddToClass( pyclass, "__imul__", (PyCFunction) THNIMul, METH_O );

   else if ( name == "TF1" )       // allow instantiation with python callable
      Utility::AddToClass( pyclass, "__init__", new TF1InitWithPyFunc );

   else if ( name == "TF2" )       // allow instantiation with python callable
      Utility::AddToClass( pyclass, "__init__", new TF2InitWithPyFunc );

   else if ( name == "TF3" )       // allow instantiation with python callable
      Utility::AddToClass( pyclass, "__init__", new TF3InitWithPyFunc );

   else if ( name == "TFunction" ) // allow direct call
      Utility::AddToClass( pyclass, "__call__", (PyCFunction) TFunctionCall );

   else if ( name == "TMinuit" )   // allow call with python callable
      Utility::AddToClass( pyclass, "SetFCN", new TMinuitSetFCN );

   else if ( name == "TFitter" )   // allow call with python callable (this is not correct)
       Utility::AddToClass( pyclass, "SetFCN", new TMinuitFitterSetFCN );

   else if ( name == "Fitter" )    // really Fit::Fitter, allow call with python callable
      Utility::AddToClass( pyclass, "FitFCN", new TFitterFitFCN );

   else if ( name == "TFile" ) {
   // TFile::Open really is a constructor, really
      PyObject* attr = PyObject_GetAttrString( pyclass, (char*)"Open" );
      if ( MethodProxy_Check( attr ) )
         ((MethodProxy*)attr)->fMethodInfo->fFlags |= TCallContext::kIsCreator;
      Py_XDECREF( attr );

   // allow member-style access to entries in file
      Utility::AddToClass( pyclass, "__getattr__", (PyCFunction) TFileGetAttr, METH_O );
   }

   else if ( name.substr(0,8) == "TVector3" ) {
      Utility::AddToClass( pyclass, "__len__", (PyCFunction) ReturnThree, METH_NOARGS );
      Utility::AddToClass( pyclass, "_getitem__unchecked", "__getitem__" );
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) CheckedGetItem, METH_O );

   }

   else if ( name.substr(0,8) == "TVectorT" ) {
      // allow proper iteration
      Utility::AddToClass( pyclass, "__len__", "GetNoElements" );
      Utility::AddToClass( pyclass, "_getitem__unchecked", "__getitem__" );
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) CheckedGetItem, METH_O );
   }

   else if ( name.substr(0,6) == "TArray" && name != "TArray" ) {    // allow proper iteration
      // __len__ is already set from GetSize()
      Utility::AddToClass( pyclass, "_getitem__unchecked", "__getitem__" );
      Utility::AddToClass( pyclass, "__getitem__", (PyCFunction) CheckedGetItem, METH_O );
   }

// Make RooFit 'using' member functions available (not supported by dictionary)
   else if ( name == "RooDataHist" )
      Utility::AddUsingToClass( pyclass, "plotOn" );

   else if ( name == "RooSimultaneous" )
      Utility::AddUsingToClass( pyclass, "plotOn" );
   else if ( name == "RooAbsPdf" ) {
      Utility::AddUsingToClass( pyclass, "createChi2" );
      Utility::AddUsingToClass( pyclass, "chi2FitTo" );
   }
   else if ( name == "RooDataSet" ) {
      Utility::AddUsingToClass( pyclass, "createHistogram" );
   }
// Same for TH2
   else if ( name == "TH2" ) {
      Utility::AddUsingToClass( pyclass, "GetBinErrorUp" );
      Utility::AddUsingToClass( pyclass, "GetBinErrorLow" );
   }

   // TODO: store these on the pythonizations module, not on gRootModule
   // TODO: externalize this code and use update handlers on the python side
   PyObject* userPythonizations = PyObject_GetAttrString( gRootModule, "UserPythonizations" );
   PyObject* pythonizationScope = PyObject_GetAttrString( gRootModule, "PythonizationScope" );

   std::vector< std::string > pythonization_scopes;
   pythonization_scopes.push_back( "__global__" );

   std::string user_scope = PyROOT_PyUnicode_AsString( pythonizationScope );
   if ( user_scope != "__global__" ) {
      if ( PyDict_Contains( userPythonizations, pythonizationScope ) ) {
          pythonization_scopes.push_back( user_scope );
      }
   }

   Bool_t pstatus = kTRUE;

   for ( auto key = pythonization_scopes.cbegin(); key != pythonization_scopes.cend(); ++key ) {
      PyObject* tmp = PyDict_GetItemString( userPythonizations, key->c_str() );
      Py_ssize_t num_pythonizations = PyList_Size( tmp );
      PyObject* arglist = nullptr;
      if ( num_pythonizations )
         arglist = Py_BuildValue( "O,s", pyclass, name.c_str() );
      for ( Py_ssize_t i = 0; i < num_pythonizations; ++i ) {
         PyObject* pythonizor = PyList_GetItem( tmp, i );
      // TODO: detail error handling for the pythonizors
         PyObject* result = PyObject_CallObject( pythonizor, arglist );
         if ( !result ) {
            pstatus = kFALSE;
            break;
         } else
            Py_DECREF( result );
      }
      Py_XDECREF( arglist );
   }

   Py_DECREF( userPythonizations );
   Py_DECREF( pythonizationScope );


// phew! all done ...
   return pstatus;
}
