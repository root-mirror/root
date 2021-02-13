// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

#ifndef PYROOT_OBJECTPROXY_H
#define PYROOT_OBJECTPROXY_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// PyROOT::ObjectProxy                                                      //
//                                                                          //
// Python-side proxy, encapsulaties a C++ object.                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


// Bindings
#include "PyRootType.h"
#include "Cppyy.h"
#include "TCallContext.h"

// ROOT
#include "DllImport.h"

// TODO: have an ObjectProxy derived or alternative type for smart pointers

namespace PyROOT {

   class ObjectProxy {
   public:
      enum EFlags { kNone = 0x0, kIsOwner = 0x0001, kIsReference = 0x0002, kIsValue = 0x0004, kIsSmartPtr = 0x0008 };

   public:
      void Set( void* address, EFlags flags = kNone )
      {
      // Initialize the proxy with the pointer value 'address.'
         fObject = address;
         fFlags  = flags;
      }

      void SetSmartPtr ( void* address, Cppyy::TCppType_t ptrType )
      {
        fFlags |= kIsSmartPtr;
        fSmartPtr = address;
        fSmartPtrType = ptrType;
      }

      void* GetObject() const
      {
      // Retrieve a pointer to the held C++ object.

      // We get the raw pointer from the smart pointer each time, in case
      // it has changed or has been freed.
         if ( fFlags & kIsSmartPtr ) {
         // TODO: this is icky and slow
            std::vector< Cppyy::TCppMethod_t > methods = Cppyy::GetMethodsFromName( fSmartPtrType, "operator->", /*bases?*/ true);
            std::vector<TParameter> args;
            return Cppyy::CallR( methods[0], fSmartPtr, &args );
         }

         if ( fObject && ( fFlags & kIsReference ) )
            return *(reinterpret_cast< void** >( const_cast< void* >( fObject ) ));
         else
            return const_cast< void* >( fObject );          // may be null
      }

      Cppyy::TCppType_t ObjectIsA() const
      {
      // Retrieve a pointer to the C++ type; may return NULL.
         return ((PyRootClass*)Py_TYPE(this))->fCppType;
      }

      void HoldOn() { fFlags |= kIsOwner; }
      void Release() { fFlags &= ~kIsOwner; }

   public:               // public, as the python C-API works with C structs
      PyObject_HEAD
      void*     fObject;
      int       fFlags;
      void*     fSmartPtr;
      Cppyy::TCppType_t fSmartPtrType;

   private:              // private, as the python C-API will handle creation
      ObjectProxy() {}
   };


//- object proxy type and type verification ----------------------------------
   R__EXTERN PyTypeObject ObjectProxy_Type;

   template< typename T >
   inline Bool_t ObjectProxy_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &ObjectProxy_Type );
   }

   template< typename T >
   inline Bool_t ObjectProxy_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &ObjectProxy_Type;
   }


//- helper for memory regulation (no PyTypeObject equiv. member in p2.2) -----
   void op_dealloc_nofree( ObjectProxy* );

} // namespace PyROOT

#endif // !PYROOT_OBJECTPROXY_H
