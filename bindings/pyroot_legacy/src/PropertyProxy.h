// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

#ifndef PYROOT_PROPERTYPROXY_H
#define PYROOT_PROPERTYPROXY_H

// Bindings
#include "Converters.h"

// ROOT
#include "DllImport.h"
#include "TClassRef.h"
class TDataMember;
class TEnumConstant;
class TGlobal;

// Standard
#include <string>


namespace PyROOT {

   class ObjectProxy;

   class PropertyProxy {
   public:
      void Set( Cppyy::TCppScope_t scope, Cppyy::TCppIndex_t idata );
      void Set( Cppyy::TCppScope_t scope, const std::string& name, void* address, TEnum* en );

      std::string GetName() { return fName; }
      void* GetAddress( ObjectProxy* pyobj /* owner */ );

   public:               // public, as the python C-API works with C structs
      PyObject_HEAD
      ptrdiff_t          fOffset;
      Long_t             fProperty;
      TConverter*        fConverter;
      Cppyy::TCppScope_t fEnclosingScope;
      std::string        fName;

   private:              // private, as the python C-API will handle creation
      PropertyProxy() {}
   };


//- property proxy type and type verification --------------------------------
   R__EXTERN PyTypeObject PropertyProxy_Type;

   template< typename T >
   inline Bool_t PropertyProxy_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &PropertyProxy_Type );
   }

   template< typename T >
   inline Bool_t PropertyProxy_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &PropertyProxy_Type;
   }

//- creation -----------------------------------------------------------------
   inline PropertyProxy* PropertyProxy_New(
      Cppyy::TCppScope_t scope, Cppyy::TCppIndex_t idata )
   {
   // Create an initialize a new property descriptor, given the C++ datum.
      PropertyProxy* pyprop =
         (PropertyProxy*)PropertyProxy_Type.tp_new( &PropertyProxy_Type, 0, 0 );
      pyprop->Set( scope, idata );
      return pyprop;
   }

   inline PropertyProxy* PropertyProxy_NewConstant(
      Cppyy::TCppScope_t scope, const std::string& name, void* address, TEnum* en )
   {
   // Create an initialize a new property descriptor, given the C++ datum.
      PropertyProxy* pyprop =
         (PropertyProxy*)PropertyProxy_Type.tp_new( &PropertyProxy_Type, 0, 0 );
      pyprop->Set( scope, name, address, en );
      return pyprop;
   }

} // namespace PyROOT

#endif // !PYROOT_PROPERTYPROXY_H
