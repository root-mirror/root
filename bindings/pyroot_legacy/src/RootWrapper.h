// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

#ifndef PYROOT_ROOTWRAPPER_H
#define PYROOT_ROOTWRAPPER_H

// ROOT
class TGlobal;

// Standard
#include <string>


namespace PyROOT {

// initialize ROOT
   void InitRoot();

// construct a Python shadow class for the named C++ class
   PyObject* GetScopeProxy( Cppyy::TCppScope_t );
   PyObject* CreateScopeProxy( Cppyy::TCppScope_t );
   PyObject* CreateScopeProxy( PyObject*, PyObject* args );
   PyObject* CreateScopeProxy(
      const std::string& scope_name, PyObject* parent = 0 );

// convenience function to retrieve global variables and enums
   PyObject* GetCppGlobal( const std::string& name );
   PyObject* GetCppGlobal( PyObject*, PyObject* args );

// clean up all objects controlled by TMemoryRegulator
   PyObject *ClearProxiedObjects();

// bind a ROOT object into a Python object
   PyObject* BindCppObjectNoCast( Cppyy::TCppObject_t object, Cppyy::TCppType_t klass,
      Bool_t isRef = kFALSE, Bool_t isValue = kFALSE );
   PyObject* BindCppObject(
      Cppyy::TCppObject_t object, Cppyy::TCppType_t klass, Bool_t isRef = kFALSE );
   inline PyObject* BindCppObject(
      Cppyy::TCppObject_t object, const std::string& clName, Bool_t isRef = kFALSE )
   {
      return BindCppObject( object, Cppyy::GetScope( clName ), isRef );
   }

   PyObject* BindCppObjectArray( Cppyy::TCppObject_t address, Cppyy::TCppType_t klass, Int_t size );
   PyObject* BindCppGlobal( TGlobal* );

} // namespace PyROOT

#endif // !PYROOT_ROOTWRAPPER_H
