#ifndef PYROOT_TTUPLEOFINSTANCES_H
#define PYROOT_TTUPLEOFINSTANCES_H

// ROOT
#include "DllImport.h"


namespace PyROOT {

/** Representation of C-style array of instances
      @author  WLAV
      @date    02/10/2014
      @version 1.0
 */

//- custom tuple type that can pass through C-style arrays -------------------
   R__EXTERN PyTypeObject TTupleOfInstances_Type;

   template< typename T >
   inline Bool_t TTupleOfInstances_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &TTupleOfInstances_Type );
   }

   template< typename T >
   inline Bool_t TTupleOfInstances_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &TTupleOfInstances_Type;
   }

   PyObject* TTupleOfInstances_New(
      Cppyy::TCppObject_t address, Cppyy::TCppType_t klass, Py_ssize_t nelems );

} // namespace PyROOT

#endif // !PYROOT_TTUPLEOFINSTANCES_H
