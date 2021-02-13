// Author: Wim Lavrijsen, Dec 2006

#ifndef PYROOT_TCUSTOMPYTYPES_H
#define PYROOT_TCUSTOMPYTYPES_H

// ROOT
#include "DllImport.h"


namespace PyROOT {

/** Custom builtins, detectable by type, for pass by ref
      @author  WLAV
      @date    12/13/2006
      @version 1.0
 */

//- custom float object type and type verification ---------------------------
   R__EXTERN PyTypeObject TCustomFloat_Type;

   template< typename T >
   inline Bool_t TCustomFloat_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &TCustomFloat_Type );
   }

   template< typename T >
   inline Bool_t TCustomFloat_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &TCustomFloat_Type;
   }

//- custom long object type and type verification ----------------------------
   R__EXTERN PyTypeObject TCustomInt_Type;

   template< typename T >
   inline Bool_t TCustomInt_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &TCustomInt_Type );
   }

   template< typename T >
   inline Bool_t TCustomInt_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &TCustomInt_Type;
   }

//- custom instance method object type and type verification -----------------
   R__EXTERN PyTypeObject TCustomInstanceMethod_Type;

   template< typename T >
   inline Bool_t TCustomInstanceMethod_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &TCustomInstanceMethod_Type );
   }

   template< typename T >
   inline Bool_t TCustomInstanceMethod_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &TCustomInstanceMethod_Type;
   }

   PyObject* TCustomInstanceMethod_New( PyObject* func, PyObject* self, PyObject* pyclass );

} // namespace PyROOT

#endif // !PYROOT_TCUSTOMPYTYPES_H
