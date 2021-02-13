// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

#ifndef PYROOT_PYROOTTYPE_H
#define PYROOT_PYROOTTYPE_H

// ROOT
#include "DllImport.h"

#if PY_MAJOR_VERSION == 2 && PY_MINOR_VERSION == 2

// In p2.2, PyHeapTypeObject is not yet part of the interface
#include "structmember.h"

typedef struct {
   PyTypeObject type;
   PyNumberMethods as_number;
   PySequenceMethods as_sequence;
   PyMappingMethods as_mapping;
   PyBufferProcs as_buffer;
   PyObject *name, *slots;
   PyMemberDef members[1];
} PyHeapTypeObject;

#endif


namespace PyROOT {

/** Type object to hold TClassRef instance (this is only semantically a presentation
    of PyRootType instances, not in a C++ sense)
      @author  WLAV
      @date    03/28/2008
      @version 1.0
 */

   class PyRootClass {
   public:
      PyHeapTypeObject fType;      // placeholder, in a single block with the TClassRef
      Cppyy::TCppType_t fCppType;

   private:
      PyRootClass() {}
   };

//- metatype type and type verification --------------------------------------
   R__EXTERN PyTypeObject PyRootType_Type;

   template< typename T >
   inline Bool_t PyRootType_Check( T* object )
   {
      return object && PyObject_TypeCheck( object, &PyRootType_Type );
   }

   template< typename T >
   inline Bool_t PyRootType_CheckExact( T* object )
   {
      return object && Py_TYPE(object) == &PyRootType_Type;
   }

} // namespace PyROOT

#endif // !PYROOT_PYROOTTYPE_H
