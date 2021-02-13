// Author: Wim Lavrijsen, Dec 2006

// Bindings
#include "PyROOT.h"
#include "TTupleOfInstances.h"
#include "RootWrapper.h"


namespace PyROOT {

//= support for C-style arrays of objects ====================================
PyObject* TTupleOfInstances_New(
      Cppyy::TCppObject_t address, Cppyy::TCppType_t klass, Py_ssize_t nelems )
{
// TODO: the extra copy is inefficient, but it appears that the only way to
// initialize a subclass of a tuple is through a sequence
   PyObject* tup = PyTuple_New( nelems );
   for ( int i = 0; i < nelems; ++i ) {
   // TODO: there's an assumption here that there is no padding, which is bound
   // to be incorrect in certain cases
      PyTuple_SetItem( tup, i,
         BindCppObject( (char*)address + i*Cppyy::SizeOf( klass ), klass, kFALSE /* isRef */ ) );
   // Note: objects are bound as pointers, yet since the pointer value stays in
   // place, updates propagate just as if they were bound by-reference
   }

   PyObject* args = PyTuple_New( 1 );
   Py_INCREF( tup ); PyTuple_SET_ITEM( args, 0, tup );
   PyObject* arr = PyTuple_Type.tp_new( &TTupleOfInstances_Type, args, NULL );
   if ( PyErr_Occurred() ) PyErr_Print();

   Py_DECREF( args );
   // tup ref eaten by SET_ITEM on args

   return arr;
}

//= PyROOT custom tuple-like array type ======================================
PyTypeObject TTupleOfInstances_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.InstancesArray",        // tp_name
   0,                         // tp_basicsize
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
   0,                         // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES |
      Py_TPFLAGS_BASETYPE,    // tp_flags
   (char*)"PyROOT long object for pass by reference",      // tp_doc
   0,                         // tp_traverse
   0,                         // tp_clear
   0,                         // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   0,                         // tp_getset
   &PyTuple_Type,             // tp_base
   0,                         // tp_dict
   0,                         // tp_descr_get
   0,                         // tp_descr_set
   0,                         // tp_dictoffset
   0,                         // tp_init
   0,                         // tp_alloc
   0,                         // tp_new
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
