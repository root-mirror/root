// Author: Wim Lavrijsen, Dec 2006

// Bindings
#include "PyROOT.h"
#include "TCustomPyTypes.h"

#if PY_VERSION_HEX >= 0x03000000
// TODO: this will break functionality
#define PyMethod_GET_CLASS( meth ) Py_None
#endif


namespace PyROOT {

//= float type allowed for reference passing =================================
PyTypeObject TCustomFloat_Type = {     // python float is a C/C++ double
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.double",      // tp_name
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
   (char*)"PyROOT float object for pass by reference",     // tp_doc
   0,                         // tp_traverse
   0,                         // tp_clear
   0,                         // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   0,                         // tp_getset
   &PyFloat_Type,             // tp_base
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

//= long type allowed for reference passing ==================================
PyTypeObject TCustomInt_Type = {       // python int is a C/C++ long
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.long",        // tp_name
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
   &PyInt_Type,               // tp_base
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

//= instancemethod object with a more efficient call function ================
static PyMethodObject* free_list;
static int numfree = 0;
#ifndef PyMethod_MAXFREELIST
#define PyMethod_MAXFREELIST 256
#endif

PyObject* TCustomInstanceMethod_New( PyObject* func, PyObject* self, PyObject*
#if PY_VERSION_HEX < 0x03000000
      pyclass
#endif
   )
{
// from instancemethod, but with custom type (at issue is that instancemethod is not
// meant to be derived from)
   PyMethodObject* im;
   if ( ! PyCallable_Check( func ) ) {
      PyErr_Format( PyExc_SystemError,
                     "%s:%d: bad argument to internal function",
                     __FILE__, __LINE__ );
      return NULL;
   }

   im = free_list;
   if ( im != NULL ) {
      free_list = (PyMethodObject*)( im->im_self );
      (void)PyObject_INIT( im, &TCustomInstanceMethod_Type );
   }
   else {
      im = PyObject_GC_New( PyMethodObject, &TCustomInstanceMethod_Type );
      if ( im == NULL )
         return NULL;
   }

   im->im_weakreflist = NULL;
   Py_INCREF( func );
   im->im_func = func;
   Py_XINCREF( self );
   im->im_self = self;
#if PY_VERSION_HEX < 0x03000000
   Py_XINCREF( pyclass );
   im->im_class = pyclass;
#endif
   PyObject_GC_Track( im );
   return (PyObject*)im;
}

////////////////////////////////////////////////////////////////////////////////
/// from instancemethod, but with custom type (at issue is that instancemethod is not
/// meant to be derived from)

static void im_dealloc( PyMethodObject* im )
{
   PyObject_GC_UnTrack( im );

   if ( im->im_weakreflist != NULL )
      PyObject_ClearWeakRefs( (PyObject*) im );

   Py_DECREF( im->im_func );
   Py_XDECREF( im->im_self );
#if PY_VERSION_HEX < 0x03000000
   Py_XDECREF( im->im_class );
#endif

   if ( numfree < PyMethod_MAXFREELIST ) {
      im->im_self = (PyObject*)free_list;
      free_list = im;
      numfree++;
   } else {
      PyObject_GC_Del(im);
   }
}

////////////////////////////////////////////////////////////////////////////////
/// The mapping from a method to a function involves reshuffling of self back
/// into the list of arguments. However, the pythonized methods will then have
/// to undo that shuffling, which is inefficient. This method is the same as
/// the one for the instancemethod object, except for the shuffling.

static PyObject* im_call( PyObject* meth, PyObject* args, PyObject* kw )
{
   PyObject* self = PyMethod_GET_SELF( meth );

   if ( ! self ) {
   // unbound methods must be called with an instance of the class (or a
   // derived class) as first argument
      Py_ssize_t argc = PyTuple_GET_SIZE( args );
      PyObject* pyclass = PyMethod_GET_CLASS( meth );
      if ( 1 <= argc && PyObject_IsInstance( PyTuple_GET_ITEM( args, 0 ), pyclass ) == 1 ) {
         self = PyTuple_GET_ITEM( args, 0 );

         PyObject* newArgs = PyTuple_New( argc - 1 );
         for ( int i = 1; i < argc; ++i ) {
            PyObject* v = PyTuple_GET_ITEM( args, i );
            Py_INCREF( v );
            PyTuple_SET_ITEM( newArgs, i-1, v );
         }

         args = newArgs;

      } else
         return PyMethod_Type.tp_call( meth, args, kw );    // will set proper error msg

   } else
      Py_INCREF( args );

   PyCFunctionObject* func = (PyCFunctionObject*)PyMethod_GET_FUNCTION( meth );

// the function is globally shared, so set and reset its "self" (ok, b/c of GIL)
   Py_INCREF( self );
   func->m_self = self;
   PyObject* result = PyCFunction_Call( (PyObject*)func, args, kw );
   func->m_self = 0;
   Py_DECREF( self );
   Py_DECREF( args );
   return result;
}

////////////////////////////////////////////////////////////////////////////////
/// from instancemethod: don't rebind an already bound method, or an unbound method
/// of a class that's not a base class of pyclass

static PyObject* im_descr_get( PyObject* meth, PyObject* obj, PyObject* pyclass )
{
   if ( PyMethod_GET_SELF( meth ) != NULL
#if PY_VERSION_HEX < 0x03000000
        || ( PyMethod_GET_CLASS( meth ) != NULL &&
          ! PyObject_IsSubclass( pyclass,  PyMethod_GET_CLASS(meth) ) )
#endif
                                          ) {
      Py_INCREF( meth );
      return meth;
   }

   if ( obj == Py_None )
      obj = NULL;

   return TCustomInstanceMethod_New( PyMethod_GET_FUNCTION( meth ), obj, pyclass );
}

//= PyROOT custom instance method type =======================================
PyTypeObject TCustomInstanceMethod_Type = {
   PyVarObject_HEAD_INIT( &PyType_Type, 0 )
   (char*)"ROOT.InstanceMethod",      // tp_name
   0,                         // tp_basicsize
   0,                         // tp_itemsize
   (destructor)im_dealloc,    // tp_dealloc
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
   im_call,                   // tp_call
   0,                         // tp_str
   0,                         // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_CHECKTYPES |
      Py_TPFLAGS_BASETYPE,    // tp_flags
   (char*)"PyROOT custom instance method (internal)",     // tp_doc
   0,                         // tp_traverse
   0,                         // tp_clear
   0,                         // tp_richcompare
   0,                         // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   0,                         // tp_getset
   &PyMethod_Type,            // tp_base
   0,                         // tp_dict
   im_descr_get,              // tp_descr_get
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
