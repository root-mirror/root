// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

#ifndef PYROOT_PYROOT_H
#define PYROOT_PYROOT_H

#ifdef _WIN32
// Disable warning C4275: non dll-interface class
#pragma warning ( disable : 4275 )
// Disable warning C4251: needs to have dll-interface to be used by clients
#pragma warning ( disable : 4251 )
// Disable warning C4800: 'int' : forcing value to bool
#pragma warning ( disable : 4800 )
// Avoid that pyconfig.h decides using a #pragma what library python library to use
//#define MS_NO_COREDLL 1
#endif

// to prevent problems with fpos_t and redefinition warnings
#if defined(linux)

#include <stdio.h>

#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif

#ifdef _FILE_OFFSET_BITS
#undef _FILE_OFFSET_BITS
#endif

#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif

#endif // linux


#include "Python.h"
#include "Rtypes.h"

// for 3.3 support
#if PY_VERSION_HEX < 0x03030000
   typedef PyDictEntry* (*dict_lookup_func) ( PyDictObject*, PyObject*, Long_t );
#else
   struct PyDictKeyEntry;
   typedef PyDictKeyEntry* (*dict_lookup_func) ( PyDictObject*, PyObject*, Py_hash_t, PyObject*** );
#define PyDictEntry PyDictKeyEntry
#endif

// for 3.3 support
#if PY_VERSION_HEX < 0x03030000
// Use PyString_Size for a version lower than 3.0
#if PY_VERSION_HEX < 0x03000000
#define PyROOT_PyUnicode_GetSize            PyString_Size
#else
#define PyROOT_PyUnicode_GetSize            PyUnicode_GetSize
#endif
#else // if greater or equal than 3.3
#define PyROOT_PyUnicode_GetSize            PyUnicode_GetLength
#endif

// for 3.0 support (backwards compatibility, really)
#if PY_VERSION_HEX < 0x03000000
#define PyBytes_Check                  PyString_Check
#define PyBytes_CheckExact             PyString_CheckExact
#define PyBytes_AS_STRING              PyString_AS_STRING
#define PyBytes_AsString               PyString_AsString
#define PyBytes_GET_SIZE               PyString_GET_SIZE
#define PyBytes_Size                   PyString_Size
#define PyBytes_FromFormat             PyString_FromFormat
#define PyBytes_FromString             PyString_FromString
#define PyBytes_FromStringAndSize      PyString_FromStringAndSize

#define PyBytes_Type    PyString_Type

#define PyROOT_PyUnicode_Check              PyString_Check
#define PyROOT_PyUnicode_CheckExact         PyString_CheckExact
#define PyROOT_PyUnicode_AsString           PyString_AS_STRING
#define PyROOT_PyUnicode_AsStringChecked    PyString_AsString
#define PyROOT_PyUnicode_GET_SIZE           PyString_GET_SIZE
#define PyROOT_PyUnicode_FromFormat         PyString_FromFormat
#define PyROOT_PyUnicode_FromString         PyString_FromString
#define PyROOT_PyUnicode_InternFromString   PyString_InternFromString
#define PyROOT_PyUnicode_Append             PyString_Concat
#define PyROOT_PyUnicode_AppendAndDel       PyString_ConcatAndDel
#define PyROOT_PyUnicode_FromStringAndSize  PyString_FromStringAndSize

#define PyROOT_PyUnicode_Type PyString_Type

static inline PyObject* PyROOT_PyCapsule_New( void* cobj, const char* /* name */, void (*destr)(void *) )
{
   return PyCObject_FromVoidPtr( cobj, destr );
}
#define PyROOT_PyCapsule_CheckExact    PyCObject_Check
static inline void* PyROOT_PyCapsule_GetPointer( PyObject* capsule, const char* /* name */ )
{
   return (void*)PyCObject_AsVoidPtr( capsule );
}

#define PYROOT__long__ "__long__"
#define PYROOT__idiv__ "__idiv__"
#define PYROOT__div__  "__div__"
#define PYROOT__next__ "next"

#endif  // ! 3.0

// Support non-ASCII strings
#if PY_VERSION_HEX >= 0x03030000
#define PyROOT_PyUnicode_AsStringAndSize    _PyUnicode_AsStringAndSize
#endif

// for 3.0 support (backwards compatibility, really)
#if PY_VERSION_HEX >= 0x03000000
#define PyROOT_PyUnicode_Check              PyUnicode_Check
#define PyROOT_PyUnicode_CheckExact         PyUnicode_CheckExact
#define PyROOT_PyUnicode_AsString           _PyUnicode_AsString
#define PyROOT_PyUnicode_AsStringChecked    _PyUnicode_AsString
#define PyROOT_PyUnicode_GET_SIZE           PyUnicode_GET_SIZE
#define PyROOT_PyUnicode_FromFormat         PyUnicode_FromFormat
#define PyROOT_PyUnicode_FromString         PyUnicode_FromString
#define PyROOT_PyUnicode_InternFromString   PyUnicode_InternFromString
#define PyROOT_PyUnicode_Append             PyUnicode_Append
#define PyROOT_PyUnicode_AppendAndDel       PyUnicode_AppendAndDel
#define PyROOT_PyUnicode_FromStringAndSize  PyUnicode_FromStringAndSize

#define PyROOT_PyUnicode_Type PyUnicode_Type

#define PyIntObject          PyLongObject
#define PyInt_Check          PyLong_Check
#define PyInt_AsLong         PyLong_AsLong
#define PyInt_AS_LONG        PyLong_AsLong
#define PyInt_AsSsize_t      PyLong_AsSsize_t
#define PyInt_CheckExact     PyLong_CheckExact
#define PyInt_FromLong       PyLong_FromLong
#define PyInt_FromSsize_t    PyLong_FromSsize_t

#define PyInt_Type      PyLong_Type

#define PyROOT_PyCapsule_New           PyCapsule_New
#define PyROOT_PyCapsule_CheckExact    PyCapsule_CheckExact
#define PyROOT_PyCapsule_GetPointer    PyCapsule_GetPointer

#define PYROOT__long__ "__int__"
#define PYROOT__idiv__ "__itruediv__"
#define PYROOT__div__  "__truediv__"
#define PYROOT__next__ "__next__"

#define Py_TPFLAGS_HAVE_RICHCOMPARE 0
#define Py_TPFLAGS_CHECKTYPES 0

#define PyClass_Check   PyType_Check

#define PyBuffer_Type   PyMemoryView_Type
#endif  // ! 3.0

#if PY_VERSION_HEX >= 0x03020000
#define PyROOT_PySliceCast   PyObject*
#else
#define PyROOT_PySliceCast   PySliceObject*
#endif  // >= 3.2

// feature of 3.0 not in 2.5 and earlier
#if PY_VERSION_HEX < 0x02060000
#define PyVarObject_HEAD_INIT(type, size)       \
    PyObject_HEAD_INIT(type) size,
#define Py_TYPE(ob)             (((PyObject*)(ob))->ob_type)
#endif

// backwards compatibility, pre python 2.5
#if PY_VERSION_HEX < 0x02050000
typedef int Py_ssize_t;
#define PyInt_AsSsize_t PyInt_AsLong
#define PyInt_FromSsize_t PyInt_FromLong
# define PY_SSIZE_T_FORMAT "%d"
# if !defined(PY_SSIZE_T_MIN)
#  define PY_SSIZE_T_MAX INT_MAX
#  define PY_SSIZE_T_MIN INT_MIN
# endif
#define ssizeobjargproc intobjargproc
#define lenfunc         inquiry
#define ssizeargfunc    intargfunc

#define PyIndex_Check(obj) \
   (PyInt_Check(obj) || PyLong_Check(obj))

inline Py_ssize_t PyNumber_AsSsize_t( PyObject* obj, PyObject* ) {
   return (Py_ssize_t)PyLong_AsLong( obj );
}

#else
# ifdef R__MACOSX
#  if SIZEOF_SIZE_T == SIZEOF_INT
#    if defined(MAC_OS_X_VERSION_10_4)
#       define PY_SSIZE_T_FORMAT "%ld"
#    else
#       define PY_SSIZE_T_FORMAT "%d"
#    endif
#  elif SIZEOF_SIZE_T == SIZEOF_LONG
#    define PY_SSIZE_T_FORMAT "%ld"
#  endif
# else
#  define PY_SSIZE_T_FORMAT "%zd"
# endif
#endif

#if PY_VERSION_HEX < 0x02020000
#define PyBool_FromLong  PyInt_FromLong
#endif

#if PY_VERSION_HEX < 0x03000000
// the following should quiet Solaris
#ifdef Py_False
#undef Py_False
#define Py_False ( (PyObject*)(void*)&_Py_ZeroStruct )
#endif

#ifdef Py_True
#undef Py_True
#define Py_True ( (PyObject*)(void*)&_Py_TrueStruct )
#endif
#endif

// C++ version of the cppyy API
#include "Cppyy.h"

#endif // !PYROOT_PYROOT_H
