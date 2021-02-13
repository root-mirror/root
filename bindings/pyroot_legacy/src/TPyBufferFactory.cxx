// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Apr 2004

// Bindings
#include "PyROOT.h"
#include "TPyBufferFactory.h"

// Standard
#include <map>

#if PY_VERSION_HEX >= 0x03000000
static PyObject* PyBuffer_FromReadWriteMemory( void* ptr, int size ) {
#if PY_VERSION_HEX > 0x03000000
   if ( !ptr ) {        // p3 will set an exception if nullptr, just rely on size == 0
      static long dummy[1];
      ptr = dummy;
      size = 0;
   }
#endif

   Py_buffer bufinfo = { ptr, NULL, size, 1, 0, 1, NULL, NULL, NULL, NULL,
#if PY_VERSION_HEX < 0x03030000
      { 0, 0 },
#endif
      NULL };
   return PyMemoryView_FromBuffer( &bufinfo );
}
#endif


//- data ---------------------------------------------------------------------
namespace {

// top of buffer (rest of buffer object has changed across python versions)
   struct PyBufferTop_t {
      PyObject_HEAD
      PyObject*  fBase;            // b_base in python
      void*      fPtr;             // b_ptr in python
      Py_ssize_t fSize;            // b_size in python
      Py_ssize_t fItemSize;        // b_itemsize in python
   };

// callable cache
   std::map< PyObject*, PyObject* > gSizeCallbacks;

// create buffer types and copy methods to be adjusted
#define PYROOT_PREPARE_PYBUFFER_TYPE( name )                                 \
   PyTypeObject      Py##name##Buffer_Type;                                  \
   PySequenceMethods Py##name##Buffer_SeqMethods = *(PyBuffer_Type.tp_as_sequence);\
   PyMappingMethods  Py##name##Buffer_MapMethods;

   PYROOT_PREPARE_PYBUFFER_TYPE( Bool )
   PYROOT_PREPARE_PYBUFFER_TYPE( Char )
   PYROOT_PREPARE_PYBUFFER_TYPE( UChar )
   PYROOT_PREPARE_PYBUFFER_TYPE( Short )
   PYROOT_PREPARE_PYBUFFER_TYPE( UShort )
   PYROOT_PREPARE_PYBUFFER_TYPE( Int )
   PYROOT_PREPARE_PYBUFFER_TYPE( UInt )
   PYROOT_PREPARE_PYBUFFER_TYPE( Long )
   PYROOT_PREPARE_PYBUFFER_TYPE( ULong )
   PYROOT_PREPARE_PYBUFFER_TYPE( Float )
   PYROOT_PREPARE_PYBUFFER_TYPE( Double )

// implement get, str, and length functions
   Py_ssize_t buffer_length( PyObject* self )
   {
   // Retrieve the (type-strided) size of the the buffer; may be a guess.
#if PY_VERSION_HEX < 0x03000000
      Py_ssize_t nlen = ((PyBufferTop_t*)self)->fSize;
      Py_ssize_t item = ((PyBufferTop_t*)self)->fItemSize;
#else
      Py_buffer* bufinfo = PyMemoryView_GET_BUFFER(self);
      Py_ssize_t nlen = bufinfo->len;
      Py_ssize_t item = bufinfo->itemsize;
#endif
      if ( nlen != INT_MAX )  // INT_MAX is the default, i.e. unknown actual length
         return nlen/item;

      std::map< PyObject*, PyObject* >::iterator iscbp = gSizeCallbacks.find( self );
      if ( iscbp != gSizeCallbacks.end() ) {
         PyObject* pylen = PyObject_CallObject( iscbp->second, NULL );
         Py_ssize_t nlen2 = PyInt_AsSsize_t( pylen );
         Py_DECREF( pylen );

         if ( nlen2 == (Py_ssize_t)-1 && PyErr_Occurred() )
            PyErr_Clear();
         else
            return nlen2;
      }

      return nlen;            // return nlen after all, since have nothing better
   }

////////////////////////////////////////////////////////////////////////////////
/// Retrieve the buffer as a linear char array.

   const char* buffer_get( PyObject* self, int idx )
   {
      if ( idx < 0 || idx >= buffer_length( self ) ) {
         PyErr_SetString( PyExc_IndexError, "buffer index out of range" );
         return 0;
      }

#if PY_VERSION_HEX < 0x02050000
      const char* buf = 0;
#else
      char* buf = 0;     // interface change in 2.5, no other way to handle it
#endif
#if PY_VERSION_HEX < 0x03000000
      (*(PyBuffer_Type.tp_as_buffer->bf_getcharbuffer))( self, 0, &buf );
#else
      Py_buffer bufinfo;
      (*(PyBuffer_Type.tp_as_buffer->bf_getbuffer))( self, &bufinfo, PyBUF_SIMPLE );
      buf = (char*)bufinfo.buf;
      (*(PyBuffer_Type.tp_as_buffer->bf_releasebuffer))(self, &bufinfo);
      Py_DECREF(bufinfo.obj);
#endif

      if ( ! buf )
         PyErr_SetString( PyExc_IndexError, "attempt to index a null-buffer" );

      return buf;
   }

////////////////////////////////////////////////////////////////////////////////

#define PYROOT_IMPLEMENT_PYBUFFER_METHODS( name, type, stype, F1, F2 )       \
   PyObject* name##_buffer_str( PyObject* self )                             \
   {                                                                         \
      Py_ssize_t l = buffer_length( self );                                  \
      return PyROOT_PyUnicode_FromFormat( "<"#type" buffer, size " PY_SSIZE_T_FORMAT ">", l );\
   }                                                                         \
                                                                             \
   PyObject* name##_buffer_item( PyObject* self, Py_ssize_t idx ) {          \
      const char* buf = buffer_get( self, idx );                             \
      if ( buf )                                                             \
         return F1( (stype)*((type*)buf + idx) );                            \
      return 0;                                                              \
   }                                                                         \
                                                                             \
   int name##_buffer_ass_item( PyObject* self, Py_ssize_t idx, PyObject* val ) {\
      const char* buf = buffer_get( self, idx );                             \
      if ( ! buf )                                                           \
         return -1;                                                          \
                                                                             \
      type value = F2( val );                                                \
      if ( value == (type)-1 && PyErr_Occurred() )                           \
         return -1;                                                          \
                                                                             \
      *((type*)buf+idx) = (type)value;                                       \
       return 0;                                                             \
   }                                                                         \
                                                                             \
   PyObject* name##_buffer_subscript( PyObject* self, PyObject* item ) {     \
      if ( PyIndex_Check( item ) ) {                                         \
          Py_ssize_t idx = PyNumber_AsSsize_t( item, PyExc_IndexError );     \
          if ( idx == -1 && PyErr_Occurred() )                               \
               return 0;                                                     \
          return name##_buffer_item( self, idx );                            \
      }                                                                      \
      return 0;                                                              \
   }

   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Bool,   Bool_t,   Long_t,   PyBool_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Char,   Char_t,   Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( UChar,  UChar_t,  Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Short,  Short_t,  Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( UShort, UShort_t, Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Int,    Int_t,    Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( UInt,   UInt_t,   Long_t,   PyInt_FromLong, PyInt_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Long,   Long_t,   Long_t,   PyLong_FromLong, PyLong_AsLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( ULong,  ULong_t,  ULong_t,  PyLong_FromUnsignedLong, PyLong_AsUnsignedLong )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Float,  Float_t,  Double_t, PyFloat_FromDouble, PyFloat_AsDouble )
   PYROOT_IMPLEMENT_PYBUFFER_METHODS( Double, Double_t, Double_t, PyFloat_FromDouble, PyFloat_AsDouble )

   int pyroot_buffer_ass_subscript( PyObject* self, PyObject* idx, PyObject* val ) {
   // Assign the given value 'val' to the item at index 'idx.'
      if ( PyIndex_Check( idx ) ) {
         Py_ssize_t i = PyNumber_AsSsize_t( idx, PyExc_IndexError );
         if ( i == -1 && PyErr_Occurred() )
            return -1;
         return Py_TYPE(self)->tp_as_sequence->sq_ass_item( self, i, val );
      } else {
         PyErr_SetString( PyExc_TypeError, "buffer indices must be integers" );
         return -1;
      }
   }


////////////////////////////////////////////////////////////////////////////////
/// Allow the user to fix up the actual (type-strided) size of the buffer.

   PyObject* buffer_setsize( PyObject* self, PyObject* pynlen )
   {
      Py_ssize_t nlen = PyInt_AsSsize_t( pynlen );
      if ( nlen == -1 && PyErr_Occurred() )
         return 0;

#if PY_VERSION_HEX < 0x03000000
      ((PyBufferTop_t*)self)->fSize = nlen * ((PyBufferTop_t*)self)->fItemSize;
#else
      PyMemoryView_GET_BUFFER(self)->len = nlen * PyMemoryView_GET_BUFFER(self)->itemsize;
#endif 

      Py_INCREF( Py_None );
      return Py_None;
   }

   PyObject *SetSize(PyObject *self, PyObject *pynlen)
   {
      if (PyErr_WarnEx(PyExc_FutureWarning,
                       "buffer.SetSize(N) is deprecated and will disappear in a future version of ROOT. "
                       "Instead, use buffer.reshape((N,))", 1) < 0) {
         return nullptr;
      }

      return buffer_setsize(self, pynlen);
   }

   PyObject *reshape(PyObject *self, PyObject *shape)
   {
      // Allow the user to fix up the actual (type-strided) size of the buffer.
      if (!PyTuple_Check(shape) || PyTuple_GET_SIZE(shape) != 1) {
         if (shape) {
            PyObject* pystr = PyObject_Str(shape);
            if (pystr) {
               PyErr_Format(PyExc_TypeError, "tuple object of length 1 expected, received %s",
                   PyROOT_PyUnicode_AsStringChecked(pystr));
               Py_DECREF(pystr);
               return nullptr;
            }
         }
         PyErr_SetString(PyExc_TypeError, "tuple object of length 1 expected");
         return nullptr;
      }

      return buffer_setsize(self, PyTuple_GET_ITEM(shape, 0));
   }

////////////////////////////////////////////////////////////////////////////////
/// return a typecode in the style of module array

   PyObject* buf_typecode( PyObject* pyobject, void* )
   {
      if ( PyObject_TypeCheck( pyobject, &PyBoolBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"b" );
      else if ( PyObject_TypeCheck( pyobject, &PyCharBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"b" );
      else if ( PyObject_TypeCheck( pyobject, &PyUCharBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"B" );
      else if ( PyObject_TypeCheck( pyobject, &PyShortBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"h" );
      else if ( PyObject_TypeCheck( pyobject, &PyUShortBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"H" );
      else if ( PyObject_TypeCheck( pyobject, &PyIntBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"i" );
      else if ( PyObject_TypeCheck( pyobject, &PyUIntBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"I" );
      else if ( PyObject_TypeCheck( pyobject, &PyLongBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"l" );
      else if ( PyObject_TypeCheck( pyobject, &PyULongBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"L" );
      else if ( PyObject_TypeCheck( pyobject, &PyFloatBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"f" );
      else if ( PyObject_TypeCheck( pyobject, &PyDoubleBuffer_Type ) )
         return PyROOT_PyUnicode_FromString( (char*)"d" );

      PyErr_SetString( PyExc_TypeError, "received unknown buffer object" );
      return 0;
   }

////////////////////////////////////////////////////////////////////////////////

   PyGetSetDef buffer_getset[] = {
      { (char*)"typecode", (getter)buf_typecode, NULL, NULL, NULL },
      { (char*)NULL, NULL, NULL, NULL, NULL }
   };

////////////////////////////////////////////////////////////////////////////////

   PyMethodDef buffer_methods[] = {
      { (char*)"SetSize", (PyCFunction)SetSize, METH_O, NULL },
      { (char*)"reshape", (PyCFunction)reshape, METH_O, NULL },
      { (char*)NULL, NULL, 0, NULL }
   };

} // unnamed namespace


//- instance handler ------------------------------------------------------------
PyROOT::TPyBufferFactory* PyROOT::TPyBufferFactory::Instance()
{
// singleton factory
   static TPyBufferFactory* fac = new TPyBufferFactory;
   return fac;
}


//- constructor/destructor ------------------------------------------------------
#define PYROOT_INSTALL_PYBUFFER_METHODS( name, type )                           \
   Py##name##Buffer_Type.tp_name            = (char*)"ROOT.Py"#name"Buffer";    \
   Py##name##Buffer_Type.tp_base            = &PyBuffer_Type;                   \
   Py##name##Buffer_Type.tp_as_buffer       = PyBuffer_Type.tp_as_buffer;       \
   Py##name##Buffer_SeqMethods.sq_item      = (ssizeargfunc)name##_buffer_item; \
   Py##name##Buffer_SeqMethods.sq_ass_item  = (ssizeobjargproc)name##_buffer_ass_item;\
   Py##name##Buffer_SeqMethods.sq_length    = (lenfunc)buffer_length;           \
   Py##name##Buffer_Type.tp_as_sequence     = &Py##name##Buffer_SeqMethods;     \
   if ( PyBuffer_Type.tp_as_mapping ) { /* p2.6 and later */                    \
      Py##name##Buffer_MapMethods.mp_length    = (lenfunc)buffer_length;        \
      Py##name##Buffer_MapMethods.mp_subscript = (binaryfunc)name##_buffer_subscript;\
      Py##name##Buffer_MapMethods.mp_ass_subscript = (objobjargproc)pyroot_buffer_ass_subscript;\
      Py##name##Buffer_Type.tp_as_mapping      = &Py##name##Buffer_MapMethods;  \
   }                                                                            \
   Py##name##Buffer_Type.tp_str             = (reprfunc)name##_buffer_str;      \
   Py##name##Buffer_Type.tp_methods         = buffer_methods;                   \
   Py##name##Buffer_Type.tp_getset          = buffer_getset;                    \
   PyType_Ready( &Py##name##Buffer_Type );

PyROOT::TPyBufferFactory::TPyBufferFactory()
{
// construct python buffer types
   PYROOT_INSTALL_PYBUFFER_METHODS( Bool,   Bool_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Char,   Char_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( UChar,  UChar_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Short,  Short_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( UShort, UShort_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Int,    Int_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( UInt,   UInt_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Long,   Long_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( ULong,  ULong_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Float,  Float_t )
   PYROOT_INSTALL_PYBUFFER_METHODS( Double, Double_t )
}

////////////////////////////////////////////////////////////////////////////////

PyROOT::TPyBufferFactory::~TPyBufferFactory()
{
}

const char* getBoolFormat()   { return "b";}
const char* getCharFormat()   { return "b";}
const char* getUCharFormat()  { return "B";}
const char* getShortFormat()  { return "h";}
const char* getUShortFormat() { return "H";}
const char* getIntFormat()    { return "i";}
const char* getUIntFormat()   { return "I";}
const char* getLongFormat()   { return "l";}
const char* getULongFormat()  { return "L";}
const char* getFloatFormat()  { return "f";}
const char* getDoubleFormat() { return "d";}

#if PY_VERSION_HEX < 0x03000000
   #define PYBUFFER_SETITEMSIZE(buf,type) ((PyBufferTop_t*)buf)->fItemSize = Py_ssize_t(sizeof(type))
   #define PYBUFFER_SETFORMAT(buf,name) 
#else
   #define PYBUFFER_SETITEMSIZE(buf,type) PyMemoryView_GET_BUFFER(buf)->itemsize = Py_ssize_t(sizeof(type))
   #define PYBUFFER_SETFORMAT(buf,name) PyMemoryView_GET_BUFFER(buf)->format = (char *)get##name##Format() 
#endif

//- public members --------------------------------------------------------------
#define PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( name, type )                     \
PyObject* PyROOT::TPyBufferFactory::PyBuffer_FromMemory( type* address, Py_ssize_t size )\
{                                                                               \
   size = size < 0 ? INT_MAX : size;                                            \
   PyObject* buf = PyBuffer_FromReadWriteMemory( (void*)address, size );        \
   if ( buf ) {                                                                 \
      Py_INCREF( (PyObject*)(void*)&Py##name##Buffer_Type );                    \
      buf->ob_type = &Py##name##Buffer_Type;                                    \
      PYBUFFER_SETITEMSIZE(buf,type);                                           \
      PYBUFFER_SETFORMAT(buf, name);                                            \
   }                                                                            \
   return buf;                                                                  \
}                                                                               \
                                                                                \
PyObject* PyROOT::TPyBufferFactory::PyBuffer_FromMemory( type* address, PyObject* scb )\
{                                                                               \
   PyObject* buf = PyBuffer_FromMemory( address, Py_ssize_t(0) );               \
   if ( buf != 0 && PyCallable_Check( scb ) ) {                                 \
      Py_INCREF( scb );                                                         \
      gSizeCallbacks[ buf ] = scb;                                              \
   }                                                                            \
   return buf;                                                                  \
}

PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Bool,   Bool_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Char,   Char_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( UChar,  UChar_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Short,  Short_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( UShort, UShort_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Int,    Int_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( UInt,   UInt_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Long,   Long_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( ULong,  ULong_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Float,  Float_t )
PYROOT_IMPLEMENT_PYBUFFER_FROM_MEMORY( Double, Double_t )
