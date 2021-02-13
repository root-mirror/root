// @(#)root/pyroot:$Id$
// Author: Wim Lavrijsen, Jan 2005

// Bindings
#include "PyROOT.h"
#include "PyStrings.h"
#include "Converters.h"
#include "TCallContext.h"
#include "ObjectProxy.h"
#include "TPyBufferFactory.h"
#include "TCustomPyTypes.h"
#include "TTupleOfInstances.h"
#include "Utility.h"
#include "RootWrapper.h"

// ROOT
#include "TClass.h"          // for checking class info existence
#include "TClassEdit.h"      // for CleanType and ShortType

// Standard
#include <limits.h>
#include <stddef.h>      // for ptrdiff_t
#include <string.h>
#include <utility>
#include <sstream>
#include <tuple>

// FIXME: Should refer to PyROOT::TParameter in the code.
#ifdef R__CXXMODULES
  #define TParameter PyROOT::TParameter
#endif

//- data ______________________________________________________________________
namespace PyROOT {

// factories
   typedef TConverter* (*ConverterFactory_t) ( Long_t size );
   typedef std::map< std::string, ConverterFactory_t > ConvFactories_t;
   ConvFactories_t gConvFactories;
   R__EXTERN PyObject* gNullPtrObject;

}

//- pretend-ctypes helpers ----------------------------------------------------
#if PY_VERSION_HEX >= 0x02050000

struct PyROOT_tagCDataObject { // non-public (but so far very stable)
    PyObject_HEAD
    char* b_ptr;
};

static inline PyTypeObject* GetCTypesType( const char* name ) {
   PyObject* ct = PyImport_ImportModule( "ctypes" );
   if ( ! ct ) return nullptr;
   PyTypeObject* ct_t = (PyTypeObject*)PyObject_GetAttrString( ct, name );
   Py_DECREF( ct );
   return ct_t;
}

#endif

//- custom helpers to check ranges --------------------------------------------

////////////////////////////////////////////////////////////////////////////////
/// range-checking python integer to C++ bool conversion
static inline Bool_t PyROOT_PyLong_AsBool( PyObject* pyobject )
{
   Long_t l = PyLong_AsLong( pyobject );
// fail to pass float -> bool; the problem is rounding (0.1 -> 0 -> False)
   if ( ! ( l == 0 || l == 1 ) || PyFloat_Check( pyobject ) ) {
      PyErr_SetString( PyExc_ValueError, "boolean value should be bool, or integer 1 or 0" );
      return (Bool_t)-1;
   }
   return (Bool_t)l;
}

////////////////////////////////////////////////////////////////////////////////
/// python string to C++ char conversion
static inline Char_t PyROOT_PyUnicode_AsChar( PyObject* pyobject ) {
   return (Char_t)PyROOT_PyUnicode_AsString( pyobject )[0];
}

////////////////////////////////////////////////////////////////////////////////
/// range-checking python integer to C++ unsigned short int conversion
static inline UShort_t PyROOT_PyLong_AsUShort( PyObject* pyobject )
{
// prevent p2.7 silent conversions and do a range check
   if ( ! (PyLong_Check( pyobject ) || PyInt_Check( pyobject )) ) {
      PyErr_SetString( PyExc_TypeError, "unsigned short conversion expects an integer object" );
      return (UShort_t)-1;
   }
   Long_t l = PyLong_AsLong( pyobject );
   if ( l < 0 || USHRT_MAX < l ) {
      PyErr_Format( PyExc_ValueError, "integer %ld out of range for unsigned short", l );
      return (UShort_t)-1;

   }
   return (UShort_t)l;
}

////////////////////////////////////////////////////////////////////////////////
/// range-checking python integer to C++ short int conversion
static inline Short_t PyROOT_PyLong_AsShort( PyObject* pyobject )
{
// prevent p2.7 silent conversions and do a range check
   if ( ! (PyLong_Check( pyobject ) || PyInt_Check( pyobject )) ) {
      PyErr_SetString( PyExc_TypeError, "short int conversion expects an integer object" );
      return (Short_t)-1;
   }
   Long_t l = PyLong_AsLong( pyobject );
   if ( l < SHRT_MIN || SHRT_MAX < l ) {
      PyErr_Format( PyExc_ValueError, "integer %ld out of range for short int", l );
      return (Short_t)-1;

   }
   return (Short_t)l;
}


////////////////////////////////////////////////////////////////////////////////
/// strict python integer to C++ integer conversion
static inline Long_t PyROOT_PyLong_AsStrictLong( PyObject* pyobject )
{
// p2.7 and later silently converts floats to long, therefore require this
// check; earlier pythons may raise a SystemError which should be avoided as
// it is confusing
   if ( ! (PyLong_Check( pyobject ) || PyInt_Check( pyobject )) ) {
      PyErr_SetString( PyExc_TypeError, "int/long conversion expects an integer object" );
      return (Long_t)-1;
   }
   return (Long_t)PyLong_AsLong( pyobject );
}


//- base converter implementation ---------------------------------------------
PyObject* PyROOT::TConverter::FromMemory( void* )
{
// could happen if no derived class override
   PyErr_SetString( PyExc_TypeError, "C++ type can not be converted from memory" );
   return 0;
}

////////////////////////////////////////////////////////////////////////////////
/// could happen if no derived class override

Bool_t PyROOT::TConverter::ToMemory( PyObject*, void* )
{
   PyErr_SetString( PyExc_TypeError, "C++ type can not be converted to memory" );
   return kFALSE;
}


//- helper macro's ------------------------------------------------------------
#define PYROOT_IMPLEMENT_BASIC_CONVERTER( name, type, stype, F1, F2, tc )     \
Bool_t PyROOT::T##name##Converter::SetArg(                                    \
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )        \
{                                                                             \
/* convert `pyobject` to C++ 'type', set arg for call */                      \
   type val = (type)F2( pyobject );                                           \
   if ( val == (type)-1 && PyErr_Occurred() )                                 \
      return kFALSE;                                                          \
   para.fValue.f##name = val;                                                 \
   para.fTypeCode = tc;                                                       \
   return kTRUE;                                                              \
}                                                                             \
                                                                              \
PyObject* PyROOT::T##name##Converter::FromMemory( void* address )             \
{                                                                             \
   return F1( (stype)*((type*)address) );                                     \
}                                                                             \
                                                                              \
Bool_t PyROOT::T##name##Converter::ToMemory( PyObject* value, void* address ) \
{                                                                             \
   type s = (type)F2( value );                                                \
   if ( s == (type)-1 && PyErr_Occurred() )                                   \
      return kFALSE;                                                          \
   *((type*)address) = (type)s;                                               \
   return kTRUE;                                                              \
}


static inline Int_t ExtractChar( PyObject* pyobject, const char* tname, Int_t low, Int_t high )
{
   Int_t lchar = -1;
   if ( PyROOT_PyUnicode_Check( pyobject ) ) {
      if ( PyROOT_PyUnicode_GET_SIZE( pyobject ) == 1 )
         lchar = (Int_t)PyROOT_PyUnicode_AsChar( pyobject );
      else
         PyErr_Format( PyExc_TypeError, "%s expected, got string of size " PY_SSIZE_T_FORMAT,
             tname, PyROOT_PyUnicode_GET_SIZE( pyobject ) );
   } else if ( ! PyFloat_Check( pyobject ) ) {    // don't allow truncating conversion
      lchar = PyLong_AsLong( pyobject );
      if ( lchar == -1 && PyErr_Occurred() )
         ; // empty, as error already set
      else if ( ! ( low <= lchar && lchar <= high ) ) {
         PyErr_Format( PyExc_ValueError,
            "integer to character: value %d not in range [%d,%d]", lchar, low, high );
         lchar = -1;
      }
   } else
      PyErr_SetString( PyExc_TypeError, "char or small int type expected" );

   return lchar;
}


#define PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( name, type, F1           )\
Bool_t PyROOT::TConst##name##RefConverter::SetArg(                            \
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )        \
{                                                                             \
   type val = (type)F1( pyobject );                                           \
   if ( val == (type)-1 && PyErr_Occurred() )                                 \
      return kFALSE;                                                          \
   para.fValue.f##name = val;                                                 \
   para.fRef = &para.fValue.f##name;                                          \
   para.fTypeCode = 'r';                                                      \
   return kTRUE;                                                              \
}

#define PYROOT_IMPLEMENT_BASIC_CONST_CHAR_REF_CONVERTER( name, type, low, high )\
Bool_t PyROOT::TConst##name##RefConverter::SetArg(                            \
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )        \
{                                                                             \
/* convert `pyobject` to C++ <<type>>, set arg for call, allow int -> char */ \
   type val = (type)ExtractChar( pyobject, #type, low, high );                \
   if ( val == (type)-1 && PyErr_Occurred() )                                 \
      return kFALSE;                                                          \
   para.fValue.fLong = val;                                                   \
   para.fTypeCode = 'l';                                                      \
   return kTRUE;                                                              \
}


////////////////////////////////////////////////////////////////////////////////

#define PYROOT_IMPLEMENT_BASIC_CHAR_CONVERTER( name, type, low, high )        \
Bool_t PyROOT::T##name##Converter::SetArg(                                    \
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )        \
{                                                                             \
/* convert `pyobject` to C++ <<type>>, set arg for call, allow int -> char */ \
   Long_t val = ExtractChar( pyobject, #type, low, high );                    \
   if ( val == -1 && PyErr_Occurred() )                                       \
      return kFALSE;                                                          \
   para.fValue.fLong = val;                                                   \
   para.fTypeCode = 'l';                                                      \
   return kTRUE;                                                              \
}                                                                             \
                                                                              \
PyObject* PyROOT::T##name##Converter::FromMemory( void* address )             \
{                                                                             \
   return PyROOT_PyUnicode_FromFormat( "%c", *((type*)address) );             \
}                                                                             \
                                                                              \
Bool_t PyROOT::T##name##Converter::ToMemory( PyObject* value, void* address ) \
{                                                                             \
   if ( PyROOT_PyUnicode_Check( value ) ) {                                   \
      const char* buf = PyROOT_PyUnicode_AsString( value );                   \
      if ( PyErr_Occurred() )                                                 \
         return kFALSE;                                                       \
      int len = PyROOT_PyUnicode_GET_SIZE( value );                           \
      if ( len != 1 ) {                                                       \
         PyErr_Format( PyExc_TypeError, #type" expected, got string of size %d", len );\
         return kFALSE;                                                       \
      }                                                                       \
      *((type*)address) = (type)buf[0];                                       \
   } else {                                                                   \
      Long_t l = PyLong_AsLong( value );                                      \
      if ( l == -1 && PyErr_Occurred() )                                      \
         return kFALSE;                                                       \
      if ( ! ( low <= l && l <= high ) ) {                                    \
         PyErr_Format( PyExc_ValueError, \
            "integer to character: value %ld not in range [%d,%d]", l, low, high );\
         return kFALSE;                                                       \
      }                                                                       \
      *((type*)address) = (type)l;                                            \
   }                                                                          \
   return kTRUE;                                                              \
}


//- converters for built-ins --------------------------------------------------
PYROOT_IMPLEMENT_BASIC_CONVERTER( Long, Long_t, Long_t, PyLong_FromLong, PyROOT_PyLong_AsStrictLong, 'l' )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ long&, set arg for call

Bool_t PyROOT::TLongRefConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
#if PY_VERSION_HEX < 0x03000000
   if ( TCustomInt_CheckExact( pyobject ) ) {
      para.fValue.fVoidp = (void*)&((PyIntObject*)pyobject)->ob_ival;
      para.fTypeCode = 'V';
      if (PyErr_WarnEx(PyExc_FutureWarning,
                       "ROOT.Long is deprecated and will disappear in a future version of ROOT. "
                       "Instead, use ctypes.c_long for pass-by-ref of longs", 1) < 0) {
         return kFALSE;
      }
      return kTRUE;
   }
#endif

#if PY_VERSION_HEX < 0x02050000
   PyErr_SetString( PyExc_TypeError, "use ROOT.Long for pass-by-ref of longs" );
   return kFALSE;
#endif

// TODO: this keeps a refcount to the type .. it should be okay to drop that
   static PyTypeObject* c_long_type = GetCTypesType( "c_long" );
   if ( Py_TYPE( pyobject ) == c_long_type ) {
      para.fValue.fVoidp = (void*)((PyROOT_tagCDataObject*)pyobject)->b_ptr;
      para.fTypeCode = 'V';
      return kTRUE;
   }

   PyErr_SetString( PyExc_TypeError, "use ctypes.c_long for pass-by-ref of longs" );
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

PYROOT_IMPLEMENT_BASIC_CONST_CHAR_REF_CONVERTER( Char,  Char_t,  CHAR_MIN,  CHAR_MAX )
PYROOT_IMPLEMENT_BASIC_CONST_CHAR_REF_CONVERTER( UChar, UChar_t,        0, UCHAR_MAX )

PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Bool,      Bool_t,    PyROOT_PyLong_AsBool )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Short,     Short_t,   PyROOT_PyLong_AsShort )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( UShort,    UShort_t,  PyROOT_PyLong_AsUShort )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Int,       Int_t,     PyROOT_PyLong_AsStrictLong )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( UInt,      UInt_t,    PyLongOrInt_AsULong )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Long,      Long_t,    PyROOT_PyLong_AsStrictLong )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( ULong,     ULong_t,   PyLongOrInt_AsULong )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( LongLong,  Long64_t,  PyLong_AsLongLong )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( ULongLong, ULong64_t, PyLongOrInt_AsULong64 )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ (pseudo)int&, set arg for call

Bool_t PyROOT::TIntRefConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
#if PY_VERSION_HEX < 0x03000000
   if ( TCustomInt_CheckExact( pyobject ) ) {
      para.fValue.fVoidp = (void*)&((PyIntObject*)pyobject)->ob_ival;
      para.fTypeCode = 'V';
      if (PyErr_WarnEx(PyExc_FutureWarning,
                       "ROOT.Long is deprecated and will disappear in a future version of ROOT. "
                       "Instead, use ctypes.c_int for pass-by-ref of ints", 1) < 0) {
         return kFALSE;
      }
      return kTRUE;
   }
#endif

#if PY_VERSION_HEX >= 0x02050000
// TODO: this keeps a refcount to the type .. it should be okay to drop that
   static PyTypeObject* c_int_type = GetCTypesType( "c_int" );
   if ( Py_TYPE( pyobject ) == c_int_type ) {
      para.fValue.fVoidp = (void*)((PyROOT_tagCDataObject*)pyobject)->b_ptr;
      para.fTypeCode = 'V';
      return kTRUE;
   }
#endif

// alternate, pass pointer from buffer
   int buflen = Utility::GetBuffer( pyobject, 'i', sizeof(int), para.fValue.fVoidp );
   if ( para.fValue.fVoidp && buflen ) {
      para.fTypeCode = 'V';
      return kTRUE;
   };

#if PY_VERSION_HEX < 0x02050000
   PyErr_SetString( PyExc_TypeError, "use ROOT.Long for pass-by-ref of ints" );
#else
   PyErr_SetString( PyExc_TypeError, "use ctypes.c_int for pass-by-ref of ints" );
#endif
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ bool, allow int/long -> bool, set arg for call

PYROOT_IMPLEMENT_BASIC_CONVERTER( Bool, Bool_t, Long_t, PyInt_FromLong, PyROOT_PyLong_AsBool, 'l' )

////////////////////////////////////////////////////////////////////////////////

PYROOT_IMPLEMENT_BASIC_CHAR_CONVERTER( Char,  Char_t,  CHAR_MIN, CHAR_MAX  )
PYROOT_IMPLEMENT_BASIC_CHAR_CONVERTER( UChar, UChar_t,        0, UCHAR_MAX )

////////////////////////////////////////////////////////////////////////////////

PYROOT_IMPLEMENT_BASIC_CONVERTER( Short,  Short_t,  Long_t, PyInt_FromLong,  PyROOT_PyLong_AsShort, 'l' )
PYROOT_IMPLEMENT_BASIC_CONVERTER( UShort, UShort_t, Long_t, PyInt_FromLong,  PyROOT_PyLong_AsUShort, 'l' )
PYROOT_IMPLEMENT_BASIC_CONVERTER( Int,    Int_t,    Long_t, PyInt_FromLong,  PyROOT_PyLong_AsStrictLong, 'l' )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ unsigned long, set arg for call

Bool_t PyROOT::TULongConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   para.fValue.fULong = PyLongOrInt_AsULong( pyobject );
   if ( PyErr_Occurred() )
      return kFALSE;
   para.fTypeCode = 'U';
   return kTRUE;
}

PyObject* PyROOT::TULongConverter::FromMemory( void* address )
{
// construct python object from C++ unsigned long read at <address>
   return PyLong_FromUnsignedLong( *((ULong_t*)address) );
}

Bool_t PyROOT::TULongConverter::ToMemory( PyObject* value, void* address )
{
// convert <value> to C++ unsigned long, write it at <address>
   ULong_t u = PyLongOrInt_AsULong( value );
   if ( PyErr_Occurred() )
      return kFALSE;
   *((ULong_t*)address) = u;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// construct python object from C++ unsigned int read at <address>

PyObject* PyROOT::TUIntConverter::FromMemory( void* address )
{
   return PyLong_FromUnsignedLong( *((UInt_t*)address) );
}

Bool_t PyROOT::TUIntConverter::ToMemory( PyObject* value, void* address )
{
// convert <value> to C++ unsigned int, write it at <address>
   ULong_t u = PyLongOrInt_AsULong( value );
   if ( PyErr_Occurred() )
      return kFALSE;

   if ( u > (ULong_t)UINT_MAX ) {
      PyErr_SetString( PyExc_OverflowError, "value too large for unsigned int" );
      return kFALSE;
   }

   *((UInt_t*)address) = (UInt_t)u;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// floating point converters

PYROOT_IMPLEMENT_BASIC_CONVERTER( Float,  Float_t,  Double_t, PyFloat_FromDouble, PyFloat_AsDouble, 'f' )
PYROOT_IMPLEMENT_BASIC_CONVERTER( Double, Double_t, Double_t, PyFloat_FromDouble, PyFloat_AsDouble, 'd' )

PYROOT_IMPLEMENT_BASIC_CONVERTER( LongDouble, LongDouble_t, LongDouble_t, PyFloat_FromDouble, PyFloat_AsDouble, 'D' )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ double&, set arg for call

Bool_t PyROOT::TDoubleRefConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( TCustomFloat_CheckExact( pyobject ) ) {
      para.fValue.fVoidp = (void*)&((PyFloatObject*)pyobject)->ob_fval;
      para.fTypeCode = 'V';
      if (PyErr_WarnEx(PyExc_FutureWarning,
                       "ROOT.Double is deprecated and will disappear in a future version of ROOT. "
                       "Instead, use ctypes.c_double for pass-by-ref of doubles", 1) < 0) {
         return kFALSE;
      }
      return kTRUE;
   }

// alternate, pass pointer from buffer
   int buflen = Utility::GetBuffer( pyobject, 'd', sizeof(double), para.fValue.fVoidp );
   if ( para.fValue.fVoidp && buflen ) {
      para.fTypeCode = 'V';
      return kTRUE;
   }

   PyErr_SetString( PyExc_TypeError, "use ctypes.c_double for pass-by-ref of doubles" );
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////

PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Float,      Float_t,      PyFloat_AsDouble )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( Double,     Double_t,     PyFloat_AsDouble )
PYROOT_IMPLEMENT_BASIC_CONST_REF_CONVERTER( LongDouble, LongDouble_t, PyFloat_AsDouble )

////////////////////////////////////////////////////////////////////////////////
/// can't happen (unless a type is mapped wrongly), but implemented for completeness

Bool_t PyROOT::TVoidConverter::SetArg( PyObject*, TParameter&, TCallContext* )
{
   PyErr_SetString( PyExc_SystemError, "void/unknown arguments can\'t be set" );
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ long long, set arg for call

Bool_t PyROOT::TLongLongConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( PyFloat_Check( pyobject ) ) {
   // special case: float implements nb_int, but allowing rounding conversions
   // interferes with overloading
      PyErr_SetString( PyExc_ValueError, "can not convert float to long long" );
      return kFALSE;
   }

   para.fValue.fLongLong = PyLong_AsLongLong( pyobject );
   if ( PyErr_Occurred() )
      return kFALSE;
   para.fTypeCode = 'k';
   return kTRUE;
}

PyObject* PyROOT::TLongLongConverter::FromMemory( void* address )
{
// construct python object from C++ long long read at <address>
   return PyLong_FromLongLong( *(Long64_t*)address );
}

Bool_t PyROOT::TLongLongConverter::ToMemory( PyObject* value, void* address )
{
// convert <value> to C++ long long, write it at <address>
   Long64_t ll = PyLong_AsLongLong( value );
   if ( ll == -1 && PyErr_Occurred() )
      return kFALSE;
   *((Long64_t*)address) = ll;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ unsigned long long, set arg for call

Bool_t PyROOT::TULongLongConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   para.fValue.fULongLong = PyLongOrInt_AsULong64( pyobject );
   if ( PyErr_Occurred() )
      return kFALSE;
   para.fTypeCode = 'K';
   return kTRUE;
}

PyObject* PyROOT::TULongLongConverter::FromMemory( void* address )
{
// construct python object from C++ unsigned long long read at <address>
   return PyLong_FromUnsignedLongLong( *(ULong64_t*)address );
}

Bool_t PyROOT::TULongLongConverter::ToMemory( PyObject* value, void* address )
{
// convert <value> to C++ unsigned long long, write it at <address>
   Long64_t ull = PyLongOrInt_AsULong64( value );
   if ( PyErr_Occurred() )
      return kFALSE;
   *((ULong64_t*)address) = ull;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////
/// construct a new string and copy it in new memory

static std::tuple<const char*,Py_ssize_t> getStringAndSizeCString(PyObject* pyobject) {
#if PY_VERSION_HEX >= 0x03030000
   // Support non-ASCII strings (get the right length in bytes)
   Py_ssize_t size = 0;
   auto charArr = PyROOT_PyUnicode_AsStringAndSize(pyobject, &size);
#else
   auto size = PyROOT_PyUnicode_GET_SIZE(pyobject);
   auto charArr = PyROOT_PyUnicode_AsStringChecked(pyobject);
#endif
   return std::tuple<const char*,Py_ssize_t>(charArr, size);
}

Bool_t PyROOT::TCStringConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if (PyROOT_PyUnicode_Check(pyobject)
#if PY_VERSION_HEX < 0x03000000
       || PyUnicode_Check(pyobject)
#endif
   ) {
      auto strAndSize = getStringAndSizeCString(pyobject);
      fBuffer = std::string(std::get<0>(strAndSize), std::get<1>(strAndSize));
   } else if (PyBytes_Check(pyobject)) {
      auto s = PyBytes_AsString(pyobject);                                    \
      auto size = PyBytes_GET_SIZE(pyobject);
      fBuffer = std::string(s, size);
   } else {
      return kFALSE;
   }

// verify (too long string will cause truncation, no crash)
   if ( fMaxSize < (UInt_t)fBuffer.size() )
      PyErr_Warn( PyExc_RuntimeWarning, (char*)"string too long for char array (truncated)" );
   else if ( fMaxSize != UINT_MAX )
      fBuffer.resize( fMaxSize, '\0' );      // padd remainder of buffer as needed

// set the value and declare success
   para.fValue.fVoidp = (void*)fBuffer.c_str();
   para.fTypeCode = 'p';
   return kTRUE;
}

PyObject* PyROOT::TCStringConverter::FromMemory( void* address )
{
// construct python object from C++ const char* read at <address>
   if ( address && *(char**)address ) {
      if ( fMaxSize != UINT_MAX ) {          // need to prevent reading beyond boundary
         std::string buf( *(char**)address, fMaxSize );     // cut on fMaxSize
         return PyROOT_PyUnicode_FromString( buf.c_str() ); // cut on \0
      }

      return PyROOT_PyUnicode_FromString( *(char**)address );
   }

// empty string in case there's no address
   Py_INCREF( PyStrings::gEmptyString );
   return PyStrings::gEmptyString;
}

Bool_t PyROOT::TCStringConverter::ToMemory( PyObject* value, void* address )
{
// convert <value> to C++ const char*, write it at <address>
   const char* s = PyROOT_PyUnicode_AsStringChecked( value );
   if ( PyErr_Occurred() )
      return kFALSE;

// verify (too long string will cause truncation, no crash)
   if ( fMaxSize < (UInt_t)PyROOT_PyUnicode_GET_SIZE( value ) )
      PyErr_Warn( PyExc_RuntimeWarning, (char*)"string too long for char array (truncated)" );

   if ( fMaxSize != UINT_MAX )
      strncpy( *(char**)address, s, fMaxSize );   // padds remainder
   else
      // coverity[secure_coding] - can't help it, it's intentional.
      strcpy( *(char**)address, s );

   return kTRUE;
}


//- pointer/array conversions -------------------------------------------------
namespace {

   using namespace PyROOT;

   inline Bool_t CArraySetArg(
      PyObject* pyobject, TParameter& para, char tc, int size )
   {
   // general case of loading a C array pointer (void* + type code) as function argument
      if ( pyobject == gNullPtrObject ) {
         para.fValue.fVoidp = NULL;
      } else {
         int buflen = Utility::GetBuffer( pyobject, tc, size, para.fValue.fVoidp );
         if ( ! para.fValue.fVoidp || buflen == 0 )
            return kFALSE;
      }
      para.fTypeCode = 'p';
      return kTRUE;
   }

} // unnamed namespace


////////////////////////////////////////////////////////////////////////////////
/// attempt base class first (i.e. passing a string), but if that fails, try a buffer

Bool_t PyROOT::TNonConstCStringConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   if ( this->TCStringConverter::SetArg( pyobject, para, ctxt ) )
      return kTRUE;

// apparently failed, try char buffer
   PyErr_Clear();
   return CArraySetArg( pyobject, para, 'c', sizeof(char) );
}

////////////////////////////////////////////////////////////////////////////////
/// assume this is a buffer access if the size is known; otherwise assume string

PyObject* PyROOT::TNonConstCStringConverter::FromMemory( void* address )
{
   if ( fMaxSize != UINT_MAX )
      return PyROOT_PyUnicode_FromStringAndSize( *(char**)address, fMaxSize );
   return this->TCStringConverter::FromMemory( address );
}

////////////////////////////////////////////////////////////////////////////////
/// attempt base class first (i.e. passing a string), but if that fails, try a buffer

Bool_t PyROOT::TNonConstUCStringConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   if ( this->TCStringConverter::SetArg( pyobject, para, ctxt ) )
      return kTRUE;

// apparently failed, try char buffer
   PyErr_Clear();
   return CArraySetArg( pyobject, para, 'B', sizeof(unsigned char) );
}

////////////////////////////////////////////////////////////////////////////////
/// (1): "null pointer" or C++11 style nullptr

Bool_t PyROOT::TVoidArrayConverter::GetAddressSpecialCase( PyObject* pyobject, void*& address )
{
   if ( pyobject == Py_None || pyobject == gNullPtrObject ) {

      if (pyobject == Py_None) {
         if (PyErr_WarnEx(PyExc_FutureWarning,
                          "The conversion from None to null pointer is deprecated "
                          "and will not be allowed anymore in a future version of ROOT. "
                          "Instead, use ROOT.nullptr or 0", 1) < 0) {
            return kFALSE;
         }
      }

      address = (void*)0;
      return kTRUE;
   }

// (2): allow integer zero to act as a null pointer, no deriveds
   if ( PyInt_CheckExact( pyobject ) || PyLong_CheckExact( pyobject ) ) {
      Long_t val = (Long_t)PyLong_AsLong( pyobject );
      if ( val == 0l ) {
         address = (void*)val;
         return kTRUE;
      }

      return kFALSE;
   }

// (3): opaque PyCapsule (CObject in older pythons) from somewhere
   if ( PyROOT_PyCapsule_CheckExact( pyobject ) ) {
      address = (void*)PyROOT_PyCapsule_GetPointer( pyobject, NULL );
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// just convert pointer if it is a ROOT object

Bool_t PyROOT::TVoidArrayConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   if ( ObjectProxy_Check( pyobject ) ) {
   // depending on memory policy, some objects are no longer owned when passed to C++
      if ( ! fKeepControl && ! UseStrictOwnership( ctxt ) )
         ((ObjectProxy*)pyobject)->Release();

   // set pointer (may be null) and declare success
      para.fValue.fVoidp = ((ObjectProxy*)pyobject)->GetObject();
      para.fTypeCode = 'p';
      return kTRUE;
   }

// handle special cases
   if ( GetAddressSpecialCase( pyobject, para.fValue.fVoidp ) ) {
      para.fTypeCode = 'p';
      return kTRUE;
   }

// final try: attempt to get buffer
   int buflen = Utility::GetBuffer( pyobject, '*', 1, para.fValue.fVoidp, kFALSE );

// ok if buffer exists (can't perform any useful size checks)
   if ( para.fValue.fVoidp && buflen != 0 ) {
      para.fTypeCode = 'p';
      return kTRUE;
   }

// give up
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// nothing sensible can be done, just return <address> as pylong

PyObject* PyROOT::TVoidArrayConverter::FromMemory( void* address )
{
   if ( ! address || *(ptrdiff_t*)address == 0 ) {
      Py_INCREF( gNullPtrObject );
      return gNullPtrObject;
   }
   return BufFac_t::Instance()->PyBuffer_FromMemory( (Long_t*)*(ptrdiff_t**)address, sizeof(void*) );
}

////////////////////////////////////////////////////////////////////////////////
/// just convert pointer if it is a ROOT object

Bool_t PyROOT::TVoidArrayConverter::ToMemory( PyObject* value, void* address )
{
   if ( ObjectProxy_Check( value ) ) {
   // depending on memory policy, some objects are no longer owned when passed to C++
      if ( ! fKeepControl && TCallContext::sMemoryPolicy != TCallContext::kUseStrict )
         ((ObjectProxy*)value)->Release();

   // set pointer (may be null) and declare success
      *(void**)address = ((ObjectProxy*)value)->GetObject();
      return kTRUE;
   }

// handle special cases
   void* ptr = 0;
   if ( GetAddressSpecialCase( value, ptr ) ) {
      *(void**)address = ptr;
      return kTRUE;
   }

// final try: attempt to get buffer
   void* buf = 0;
   int buflen = Utility::GetBuffer( value, '*', 1, buf, kFALSE );
   if ( ! buf || buflen == 0 )
      return kFALSE;

   *(void**)address = buf;
   return kTRUE;
}

////////////////////////////////////////////////////////////////////////////////

#define PYROOT_IMPLEMENT_ARRAY_CONVERTER( name, type, code )                 \
Bool_t PyROOT::T##name##ArrayConverter::SetArg(                              \
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )       \
{                                                                            \
   return CArraySetArg( pyobject, para, code, sizeof(type) );                \
}                                                                            \
                                                                             \
Bool_t PyROOT::T##name##ArrayRefConverter::SetArg(                           \
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )             \
{                                                                            \
   Bool_t result = T##name##ArrayConverter::SetArg( pyobject, para, ctxt );  \
   para.fTypeCode = 'V';                                                     \
   return result;                                                            \
}                                                                            \
                                                                             \
PyObject* PyROOT::T##name##ArrayConverter::FromMemory( void* address )       \
{                                                                            \
   return BufFac_t::Instance()->PyBuffer_FromMemory( *(type**)address, fSize * sizeof(type) );\
}                                                                            \
                                                                             \
Bool_t PyROOT::T##name##ArrayConverter::ToMemory( PyObject* value, void* address )\
{                                                                            \
   void* buf = 0;                                                            \
   int buflen = Utility::GetBuffer( value, code, sizeof(type), buf );        \
   if ( ! buf || buflen == 0 )                                               \
      return kFALSE;                                                         \
   if ( 0 <= fSize ) {                                                       \
      if ( fSize < buflen/(int)sizeof(type) ) {                              \
         PyErr_SetString( PyExc_ValueError, "buffer too large for value" );  \
         return kFALSE;                                                      \
      }                                                                      \
      memcpy( *(type**)address, buf, 0 < buflen ? ((size_t) buflen) : sizeof(type) );\
   } else                                                                    \
      *(type**)address = (type*)buf;                                         \
   return kTRUE;                                                             \
}

////////////////////////////////////////////////////////////////////////////////

PYROOT_IMPLEMENT_ARRAY_CONVERTER( Bool,   Bool_t,   'b' )   // signed char
PYROOT_IMPLEMENT_ARRAY_CONVERTER( Short,  Short_t,  'h' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( UShort, UShort_t, 'H' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( Int,    Int_t,    'i' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( UInt,   UInt_t,   'I' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( Long,   Long_t,   'l' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( ULong,  ULong_t,  'L' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( Float,  Float_t,  'f' )
PYROOT_IMPLEMENT_ARRAY_CONVERTER( Double, Double_t, 'd' )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ long long*, set arg for call

Bool_t PyROOT::TLongLongArrayConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   PyObject* pytc = PyObject_GetAttr( pyobject, PyStrings::gTypeCode );
   if ( pytc != 0 ) {              // iow, this array has a known type, but there's no
      Py_DECREF( pytc );           // such thing for long long in module array
      return kFALSE;
   }

   return TVoidArrayConverter::SetArg( pyobject, para, ctxt );
}


//- converters for special cases ----------------------------------------------

static std::tuple<const char*,Py_ssize_t> getStringAndSizeSTLString(PyObject* pyobject) {
#if PY_VERSION_HEX >= 0x03030000
   // Support non-ASCII strings (get the right length in bytes)
   Py_ssize_t size = 0;
   auto charArr = PyROOT_PyUnicode_AsStringAndSize(pyobject, &size);
#else
   auto size = PyROOT_PyUnicode_GET_SIZE(pyobject);
   auto charArr = PyROOT_PyUnicode_AsString(pyobject);
#endif
   return std::tuple<const char*,Py_ssize_t>(charArr, size);
}

#define PYROOT_IMPLEMENT_STRING_AS_PRIMITIVE_CONVERTER( name, type, F1, F2 )  \
PyROOT::T##name##Converter::T##name##Converter( Bool_t keepControl ) :        \
      TCppObjectConverter( Cppyy::GetScope( #type ), keepControl ) {}         \
                                                                              \
Bool_t PyROOT::T##name##Converter::SetArg(                                    \
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )              \
{                                                                             \
   if ( PyROOT_PyUnicode_Check( pyobject ) ) {                                \
      auto strAndSize = getStringAndSizeSTLString(pyobject);                  \
      fBuffer = type(std::get<0>(strAndSize), std::get<1>(strAndSize));       \
      para.fValue.fVoidp = &fBuffer;                                          \
      para.fTypeCode = 'V';                                                   \
      return kTRUE;                                                           \
   }                                                                          \
                                                                              \
   if (PyBytes_Check(pyobject)) {                                             \
      auto s = PyBytes_AsString(pyobject);                                    \
      auto size = PyBytes_GET_SIZE(pyobject);                                 \
      fBuffer = type(s, size);                                                \
      para.fValue.fVoidp = &fBuffer;                                          \
      para.fTypeCode = 'V';                                                   \
      return kTRUE;                                                           \
   }                                                                          \
                                                                              \
   if ( ! ( PyInt_Check( pyobject ) || PyLong_Check( pyobject ) ) ) {         \
      Bool_t result = TCppObjectConverter::SetArg( pyobject, para, ctxt );    \
      para.fTypeCode = 'V';                                                   \
      return result;                                                          \
   }                                                                          \
   return kFALSE;                                                             \
}                                                                             \
                                                                              \
PyObject* PyROOT::T##name##Converter::FromMemory( void* address )             \
{                                                                             \
   if ( address )                                                             \
      return PyROOT_PyUnicode_FromStringAndSize( ((type*)address)->F1(), ((type*)address)->F2() );\
   Py_INCREF( PyStrings::gEmptyString );                                      \
   return PyStrings::gEmptyString;                                            \
}                                                                             \
                                                                              \
Bool_t PyROOT::T##name##Converter::ToMemory( PyObject* value, void* address ) \
{                                                                             \
   if ( PyROOT_PyUnicode_Check( value ) ) {                                   \
      *((type*)address) = PyROOT_PyUnicode_AsString( value );                 \
      return kTRUE;                                                           \
   }                                                                          \
                                                                              \
   return TCppObjectConverter::ToMemory( value, address );                    \
}

PYROOT_IMPLEMENT_STRING_AS_PRIMITIVE_CONVERTER( TString,   TString,     Data, Length )
PYROOT_IMPLEMENT_STRING_AS_PRIMITIVE_CONVERTER( STLString, std::string, c_str, size )
PYROOT_IMPLEMENT_STRING_AS_PRIMITIVE_CONVERTER( STLStringView, std::string_view, data, size )

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ instance*, set arg for call

Bool_t PyROOT::TCppObjectConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   if ( ! ObjectProxy_Check( pyobject ) ) {
      if ( GetAddressSpecialCase( pyobject, para.fValue.fVoidp ) ) {
         para.fTypeCode = 'p';      // allow special cases such as NULL
         return kTRUE;
      }

   // not a PyROOT object (TODO: handle SWIG etc.)
      return kFALSE;
   }

   ObjectProxy* pyobj = (ObjectProxy*)pyobject;
   if ( pyobj->ObjectIsA() && Cppyy::IsSubtype( pyobj->ObjectIsA(), fClass ) ) {
   // depending on memory policy, some objects need releasing when passed into functions
      if ( ! KeepControl() && ! UseStrictOwnership( ctxt ) )
         ((ObjectProxy*)pyobject)->Release();

   // calculate offset between formal and actual arguments
      para.fValue.fVoidp = pyobj->GetObject();
      if ( pyobj->ObjectIsA() != fClass ) {
         para.fValue.fLong += Cppyy::GetBaseOffset(
            pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */ );
      }

   // set pointer (may be null) and declare success
      para.fTypeCode = 'p';
      return kTRUE;

   } else if ( ! TClass::GetClass( Cppyy::GetFinalName( fClass ).c_str() )->GetClassInfo() ) {
   // assume "user knows best" to allow anonymous pointer passing
      para.fValue.fVoidp = pyobj->GetObject();
      para.fTypeCode = 'p';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// construct python object from C++ instance read at <address>

PyObject* PyROOT::TCppObjectConverter::FromMemory( void* address )
{
   return BindCppObject( address, fClass, kFALSE );
}

////////////////////////////////////////////////////////////////////////////////
/// convert <value> to C++ instance, write it at <address>

Bool_t PyROOT::TCppObjectConverter::ToMemory( PyObject* value, void* address )
{
   if ( ! ObjectProxy_Check( value ) ) {
      void* ptr = 0;
      if ( GetAddressSpecialCase( value, ptr ) ) {
         *(void**)address = ptr;             // allow special cases such as NULL
         return kTRUE;
      }

   // not a PyROOT object (TODO: handle SWIG etc.)
      return kFALSE;
   }

   if ( Cppyy::IsSubtype( ((ObjectProxy*)value)->ObjectIsA(), fClass ) ) {
   // depending on memory policy, some objects need releasing when passed into functions
      if ( ! KeepControl() && TCallContext::sMemoryPolicy != TCallContext::kUseStrict )
         ((ObjectProxy*)value)->Release();

   // call assignment operator through a temporarily wrapped object proxy
      PyObject* pyobj = BindCppObjectNoCast( address, fClass );
      ((ObjectProxy*)pyobj)->Release();     // TODO: might be recycled (?)
      PyObject* result = PyObject_CallMethod( pyobj, (char*)"__assign__", (char*)"O", value );
      Py_DECREF( pyobj );
      if ( result ) {
         Py_DECREF( result );
         return kTRUE;
      }
   }

   return kFALSE;
}

// TODO: CONSOLIDATE ValueCpp, RefCpp, and CppObject ...

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ instance, set arg for call

Bool_t PyROOT::TValueCppObjectConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( ObjectProxy_Check( pyobject ) ) {
      ObjectProxy* pyobj = (ObjectProxy*)pyobject;
      if ( pyobj->ObjectIsA() && Cppyy::IsSubtype( pyobj->ObjectIsA(), fClass ) ) {
         // calculate offset between formal and actual arguments
         para.fValue.fVoidp = pyobj->GetObject();
         if ( ! para.fValue.fVoidp )
            return kFALSE;

         if ( pyobj->ObjectIsA() != fClass ) {
            para.fValue.fLong += Cppyy::GetBaseOffset(
            pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */ );
         }

         para.fTypeCode = 'V';
         return kTRUE;
      }
   }
   else if ( PyTuple_Check( pyobject ) ){  // It is a Python tuple (equivalent to a C++ initializer list)

      // instantiate an object proxy of this class
      if( ! fObjProxy ) {
         // retrieve the python class from which we will create an instance
         PyObject* pyclass = CreateScopeProxy( fClass );
         if ( ! pyclass ) return kFALSE;                  // error has been set in CreateScopeProxy
         fObjProxy = (ObjectProxy*)((PyTypeObject*)pyclass)->tp_new( (PyTypeObject*)pyclass, NULL, NULL );
         Py_DECREF( pyclass );
      }

      if( fObjProxy->GetObject() ) {
         // the actual C++ object was already created (in a previous call), so we need to destroy it
         Cppyy::CallDestructor( fObjProxy->ObjectIsA(), fObjProxy->GetObject() );
         Cppyy::Deallocate( fObjProxy->ObjectIsA(), fObjProxy->GetObject() );
         fObjProxy->Set(nullptr);
      }

      // get the constructor (i.e. __init__)
      PyObject* constructor = PyObject_GetAttr( (PyObject*)fObjProxy, PyStrings::gInit );
      if( ! constructor ) return kFALSE;

      // call the constructor with the arguments in the given tuple
      PyObject* obj = PyObject_CallObject( constructor, pyobject );
      Py_DECREF( constructor );
      if ( ! obj ) return kFALSE;
      Py_DECREF( obj );

      para.fValue.fVoidp = fObjProxy->GetObject();
      para.fTypeCode = 'V';
      return kTRUE;

   }
   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ instance&, set arg for call

Bool_t PyROOT::TRefCppObjectConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( ObjectProxy_Check( pyobject ) ) {   // It is PyROOT object
      ObjectProxy* pyobj = (ObjectProxy*)pyobject;
      if ( pyobj->ObjectIsA() && Cppyy::IsSubtype( pyobj->ObjectIsA(), fClass ) ) {
      // calculate offset between formal and actual arguments
            para.fValue.fVoidp = pyobj->GetObject();
            if ( pyobj->ObjectIsA() != fClass ) {
            para.fValue.fLong += Cppyy::GetBaseOffset(
                  pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */ );
            }

            para.fTypeCode = 'V';
            return kTRUE;
      } else if ( ! TClass::GetClass( Cppyy::GetFinalName( fClass ).c_str() )->GetClassInfo() ) {
      // assume "user knows best" to allow anonymous reference passing
            para.fValue.fVoidp = pyobj->GetObject();
            para.fTypeCode = 'V';
            return kTRUE;
      }
   }
   else if ( PyTuple_Check( pyobject ) ){  // It is a Python tuple (equivalent to a C++ initializer list)

      // instantiate an object proxy of this class
      if( ! fObjProxy ) {
         // retrieve the python class from which we will create an instance
         PyObject* pyclass = CreateScopeProxy( fClass );
         if ( ! pyclass ) return kFALSE;                  // error has been set in CreateScopeProxy
         fObjProxy = (ObjectProxy*)((PyTypeObject*)pyclass)->tp_new( (PyTypeObject*)pyclass, NULL, NULL );
         Py_DECREF( pyclass );
      }

      if( fObjProxy->GetObject() ) {
         // the actual C++ object was already created (in a previous call), so we need to destroy it
         Cppyy::CallDestructor( fObjProxy->ObjectIsA(), fObjProxy->GetObject() );
         Cppyy::Deallocate( fObjProxy->ObjectIsA(), fObjProxy->GetObject() );
         fObjProxy->Set(nullptr);
      }

      // get the constructor (i.e. __init__)
      PyObject* constructor = PyObject_GetAttr( (PyObject*)fObjProxy, PyStrings::gInit );
      if( ! constructor ) return kFALSE;

      // call the constructor with the arguments in the given tuple
      PyObject* obj = PyObject_CallObject( constructor, pyobject );
      Py_DECREF( constructor );
      if ( ! obj ) return kFALSE;
      Py_DECREF( obj );

      para.fValue.fVoidp = fObjProxy->GetObject();
      para.fTypeCode = 'V';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ instance**, set arg for call

template <bool ISREFERENCE>
Bool_t PyROOT::TCppObjectPtrConverter<ISREFERENCE>::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   if ( ! ObjectProxy_Check( pyobject ) )
      return kFALSE;              // not a PyROOT object (TODO: handle SWIG etc.)

   if ( Cppyy::IsSubtype( ((ObjectProxy*)pyobject)->ObjectIsA(), fClass ) ) {
   // depending on memory policy, some objects need releasing when passed into functions
      if ( ! KeepControl() && ! UseStrictOwnership( ctxt ) )
         ((ObjectProxy*)pyobject)->Release();

   // set pointer (may be null) and declare success
      if( ((ObjectProxy*)pyobject)->fFlags & ObjectProxy::kIsReference)
        // If given object is already a reference (aka pointer) then we should not take the address of it
        para.fValue.fVoidp = ((ObjectProxy*)pyobject)->fObject;
      else
        para.fValue.fVoidp = &((ObjectProxy*)pyobject)->fObject;
      para.fTypeCode = ISREFERENCE ? 'V' : 'p';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// construct python object from C++ instance* read at <address>

template <bool ISREFERENCE>
PyObject* PyROOT::TCppObjectPtrConverter<ISREFERENCE>::FromMemory( void* address )
{
   return BindCppObject( address, fClass, kTRUE );
}

////////////////////////////////////////////////////////////////////////////////
/// convert <value> to C++ instance*, write it at <address>

template <bool ISREFERENCE>
Bool_t PyROOT::TCppObjectPtrConverter<ISREFERENCE>::ToMemory( PyObject* value, void* address )
{
   if ( ! ObjectProxy_Check( value ) )
      return kFALSE;              // not a PyROOT object (TODO: handle SWIG etc.)

   if ( Cppyy::IsSubtype( ((ObjectProxy*)value)->ObjectIsA(), fClass ) ) {
   // depending on memory policy, some objects need releasing when passed into functions
      if ( ! KeepControl() && TCallContext::sMemoryPolicy != TCallContext::kUseStrict )
         ((ObjectProxy*)value)->Release();

   // set pointer (may be null) and declare success
      *(void**)address = ((ObjectProxy*)value)->GetObject();
      return kTRUE;
   }

   return kFALSE;
}


namespace PyROOT {
////////////////////////////////////////////////////////////////////////////////
/// Instantiate the templates

template class TCppObjectPtrConverter<true>;
template class TCppObjectPtrConverter<false>;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ instance**, set arg for call

Bool_t PyROOT::TCppObjectArrayConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* txt */ )
{
   if ( ! TTupleOfInstances_CheckExact( pyobject ) )
      return kFALSE;              // no guarantee that the tuple is okay

// treat the first instance of the tuple as the start of the array, and pass it
// by pointer (TODO: store and check sizes)
   if ( PyTuple_Size( pyobject ) < 1 )
      return kFALSE;

   PyObject* first = PyTuple_GetItem( pyobject, 0 );
   if ( ! ObjectProxy_Check( first ) )
      return kFALSE;              // should not happen

   if ( Cppyy::IsSubtype( ((ObjectProxy*)first)->ObjectIsA(), fClass ) ) {
   // no memory policies supported; set pointer (may be null) and declare success
      para.fValue.fVoidp = ((ObjectProxy*)first)->fObject;
      para.fTypeCode = 'p';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// construct python tuple of instances from C++ array read at <address>

PyObject* PyROOT::TCppObjectArrayConverter::FromMemory( void* address )
{
   if ( m_size <= 0 )   // if size unknown, just hand out the first object
      return BindCppObjectNoCast( address, fClass );

   return BindCppObjectArray( address, fClass, m_size );
}

////////////////////////////////////////////////////////////////////////////////
/// convert <value> to C++ array of instances, write it at <address>

Bool_t PyROOT::TCppObjectArrayConverter::ToMemory( PyObject* /* value */, void* /* address */ )
{
// TODO: need to have size both for the array and from the input
   PyErr_SetString( PyExc_NotImplementedError,
      "access to C-arrays of objects not yet implemented!" );
   return kFALSE;
}

//____________________________________________________________________________
// CLING WORKAROUND -- classes for STL iterators are completely undefined in that
// they come in a bazillion different guises, so just do whatever
Bool_t PyROOT::TSTLIteratorConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( ! ObjectProxy_Check( pyobject ) )
      return kFALSE;

// just set the pointer value, no check
   ObjectProxy* pyobj = (ObjectProxy*)pyobject;
   para.fValue.fVoidp = pyobj->GetObject();
   para.fTypeCode = 'V';
   return kTRUE;
}
// -- END CLING WORKAROUND

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ void*&, set arg for call

Bool_t PyROOT::TVoidPtrRefConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( ObjectProxy_Check( pyobject ) ) {
      para.fValue.fVoidp = &((ObjectProxy*)pyobject)->fObject;
      para.fTypeCode = 'V';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// convert `pyobject` to C++ void**, set arg for call

Bool_t PyROOT::TVoidPtrPtrConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   if ( ObjectProxy_Check( pyobject ) ) {
   // this is a ROOT object, take and set its address
      para.fValue.fVoidp = &((ObjectProxy*)pyobject)->fObject;
      para.fTypeCode = 'p';
      return kTRUE;
   }

// buffer objects are allowed under "user knows best"
   int buflen = Utility::GetBuffer( pyobject, '*', 1, para.fValue.fVoidp, kFALSE );

// ok if buffer exists (can't perform any useful size checks)
   if ( para.fValue.fVoidp && buflen != 0 ) {
      para.fTypeCode = 'p';
      return kTRUE;
   }

   return kFALSE;
}

////////////////////////////////////////////////////////////////////////////////
/// read a void** from address; since this is unknown, long is used (user can cast)

PyObject* PyROOT::TVoidPtrPtrConverter::FromMemory( void* address )
{
   if ( ! address || *(ptrdiff_t*)address == 0 ) {
      Py_INCREF( gNullPtrObject );
      return gNullPtrObject;
   }
   return BufFac_t::Instance()->PyBuffer_FromMemory( (Long_t*)*(ptrdiff_t**)address, sizeof(void*) );
}

////////////////////////////////////////////////////////////////////////////////
/// by definition: set and declare success

Bool_t PyROOT::TPyObjectConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* /* ctxt */ )
{
   para.fValue.fVoidp = pyobject;
   para.fTypeCode = 'p';
   return kTRUE;
}

PyObject* PyROOT::TPyObjectConverter::FromMemory( void* address )
{
// construct python object from C++ PyObject* read at <address>
   PyObject* pyobject = *((PyObject**)address);

   if ( ! pyobject ) {
      Py_INCREF( Py_None );
      return Py_None;
   }

   Py_INCREF( pyobject );
   return pyobject;
}

Bool_t PyROOT::TPyObjectConverter::ToMemory( PyObject* value, void* address )
{
// no conversion needed, write <value> at <address>
   Py_INCREF( value );
   *((PyObject**)address) = value;
   return kTRUE;
}


////////////////////////////////////////////////////////////////////////////////
/// smart pointer converter

Bool_t PyROOT::TSmartPtrCppObjectConverter::SetArg(
      PyObject* pyobject, TParameter& para, TCallContext* ctxt )
{
   char typeCode = fHandlePtr ? 'p' : 'V';

   if ( ! ObjectProxy_Check( pyobject ) ) {
      if ( fHandlePtr && GetAddressSpecialCase( pyobject, para.fValue.fVoidp ) ) {
         para.fTypeCode = typeCode;      // allow special cases such as NULL
         return kTRUE;
      }

      return kFALSE;
   }

   ObjectProxy* pyobj = (ObjectProxy*)pyobject;

// for the case where we have a 'hidden' smart pointer:
   if ( pyobj->fFlags & ObjectProxy::kIsSmartPtr && Cppyy::IsSubtype( pyobj->fSmartPtrType, fClass ) ) {
   // depending on memory policy, some objects need releasing when passed into functions
      if ( fKeepControl && ! UseStrictOwnership( ctxt ) )
         ((ObjectProxy*)pyobject)->Release();

   // calculate offset between formal and actual arguments
      para.fValue.fVoidp = pyobj->fSmartPtr;
      if ( pyobj->fSmartPtrType != fClass ) {
         para.fValue.fLong += Cppyy::GetBaseOffset(
            pyobj->fSmartPtrType, fClass, para.fValue.fVoidp, 1 /* up-cast */ );
      }

   // set pointer (may be null) and declare success
      para.fTypeCode = typeCode;
      return kTRUE;
   }

// for the case where we have an 'exposed' smart pointer:
   if ( pyobj->ObjectIsA() && Cppyy::IsSubtype( pyobj->ObjectIsA(), fClass ) ) {
   // calculate offset between formal and actual arguments
      para.fValue.fVoidp = pyobj->GetObject();
      if ( pyobj->ObjectIsA() != fClass ) {
         para.fValue.fLong += Cppyy::GetBaseOffset(
            pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */ );
      }

   // set pointer (may be null) and declare success
      para.fTypeCode = typeCode;
      return kTRUE;
   }

   return kFALSE;
}

PyObject* PyROOT::TSmartPtrCppObjectConverter::FromMemory( void* address )
{
   if ( !address || !fClass )
      return nullptr;

// obtain raw pointer
   std::vector<TParameter> args;
   ObjectProxy* pyobj = (ObjectProxy*) BindCppObject(
      Cppyy::CallR( (Cppyy::TCppMethod_t)fDereferencer, address, &args ), fRawPtrType );
   if ( pyobj )
      pyobj->SetSmartPtr( (void*)address, fClass );

   return (PyObject*)pyobj;
}


////////////////////////////////////////////////////////////////////////////////
/// raise a NotImplemented exception to take a method out of overload resolution

Bool_t PyROOT::TNotImplementedConverter::SetArg( PyObject*, TParameter&, TCallContext* )
{
   PyErr_SetString( PyExc_NotImplementedError, "this method can not (yet) be called" );
   return kFALSE;
}


//- factories -----------------------------------------------------------------
PyROOT::TConverter* PyROOT::CreateConverter( const std::string& fullType, Long_t size )
{
// The matching of the fulltype to a converter factory goes through up to five levels:
//   1) full, exact match
//   2) match of decorated, unqualified type
//   3) accept const ref as by value
//   4) accept ref as pointer
//   5) generalized cases (covers basically all ROOT classes)
//
// If all fails, void is used, which will generate a run-time warning when used.

// an exactly matching converter is best
   ConvFactories_t::iterator h = gConvFactories.find( fullType );
   if ( h != gConvFactories.end() )
      return (h->second)( size );

// resolve typedefs etc.
   std::string resolvedType = Cppyy::ResolveName( fullType );

// a full, qualified matching converter is preferred
   h = gConvFactories.find( resolvedType );
   if ( h != gConvFactories.end() )
      return (h->second)( size );

//-- nothing? ok, collect information about the type and possible qualifiers/decorators
   const std::string& cpd = Utility::Compound( resolvedType );
   std::string realType   = TClassEdit::ShortType( resolvedType.c_str(), 1 );

// accept unqualified type (as python does not know about qualifiers)
   h = gConvFactories.find( realType + cpd );
   if ( h != gConvFactories.end() )
      return (h->second)( size );

// CLING WORKAROUND -- if the type is a fixed-size array, it will have a funky
// resolved type like MyClass(&)[N], which TClass::GetClass() fails on. So, strip
// it down:
   if ( cpd == "[]" )
      realType = TClassEdit::CleanType( realType.substr( 0, realType.rfind("(") ).c_str(), 1 );
// -- CLING WORKAROUND

//-- still nothing? try pointer instead of array (for builtins)
   if ( cpd == "[]" ) {
      h = gConvFactories.find( realType + "*" );
      if ( h != gConvFactories.end() )
         return (h->second)( size );
   }

//-- still nothing? use a generalized converter
   Bool_t isConst = resolvedType.substr(0, 5) == "const";
   Bool_t control = cpd == "&" || isConst;

// converters for known/ROOT classes and default (void*)
   TConverter* result = 0;
   if ( Cppyy::TCppScope_t klass = Cppyy::GetScope( realType ) ) {
      if ( Cppyy::IsSmartPtr( realType ) ) {
	const std::vector< Cppyy::TCppMethod_t > methods = Cppyy::GetMethodsFromName( klass, "operator->", /*bases?*/ true );
         if ( ! methods.empty() ) {
            Cppyy::TCppType_t rawPtrType = Cppyy::GetScope(
               TClassEdit::ShortType( Cppyy::GetMethodResultType( methods[0] ).c_str(), 1 ) );
            if ( rawPtrType ) {
               if ( cpd == "" ) {
                  result = new TSmartPtrCppObjectConverter( klass, rawPtrType, methods[0], control );
               } else if ( cpd == "&" ) {
                  result = new TSmartPtrCppObjectConverter( klass, rawPtrType, methods[0] );
               } else if ( cpd == "*" && size <= 0 ) {
                  result = new TSmartPtrCppObjectConverter( klass, rawPtrType, methods[0], control, kTRUE );
               } /* else if ( cpd == "**" || cpd == "*&" || cpd == "&*" ) {
               } else if ( cpd == "[]" || size > 0 ) {
               } else {
               } */
            }
         }
      }

      if ( ! result ) {
        // CLING WORKAROUND -- special case for STL iterators
        if ( realType.find( "__gnu_cxx::__normal_iterator", 0 ) /* vector */ == 0 )
          result = new TSTLIteratorConverter();
        else
          // -- CLING WORKAROUND
        if ( cpd == "**" || cpd == "&*" )
          result = new TCppObjectPtrConverter<false>( klass, control);
        else if ( cpd == "*&" )
          result = new TCppObjectPtrConverter<true>( klass, control);
        else if ( cpd == "*" && size <= 0 )
          result = new TCppObjectConverter( klass, control );
        else if ( cpd == "&" )
          result = new TRefCppObjectConverter( klass );
        else if ( cpd == "[]" || size > 0 )
          result = new TCppObjectArrayConverter( klass, size, kFALSE );
        else if ( cpd == "" )               // by value
          result = new TValueCppObjectConverter( klass, kTRUE );
      }
   } else if ( Cppyy::IsEnum( realType ) ) {
      // Get underlying type of enum
      std::string et(TClassEdit::ResolveTypedef(Cppyy::ResolveEnum(realType).c_str()));
      if (cpd == "&") {
         auto reft = et + "&";
         h = isConst ? gConvFactories.find("const " + reft) : gConvFactories.find(reft);
      } else
         h = gConvFactories.find(et);
   } else if ( realType.find( "(*)" ) != std::string::npos ||
             ( realType.find( "::*)" ) != std::string::npos ) ) {
   // this is a function function pointer
   // TODO: find better way of finding the type
   // TODO: a converter that generates wrappers as appropriate
      h = gConvFactories.find( "void*" );
   }

   if ( ! result && cpd == "&&" )                  // moves
      result = new TNotImplementedConverter();

   if ( ! result && h != gConvFactories.end() )
   // converter factory available, use it to create converter
      result = (h->second)( size );
   else if ( ! result ) {
      if ( cpd != "" ) {
         std::stringstream s;
         s << "creating converter for unknown type \"" << fullType << "\"" << std::ends;
         PyErr_Warn( PyExc_RuntimeWarning, (char*)s.str().c_str() );
         result = new TVoidArrayConverter();       // "user knows best"
      } else
         result = new TVoidConverter();            // fails on use
   }

   return result;
}

////////////////////////////////////////////////////////////////////////////////

#define PYROOT_BASIC_CONVERTER_FACTORY( name )                               \
TConverter* Create##name##Converter( Long_t )                                \
{                                                                            \
   return new T##name##Converter();                                          \
}

#define PYROOT_ARRAY_CONVERTER_FACTORY( name )                               \
TConverter* Create##name##Converter( Long_t size )                           \
{                                                                            \
   return new T##name##Converter( size );                                    \
}

////////////////////////////////////////////////////////////////////////////////

namespace {
   using namespace PyROOT;

// use macro rather than template for portability ...
   PYROOT_BASIC_CONVERTER_FACTORY( Bool )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstBoolRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Char )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstCharRef )
   PYROOT_BASIC_CONVERTER_FACTORY( UChar )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstUCharRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Short )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstShortRef )
   PYROOT_BASIC_CONVERTER_FACTORY( UShort )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstUShortRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Int )
   PYROOT_BASIC_CONVERTER_FACTORY( IntRef )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstIntRef )
   PYROOT_BASIC_CONVERTER_FACTORY( UInt )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstUIntRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Long )
   PYROOT_BASIC_CONVERTER_FACTORY( LongRef )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstLongRef )
   PYROOT_BASIC_CONVERTER_FACTORY( ULong )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstULongRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Float )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstFloatRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Double )
   PYROOT_BASIC_CONVERTER_FACTORY( DoubleRef )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstDoubleRef )
   PYROOT_BASIC_CONVERTER_FACTORY( LongDouble )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstLongDoubleRef )
   PYROOT_BASIC_CONVERTER_FACTORY( Void )
   PYROOT_BASIC_CONVERTER_FACTORY( LongLong )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstLongLongRef )
   PYROOT_BASIC_CONVERTER_FACTORY( ULongLong )
   PYROOT_BASIC_CONVERTER_FACTORY( ConstULongLongRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( CString )
   PYROOT_ARRAY_CONVERTER_FACTORY( NonConstCString )
   PYROOT_ARRAY_CONVERTER_FACTORY( NonConstUCString )
   PYROOT_ARRAY_CONVERTER_FACTORY( BoolArray )
   PYROOT_BASIC_CONVERTER_FACTORY( BoolArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( ShortArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( ShortArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( UShortArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( UShortArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( IntArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( UIntArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( UIntArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( LongArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( ULongArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( ULongArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( FloatArray )
   PYROOT_ARRAY_CONVERTER_FACTORY( FloatArrayRef )
   PYROOT_ARRAY_CONVERTER_FACTORY( DoubleArray )
   PYROOT_BASIC_CONVERTER_FACTORY( VoidArray )
   PYROOT_BASIC_CONVERTER_FACTORY( LongLongArray )
   PYROOT_BASIC_CONVERTER_FACTORY( TString )
   PYROOT_BASIC_CONVERTER_FACTORY( STLString )
   PYROOT_BASIC_CONVERTER_FACTORY( STLStringView )
   PYROOT_BASIC_CONVERTER_FACTORY( VoidPtrRef )
   PYROOT_BASIC_CONVERTER_FACTORY( VoidPtrPtr )
   PYROOT_BASIC_CONVERTER_FACTORY( PyObject )

// converter factories for ROOT types
   typedef std::pair< const char*, ConverterFactory_t > NFp_t;

   // clang-format off
   NFp_t factories_[] = {
   // factories for built-ins
      NFp_t( "bool",                      &CreateBoolConverter               ),
      NFp_t( "const bool&",               &CreateConstBoolRefConverter       ),
      NFp_t( "char",                      &CreateCharConverter               ),
      NFp_t( "const char&",               &CreateConstCharRefConverter       ),
      NFp_t( "signed char",               &CreateCharConverter               ),
      NFp_t( "const signed char&",        &CreateConstCharRefConverter       ),
      NFp_t( "unsigned char",             &CreateUCharConverter              ),
      NFp_t( "const unsigned char&",      &CreateConstUCharRefConverter      ),
      NFp_t( "short",                     &CreateShortConverter              ),
      NFp_t( "const short&",              &CreateConstShortRefConverter      ),
      NFp_t( "unsigned short",            &CreateUShortConverter             ),
      NFp_t( "const unsigned short&",     &CreateConstUShortRefConverter     ),
      NFp_t( "int",                       &CreateIntConverter                ),
      NFp_t( "int&",                      &CreateIntRefConverter             ),
      NFp_t( "const int&",                &CreateConstIntRefConverter        ),
      NFp_t( "unsigned int",              &CreateUIntConverter               ),
      NFp_t( "const unsigned int&",       &CreateConstUIntRefConverter       ),
      NFp_t( "long",                      &CreateLongConverter               ),
      NFp_t( "long&",                     &CreateLongRefConverter            ),
      NFp_t( "const long&",               &CreateConstLongRefConverter       ),
      NFp_t( "unsigned long",             &CreateULongConverter              ),
      NFp_t( "const unsigned long&",      &CreateConstULongRefConverter      ),
      NFp_t( "long long",                 &CreateLongLongConverter           ),
      NFp_t( "const long long&",          &CreateConstLongLongRefConverter   ),
      NFp_t( "Long64_t",                  &CreateLongLongConverter           ),
      NFp_t( "const Long64_t&",           &CreateConstLongLongRefConverter   ),
      NFp_t( "unsigned long long",        &CreateULongLongConverter          ),
      NFp_t( "const unsigned long long&", &CreateConstULongLongRefConverter  ),
      NFp_t( "ULong64_t",                 &CreateULongLongConverter          ),
      NFp_t( "const ULong64_t&",          &CreateConstULongLongRefConverter  ),

      NFp_t( "float",                     &CreateFloatConverter              ),
      NFp_t( "const float&",              &CreateConstFloatRefConverter      ),
      NFp_t( "double",                    &CreateDoubleConverter             ),
      NFp_t( "double&",                   &CreateDoubleRefConverter          ),
      NFp_t( "const double&",             &CreateConstDoubleRefConverter     ),
      NFp_t( "long double",               &CreateLongDoubleConverter         ),
      NFp_t( "const long double&",        &CreateConstLongDoubleRefConverter ),
      NFp_t( "void",                      &CreateVoidConverter               ),

   // pointer/array factories
      NFp_t( "bool*",                     &CreateBoolArrayConverter          ),
      NFp_t( "bool&",                     &CreateBoolArrayRefConverter       ),
      NFp_t( "const unsigned char*",      &CreateCStringConverter            ),
      NFp_t( "unsigned char*",            &CreateNonConstUCStringConverter   ),
      NFp_t( "short*",                    &CreateShortArrayConverter         ),
      NFp_t( "short&",                    &CreateShortArrayRefConverter      ),
      NFp_t( "unsigned short*",           &CreateUShortArrayConverter        ),
      NFp_t( "unsigned short&",           &CreateUShortArrayRefConverter     ),
      NFp_t( "int*",                      &CreateIntArrayConverter           ),
      NFp_t( "unsigned int*",             &CreateUIntArrayConverter          ),
      NFp_t( "unsigned int&",             &CreateUIntArrayRefConverter       ),
      NFp_t( "long*",                     &CreateLongArrayConverter          ),
      NFp_t( "unsigned long*",            &CreateULongArrayConverter         ),
      NFp_t( "unsigned long&",            &CreateULongArrayRefConverter      ),
      NFp_t( "float*",                    &CreateFloatArrayConverter         ),
      NFp_t( "float&",                    &CreateFloatArrayRefConverter      ),
      NFp_t( "double*",                   &CreateDoubleArrayConverter        ),
      NFp_t( "long long*",                &CreateLongLongArrayConverter      ),
      NFp_t( "Long64_t*",                 &CreateLongLongArrayConverter      ),
      NFp_t( "unsigned long long*",       &CreateLongLongArrayConverter      ),  // TODO: ULongLong
      NFp_t( "ULong64_t*",                &CreateLongLongArrayConverter      ),  // TODO: ULongLong
      NFp_t( "void*",                     &CreateVoidArrayConverter          ),

   // factories for special cases
      NFp_t( "const char*",               &CreateCStringConverter            ),
      NFp_t( "char*",                     &CreateNonConstCStringConverter    ),
      NFp_t( "TString",                   &CreateTStringConverter            ),
      NFp_t( "const TString&",            &CreateTStringConverter            ),
      NFp_t( "std::string",               &CreateSTLStringConverter          ),
      NFp_t( "string",                    &CreateSTLStringConverter          ),
      NFp_t( "const std::string&",        &CreateSTLStringConverter          ),
      NFp_t( "const string&",             &CreateSTLStringConverter          ),
      NFp_t( "std::string_view",          &CreateSTLStringViewConverter      ),
      NFp_t( "string_view",               &CreateSTLStringViewConverter      ),
      NFp_t( "experimental::basic_string_view<char,char_traits<char> >",&CreateSTLStringViewConverter),
      NFp_t( "basic_string_view<char,char_traits<char> >",&CreateSTLStringViewConverter),
      NFp_t( "void*&",                    &CreateVoidPtrRefConverter         ),
      NFp_t( "void**",                    &CreateVoidPtrPtrConverter         ),
      NFp_t( "PyObject*",                 &CreatePyObjectConverter           ),
      NFp_t( "_object*",                  &CreatePyObjectConverter           ),
      NFp_t( "FILE*",                     &CreateVoidArrayConverter          ),
      NFp_t( "Float16_t",                 &CreateFloatConverter              ),
      NFp_t( "const Float16_t&",          &CreateConstFloatRefConverter      ),
      NFp_t( "Double32_t",                &CreateDoubleConverter             ),
      NFp_t( "Double32_t&",               &CreateDoubleRefConverter          ),
      NFp_t( "const Double32_t&",         &CreateConstDoubleRefConverter     )
   };
   // clang-format on

   struct InitConvFactories_t {
   public:
      InitConvFactories_t()
      {
      // load all converter factories in the global map 'gConvFactories'
         int nf = sizeof( factories_ ) / sizeof( factories_[ 0 ] );
         for ( int i = 0; i < nf; ++i ) {
            gConvFactories[ factories_[ i ].first ] = factories_[ i ].second;
         }
      }
   } initConvFactories_;

} // unnamed namespace
