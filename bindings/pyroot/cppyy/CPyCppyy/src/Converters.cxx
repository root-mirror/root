// Bindings
#include "CPyCppyy.h"
#include "DeclareConverters.h"
#include "CallContext.h"
#include "CPPExcInstance.h"
#include "CPPInstance.h"
#include "CPPOverload.h"
#include "CustomPyTypes.h"
#include "LowLevelViews.h"
#include "MemoryRegulator.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"
#include "TemplateProxy.h"
#include "TupleOfInstances.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <limits.h>
#include <stddef.h>      // for ptrdiff_t
#include <string.h>
#include <array>
#include <utility>
#include <sstream>
#if __cplusplus > 201402L
#include <cstddef>
#endif
#include "ROOT/RStringView.hxx"

#define UNKNOWN_SIZE         -1
#define UNKNOWN_ARRAY_SIZE   -2


//- data _____________________________________________________________________
namespace CPyCppyy {

// factories
    typedef std::map<std::string, cf_t> ConvFactories_t;
    static ConvFactories_t gConvFactories;
    extern PyObject* gNullPtrObject;

}

#if PY_VERSION_HEX < 0x03000000
const size_t MOVE_REFCOUNT_CUTOFF = 1;
#else
// p3 has at least 2 ref-counts, as contrary to p2, it will create a descriptor
// copy for the method holding self in the case of __init__; but there can also
// be a reference held by the frame object, which is indistinguishable from a
// local variable reference, so the cut-off has to remain 2.
const size_t MOVE_REFCOUNT_CUTOFF = 2;
#endif

//- pretend-ctypes helpers ---------------------------------------------------
struct CPyCppyy_tagCDataObject {       // non-public (but stable)
    PyObject_HEAD
    char* b_ptr;
    int   b_needsfree;
};

struct CPyCppyy_tagPyCArgObject {      // not public (but stable; note that older
    PyObject_HEAD                      // Pythons protect 'D' with HAVE_LONG_LONG)
    void* pffi_type;
    char tag;
    union {                            // for convenience, kept only relevant vals
        long long q;
        long double D;
        void *p;
    } value;
    PyObject* obj;
};

// indices of ctypes types into the array caches (not that c_complex does not
// exist as a type in ctypes)
#define ct_c_bool        0
#define ct_c_char        1
#define ct_c_shar        1
#define ct_c_wchar       2
#define ct_c_byte        3
#define ct_c_int8        3
#define ct_c_ubyte       4
#define ct_c_uchar       4
#define ct_c_uint8       4
#define ct_c_short       5
#define ct_c_ushort      6
#define ct_c_uint16      7
#define ct_c_int         8
#define ct_c_uint        9
#define ct_c_uint32     10
#define ct_c_long       11
#define ct_c_ulong      12
#define ct_c_longlong   13
#define ct_c_ulonglong  14
#define ct_c_float      15
#define ct_c_double     16
#define ct_c_longdouble 17
#define ct_c_char_p     18
#define ct_c_wchar_p    19
#define ct_c_void_p     20
#define ct_c_complex    21
#define NTYPES          22

static std::array<const char*, NTYPES> gCTypesNames = {
    "c_bool", "c_char", "c_wchar", "c_byte", "c_ubyte", "c_short", "c_ushort", "c_uint16",
    "c_int", "c_uint", "c_uint32", "c_long", "c_ulong", "c_longlong", "c_ulonglong",
    "c_float", "c_double", "c_longdouble",
    "c_char_p", "c_wchar_p", "c_void_p", "c_complex" };
static std::array<PyTypeObject*, NTYPES> gCTypesTypes;
static std::array<PyTypeObject*, NTYPES> gCTypesPtrTypes;

// Both GetCTypesType and GetCTypesPtrType, rely on the ctypes module itself
// caching the types (thus also making them unique), so no ref-count is needed.
// Further, by keeping a ref-count on the module, it won't be off-loaded until
// the 2nd cleanup cycle.
static PyTypeObject* GetCTypesType(int nidx)
{
    static PyObject* ctmod = PyImport_ImportModule("ctypes");   // ref-count kept
    if (!ctmod) {
        PyErr_Clear();
        return nullptr;
    }
    PyTypeObject* ct_t = gCTypesTypes[nidx];
    if (!ct_t) {
        ct_t = (PyTypeObject*)PyObject_GetAttrString(ctmod, gCTypesNames[nidx]);
        if (!ct_t) PyErr_Clear();
        else {
            gCTypesTypes[nidx] = ct_t;
            Py_DECREF(ct_t);
        }
    }
    return ct_t;
}

static PyTypeObject* GetCTypesPtrType(int nidx)
{
    static PyObject* ctmod = PyImport_ImportModule("ctypes");   // ref-count kept
    if (!ctmod) {
        PyErr_Clear();
        return nullptr;
    }
    PyTypeObject* cpt_t = gCTypesPtrTypes[nidx];
    if (!cpt_t) {
        if (strcmp(gCTypesNames[nidx], "c_char") == 0)  {
            cpt_t = (PyTypeObject*)PyObject_GetAttrString(ctmod, "c_char_p");
        } else {
            PyObject* ct_t = (PyObject*)GetCTypesType(nidx);
            if (ct_t) {
                PyObject* ptrcreat = PyObject_GetAttrString(ctmod, "POINTER");
                cpt_t = (PyTypeObject*)PyObject_CallFunctionObjArgs(ptrcreat, ct_t, NULL);
                Py_DECREF(ptrcreat);
            }
        }
        if (cpt_t) {
            gCTypesPtrTypes[nidx] = cpt_t;
            Py_DECREF(cpt_t);
        }
    }
    return cpt_t;
}

static bool IsPyCArgObject(PyObject* pyobject)
{
    static PyTypeObject* pycarg_type = nullptr;
    if (!pycarg_type) {
        PyObject* ctmod = PyImport_ImportModule("ctypes");
        if (!ctmod) PyErr_Clear();
        else {
             PyTypeObject* ct_t = (PyTypeObject*)PyObject_GetAttrString(ctmod, "c_int");
             PyObject* cobj = ct_t->tp_new(ct_t, nullptr, nullptr);
             PyObject* byref = PyObject_GetAttrString(ctmod, "byref");
             PyObject* pyptr = PyObject_CallFunctionObjArgs(byref, cobj, NULL);
             Py_DECREF(byref); Py_DECREF(cobj); Py_DECREF(ct_t);
             pycarg_type = Py_TYPE(pyptr);  // static, no ref-count needed
             Py_DECREF(pyptr);
        }
        Py_DECREF(ctmod);
    }
    return Py_TYPE(pyobject) == pycarg_type;
}

static bool IsCTypesArrayOrPointer(PyObject* pyobject)
{
    static PyTypeObject* cstgdict_type = nullptr;
    if (!cstgdict_type) {
    // get any pointer type to initialize the extended dictionary type
        PyTypeObject* ct_int = GetCTypesType(ct_c_int);
        if (ct_int && ct_int->tp_dict) {
            cstgdict_type = Py_TYPE(ct_int->tp_dict);
        }
    }

    PyTypeObject* pytype = Py_TYPE(pyobject);
    if (pytype->tp_dict && Py_TYPE(pytype->tp_dict) == cstgdict_type)
        return true;
    return false;
}


//- helper to work with both CPPInstance and CPPExcInstance ------------------
static inline CPyCppyy::CPPInstance* GetCppInstance(PyObject* pyobject)
{
    using namespace CPyCppyy;
    if (CPPInstance_Check(pyobject))
        return (CPPInstance*)pyobject;
    if (CPPExcInstance_Check(pyobject))
        return (CPPInstance*)((CPPExcInstance*)pyobject)->fCppInstance;
    return nullptr;
}


//- custom helpers to check ranges -------------------------------------------
static inline bool CPyCppyy_PyLong_AsBool(PyObject* pyobject)
{
// range-checking python integer to C++ bool conversion
    long l = PyLong_AsLong(pyobject);
// fail to pass float -> bool; the problem is rounding (0.1 -> 0 -> False)
    if (!(l == 0|| l == 1) || PyFloat_Check(pyobject)) {
        PyErr_SetString(PyExc_ValueError, "boolean value should be bool, or integer 1 or 0");
        return (bool)-1;
    }
    return (bool)l;
}

static inline char CPyCppyy_PyText_AsChar(PyObject* pyobject) {
// python string to C++ char conversion
    return (char)CPyCppyy_PyText_AsString(pyobject)[0];
}

static inline uint8_t CPyCppyy_PyLong_AsUInt8(PyObject* pyobject)
{
// range-checking python integer to C++ uint8_t conversion (typically, an unsigned char)
// prevent p2.7 silent conversions and do a range check
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "short int conversion expects an integer object");
        return (uint8_t)-1;
    }
    long l = PyLong_AsLong(pyobject);
    if (l < 0 || UCHAR_MAX < l) {
        PyErr_Format(PyExc_ValueError, "integer %ld out of range for uint8_t", l);
        return (uint8_t)-1;

    }
    return (uint8_t)l;
}

static inline int8_t CPyCppyy_PyLong_AsInt8(PyObject* pyobject)
{
// range-checking python integer to C++ int8_t conversion (typically, an signed char)
// prevent p2.7 silent conversions and do a range check
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "short int conversion expects an integer object");
        return (int8_t)-1;
    }
    long l = PyLong_AsLong(pyobject);
    if (l < SCHAR_MIN || SCHAR_MAX < l) {
        PyErr_Format(PyExc_ValueError, "integer %ld out of range for int8_t", l);
        return (int8_t)-1;

    }
    return (int8_t)l;
}

static inline unsigned short CPyCppyy_PyLong_AsUShort(PyObject* pyobject)
{
// range-checking python integer to C++ unsigend short int conversion

// prevent p2.7 silent conversions and do a range check
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "unsigned short conversion expects an integer object");
        return (unsigned short)-1;
    }
    long l = PyLong_AsLong(pyobject);
    if (l < 0 || USHRT_MAX < l) {
        PyErr_Format(PyExc_ValueError, "integer %ld out of range for unsigned short", l);
        return (unsigned short)-1;

    }
    return (unsigned short)l;
}

static inline short CPyCppyy_PyLong_AsShort(PyObject* pyobject)
{
// range-checking python integer to C++ short int conversion
// prevent p2.7 silent conversions and do a range check
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "short int conversion expects an integer object");
        return (short)-1;
    }
    long l = PyLong_AsLong(pyobject);
    if (l < SHRT_MIN || SHRT_MAX < l) {
        PyErr_Format(PyExc_ValueError, "integer %ld out of range for short int", l);
        return (short)-1;

    }
    return (short)l;
}

static inline int CPyCppyy_PyLong_AsStrictInt(PyObject* pyobject)
{
// strict python integer to C++ integer conversion

// p2.7 and later silently converts floats to long, therefore require this
// check; earlier pythons may raise a SystemError which should be avoided as
// it is confusing
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "int/long conversion expects an integer object");
        return -1;
    }
    long l = PyLong_AsLong(pyobject);
    if (l < INT_MIN || INT_MAX < l) {
        PyErr_Format(PyExc_ValueError, "integer %ld out of range for int", l);
        return (int)-1;

    }
    return (int)l;
}

static inline long CPyCppyy_PyLong_AsStrictLong(PyObject* pyobject)
{
// strict python integer to C++ long integer conversion

// prevent float -> long (see CPyCppyy_PyLong_AsStrictInt)
    if (!(PyLong_Check(pyobject) || PyInt_Check(pyobject))) {
        PyErr_SetString(PyExc_TypeError, "int/long conversion expects an integer object");
        return (long)-1;
    }
    return (long)PyLong_AsLong(pyobject);
}


//- helper for pointer/array/reference conversions ---------------------------
static inline bool CArraySetArg(PyObject* pyobject, CPyCppyy::Parameter& para, char tc, int size)
{
// general case of loading a C array pointer (void* + type code) as function argument
    if (pyobject == CPyCppyy::gNullPtrObject)
        para.fValue.fVoidp = nullptr;
    else {
        Py_ssize_t buflen = CPyCppyy::Utility::GetBuffer(pyobject, tc, size, para.fValue.fVoidp);
        if (!buflen) {
        // stuck here as it's the least common
            if (CPyCppyy_PyLong_AsStrictInt(pyobject) == 0)
                para.fValue.fVoidp = nullptr;
            else {
                PyErr_Format(PyExc_TypeError,     // ValueError?
                   "could not convert argument to buffer or nullptr");
                return false;
            }
        }
    }
    para.fTypeCode = 'p';
    return true;
}


//- helper for implicit conversions ------------------------------------------
static inline bool ConvertImplicit(Cppyy::TCppType_t klass,
    PyObject* pyobject, CPyCppyy::Parameter& para, CPyCppyy::CallContext* ctxt)
{
    using namespace CPyCppyy;

// filter out copy and move constructors
    if (IsConstructor(ctxt->fFlags) && klass == ctxt->fCurScope && ctxt->GetSize() == 1)
        return false;

// only proceed if implicit conversions are allowed (in "round 2") or if the
// argument is exactly a tuple or list, as these are the equivalent of
// initializer lists and thus "syntax" not a conversion
    if (!AllowImplicit(ctxt)) {
        PyTypeObject* pytype = (PyTypeObject*)Py_TYPE(pyobject);
        if (!(pytype == &PyList_Type || pytype == &PyTuple_Type)) {
            if (!NoImplicit(ctxt)) ctxt->fFlags |= CallContext::kHaveImplicit;
            return false;
        }
    }

// exercise implicit conversion
    PyObject* pyscope = CreateScopeProxy(klass);
    if (!CPPScope_Check(pyscope)) {
        Py_XDECREF(pyscope);
        return false;
    }

// add a pseudo-keyword argument to prevent recursion
    PyObject* kwds = PyDict_New();
    PyDict_SetItem(kwds, PyStrings::gNoImplicit, Py_True);
    PyObject* args = PyTuple_New(1);
    Py_INCREF(pyobject); PyTuple_SET_ITEM(args, 0, pyobject);

// call constructor of argument type to attempt implicit conversion
    CPPInstance* pytmp = (CPPInstance*)PyObject_Call(pyscope, args, kwds);
    if (!pytmp && PyTuple_CheckExact(pyobject)) {
    // special case: allow implicit conversion from given set of arguments in tuple
        PyErr_Clear();
        PyDict_SetItem(kwds, PyStrings::gNoImplicit, Py_True); // was deleted
        pytmp = (CPPInstance*)PyObject_Call(pyscope, pyobject, kwds);
    }

    Py_DECREF(args);
    Py_DECREF(kwds);
    Py_DECREF(pyscope);

    if (pytmp) {
    // implicit conversion succeeded!
        ctxt->AddTemporary((PyObject*)pytmp);
        para.fValue.fVoidp = pytmp->GetObject();
        para.fTypeCode = 'V';
        return true;
    }

    PyErr_Clear();
    return false;
}


//- base converter implementation --------------------------------------------
CPyCppyy::Converter::~Converter()
{
    /* empty */
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::Converter::FromMemory(void*)
{
// could happen if no derived class override
    PyErr_SetString(PyExc_TypeError, "C++ type cannot be converted from memory");
    return nullptr;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Converter::ToMemory(PyObject*, void*)
{
// could happen if no derived class override
    PyErr_SetString(PyExc_TypeError, "C++ type cannot be converted to memory");
    return false;
}


//- helper macro's -----------------------------------------------------------
#define CPPYY_IMPL_BASIC_CONVERTER(name, type, stype, ctype, F1, F2, tc)     \
bool CPyCppyy::name##Converter::SetArg(                                      \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
/* convert <pyobject> to C++ 'type', set arg for call */                     \
    type val = (type)F2(pyobject);                                           \
    if (val == (type)-1 && PyErr_Occurred()) {                               \
        static PyTypeObject* ctypes_type = nullptr;                          \
        if (!ctypes_type) {                                                  \
            PyObject* pytype = 0, *pyvalue = 0, *pytrace = 0;                \
            PyErr_Fetch(&pytype, &pyvalue, &pytrace);                        \
            ctypes_type = GetCTypesType(ct_##ctype);                         \
            PyErr_Restore(pytype, pyvalue, pytrace);                         \
        }                                                                    \
        if (Py_TYPE(pyobject) == ctypes_type) {                              \
            PyErr_Clear();                                                   \
            val = *((type*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr);     \
        } else                                                               \
            return false;                                                    \
    }                                                                        \
    para.fValue.f##name = val;                                               \
    para.fTypeCode = tc;                                                     \
    return true;                                                             \
}                                                                            \
                                                                             \
PyObject* CPyCppyy::name##Converter::FromMemory(void* address)               \
{                                                                            \
    return F1((stype)*((type*)address));                                     \
}                                                                            \
                                                                             \
bool CPyCppyy::name##Converter::ToMemory(PyObject* value, void* address)     \
{                                                                            \
    type s = (type)F2(value);                                                \
    if (s == (type)-1 && PyErr_Occurred())                                   \
        return false;                                                        \
    *((type*)address) = (type)s;                                             \
    return true;                                                             \
}

//----------------------------------------------------------------------------
static inline int ExtractChar(PyObject* pyobject, const char* tname, int low, int high)
{
    int lchar = -1;
    if (CPyCppyy_PyText_Check(pyobject)) {
        if (CPyCppyy_PyText_GET_SIZE(pyobject) == 1)
            lchar = (int)CPyCppyy_PyText_AsChar(pyobject);
        else
            PyErr_Format(PyExc_ValueError, "%s expected, got string of size " PY_SSIZE_T_FORMAT,
                tname, CPyCppyy_PyText_GET_SIZE(pyobject));
    } else if (!PyFloat_Check(pyobject)) {   // don't allow truncating conversion
        lchar = (int)PyLong_AsLong(pyobject);
        if (lchar == -1 && PyErr_Occurred())
            ; // empty, as error already set
        else if (!(low <= lchar && lchar <= high)) {
            PyErr_Format(PyExc_ValueError,
                "integer to character: value %d not in range [%d,%d]", lchar, low, high);
            lchar = -1;
        }
    } else
        PyErr_SetString(PyExc_TypeError, "char or small int type expected");

    return lchar;
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(name, ctype)                     \
PyObject* CPyCppyy::name##RefConverter::FromMemory(void* ptr)                \
{                                                                            \
/* convert a reference to int to Python through ctypes pointer object */     \
    PyTypeObject* ctypes_type = GetCTypesType(ct_##ctype);                   \
    if (!ctypes_type) {                                                      \
        PyErr_SetString(PyExc_RuntimeError, "no ctypes available");          \
        return nullptr;                                                      \
    }                                                                        \
    PyObject* ref = ctypes_type->tp_new(ctypes_type, nullptr, nullptr);      \
    ((CPyCppyy_tagCDataObject*)ref)->b_ptr = (char*)ptr;                     \
    ((CPyCppyy_tagCDataObject*)ref)->b_needsfree = 0;                        \
    return ref;                                                              \
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_BASIC_CONST_REFCONVERTER(name, type, ctype, F1)           \
bool CPyCppyy::Const##name##RefConverter::SetArg(                            \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
    type val = (type)F1(pyobject);                                           \
    if (val == (type)-1 && PyErr_Occurred())                                 \
        return false;                                                        \
    para.fValue.f##name = val;                                               \
    para.fRef = &para.fValue.f##name;                                        \
    para.fTypeCode = 'r';                                                    \
    return true;                                                             \
}                                                                            \
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(Const##name, ctype)

//----------------------------------------------------------------------------
#define CPPYY_IMPL_BASIC_CONST_CHAR_REFCONVERTER(name, type, ctype, low, high)\
bool CPyCppyy::Const##name##RefConverter::SetArg(                            \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
/* convert <pyobject> to C++ <<type>>, set arg for call, allow int -> char */\
    type val = (type)ExtractChar(pyobject, #type, low, high);                \
    if (val == (type)-1 && PyErr_Occurred())                                 \
        return false;                                                        \
    para.fValue.fLong = val;                                                 \
    para.fTypeCode = 'l';                                                    \
    return true;                                                             \
}                                                                            \
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(Const##name, ctype)


//----------------------------------------------------------------------------
#define CPPYY_IMPL_BASIC_CHAR_CONVERTER(name, type, low, high)               \
bool CPyCppyy::name##Converter::SetArg(                                      \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
/* convert <pyobject> to C++ <<type>>, set arg for call, allow int -> char */\
    long val = ExtractChar(pyobject, #type, low, high);                      \
    if (val == -1 && PyErr_Occurred())                                       \
        return false;                                                        \
    para.fValue.fLong = val;                                                 \
    para.fTypeCode = 'l';                                                    \
    return true;                                                             \
}                                                                            \
                                                                             \
PyObject* CPyCppyy::name##Converter::FromMemory(void* address)               \
{                                                                            \
    return CPyCppyy_PyText_FromFormat("%c", *((type*)address));              \
}                                                                            \
                                                                             \
bool CPyCppyy::name##Converter::ToMemory(PyObject* value, void* address)     \
{                                                                            \
    Py_ssize_t len;                                                          \
    const char* cstr = CPyCppyy_PyText_AsStringAndSize(value, &len);         \
    if (cstr) {                                                              \
        if (len != 1) {                                                      \
            PyErr_Format(PyExc_TypeError, #type" expected, got string of size %zd", len);\
            return false;                                                    \
        }                                                                    \
        *((type*)address) = (type)cstr[0];                                   \
    } else {                                                                 \
        PyErr_Clear();                                                       \
        long l = PyLong_AsLong(value);                                       \
        if (l == -1 && PyErr_Occurred())                                     \
            return false;                                                    \
        if (!(low <= l && l <= high)) {                                      \
            PyErr_Format(PyExc_ValueError,                                   \
                "integer to character: value %ld not in range [%d,%d]", l, low, high);\
            return false;                                                    \
        }                                                                    \
        *((type*)address) = (type)l;                                         \
    }                                                                        \
    return true;                                                             \
}


//- converters for built-ins -------------------------------------------------
CPPYY_IMPL_BASIC_CONVERTER(Long, long, long, c_long, PyLong_FromLong, CPyCppyy_PyLong_AsStrictLong, 'l')

//----------------------------------------------------------------------------
bool CPyCppyy::LongRefConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ long&, set arg for call
#if PY_VERSION_HEX < 0x03000000
    if (RefInt_CheckExact(pyobject)) {
        para.fValue.fVoidp = (void*)&((PyIntObject*)pyobject)->ob_ival;
        para.fTypeCode = 'V';
        return true;
    }
#endif

    if (Py_TYPE(pyobject) == GetCTypesType(ct_c_long)) {
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;
        para.fTypeCode = 'V';
        return true;
    }

    if (CArraySetArg(pyobject, para, 'l', sizeof(long))) {
        para.fTypeCode = 'V';
        return true;
    }

    PyErr_SetString(PyExc_TypeError, "use ctypes.c_long for pass-by-ref of longs");
    return false;
}

//----------------------------------------------------------------------------
CPPYY_IMPL_BASIC_CONST_CHAR_REFCONVERTER(Char,  char,          c_char,  CHAR_MIN,  CHAR_MAX)
CPPYY_IMPL_BASIC_CONST_CHAR_REFCONVERTER(UChar, unsigned char, c_uchar,        0, UCHAR_MAX)

CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Bool,   bool,           c_bool,      CPyCppyy_PyLong_AsBool)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Int8,   int8_t,         c_int8,      CPyCppyy_PyLong_AsInt8)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(UInt8,  uint8_t,        c_uint8,     CPyCppyy_PyLong_AsUInt8)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Short,  short,          c_short,     CPyCppyy_PyLong_AsShort)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(UShort, unsigned short, c_ushort,    CPyCppyy_PyLong_AsUShort)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Int,    int,            c_int,       CPyCppyy_PyLong_AsStrictInt)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(UInt,   unsigned int,   c_uint,      PyLongOrInt_AsULong)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Long,   long,           c_long,      CPyCppyy_PyLong_AsStrictLong)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(ULong,  unsigned long,  c_ulong,     PyLongOrInt_AsULong)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(LLong,  Long64_t,       c_longlong,  PyLong_AsLongLong)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(ULLong, ULong64_t,      c_ulonglong, PyLongOrInt_AsULong64)

//----------------------------------------------------------------------------
bool CPyCppyy::IntRefConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ (pseudo)int&, set arg for call
#if PY_VERSION_HEX < 0x03000000
    if (RefInt_CheckExact(pyobject)) {
        para.fValue.fVoidp = (void*)&((PyIntObject*)pyobject)->ob_ival;
        para.fTypeCode = 'V';
        return true;
    }
#endif

#if PY_VERSION_HEX >= 0x02050000
    if (Py_TYPE(pyobject) == GetCTypesType(ct_c_int)) {
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;
        para.fTypeCode = 'V';
        return true;
    }
#endif

// alternate, pass pointer from buffer
    Py_ssize_t buflen = Utility::GetBuffer(pyobject, 'i', sizeof(int), para.fValue.fVoidp);
    if (para.fValue.fVoidp && buflen) {
        para.fTypeCode = 'V';
        return true;
    };

#if PY_VERSION_HEX < 0x02050000
    PyErr_SetString(PyExc_TypeError, "use cppyy.Long for pass-by-ref of ints");
#else
    PyErr_SetString(PyExc_TypeError, "use ctypes.c_int for pass-by-ref of ints");
#endif
    return false;
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_REFCONVERTER(name, ctype, type, code)                     \
bool CPyCppyy::name##RefConverter::SetArg(                                   \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
/* convert a reference to int to Python through ctypes pointer object */     \
    if (Py_TYPE(pyobject) == GetCTypesType(ct_##ctype)) {                    \
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;\
        para.fTypeCode = 'V';                                                \
        return true;                                                         \
    }                                                                        \
    bool res = CArraySetArg(pyobject, para, code, sizeof(type));             \
    if (!res) {                                                              \
        PyErr_SetString(PyExc_TypeError, "use ctypes."#ctype" for pass-by-ref of "#type);\
        return false;                                                        \
    }                                                                        \
    para.fTypeCode = 'V';                                                    \
    return res;                                                              \
}                                                                            \
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(name, ctype)

CPPYY_IMPL_REFCONVERTER(Bool,    c_bool,       bool,               '?');
CPPYY_IMPL_REFCONVERTER(Char,    c_char,       char,               'b');
CPPYY_IMPL_REFCONVERTER(WChar,   c_wchar,      wchar_t,            'u');
CPPYY_IMPL_REFCONVERTER(Char16,  c_uint16,     char16_t,           'H');
CPPYY_IMPL_REFCONVERTER(Char32,  c_uint32,     char32_t,           'I');
CPPYY_IMPL_REFCONVERTER(SChar,   c_byte,       signed char,        'b');
CPPYY_IMPL_REFCONVERTER(UChar,   c_ubyte,      unsigned char,      'B');
CPPYY_IMPL_REFCONVERTER(Int8,    c_int8,       int8_t,             'b');
CPPYY_IMPL_REFCONVERTER(UInt8,   c_uint8,      uint8_t,            'B');
CPPYY_IMPL_REFCONVERTER(Short,   c_short,      short,              'h');
CPPYY_IMPL_REFCONVERTER(UShort,  c_ushort,     unsigned short,     'H');
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(Int, c_int);
CPPYY_IMPL_REFCONVERTER(UInt,    c_uint,       unsigned int,       'I');
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(Long, c_long);
CPPYY_IMPL_REFCONVERTER(ULong,   c_ulong,      unsigned long,      'L');
CPPYY_IMPL_REFCONVERTER(LLong,   c_longlong,   long long,          'q');
CPPYY_IMPL_REFCONVERTER(ULLong,  c_ulonglong,  unsigned long long, 'Q');
CPPYY_IMPL_REFCONVERTER(Float,   c_float,      float,              'f');
CPPYY_IMPL_REFCONVERTER_FROM_MEMORY(Double, c_double);
CPPYY_IMPL_REFCONVERTER(LDouble, c_longdouble, LongDouble_t,       'D');


//----------------------------------------------------------------------------
// convert <pyobject> to C++ bool, allow int/long -> bool, set arg for call
CPPYY_IMPL_BASIC_CONVERTER(
    Bool, bool, long, c_bool, PyInt_FromLong, CPyCppyy_PyLong_AsBool, 'l')

//----------------------------------------------------------------------------
CPPYY_IMPL_BASIC_CHAR_CONVERTER(Char,  char,          CHAR_MIN,  CHAR_MAX)
CPPYY_IMPL_BASIC_CHAR_CONVERTER(UChar, unsigned char,        0, UCHAR_MAX)

PyObject* CPyCppyy::UCharAsIntConverter::FromMemory(void* address)
{
// special case to be used with arrays: return a Python int instead of str
// (following the same convention as module array.array)
    return PyInt_FromLong((long)*((unsigned char*)address));
}

//----------------------------------------------------------------------------
bool CPyCppyy::WCharConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ <wchar_t>, set arg for call
    if (!PyUnicode_Check(pyobject) || PyUnicode_GET_SIZE(pyobject) != 1) {
        PyErr_SetString(PyExc_ValueError, "single wchar_t character expected");
        return false;
    }
    wchar_t val;
    Py_ssize_t res = CPyCppyy_PyUnicode_AsWideChar(pyobject, &val, 1);
    if (res == -1)
        return false;
    para.fValue.fLong = (long)val;
    para.fTypeCode = 'U';
    return true;
}

PyObject* CPyCppyy::WCharConverter::FromMemory(void* address)
{
    return PyUnicode_FromWideChar((const wchar_t*)address, 1);
}

bool CPyCppyy::WCharConverter::ToMemory(PyObject* value, void* address)
{
    if (!PyUnicode_Check(value) || PyUnicode_GET_SIZE(value) != 1) {
        PyErr_SetString(PyExc_ValueError, "single wchar_t character expected");
        return false;
    }
    wchar_t val;
    Py_ssize_t res = CPyCppyy_PyUnicode_AsWideChar(value, &val, 1);
    if (res == -1)
        return false;
    *((wchar_t*)address) = val;
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Char16Converter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ <char16_t>, set arg for call
    if (!PyUnicode_Check(pyobject) || PyUnicode_GET_SIZE(pyobject) != 1) {
        PyErr_SetString(PyExc_ValueError, "single char16_t character expected");
        return false;
    }

    PyObject* bstr = PyUnicode_AsUTF16String(pyobject);
    if (!bstr) return false;

    char16_t val = *(char16_t*)(PyBytes_AS_STRING(bstr) + sizeof(char16_t) /*BOM*/);
    Py_DECREF(bstr);
    para.fValue.fLong = (long)val;
    para.fTypeCode = 'U';
    return true;
}

PyObject* CPyCppyy::Char16Converter::FromMemory(void* address)
{
    return PyUnicode_DecodeUTF16((const char*)address, sizeof(char16_t), nullptr, nullptr);
}

bool CPyCppyy::Char16Converter::ToMemory(PyObject* value, void* address)
{
    if (!PyUnicode_Check(value) || PyUnicode_GET_SIZE(value) != 1) {
        PyErr_SetString(PyExc_ValueError, "single char16_t character expected");
        return false;
    }

    PyObject* bstr = PyUnicode_AsUTF16String(value);
    if (!bstr) return false;

    *((char16_t*)address) = *(char16_t*)(PyBytes_AS_STRING(bstr) + sizeof(char16_t) /*BOM*/);
    Py_DECREF(bstr);
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Char32Converter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ <char32_t>, set arg for call
    if (!PyUnicode_Check(pyobject) || 2 < PyUnicode_GET_SIZE(pyobject)) {
        PyErr_SetString(PyExc_ValueError, "single char32_t character expected");
        return false;
    }

    PyObject* bstr = PyUnicode_AsUTF32String(pyobject);
    if (!bstr) return false;

    char32_t val = *(char32_t*)(PyBytes_AS_STRING(bstr) + sizeof(char32_t) /*BOM*/);
    Py_DECREF(bstr);
    para.fValue.fLong = (long)val;
    para.fTypeCode = 'U';
    return true;
}

PyObject* CPyCppyy::Char32Converter::FromMemory(void* address)
{
    return PyUnicode_DecodeUTF32((const char*)address, sizeof(char32_t), nullptr, nullptr);
}

bool CPyCppyy::Char32Converter::ToMemory(PyObject* value, void* address)
{
    if (!PyUnicode_Check(value) || 2 < PyUnicode_GET_SIZE(value)) {
        PyErr_SetString(PyExc_ValueError, "single char32_t character expected");
        return false;
    }

    PyObject* bstr = PyUnicode_AsUTF32String(value);
    if (!bstr) return false;

    *((char32_t*)address) = *(char32_t*)(PyBytes_AS_STRING(bstr) + sizeof(char32_t) /*BOM*/);
    Py_DECREF(bstr);
    return true;
}

//----------------------------------------------------------------------------
CPPYY_IMPL_BASIC_CONVERTER(
    Int8,  int8_t,  long, c_int8, PyInt_FromLong, CPyCppyy_PyLong_AsInt8,  'l')
CPPYY_IMPL_BASIC_CONVERTER(
    UInt8, uint8_t, long, c_uint8, PyInt_FromLong, CPyCppyy_PyLong_AsUInt8, 'l')
CPPYY_IMPL_BASIC_CONVERTER(
    Short, short, long, c_short, PyInt_FromLong, CPyCppyy_PyLong_AsShort, 'l')
CPPYY_IMPL_BASIC_CONVERTER(
    UShort, unsigned short, long, c_ushort, PyInt_FromLong, CPyCppyy_PyLong_AsUShort, 'l')
CPPYY_IMPL_BASIC_CONVERTER(
    Int, int, long, c_uint, PyInt_FromLong, CPyCppyy_PyLong_AsStrictInt, 'l')

//----------------------------------------------------------------------------
bool CPyCppyy::ULongConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ unsigned long, set arg for call
    para.fValue.fULong = PyLongOrInt_AsULong(pyobject);
    if (para.fValue.fULong == (unsigned long)-1 && PyErr_Occurred())
        return false;
    para.fTypeCode = 'L';
    return true;
}

PyObject* CPyCppyy::ULongConverter::FromMemory(void* address)
{
// construct python object from C++ unsigned long read at <address>
    return PyLong_FromUnsignedLong(*((unsigned long*)address));
}

bool CPyCppyy::ULongConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ unsigned long, write it at <address>
    unsigned long u = PyLongOrInt_AsULong(value);
    if (u == (unsigned long)-1 && PyErr_Occurred())
        return false;
    *((unsigned long*)address) = u;
    return true;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::UIntConverter::FromMemory(void* address)
{
// construct python object from C++ unsigned int read at <address>
    return PyLong_FromUnsignedLong(*((UInt_t*)address));
}

bool CPyCppyy::UIntConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ unsigned int, write it at <address>
    ULong_t u = PyLongOrInt_AsULong(value);
    if (u == (unsigned long)-1 && PyErr_Occurred())
        return false;

    if (u > (ULong_t)UINT_MAX) {
        PyErr_SetString(PyExc_OverflowError, "value too large for unsigned int");
        return false;
    }

    *((UInt_t*)address) = (UInt_t)u;
    return true;
}

//- floating point converters ------------------------------------------------
CPPYY_IMPL_BASIC_CONVERTER(
    Float,  float,  double, c_float, PyFloat_FromDouble, PyFloat_AsDouble, 'f')
CPPYY_IMPL_BASIC_CONVERTER(
    Double, double, double, c_double, PyFloat_FromDouble, PyFloat_AsDouble, 'd')

CPPYY_IMPL_BASIC_CONVERTER(
    LDouble, LongDouble_t, LongDouble_t, c_longdouble, PyFloat_FromDouble, PyFloat_AsDouble, 'g')

CPyCppyy::ComplexDConverter::ComplexDConverter(bool keepControl) :
    InstanceConverter(Cppyy::GetScope("std::complex<double>"), keepControl) {}

// special case for std::complex<double>, maps it to/from Python's complex
bool CPyCppyy::ComplexDConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
    const Py_complex& pc = PyComplex_AsCComplex(pyobject);
    if (pc.real != -1.0 || !PyErr_Occurred()) {
        fBuffer.real(pc.real);
        fBuffer.imag(pc.imag);
        para.fValue.fVoidp = &fBuffer;
        para.fTypeCode = 'V';
        return true;
    }

    return this->InstanceConverter::SetArg(pyobject, para, ctxt);
}                                                                            \
                                                                             \
PyObject* CPyCppyy::ComplexDConverter::FromMemory(void* address)
{
    std::complex<double>* dc = (std::complex<double>*)address;
    return PyComplex_FromDoubles(dc->real(), dc->imag());
}

bool CPyCppyy::ComplexDConverter::ToMemory(PyObject* value, void* address)
{
    const Py_complex& pc = PyComplex_AsCComplex(value);
    if (pc.real != -1.0 || !PyErr_Occurred()) {
         std::complex<double>* dc = (std::complex<double>*)address;
         dc->real(pc.real);
         dc->imag(pc.imag);
         return true;
    }
    return this->InstanceConverter::ToMemory(value, address);
}

//----------------------------------------------------------------------------
bool CPyCppyy::DoubleRefConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ double&, set arg for call
    if (RefFloat_CheckExact(pyobject)) {
        para.fValue.fVoidp = (void*)&((PyFloatObject*)pyobject)->ob_fval;
        para.fTypeCode = 'V';
        return true;
    }

// alternate, pass pointer from buffer
    Py_ssize_t buflen = Utility::GetBuffer(pyobject, 'd', sizeof(double), para.fValue.fVoidp);
    if (para.fValue.fVoidp && buflen) {
        para.fTypeCode = 'V';
        return true;
    }

#if PY_VERSION_HEX < 0x02050000
    PyErr_SetString(PyExc_TypeError, "use cppyy.Double for pass-by-ref of doubles");
#else
    PyErr_SetString(PyExc_TypeError, "use ctypes.c_double for pass-by-ref of doubles");
#endif
    return false;
}

//----------------------------------------------------------------------------
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Float,   float,        c_float,      PyFloat_AsDouble)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(Double,  double,       c_double,     PyFloat_AsDouble)
CPPYY_IMPL_BASIC_CONST_REFCONVERTER(LDouble, LongDouble_t, c_longdouble, PyFloat_AsDouble)

//----------------------------------------------------------------------------
bool CPyCppyy::VoidConverter::SetArg(PyObject*, Parameter&, CallContext*)
{
// can't happen (unless a type is mapped wrongly), but implemented for completeness
    PyErr_SetString(PyExc_SystemError, "void/unknown arguments can\'t be set");
    return false;
}

//----------------------------------------------------------------------------
bool CPyCppyy::LLongConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ long long, set arg for call
    if (PyFloat_Check(pyobject)) {
    // special case: float implements nb_int, but allowing rounding conversions
    // interferes with overloading
        PyErr_SetString(PyExc_ValueError, "cannot convert float to long long");
        return false;
    }

    para.fValue.fLLong = PyLong_AsLongLong(pyobject);
    if (PyErr_Occurred())
        return false;
    para.fTypeCode = 'q';
    return true;
}

PyObject* CPyCppyy::LLongConverter::FromMemory(void* address)
{
// construct python object from C++ long long read at <address>
    return PyLong_FromLongLong(*(Long64_t*)address);
}

bool CPyCppyy::LLongConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ long long, write it at <address>
    Long64_t ll = PyLong_AsLongLong(value);
    if (ll == -1 && PyErr_Occurred())
        return false;
    *((Long64_t*)address) = ll;
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::ULLongConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ unsigned long long, set arg for call
    para.fValue.fULLong = PyLongOrInt_AsULong64(pyobject);
    if (PyErr_Occurred())
        return false;
    para.fTypeCode = 'Q';
    return true;
}

PyObject* CPyCppyy::ULLongConverter::FromMemory(void* address)
{
// construct python object from C++ unsigned long long read at <address>
    return PyLong_FromUnsignedLongLong(*(ULong64_t*)address);
}

bool CPyCppyy::ULLongConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ unsigned long long, write it at <address>
    Long64_t ull = PyLongOrInt_AsULong64(value);
    if (PyErr_Occurred())
        return false;
    *((ULong64_t*)address) = ull;
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CStringConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// construct a new string and copy it in new memory
    Py_ssize_t len;
    const char* cstr = CPyCppyy_PyText_AsStringAndSize(pyobject, &len);
    if (!cstr) {
    // special case: allow ctypes c_char_p
        PyObject* pytype = 0, *pyvalue = 0, *pytrace = 0;
        PyErr_Fetch(&pytype, &pyvalue, &pytrace);
        if (Py_TYPE(pyobject) == GetCTypesType(ct_c_char_p)) {
            para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;
            para.fTypeCode = 'V';
            Py_XDECREF(pytype); Py_XDECREF(pyvalue); Py_XDECREF(pytrace);
            return true;
        }
        PyErr_Restore(pytype, pyvalue, pytrace);
        return false;
    }

    fBuffer = std::string(cstr, len);

// verify (too long string will cause truncation, no crash)
    if (fMaxSize != -1 && fMaxSize < (long)fBuffer.size())
        PyErr_Warn(PyExc_RuntimeWarning, (char*)"string too long for char array (truncated)");
    else if (fMaxSize != -1)
        fBuffer.resize(fMaxSize, '\0');      // padd remainder of buffer as needed

// set the value and declare success
    para.fValue.fVoidp = (void*)fBuffer.c_str();
    para.fTypeCode = 'p';
    return true;
}

PyObject* CPyCppyy::CStringConverter::FromMemory(void* address)
{
// construct python object from C++ const char* read at <address>
    if (address && *(char**)address) {
        if (fMaxSize != -1) {      // need to prevent reading beyond boundary
            std::string buf(*(char**)address, fMaxSize);    // cut on fMaxSize
            return CPyCppyy_PyText_FromString(buf.c_str());   // cut on \0
        }

        return CPyCppyy_PyText_FromString(*(char**)address);
    }

// empty string in case there's no address
    Py_INCREF(PyStrings::gEmptyString);
    return PyStrings::gEmptyString;
}

bool CPyCppyy::CStringConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ const char*, write it at <address>
    Py_ssize_t len;
    const char* cstr = CPyCppyy_PyText_AsStringAndSize(value, &len);
    if (!cstr) return false;

// verify (too long string will cause truncation, no crash)
    if (fMaxSize != -1 && fMaxSize < (long)len)
        PyErr_Warn(PyExc_RuntimeWarning, (char*)"string too long for char array (truncated)");

    if (fMaxSize != -1)
        strncpy(*(char**)address, cstr, fMaxSize);    // padds remainder
    else
    // coverity[secure_coding] - can't help it, it's intentional.
        strcpy(*(char**)address, cstr);

   return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::WCStringConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// construct a new string and copy it in new memory
    Py_ssize_t len = PyUnicode_GetSize(pyobject);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

    fBuffer = (wchar_t*)realloc(fBuffer, sizeof(wchar_t)*(len+1));
    Py_ssize_t res = CPyCppyy_PyUnicode_AsWideChar(pyobject, fBuffer, len);
    if (res == -1)
        return false;   // could free the buffer here

// set the value and declare success
    fBuffer[len] = L'\0';
    para.fValue.fVoidp = (void*)fBuffer;
    para.fTypeCode = 'p';
    return true;
}

PyObject* CPyCppyy::WCStringConverter::FromMemory(void* address)
{
// construct python object from C++ wchar_t* read at <address>
    if (address && *(wchar_t**)address) {
        if (fMaxSize != -1)        // need to prevent reading beyond boundary
            return PyUnicode_FromWideChar(*(wchar_t**)address, fMaxSize);
    // with unknown size
        return PyUnicode_FromWideChar(*(wchar_t**)address, wcslen(*(wchar_t**)address));
    }

// empty string in case there's no valid address
    wchar_t w = L'\0';
    return PyUnicode_FromWideChar(&w, 0);
}

bool CPyCppyy::WCStringConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ wchar_t*, write it at <address>
    Py_ssize_t len = PyUnicode_GetSize(value);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

// verify (too long string will cause truncation, no crash)
    if (fMaxSize != -1 && fMaxSize < len)
        PyErr_Warn(PyExc_RuntimeWarning, (char*)"string too long for wchar_t array (truncated)");

    Py_ssize_t res = -1;
    if (fMaxSize != -1)
        res = CPyCppyy_PyUnicode_AsWideChar(value, *(wchar_t**)address, fMaxSize);
    else
    // coverity[secure_coding] - can't help it, it's intentional.
        res = CPyCppyy_PyUnicode_AsWideChar(value, *(wchar_t**)address, len);

    if (res == -1) return false;
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CString16Converter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// construct a new string and copy it in new memory
    Py_ssize_t len = PyUnicode_GetSize(pyobject);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

    PyObject* bstr = PyUnicode_AsUTF16String(pyobject);
    if (!bstr) return false;

    fBuffer = (char16_t*)realloc(fBuffer, sizeof(char16_t)*(len+1));
    memcpy(fBuffer, PyBytes_AS_STRING(bstr) + sizeof(char16_t) /*BOM*/, len*sizeof(char16_t));
    Py_DECREF(bstr);

// set the value and declare success
    fBuffer[len] = u'\0';
    para.fValue.fVoidp = (void*)fBuffer;
    para.fTypeCode = 'p';
    return true;
}

PyObject* CPyCppyy::CString16Converter::FromMemory(void* address)
{
// construct python object from C++ char16_t* read at <address>
    if (address && *(char16_t**)address) {
        if (fMaxSize != -1)        // need to prevent reading beyond boundary
            return PyUnicode_DecodeUTF16(*(const char**)address, fMaxSize, nullptr, nullptr);
    // with unknown size
        return PyUnicode_DecodeUTF16(*(const char**)address,
            std::char_traits<char16_t>::length(*(char16_t**)address)*sizeof(char16_t), nullptr, nullptr);
    }

// empty string in case there's no valid address
    char16_t w = u'\0';
    return PyUnicode_DecodeUTF16((const char*)&w, 0, nullptr, nullptr);
}

bool CPyCppyy::CString16Converter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ char16_t*, write it at <address>
    Py_ssize_t len = PyUnicode_GetSize(value);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

// verify (too long string will cause truncation, no crash)
    if (fMaxSize != -1 && fMaxSize < len) {
        PyErr_Warn(PyExc_RuntimeWarning, (char*)"string too long for char16_t array (truncated)");
        len = fMaxSize-1;
    }

    PyObject* bstr = PyUnicode_AsUTF16String(value);
    if (!bstr) return false;

    memcpy(*((void**)address), PyBytes_AS_STRING(bstr) + sizeof(char16_t) /*BOM*/, len*sizeof(char16_t));
    Py_DECREF(bstr);
    *((char16_t**)address)[len] = u'\0';
    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CString32Converter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// construct a new string and copy it in new memory
    Py_ssize_t len = PyUnicode_GetSize(pyobject);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

    PyObject* bstr = PyUnicode_AsUTF32String(pyobject);
    if (!bstr) return false;

    fBuffer = (char32_t*)realloc(fBuffer, sizeof(char32_t)*(len+1));
    memcpy(fBuffer, PyBytes_AS_STRING(bstr) + sizeof(char32_t) /*BOM*/, len*sizeof(char32_t));
    Py_DECREF(bstr);

// set the value and declare success
    fBuffer[len] = U'\0';
    para.fValue.fVoidp = (void*)fBuffer;
    para.fTypeCode = 'p';
    return true;
}

PyObject* CPyCppyy::CString32Converter::FromMemory(void* address)
{
// construct python object from C++ char32_t* read at <address>
    if (address && *(char32_t**)address) {
        if (fMaxSize != -1)        // need to prevent reading beyond boundary
            return PyUnicode_DecodeUTF32(*(const char**)address, fMaxSize, nullptr, nullptr);
    // with unknown size
        return PyUnicode_DecodeUTF32(*(const char**)address,
            std::char_traits<char32_t>::length(*(char32_t**)address)*sizeof(char32_t), nullptr, nullptr);
    }

// empty string in case there's no valid address
    char32_t w = U'\0';
    return PyUnicode_DecodeUTF32((const char*)&w, 0, nullptr, nullptr);
}

bool CPyCppyy::CString32Converter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ char32_t*, write it at <address>
    Py_ssize_t len = PyUnicode_GetSize(value);
    if (len == (Py_ssize_t)-1 && PyErr_Occurred())
        return false;

// verify (too long string will cause truncation, no crash)
    if (fMaxSize != -1 && fMaxSize < len) {
        PyErr_Warn(PyExc_RuntimeWarning, (char*)"string too long for char32_t array (truncated)");
        len = fMaxSize-1;
    }

    PyObject* bstr = PyUnicode_AsUTF32String(value);
    if (!bstr) return false;

    memcpy(*((void**)address), PyBytes_AS_STRING(bstr) + sizeof(char32_t) /*BOM*/, len*sizeof(char32_t));
    Py_DECREF(bstr);
    *((char32_t**)address)[len] = U'\0';
    return true;
}


//----------------------------------------------------------------------------
bool CPyCppyy::NonConstCStringConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// attempt base class first (i.e. passing a string), but if that fails, try a buffer
    if (this->CStringConverter::SetArg(pyobject, para, ctxt))
        return true;

// apparently failed, try char buffer
    PyErr_Clear();
    return CArraySetArg(pyobject, para, 'c', sizeof(char));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::NonConstCStringConverter::FromMemory(void* address)
{
// assume this is a buffer access if the size is known; otherwise assume string
    if (fMaxSize != -1)
        return CPyCppyy_PyText_FromStringAndSize(*(char**)address, fMaxSize);
    return this->CStringConverter::FromMemory(address);
}

//----------------------------------------------------------------------------
bool CPyCppyy::VoidArrayConverter::GetAddressSpecialCase(PyObject* pyobject, void*& address)
{
// (1): C++11 style "null pointer"
    if (pyobject == gNullPtrObject) {
        address = nullptr;
        return true;
    }

// (2): allow integer zero to act as a null pointer (C NULL), no deriveds
    if (PyInt_CheckExact(pyobject) || PyLong_CheckExact(pyobject)) {
        intptr_t val = (intptr_t)PyLong_AsLongLong(pyobject);
        if (val == 0l) {
            address = (void*)val;
            return true;
        }

        return false;
    }

// (3): opaque PyCapsule (CObject in older pythons) from somewhere
    if (CPyCppyy_PyCapsule_CheckExact(pyobject)) {
        address = (void*)CPyCppyy_PyCapsule_GetPointer(pyobject, nullptr);
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
bool CPyCppyy::VoidArrayConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// just convert pointer if it is a C++ object
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (pyobj) {
    // depending on memory policy, some objects are no longer owned when passed to C++
        if (!fKeepControl && !UseStrictOwnership(ctxt))
            pyobj->CppOwns();

   // set pointer (may be null) and declare success
        para.fValue.fVoidp = pyobj->GetObject();
        para.fTypeCode = 'p';
        return true;
    }

// handle special cases
    if (GetAddressSpecialCase(pyobject, para.fValue.fVoidp)) {
        para.fTypeCode = 'p';
        return true;
    }

// allow ctypes voidp (which if got as a buffer will return void**, not void*); use
// isintance instead of an exact check, b/c c_void_p is the type mapper for typedefs
// of void* (typically opaque handles)
    if (PyObject_IsInstance(pyobject, (PyObject*)GetCTypesType(ct_c_void_p))) {
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;
        para.fTypeCode = 'V';
        return true;
    }

// allow any other ctypes pointer type
    if (IsCTypesArrayOrPointer(pyobject)) {
        void** payload = (void**)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;
        if (payload) {
            para.fValue.fVoidp = *payload;
            para.fTypeCode = 'p';
            return true;
        }
    }

// final try: attempt to get buffer
    Py_ssize_t buflen = Utility::GetBuffer(pyobject, '*', 1, para.fValue.fVoidp, false);

// ok if buffer exists (can't perform any useful size checks)
    if (para.fValue.fVoidp && buflen != 0) {
        para.fTypeCode = 'p';
        return true;
    }

// give up
    return false;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::VoidArrayConverter::FromMemory(void* address)
{
// nothing sensible can be done, just return <address> as pylong
    if (!address || *(ptrdiff_t*)address == 0) {
        Py_INCREF(gNullPtrObject);
        return gNullPtrObject;
    }
    return CreatePointerView(*(ptrdiff_t**)address);
}

//----------------------------------------------------------------------------
bool CPyCppyy::VoidArrayConverter::ToMemory(PyObject* value, void* address)
{
// just convert pointer if it is a C++ object
    CPPInstance* pyobj = GetCppInstance(value);
    if (pyobj) {
    // depending on memory policy, some objects are no longer owned when passed to C++
        if (!fKeepControl && CallContext::sMemoryPolicy != CallContext::kUseStrict)
            pyobj->CppOwns();

    // set pointer (may be null) and declare success
        *(void**)address = pyobj->GetObject();
        return true;
    }

// handle special cases
    void* ptr = nullptr;
    if (GetAddressSpecialCase(value, ptr)) {
        *(void**)address = ptr;
        return true;
    }

// final try: attempt to get buffer
    void* buf = nullptr;
    Py_ssize_t buflen = Utility::GetBuffer(value, '*', 1, buf, false);
    if (!buf || buflen == 0)
        return false;

    *(void**)address = buf;
    return true;
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_ARRAY_CONVERTER(name, ctype, type, code)                  \
CPyCppyy::name##ArrayConverter::name##ArrayConverter(dims_t dims) {          \
    int nalloc = (dims && 0 < dims[0]) ? (int)dims[0]+1: 2;                  \
    fShape = new Py_ssize_t[nalloc];                                         \
    if (dims) {                                                              \
        for (int i = 0; i < nalloc; ++i) fShape[i] = (Py_ssize_t)dims[i];    \
    } else {                                                                 \
        fShape[0] = 1; fShape[1] = -1;                                       \
    }                                                                        \
}                                                                            \
                                                                             \
bool CPyCppyy::name##ArrayConverter::SetArg(                                 \
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)            \
{                                                                            \
    /* filter ctypes first b/c their buffer conversion will be wrong */      \
    PyTypeObject* ctypes_type = GetCTypesType(ct_##ctype);                   \
    if (Py_TYPE(pyobject) == ctypes_type) {                                  \
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;\
        para.fTypeCode = 'p';                                                \
        return true;                                                         \
    } else if (Py_TYPE(pyobject) == GetCTypesPtrType(ct_##ctype)) {          \
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;\
        para.fTypeCode = 'V';                                                \
        return true;                                                         \
    } else if (IsPyCArgObject(pyobject)) {                                   \
        CPyCppyy_tagPyCArgObject* carg = (CPyCppyy_tagPyCArgObject*)pyobject;\
        if (carg->obj && Py_TYPE(carg->obj) == ctypes_type) {                \
            para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)carg->obj)->b_ptr;\
            para.fTypeCode = 'p';                                            \
            return true;                                                     \
        }                                                                    \
    }                                                                        \
    return CArraySetArg(pyobject, para, code, sizeof(type));                 \
}                                                                            \
                                                                             \
PyObject* CPyCppyy::name##ArrayConverter::FromMemory(void* address)          \
{                                                                            \
    if (fShape[1] == UNKNOWN_SIZE)                                           \
        return CreateLowLevelView((type**)address, fShape);                  \
    return CreateLowLevelView(*(type**)address, fShape);                     \
}                                                                            \
                                                                             \
bool CPyCppyy::name##ArrayConverter::ToMemory(PyObject* value, void* address)\
{                                                                            \
    if (fShape[0] != 1) {                                                    \
        PyErr_SetString(PyExc_ValueError, "only 1-dim arrays supported");    \
        return false;                                                        \
    }                                                                        \
    void* buf = nullptr;                                                     \
    Py_ssize_t buflen = Utility::GetBuffer(value, code, sizeof(type), buf);  \
    if (buflen == 0)                                                         \
        return false;                                                        \
    if (0 <= fShape[1]) {                                                    \
        if (fShape[1] < buflen) {                                            \
            PyErr_SetString(PyExc_ValueError, "buffer too large for value"); \
            return false;                                                    \
        }                                                                    \
        memcpy(*(type**)address, buf, (0 < buflen ? buflen : 1)*sizeof(type));\
    } else                                                                   \
        *(type**)address = (type*)buf;                                       \
    return true;                                                             \
}                                                                            \
                                                                             \
bool CPyCppyy::name##ArrayPtrConverter::SetArg(                              \
    PyObject* pyobject, Parameter& para, CallContext* ctxt )                 \
{                                                                            \
    if (Py_TYPE(pyobject) == GetCTypesPtrType(ct_##ctype)) {                 \
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;\
        para.fTypeCode = 'p';                                                \
        return true;                                                         \
    } else if (Py_TYPE(pyobject) == GetCTypesType(ct_c_void_p)) {            \
    /* special case: pass address of c_void_p buffer to return the address */\
        para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)pyobject)->b_ptr;\
        para.fTypeCode = 'p';                                                \
        return true;                                                         \
    }                                                                        \
    bool res = name##ArrayConverter::SetArg(pyobject, para, ctxt);           \
    if (res && para.fTypeCode == 'p') {                                      \
        para.fRef = para.fValue.fVoidp;                                      \
        para.fValue.fVoidp = &para.fRef;                                     \
        return true;                                                         \
    }                                                                        \
    return false;                                                            \
}


//----------------------------------------------------------------------------
CPPYY_IMPL_ARRAY_CONVERTER(Bool,     c_bool,       bool,                 '?')
CPPYY_IMPL_ARRAY_CONVERTER(SChar,    c_char,       signed char,          'b')
CPPYY_IMPL_ARRAY_CONVERTER(UChar,    c_ubyte,      unsigned char,        'B')
#if __cplusplus > 201402L
CPPYY_IMPL_ARRAY_CONVERTER(Byte,     c_ubyte,      std::byte,            'B')
#endif
CPPYY_IMPL_ARRAY_CONVERTER(Short,    c_short,      short,                'h')
CPPYY_IMPL_ARRAY_CONVERTER(UShort,   c_ushort,     unsigned short,       'H')
CPPYY_IMPL_ARRAY_CONVERTER(Int,      c_int,        int,                  'i')
CPPYY_IMPL_ARRAY_CONVERTER(UInt,     c_uint,       unsigned int,         'I')
CPPYY_IMPL_ARRAY_CONVERTER(Long,     c_long,       long,                 'l')
CPPYY_IMPL_ARRAY_CONVERTER(ULong,    c_ulong,      unsigned long,        'L')
CPPYY_IMPL_ARRAY_CONVERTER(LLong,    c_longlong,   long long,            'q')
CPPYY_IMPL_ARRAY_CONVERTER(ULLong,   c_ulonglong,  unsigned long long,   'Q')
CPPYY_IMPL_ARRAY_CONVERTER(Float,    c_float,      float,                'f')
CPPYY_IMPL_ARRAY_CONVERTER(Double,   c_double,     double,               'd')
CPPYY_IMPL_ARRAY_CONVERTER(LDouble,  c_longdouble, long double,          'D')
CPPYY_IMPL_ARRAY_CONVERTER(ComplexD, c_complex,    std::complex<double>, 'Z')


//----------------------------------------------------------------------------
PyObject* CPyCppyy::CStringArrayConverter::FromMemory(void* address)
{
    if (fShape[1] == UNKNOWN_SIZE)
        return CreateLowLevelView((const char**)address, fShape);
    return CreateLowLevelView(*(const char***)address, fShape);
}


//- converters for special cases ---------------------------------------------
bool CPyCppyy::NullptrConverter::SetArg(PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// Only allow C++11 style nullptr to pass
    if (pyobject == gNullPtrObject) {
        para.fValue.fVoidp = nullptr;
        para.fTypeCode = 'p';
        return true;
    }
    return false;
}


//----------------------------------------------------------------------------
#define CPPYY_IMPL_STRING_AS_PRIMITIVE_CONVERTER(name, type, F1, F2)         \
CPyCppyy::name##Converter::name##Converter(bool keepControl) :               \
    InstanceConverter(Cppyy::GetScope(#type), keepControl) {}                \
                                                                             \
bool CPyCppyy::name##Converter::SetArg(                                      \
    PyObject* pyobject, Parameter& para, CallContext* ctxt)                  \
{                                                                            \
    Py_ssize_t len;                                                          \
    const char* cstr = CPyCppyy_PyText_AsStringAndSize(pyobject, &len);      \
    if (cstr) {                                                              \
        fBuffer = type(cstr, len);                                           \
        para.fValue.fVoidp = &fBuffer;                                       \
        para.fTypeCode = 'V';                                                \
        return true;                                                         \
    }                                                                        \
                                                                             \
    PyErr_Clear();                                                           \
    if (!(PyInt_Check(pyobject) || PyLong_Check(pyobject))) {                \
        bool result = InstanceConverter::SetArg(pyobject, para, ctxt);       \
        para.fTypeCode = 'V';                                                \
        return result;                                                       \
    }                                                                        \
    return false;                                                            \
}                                                                            \
                                                                             \
PyObject* CPyCppyy::name##Converter::FromMemory(void* address)               \
{                                                                            \
    if (address)                                                             \
        return CPyCppyy_PyText_FromStringAndSize(((type*)address)->F1(), ((type*)address)->F2()); \
    Py_INCREF(PyStrings::gEmptyString);                                      \
    return PyStrings::gEmptyString;                                          \
}                                                                            \
                                                                             \
bool CPyCppyy::name##Converter::ToMemory(PyObject* value, void* address)     \
{                                                                            \
    if (CPyCppyy_PyText_Check(value)) {                                      \
        *((type*)address) = CPyCppyy_PyText_AsString(value);                 \
        return true;                                                         \
    }                                                                        \
                                                                             \
    return InstanceConverter::ToMemory(value, address);                      \
}

CPPYY_IMPL_STRING_AS_PRIMITIVE_CONVERTER(TString, TString, Data, Length)
CPPYY_IMPL_STRING_AS_PRIMITIVE_CONVERTER(STLString, std::string, c_str, size)
CPPYY_IMPL_STRING_AS_PRIMITIVE_CONVERTER(STLStringViewBase, std::string_view, data, size)
bool CPyCppyy::STLStringViewConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
    if (this->STLStringViewBaseConverter::SetArg(pyobject, para, ctxt))
        return true;

    if (!CPPInstance_Check(pyobject))
        return false;

    static Cppyy::TCppScope_t sStringID = Cppyy::GetScope("std::string");
    CPPInstance* pyobj = (CPPInstance*)pyobject;
    if (pyobj->ObjectIsA() == sStringID) {
        void* ptr = pyobj->GetObject();
        if (!ptr)
            return false;

        fBuffer = *((std::string*)ptr);
        para.fValue.fVoidp = &fBuffer;
        para.fTypeCode = 'V';
        return true;
    }

    return false;
}

CPyCppyy::STLWStringConverter::STLWStringConverter(bool keepControl) :
    InstanceConverter(Cppyy::GetScope("std::wstring"), keepControl) {}

bool CPyCppyy::STLWStringConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
    if (PyUnicode_Check(pyobject)) {
        Py_ssize_t len = PyUnicode_GET_SIZE(pyobject);
        fBuffer.resize(len);
        CPyCppyy_PyUnicode_AsWideChar(pyobject, &fBuffer[0], len);
        para.fValue.fVoidp = &fBuffer;
        para.fTypeCode = 'V';
        return true;
    }

    if (!(PyInt_Check(pyobject) || PyLong_Check(pyobject))) {
        bool result = InstancePtrConverter::SetArg(pyobject, para, ctxt);
        para.fTypeCode = 'V';
        return result;
    }

    return false;
}

PyObject* CPyCppyy::STLWStringConverter::FromMemory(void* address)
{
    if (address)
        return PyUnicode_FromWideChar(((std::wstring*)address)->c_str(), ((std::wstring*)address)->size());
    wchar_t w = L'\0';
    return PyUnicode_FromWideChar(&w, 0);
}

bool CPyCppyy::STLWStringConverter::ToMemory(PyObject* value, void* address)
{
    if (PyUnicode_Check(value)) {
        Py_ssize_t len = PyUnicode_GET_SIZE(value);
        wchar_t* buf = new wchar_t[len+1];
        CPyCppyy_PyUnicode_AsWideChar(value, buf, len);
        *((std::wstring*)address) = std::wstring(buf, len);
        delete[] buf;
        return true;
    }
    return InstanceConverter::ToMemory(value, address);
}


bool CPyCppyy::STLStringMoveConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ std::string&&, set arg for call
    int moveit_reason = 3;   // move on temporary fBuffer
    if (CPPInstance_Check(pyobject)) {
        CPPInstance* pyobj = (CPPInstance*)pyobject;
        if (pyobj->fFlags & CPPInstance::kIsRValue) {
            pyobj->fFlags &= ~CPPInstance::kIsRValue;
            moveit_reason = 2;
        } else if (pyobject->ob_refcnt == MOVE_REFCOUNT_CUTOFF) {
            moveit_reason = 1;
        } else
            moveit_reason = 0;
    }

    if (moveit_reason) {
        bool result = this->STLStringConverter::SetArg(pyobject, para, ctxt);
        if (!result && moveit_reason == 2)       // restore the movability flag?
            ((CPPInstance*)pyobject)->fFlags |= CPPInstance::kIsRValue;
        return result;
    }

    PyErr_SetString(PyExc_ValueError, "object is not an rvalue");
    return false;      // not a temporary or movable object
}


//----------------------------------------------------------------------------
bool CPyCppyy::InstancePtrConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ instance*, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (!pyobj) {
        if (GetAddressSpecialCase(pyobject, para.fValue.fVoidp)) {
            para.fTypeCode = 'p';      // allow special cases such as nullptr
            return true;
        }

   // not a cppyy object (TODO: handle SWIG etc.)
        return false;
    }

    if (pyobj->ObjectIsA() && Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
    // depending on memory policy, some objects need releasing when passed into functions
        if (!KeepControl() && !UseStrictOwnership(ctxt))
            pyobj->CppOwns();

    // calculate offset between formal and actual arguments
        para.fValue.fVoidp = pyobj->GetObject();
        if (pyobj->ObjectIsA() != fClass) {
            para.fValue.fIntPtr += Cppyy::GetBaseOffset(
                pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */);
        }

   // set pointer (may be null) and declare success
        para.fTypeCode = 'p';
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstancePtrConverter::FromMemory(void* address)
{
// construct python object from C++ instance read at <address>
    return BindCppObject(address, fClass, CPPInstance::kIsReference);
}

//----------------------------------------------------------------------------
bool CPyCppyy::InstancePtrConverter::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ instance, write it at <address>
    CPPInstance* pyobj = GetCppInstance(value);
    if (!pyobj) {
        void* ptr = nullptr;
        if (GetAddressSpecialCase(value, ptr)) {
            *(void**)address = ptr;          // allow special cases such as nullptr
            return true;
        }

    // not a cppyy object (TODO: handle SWIG etc.)
        return false;
    }

    if (Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
    // depending on memory policy, some objects need releasing when passed into functions
        if (!KeepControl() && CallContext::sMemoryPolicy != CallContext::kUseStrict)
            ((CPPInstance*)value)->CppOwns();

        *(void**)address = pyobj->GetObject();
        return true;
    }

    return false;
}

// TODO: CONSOLIDATE Instance, InstanceRef, InstancePtr ...

//----------------------------------------------------------------------------
bool CPyCppyy::InstanceConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ instance, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (pyobj) {
        if (pyobj->ObjectIsA() && Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
        // calculate offset between formal and actual arguments
            para.fValue.fVoidp = pyobj->GetObject();
            if (!para.fValue.fVoidp)
                return false;

            if (pyobj->ObjectIsA() != fClass) {
                para.fValue.fIntPtr += Cppyy::GetBaseOffset(
                    pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */);
            }

            para.fTypeCode = 'V';
            return true;
        }
    }

    return ConvertImplicit(fClass, pyobject, para, ctxt);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceConverter::FromMemory(void* address)
{
    return BindCppObjectNoCast((Cppyy::TCppObject_t)address, fClass);
}

//----------------------------------------------------------------------------
bool CPyCppyy::InstanceConverter::ToMemory(PyObject* value, void* address)
{
// assign value to C++ instance living at <address> through assignment operator
    PyObject* pyobj = BindCppObjectNoCast(address, fClass);
    PyObject* result = PyObject_CallMethod(pyobj, (char*)"__assign__", (char*)"O", value);
    Py_DECREF(pyobj);

    if (result) {
        Py_DECREF(result);
        return true;
    }
    return false;
}


//----------------------------------------------------------------------------
bool CPyCppyy::InstanceRefConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ instance&, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (pyobj) {

    // reject moves
        if (pyobj->fFlags & CPPInstance::kIsRValue)
            return false;

        if (pyobj->ObjectIsA() && Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
        // calculate offset between formal and actual arguments
            para.fValue.fVoidp = pyobj->GetObject();
            if (pyobj->ObjectIsA() != fClass) {
                para.fValue.fIntPtr += Cppyy::GetBaseOffset(
                    pyobj->ObjectIsA(), fClass, para.fValue.fVoidp, 1 /* up-cast */);
            }

            para.fTypeCode = 'V';
            return true;
        }
    }

    if (!fIsConst)      // no implicit conversion possible
        return false;

    return ConvertImplicit(fClass, pyobject, para, ctxt);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceRefConverter::FromMemory(void* address)
{
    return BindCppObjectNoCast((Cppyy::TCppObject_t)address, fClass);
}

//----------------------------------------------------------------------------
bool CPyCppyy::InstanceMoveConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ instance&&, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (!pyobj) {
    // implicit conversion is fine as it the temporary by definition is moveable
        return ConvertImplicit(fClass, pyobject, para, ctxt);
    }

// moving is same as by-ref, but have to check that move is allowed
    int moveit_reason = 0;
    if (pyobj->fFlags & CPPInstance::kIsRValue) {
        pyobj->fFlags &= ~CPPInstance::kIsRValue;
        moveit_reason = 2;
    } else if (pyobject->ob_refcnt == MOVE_REFCOUNT_CUTOFF) {
        moveit_reason = 1;
    }

    if (moveit_reason) {
        bool result = this->InstanceRefConverter::SetArg(pyobject, para, ctxt);
        if (!result && moveit_reason == 2)       // restore the movability flag?
            ((CPPInstance*)pyobject)->fFlags |= CPPInstance::kIsRValue;
        return result;
    }

    PyErr_SetString(PyExc_ValueError, "object is not an rvalue");
    return false;      // not a temporary or movable object
}

//----------------------------------------------------------------------------
template <bool ISREFERENCE>
bool CPyCppyy::InstancePtrPtrConverter<ISREFERENCE>::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// convert <pyobject> to C++ instance**, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (!pyobj)
        return false;              // not a cppyy object (TODO: handle SWIG etc.)

    if (Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
    // depending on memory policy, some objects need releasing when passed into functions
        if (!KeepControl() && !UseStrictOwnership(ctxt))
            pyobj->CppOwns();

    // set pointer (may be null) and declare success
        if (pyobj->fFlags & CPPInstance::kIsReference) // already a ptr to object?
            para.fValue.fVoidp = pyobj->GetObjectRaw();
        else
            para.fValue.fVoidp = &pyobj->GetObjectRaw();
        para.fTypeCode = ISREFERENCE ? 'V' : 'p';
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
template <bool ISREFERENCE>
PyObject* CPyCppyy::InstancePtrPtrConverter<ISREFERENCE>::FromMemory(void* address)
{
// construct python object from C++ instance* read at <address>
    return BindCppObject(address, fClass, CPPInstance::kIsReference);
}

//----------------------------------------------------------------------------
template <bool ISREFERENCE>
bool CPyCppyy::InstancePtrPtrConverter<ISREFERENCE>::ToMemory(PyObject* value, void* address)
{
// convert <value> to C++ instance*, write it at <address>
    CPPInstance* pyobj = GetCppInstance(value);
    if (!pyobj)
        return false;              // not a cppyy object (TODO: handle SWIG etc.)

    if (Cppyy::IsSubtype(pyobj->ObjectIsA(), fClass)) {
    // depending on memory policy, some objects need releasing when passed into functions
        if (!KeepControl() && CallContext::sMemoryPolicy != CallContext::kUseStrict)
            pyobj->CppOwns();

    // register the value for potential recycling
        MemoryRegulator::RegisterPyObject(pyobj, pyobj->GetObject());

    // set pointer (may be null) and declare success
        *(void**)address = pyobj->GetObject();
        return true;
    }

    return false;
}


namespace CPyCppyy {
// Instantiate the templates
    template class CPyCppyy::InstancePtrPtrConverter<true>;
    template class CPyCppyy::InstancePtrPtrConverter<false>;
}

//----------------------------------------------------------------------------
bool CPyCppyy::InstanceArrayConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* txt */)
{
// convert <pyobject> to C++ instance**, set arg for call
    if (!TupleOfInstances_CheckExact(pyobject))
        return false;              // no guarantee that the tuple is okay

// treat the first instance of the tuple as the start of the array, and pass it
// by pointer (TODO: store and check sizes)
    if (PyTuple_Size(pyobject) < 1)
        return false;

    PyObject* first = PyTuple_GetItem(pyobject, 0);
    if (!CPPInstance_Check(first))
        return false;              // should not happen

    if (Cppyy::IsSubtype(((CPPInstance*)first)->ObjectIsA(), fClass)) {
    // no memory policies supported; set pointer (may be null) and declare success
        para.fValue.fVoidp = ((CPPInstance*)first)->GetObject();
        para.fTypeCode = 'p';
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceArrayConverter::FromMemory(void* address)
{
// construct python tuple of instances from C++ array read at <address>
    return BindCppObjectArray(*(char**)address, fClass, m_dims);
}

//----------------------------------------------------------------------------
bool CPyCppyy::InstanceArrayConverter::ToMemory(PyObject* /* value */, void* /* address */)
{
// convert <value> to C++ array of instances, write it at <address>

// TODO: need to have size both for the array and from the input
    PyErr_SetString(PyExc_NotImplementedError,
        "access to C-arrays of objects not yet implemented!");
    return false;
}

//___________________________________________________________________________
// CLING WORKAROUND -- classes for STL iterators are completely undefined in that
// they come in a bazillion different guises, so just do whatever
bool CPyCppyy::STLIteratorConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
    if (!CPPInstance_Check(pyobject))
        return false;

// just set the pointer value, no check
    CPPInstance* pyobj = (CPPInstance*)pyobject;
    para.fValue.fVoidp = pyobj->GetObject();
    para.fTypeCode = 'V';
    return true;
}
// -- END CLING WORKAROUND

//----------------------------------------------------------------------------
bool CPyCppyy::VoidPtrRefConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ void*&, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (pyobj) {
        para.fValue.fVoidp = &pyobj->GetObjectRaw();
        para.fTypeCode = 'V';
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
bool CPyCppyy::VoidPtrPtrConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// convert <pyobject> to C++ void**, set arg for call
    CPPInstance* pyobj = GetCppInstance(pyobject);
    if (pyobj) {
    // this is a C++ object, take and set its address
        para.fValue.fVoidp = &pyobj->GetObjectRaw();
        para.fTypeCode = 'p';
        return true;
    } else if (IsPyCArgObject(pyobject)) {
        CPyCppyy_tagPyCArgObject* carg = (CPyCppyy_tagPyCArgObject*)pyobject;
        if (carg->obj) {
            para.fValue.fVoidp = (void*)((CPyCppyy_tagCDataObject*)carg->obj)->b_ptr;
            para.fTypeCode = 'p';
            return true;
        }
    }

// buffer objects are allowed under "user knows best" (this includes the buffer
// interface to ctypes.c_void_p, which results in a void**)
    Py_ssize_t buflen = Utility::GetBuffer(pyobject, '*', 1, para.fValue.fVoidp, false);

// ok if buffer exists (can't perform any useful size checks)
    if (para.fValue.fVoidp && buflen != 0) {
        para.fTypeCode = 'p';
        return true;
    }

    return false;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::VoidPtrPtrConverter::FromMemory(void* address)
{
// read a void** from address; since this is unknown, long is used (user can cast)
    if (!address || *(ptrdiff_t*)address == 0) {
        Py_INCREF(gNullPtrObject);
        return gNullPtrObject;
    }
    return CreatePointerView(*(ptrdiff_t**)address, fSize);
}

//----------------------------------------------------------------------------
bool CPyCppyy::PyObjectConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /* ctxt */)
{
// by definition: set and declare success
    para.fValue.fVoidp = pyobject;
    para.fTypeCode = 'p';
    return true;
}

PyObject* CPyCppyy::PyObjectConverter::FromMemory(void* address)
{
// construct python object from C++ PyObject* read at <address>
    PyObject* pyobject = *((PyObject**)address);

    if (!pyobject) {
        Py_RETURN_NONE;
    }

    Py_INCREF(pyobject);
    return pyobject;
}

bool CPyCppyy::PyObjectConverter::ToMemory(PyObject* value, void* address)
{
// no conversion needed, write <value> at <address>
    Py_INCREF(value);
    Py_XDECREF(*((PyObject**)address));
    *((PyObject**)address) = value;
    return true;
}


//- function pointer converter -----------------------------------------------
static unsigned int sWrapperCounter = 0;
// cache mapping signature/return type to python callable and corresponding wrapper
typedef std::pair<std::string, std::string> RetSigKey_t;
static std::map<RetSigKey_t, std::vector<void*>> sWrapperFree;
static std::map<RetSigKey_t, std::map<PyObject*, void*>> sWrapperLookup;
static std::map<PyObject*, std::pair<void*, RetSigKey_t>> sWrapperWeakRefs;
static std::map<void*, PyObject**> sWrapperReference;

static PyObject* WrapperCacheEraser(PyObject*, PyObject* pyref)
{
    auto ipos = sWrapperWeakRefs.find(pyref);
    if (ipos != sWrapperWeakRefs.end()) {
    // disable this callback and store for possible re-use
        void* wpraddress = ipos->second.first;
        *sWrapperReference[wpraddress] = nullptr;
        sWrapperFree[ipos->second.second].push_back(wpraddress);
    }

    Py_RETURN_NONE;
}
static PyMethodDef gWrapperCacheEraserMethodDef = {
    const_cast<char*>("interal_WrapperCacheEraser"),
    (PyCFunction)WrapperCacheEraser,
    METH_O, nullptr
};

static void* PyFunction_AsCPointer(PyObject* pyobject,
    const std::string& rettype, const std::string& signature)
{
// Convert a bound C++ function pointer or callable python object to a C-style
// function pointer. The former is direct, the latter involves a JIT-ed wrapper.
    static PyObject* sWrapperCacheEraser = PyCFunction_New(&gWrapperCacheEraserMethodDef, nullptr);

    using namespace CPyCppyy;

    if (CPPOverload_Check(pyobject)) {
        CPPOverload* ol = (CPPOverload*)pyobject;
        if (!ol->fMethodInfo || ol->fMethodInfo->fMethods.empty())
            return nullptr;

    // find the overload with matching signature
        for (auto& m : ol->fMethodInfo->fMethods) {
            PyObject* sig = m->GetSignature(false);
            bool found = signature == CPyCppyy_PyText_AsString(sig);
            Py_DECREF(sig);
            if (found) {
                void* fptr = (void*)m->GetFunctionAddress();
                if (fptr) return fptr;
                break;  // fall-through, with calling through Python
            }
        }
    }

    if (TemplateProxy_Check(pyobject)) {
    // get the actual underlying template matching the signature
        TemplateProxy* pytmpl = (TemplateProxy*)pyobject;
        std::string fullname = CPyCppyy_PyText_AsString(pytmpl->fTI->fCppName);
        if (pytmpl->fTemplateArgs)
            fullname += CPyCppyy_PyText_AsString(pytmpl->fTemplateArgs);
        Cppyy::TCppScope_t scope = ((CPPClass*)pytmpl->fTI->fPyClass)->fCppType;
        Cppyy::TCppMethod_t cppmeth = Cppyy::GetMethodTemplate(scope, fullname, signature);
        if (cppmeth) {
            void* fptr = (void*)Cppyy::GetFunctionAddress(cppmeth, false);
            if (fptr) return fptr;
        }
        // fall-through, with calling through Python
    }

    if (PyCallable_Check(pyobject)) {
    // generic python callable: create a C++ wrapper function
        void* wpraddress = nullptr;

    // re-use existing wrapper if possible
        auto key = std::make_pair(rettype, signature);
        const auto& lookup = sWrapperLookup.find(key);
        if (lookup != sWrapperLookup.end()) {
            const auto& existing = lookup->second.find(pyobject);
            if (existing != lookup->second.end() && *sWrapperReference[existing->second] == pyobject)
                wpraddress = existing->second;
        }

     // check for a pre-existing, unused, wrapper if not found
        if (!wpraddress) {
           const auto& freewrap = sWrapperFree.find(key);
           if (freewrap != sWrapperFree.end() && !freewrap->second.empty()) {
               wpraddress = freewrap->second.back();
               freewrap->second.pop_back();
               *sWrapperReference[wpraddress] = pyobject;
               PyObject* wref = PyWeakref_NewRef(pyobject, sWrapperCacheEraser);
               if (wref) sWrapperWeakRefs[wref] = std::make_pair(wpraddress, key);
               else PyErr_Clear();     // happens for builtins which don't need this
           }
        }

     // create wrapper if no re-use possible
        if (!wpraddress) {
            if (!Utility::IncludePython())
                return nullptr;

        // extract argument types
            const std::vector<std::string>& argtypes = TypeManip::extract_arg_types(signature);
            int nArgs = (int)argtypes.size();

        // wrapper name
            std::ostringstream wname;
            wname << "fptr_wrap" << ++sWrapperCounter;

       // build wrapper function code
            std::ostringstream code;
            code << "namespace __cppyy_internal {\n  "
                 << rettype << " " << wname.str() << "(";
            for (int i = 0; i < nArgs; ++i) {
                code << argtypes[i] << " arg" << i;
                if (i != nArgs-1) code << ", ";
            }
            code << ") {\n";

        // start function body
            Utility::ConstructCallbackPreamble(rettype, argtypes, code);

        // create a referencable pointer
            PyObject** ref = new PyObject*{pyobject};

        // function call itself and cleanup
            code << "    PyObject** ref = (PyObject**)" << (intptr_t)ref << ";\n"
                    "    PyObject* pyresult = nullptr;\n"
                    "    if (*ref) pyresult = PyObject_CallFunctionObjArgs(*ref";
            for (int i = 0; i < nArgs; ++i)
                code << ", pyargs[" << i << "]";
            code << ", NULL);\n"
                    "    else PyErr_SetString(PyExc_TypeError, \"callable was deleted\");\n";

        // close
            Utility::ConstructCallbackReturn(rettype, nArgs, code);

        // end of namespace
            code << "}";

        // finally, compile the code
            if (!Cppyy::Compile(code.str()))
                return nullptr;

        // TODO: is there no easier way?
            static Cppyy::TCppScope_t scope = Cppyy::GetScope("__cppyy_internal");
            const auto& idx = Cppyy::GetMethodIndicesFromName(scope, wname.str());
            wpraddress = Cppyy::GetFunctionAddress(Cppyy::GetMethod(scope, idx[0]), false);
            sWrapperReference[wpraddress] = ref;

        // cache the new wrapper
            sWrapperLookup[key][pyobject] = wpraddress;
            PyObject* wref = PyWeakref_NewRef(pyobject, sWrapperCacheEraser);
            if (wref) sWrapperWeakRefs[wref] = std::make_pair(wpraddress, key);
            else PyErr_Clear();     // happens for builtins which don't need this
        }

    // now pass the pointer to the wrapper function (may be null)
        return wpraddress;
    }

    return nullptr;
}

bool CPyCppyy::FunctionPointerConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /*ctxt*/)
{
// special case: allow nullptr singleton:
    if (gNullPtrObject == pyobject) {
        para.fValue.fVoidp = nullptr;
        para.fTypeCode = 'p';
        return true;
    }

// normal case, get a function pointer
    void* fptr = PyFunction_AsCPointer(pyobject, fRetType, fSignature);
    if (fptr) {
        para.fValue.fVoidp = fptr;
        para.fTypeCode = 'p';
        return true;
    }

    return false;
}

static std::map<void*, std::string> sFuncWrapperLookup;
static const char* FPCFM_ERRMSG = "conversion to std::function failed";
PyObject* CPyCppyy::FunctionPointerConverter::FromMemory(void* address)
{
// A function pointer in clang is represented by a Type, not a FunctionDecl and it's
// not possible to get the latter from the former: the backend will need to support
// both. Since that is far in the future, we'll use a std::function instead.
    static int func_count = 0;

    if (!(address && *(void**)address)) {
        PyErr_SetString(PyExc_TypeError, FPCFM_ERRMSG);
        return nullptr;
    }

    void* faddr = *(void**)address;
    auto cached = sFuncWrapperLookup.find(faddr);
    if (cached == sFuncWrapperLookup.end()) {
        std::ostringstream fname;
        fname << "ptr2func" << ++func_count;

        std::ostringstream code;
        code << "namespace __cppyy_internal {\n  std::function<"
             << fRetType << fSignature << "> " << fname.str()
             << " = (" << fRetType << "(*)" << fSignature << ")" << (intptr_t)faddr
             << ";\n}";

        if (!Cppyy::Compile(code.str())) {
            PyErr_SetString(PyExc_TypeError, FPCFM_ERRMSG);
            return nullptr;
        }

     // cache the new wrapper (TODO: does it make sense to use weakrefs on the data
     // member?)
        sFuncWrapperLookup[faddr] = fname.str();
        cached = sFuncWrapperLookup.find(faddr);
    }

    static Cppyy::TCppScope_t scope = Cppyy::GetScope("__cppyy_internal");
    PyObject* pyscope = CreateScopeProxy(scope);
    PyObject* func = PyObject_GetAttrString(pyscope, cached->second.c_str());
    Py_DECREF(pyscope);

    return func;
}

bool CPyCppyy::FunctionPointerConverter::ToMemory(PyObject* pyobject, void* address)
{
// special case: allow nullptr singleton:
    if (gNullPtrObject == pyobject) {
        *((void**)address) = nullptr;
        return true;
    }

// normal case, get a function pointer
    void* fptr = PyFunction_AsCPointer(pyobject, fRetType, fSignature);
    if (fptr) {
        *((void**)address) = fptr;
        return true;
    }

    return false;
}


//- std::function converter --------------------------------------------------
bool CPyCppyy::StdFunctionConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
// prefer normal "object" conversion
    bool rf = ctxt->fFlags & CallContext::kNoImplicit;
    ctxt->fFlags |= CallContext::kNoImplicit;
    if (fConverter->SetArg(pyobject, para, ctxt)) {
        if (!rf) ctxt->fFlags &= ~CallContext::kNoImplicit;
        return true;
    }

    PyErr_Clear();

// else create a wrapper function
    if (this->FunctionPointerConverter::SetArg(pyobject, para, ctxt)) {
    // retrieve the wrapper pointer and capture it in a temporary std::function,
    // then try normal conversion a second time
        PyObject* func = this->FunctionPointerConverter::FromMemory(&para.fValue.fVoidp);
        if (func) {
            Py_XDECREF(fFuncWrap); fFuncWrap = func;
            bool result = fConverter->SetArg(fFuncWrap, para, ctxt);
            if (!rf) ctxt->fFlags &= ~CallContext::kNoImplicit;
            return result;
        }
    }

    if (!rf) ctxt->fFlags &= ~CallContext::kNoImplicit;
    return false;
}

PyObject* CPyCppyy::StdFunctionConverter::FromMemory(void* address)
{
    return fConverter->FromMemory(address);
}

bool CPyCppyy::StdFunctionConverter::ToMemory(PyObject* value, void* address)
{
    return fConverter->ToMemory(value, address);
}


//- smart pointer converters -------------------------------------------------
bool CPyCppyy::SmartPtrConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* ctxt)
{
    char typeCode = fIsRef ? 'p' : 'V';

    if (!CPPInstance_Check(pyobject)) {
        // TODO: not sure how this is correct for pass-by-ref nor does it help with
        // implicit conversions for pass-by-value
        if (fIsRef && GetAddressSpecialCase(pyobject, para.fValue.fVoidp)) {
            para.fTypeCode = typeCode;      // allow special cases such as nullptr
            return true;
        }

        return false;
    }

    CPPInstance* pyobj = (CPPInstance*)pyobject;

// for the case where we have a 'hidden' smart pointer:
    if (Cppyy::TCppType_t tsmart = pyobj->GetSmartIsA()) {
        if (Cppyy::IsSubtype(tsmart, fSmartPtrType)) {
        // depending on memory policy, some objects need releasing when passed into functions
            if (fKeepControl && !UseStrictOwnership(ctxt))
                ((CPPInstance*)pyobject)->CppOwns();

        // calculate offset between formal and actual arguments
            para.fValue.fVoidp = pyobj->GetSmartObject();
            if (tsmart != fSmartPtrType) {
                para.fValue.fIntPtr += Cppyy::GetBaseOffset(
                    tsmart, fSmartPtrType, para.fValue.fVoidp, 1 /* up-cast */);
            }

        // set pointer (may be null) and declare success
            para.fTypeCode = typeCode;
            return true;
        }
    }

// for the case where we have an 'exposed' smart pointer:
    if (!pyobj->IsSmart() && Cppyy::IsSubtype(pyobj->ObjectIsA(), fSmartPtrType)) {
    // calculate offset between formal and actual arguments
        para.fValue.fVoidp = pyobj->GetObject();
        if (pyobj->ObjectIsA() != fSmartPtrType) {
            para.fValue.fIntPtr += Cppyy::GetBaseOffset(
                pyobj->ObjectIsA(), fSmartPtrType, para.fValue.fVoidp, 1 /* up-cast */);
        }

    // set pointer (may be null) and declare success
        para.fTypeCode = typeCode;
        return true;
    }

// final option, try mapping pointer types held (TODO: do not allow for non-const ref)
    if (pyobj->IsSmart() && Cppyy::IsSubtype(pyobj->ObjectIsA(), fUnderlyingType)) {
        para.fValue.fVoidp = ((CPPInstance*)pyobject)->GetSmartObject();
        para.fTypeCode = 'V';
        return true;
    }

    return false;
}

PyObject* CPyCppyy::SmartPtrConverter::FromMemory(void* address)
{
    if (!address || !fSmartPtrType)
        return nullptr;

    return BindCppObjectNoCast(address, fSmartPtrType);
}


//----------------------------------------------------------------------------
namespace {

// clang libcpp and gcc use the same structure (ptr, size)
#if defined (_LIBCPP_INITIALIZER_LIST) || defined(__GNUC__)
struct faux_initlist
{
    typedef size_t size_type;
    typedef void*  iterator;
    iterator  _M_array;
    size_type _M_len;
};
#elif defined (_MSC_VER)
struct faux_initlist
{
     typedef char* iterator;
     iterator _M_array; // ie. _First;
     iterator _Last;
};
#else
#define NO_KNOWN_INITIALIZER_LIST 1
#endif

} // unnamed namespace

CPyCppyy::InitializerListConverter::~InitializerListConverter()
{
    if (fConverter && fConverter->HasState()) delete fConverter;
}

bool CPyCppyy::InitializerListConverter::SetArg(
    PyObject* pyobject, Parameter& para, CallContext* /*ctxt*/)
{
#ifdef NO_KNOWN_INITIALIZER_LIST
    return false;
#else
// convert the given argument to an initializer list temporary; this is purely meant
// to be a syntactic thing, so only _python_ sequences are allowed; bound C++ proxies
// are therefore rejected (should go through eg. a copy constructor etc.)
    if (CPPInstance_Check(pyobject) || !PySequence_Check(pyobject) || CPyCppyy_PyText_Check(pyobject)
#if PY_VERSION_HEX >= 0x03000000
        || PyBytes_Check(pyobject)
#endif
        )
        return false;

    void* buf;
    Py_ssize_t buflen = Utility::GetBuffer(pyobject, '*', (int)fValueSize, buf, true);
    faux_initlist* fake = nullptr;
    if (buf && buflen) {
    // dealing with an array here, pass on whole-sale
        fake = (faux_initlist*)malloc(sizeof(faux_initlist));
	fake->_M_array = (faux_initlist::iterator)buf;
#if defined (_LIBCPP_INITIALIZER_LIST) || defined(__GNUC__)
        fake->_M_len = (faux_initlist::size_type)buflen;
#elif defined (_MSC_VER)
        fake->_Last = fake->_M_array+buflen*fValueSize;
#endif
    } else {
    // can only construct empty lists, so use a fake initializer list
        size_t len = (size_t)PySequence_Size(pyobject);
        fake = (faux_initlist*)malloc(sizeof(faux_initlist)+fValueSize*len);
        fake->_M_array = (faux_initlist::iterator)((char*)fake+sizeof(faux_initlist));
#if defined (_LIBCPP_INITIALIZER_LIST) || defined(__GNUC__)
        fake->_M_len = (faux_initlist::size_type)len;
        for (faux_initlist::size_type i = 0; i < fake->_M_len; ++i) {
#elif defined (_MSC_VER)
        fake->_Last = fake->_M_array+len*fValueSize;
        for (size_t i = 0; (fake->_M_array+i*fValueSize) != fake->_Last; ++i) {
#endif
            PyObject* item = PySequence_GetItem(pyobject, i);
            bool convert_ok = false;
            if (item) {
                if (!fConverter) {
                    if (CPPInstance_Check(item)) {
                    // by convention, use byte copy
                        memcpy((char*)fake->_M_array + i*fValueSize,
                               ((CPPInstance*)item)->GetObject(), fValueSize);
                        convert_ok = true;
                    }
                } else
                    convert_ok = fConverter->ToMemory(item, (char*)fake->_M_array + i*fValueSize);

                Py_DECREF(item);
            } else
                PyErr_Format(PyExc_TypeError, "failed to get item %d from sequence", (int)i);

            if (!convert_ok) {
                free((void*)fake);
                return false;
            }
        }
    }

    para.fValue.fVoidp = (void*)fake;
    para.fTypeCode = 'X';     // means ptr that backend has to free after call
    return true;
#endif
}

//----------------------------------------------------------------------------
bool CPyCppyy::NotImplementedConverter::SetArg(PyObject*, Parameter&, CallContext*)
{
// raise a NotImplemented exception to take a method out of overload resolution
    PyErr_SetString(PyExc_NotImplementedError, "this method cannot (yet) be called");
    return false;
}


//- helper to refactor some code from CreateConverter ------------------------
static inline CPyCppyy::Converter* selectInstanceCnv(
    Cppyy::TCppScope_t klass, const std::string& cpd, long size, dims_t dims, bool isConst, bool control)
{
    using namespace CPyCppyy;
    Converter* result = nullptr;

    if (cpd == "**" || cpd == "*[]" || cpd == "&*")
        result = new InstancePtrPtrConverter<false>(klass, control);
    else if (cpd == "*&")
        result = new InstancePtrPtrConverter<true>(klass, control);
    else if (cpd == "*" && size <= 0)
        result = new InstancePtrConverter(klass, control);
    else if (cpd == "&")
        result = new InstanceRefConverter(klass, isConst);
    else if (cpd == "&&")
        result = new InstanceMoveConverter(klass);
    else if (cpd == "[]" || size > 0)
        result = new InstanceArrayConverter(klass, dims, false);
    else if (cpd == "")             // by value
        result = new InstanceConverter(klass, true);

    return result;
}

//- factories ----------------------------------------------------------------
CPYCPPYY_EXPORT
CPyCppyy::Converter* CPyCppyy::CreateConverter(const std::string& fullType, dims_t dims)
{
// The matching of the fulltype to a converter factory goes through up to five levels:
//   1) full, exact match
//   2) match of decorated, unqualified type
//   3) accept const ref as by value
//   4) accept ref as pointer
//   5) generalized cases (covers basically all C++ classes)
//
// If all fails, void is used, which will generate a run-time warning when used.

    dim_t size = (dims && dims[0] != -1) ? dims[1] : -1;

// an exactly matching converter is best
    ConvFactories_t::iterator h = gConvFactories.find(fullType);
    if (h != gConvFactories.end())
        return (h->second)(dims);

// resolve typedefs etc.
    const std::string& resolvedType = Cppyy::ResolveName(fullType);

// a full, qualified matching converter is preferred
    if (resolvedType != fullType) {
        h = gConvFactories.find(resolvedType);
        if (h != gConvFactories.end())
            return (h->second)(dims);
    }

//-- nothing? ok, collect information about the type and possible qualifiers/decorators
    bool isConst = strncmp(resolvedType.c_str(), "const", 5) == 0;
    const std::string& cpd = Utility::Compound(resolvedType);
    std::string realType   = TypeManip::clean_type(resolvedType, false, true);

// accept unqualified type (as python does not know about qualifiers)
    h = gConvFactories.find(realType + cpd);
    if (h != gConvFactories.end())
        return (h->second)(dims);

// drop const, as that is mostly meaningless to python (with the exception
// of c-strings, but those are specialized in the converter map)
    if (isConst) {
        realType = TypeManip::remove_const(realType);
        h = gConvFactories.find(realType + cpd);
        if (h != gConvFactories.end())
            return (h->second)(dims);
    }

//-- still nothing? try pointer instead of array (for builtins)
    if (cpd == "[]") {
    // simple array
        h = gConvFactories.find(realType + "*");
        if (h != gConvFactories.end()) {
            if (dims && dims[1] == UNKNOWN_SIZE) dims[1] = UNKNOWN_ARRAY_SIZE;
            return (h->second)(dims);
        }
    } else if (cpd == "*[]") {
    // array of pointers
        h = gConvFactories.find(realType + "*");
        if (h != gConvFactories.end()) {
        // upstream treats the pointer type as the array element type, but that pointer is
        // treated as a low-level view as well, so adjust the dims
            dim_t newdim = (dims && 0 < dims[0]) ? dims[0]+1 : 2;
            dims_t newdims = new dim_t[newdim+1];
            newdims[0] = newdim;
            newdims[1] = (0 < size ? size : UNKNOWN_ARRAY_SIZE);      // the array
            newdims[2] = UNKNOWN_SIZE;                                // the pointer
            if (dims && 2 < newdim) {
                for (int i = 2; i < (newdim-1); ++i)
                    newdims[i+1] = dims[i];
            }
            Converter* cnv = (h->second)(newdims);
            delete [] newdims;
            return cnv;
        }
    }

//-- special case: initializer list
    if (realType.compare(0, 16, "initializer_list") == 0) {
    // get the type of the list and create a converter (TODO: get hold of value_type?)
        auto pos = realType.find('<');
        std::string value_type = realType.substr(pos+1, realType.size()-pos-2);
        Converter* cnv = nullptr; bool use_byte_cnv = false;
        if (cpd == "" && Cppyy::GetScope(value_type)) {
        // initializer list of object values does not work as the target is raw
        // memory; simply use byte copies

        // by convention, leave cnv as nullptr
            use_byte_cnv = true;
        } else
            cnv = CreateConverter(value_type);
        if (cnv || use_byte_cnv)
            return new InitializerListConverter(cnv, Cppyy::SizeOf(value_type));
    }

//-- still nothing? use a generalized converter
    bool control = cpd == "&" || isConst;

//-- special case: std::function
    auto pos = resolvedType.find("function<");
    if (pos == 0 /* no std:: */ || pos == 5 /* with std:: */ ||
        pos == 6 /* const no std:: */ || pos == 11 /* const with std:: */ ) {

    // get actual converter for normal passing
        Converter* cnv = selectInstanceCnv(
            Cppyy::GetScope(realType), cpd, size, dims, isConst, control);

        if (cnv) {
        // get the type of the underlying (TODO: use target_type?)
            auto pos1 = resolvedType.find("(", pos+9);
            auto pos2 = resolvedType.rfind(")");
            if (pos1 != std::string::npos && pos2 != std::string::npos) {
                auto sz1 = pos1-pos-9;
                if (resolvedType[pos+9+sz1-1] == ' ') sz1 -= 1;

                return new StdFunctionConverter(cnv,
                    resolvedType.substr(pos+9, sz1), resolvedType.substr(pos1, pos2-pos1+1));
            }
        }
    }

// converters for known C++ classes and default (void*)
    Converter* result = nullptr;
    if (Cppyy::TCppScope_t klass = Cppyy::GetScope(realType)) {
        Cppyy::TCppType_t raw{0};
        if (Cppyy::GetSmartPtrInfo(realType, &raw, nullptr)) {
            if (cpd == "") {
                result = new SmartPtrConverter(klass, raw, control);
            } else if (cpd == "&") {
                result = new SmartPtrConverter(klass, raw);
            } else if (cpd == "*" && size <= 0) {
                result = new SmartPtrConverter(klass, raw, control, true);
            }
        }

        if (!result) {
        // CLING WORKAROUND -- special case for STL iterators
            if (realType.rfind("__gnu_cxx::__normal_iterator", 0) /* vector */ == 0
#ifdef __APPLE__
                || realType.rfind("__wrap_iter", 0) == 0
#endif
                // TODO: Windows?
               ) {
                static STLIteratorConverter c;
                result = &c;
            } else
       // -- CLING WORKAROUND
                result = selectInstanceCnv(klass, cpd, size, dims, isConst, control);
        }
    } else if (resolvedType.find("(*)") != std::string::npos ||
               (resolvedType.find("::*)") != std::string::npos)) {
    // this is a function function pointer
    // TODO: find better way of finding the type
        auto pos1 = resolvedType.find('(');
        auto pos2 = resolvedType.find("*)");
        auto pos3 = resolvedType.rfind(')');
        result = new FunctionPointerConverter(
            resolvedType.substr(0, pos1), resolvedType.substr(pos2+2, pos3-pos2-1));
    }

    if (!result && cpd == "&&") {
    // for builtin, can use const-ref for r-ref
        h = gConvFactories.find("const " + realType + "&");
        if (h != gConvFactories.end())
            return (h->second)(dims);
    // else, unhandled moves
        result = new NotImplementedConverter();
    }

    if (!result && h != gConvFactories.end())
    // converter factory available, use it to create converter
        result = (h->second)(dims);
    else if (!result) {
    // default to something reasonable, assuming "user knows best"
        if (cpd.size() == 2 && cpd != "&&") // "**", "*[]", "*&"
            result = new VoidPtrPtrConverter(size);
        else if (!cpd.empty())
            result = new VoidArrayConverter();        // "user knows best"
        else
            result = new NotImplementedConverter();   // fails on use
    }

    return result;
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
void CPyCppyy::DestroyConverter(Converter* p)
{
    if (p && p->HasState())
        delete p;  // state-less converters are always shared
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
bool CPyCppyy::RegisterConverter(const std::string& name, cf_t fac)
{
// register a custom converter
    auto f = gConvFactories.find(name);
    if (f != gConvFactories.end())
        return false;

    gConvFactories[name] = fac;
    return true;
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
bool CPyCppyy::UnregisterConverter(const std::string& name)
{
// remove a custom converter
    auto f = gConvFactories.find(name);
    if (f != gConvFactories.end()) {
        gConvFactories.erase(f);
        return true;
    }
    return false;
}


//----------------------------------------------------------------------------
namespace {

using namespace CPyCppyy;

#define STRINGVIEW "basic_string_view<char,char_traits<char> >"
#define WSTRING "basic_string<wchar_t,char_traits<wchar_t>,allocator<wchar_t> >"

static struct InitConvFactories_t {
public:
    InitConvFactories_t() {
    // load all converter factories in the global map 'gConvFactories'
        CPyCppyy::ConvFactories_t& gf = gConvFactories;

    // factories for built-ins
        gf["bool"] =                        (cf_t)+[](dims_t) { static BoolConverter c{};           return &c; };
        gf["const bool&"] =                 (cf_t)+[](dims_t) { static ConstBoolRefConverter c{};   return &c; };
        gf["bool&"] =                       (cf_t)+[](dims_t) { static BoolRefConverter c{};        return &c; };
        gf["char"] =                        (cf_t)+[](dims_t) { static CharConverter c{};           return &c; };
        gf["const char&"] =                 (cf_t)+[](dims_t) { static ConstCharRefConverter c{};   return &c; };
        gf["char&"] =                       (cf_t)+[](dims_t) { static CharRefConverter c{};        return &c; };
        gf["signed char&"] =                (cf_t)+[](dims_t) { static SCharRefConverter c{};       return &c; };
        gf["unsigned char"] =               (cf_t)+[](dims_t) { static UCharConverter c{};          return &c; };
        gf["const unsigned char&"] =        (cf_t)+[](dims_t) { static ConstUCharRefConverter c{};  return &c; };
        gf["unsigned char&"] =              (cf_t)+[](dims_t) { static UCharRefConverter c{};       return &c; };
        gf["UCharAsInt"] =                  (cf_t)+[](dims_t) { static UCharAsIntConverter c{};     return &c; };
        gf["wchar_t"] =                     (cf_t)+[](dims_t) { static WCharConverter c{};          return &c; };
        gf["char16_t"] =                    (cf_t)+[](dims_t) { static Char16Converter c{};         return &c; };
        gf["char32_t"] =                    (cf_t)+[](dims_t) { static Char32Converter c{};         return &c; };
        gf["wchar_t&"] =                    (cf_t)+[](dims_t) { static WCharRefConverter c{};       return &c; };
        gf["char16_t&"] =                   (cf_t)+[](dims_t) { static Char16RefConverter c{};      return &c; };
        gf["char32_t&"] =                   (cf_t)+[](dims_t) { static Char32RefConverter c{};      return &c; };
        gf["int8_t"] =                      (cf_t)+[](dims_t) { static Int8Converter c{};           return &c; };
        gf["int8_t&"] =                     (cf_t)+[](dims_t) { static Int8RefConverter c{};        return &c; };
        gf["const int8_t&"] =               (cf_t)+[](dims_t) { static ConstInt8RefConverter c{};   return &c; };
        gf["uint8_t"] =                     (cf_t)+[](dims_t) { static UInt8Converter c{};          return &c; };
        gf["const uint8_t&"] =              (cf_t)+[](dims_t) { static ConstUInt8RefConverter c{};  return &c; };
        gf["uint8_t&"] =                    (cf_t)+[](dims_t) { static UInt8RefConverter c{};       return &c; };
        gf["short"] =                       (cf_t)+[](dims_t) { static ShortConverter c{};          return &c; };
        gf["const short&"] =                (cf_t)+[](dims_t) { static ConstShortRefConverter c{};  return &c; };
        gf["short&"] =                      (cf_t)+[](dims_t) { static ShortRefConverter c{};       return &c; };
        gf["unsigned short"] =              (cf_t)+[](dims_t) { static UShortConverter c{};         return &c; };
        gf["const unsigned short&"] =       (cf_t)+[](dims_t) { static ConstUShortRefConverter c{}; return &c; };
        gf["unsigned short&"] =             (cf_t)+[](dims_t) { static UShortRefConverter c{};      return &c; };
        gf["int"] =                         (cf_t)+[](dims_t) { static IntConverter c{};            return &c; };
        gf["int&"] =                        (cf_t)+[](dims_t) { static IntRefConverter c{};         return &c; };
        gf["const int&"] =                  (cf_t)+[](dims_t) { static ConstIntRefConverter c{};    return &c; };
        gf["unsigned int"] =                (cf_t)+[](dims_t) { static UIntConverter c{};           return &c; };
        gf["const unsigned int&"] =         (cf_t)+[](dims_t) { static ConstUIntRefConverter c{};   return &c; };
        gf["unsigned int&"] =               (cf_t)+[](dims_t) { static UIntRefConverter c{};        return &c; };
        gf["long"] =                        (cf_t)+[](dims_t) { static LongConverter c{};           return &c; };
        gf["long&"] =                       (cf_t)+[](dims_t) { static LongRefConverter c{};        return &c; };
        gf["const long&"] =                 (cf_t)+[](dims_t) { static ConstLongRefConverter c{};   return &c; };
        gf["unsigned long"] =               (cf_t)+[](dims_t) { static ULongConverter c{};          return &c; };
        gf["const unsigned long&"] =        (cf_t)+[](dims_t) { static ConstULongRefConverter c{};  return &c; };
        gf["unsigned long&"] =              (cf_t)+[](dims_t) { static ULongRefConverter c{};       return &c; };
        gf["long long"] =                   (cf_t)+[](dims_t) { static LLongConverter c{};          return &c; };
        gf["const long long&"] =            (cf_t)+[](dims_t) { static ConstLLongRefConverter c{};  return &c; };
        gf["long long&"] =                  (cf_t)+[](dims_t) { static LLongRefConverter c{};       return &c; };
        gf["unsigned long long"] =          (cf_t)+[](dims_t) { static ULLongConverter c{};         return &c; };
        gf["const unsigned long long&"] =   (cf_t)+[](dims_t) { static ConstULLongRefConverter c{}; return &c; };
        gf["unsigned long long&"] =         (cf_t)+[](dims_t) { static ULLongRefConverter c{};      return &c; };

        gf["float"] =                       (cf_t)+[](dims_t) { static FloatConverter c{};           return &c; };
        gf["const float&"] =                (cf_t)+[](dims_t) { static ConstFloatRefConverter c{};   return &c; };
        gf["float&"] =                      (cf_t)+[](dims_t) { static FloatRefConverter c{};        return &c; };
        gf["double"] =                      (cf_t)+[](dims_t) { static DoubleConverter c{};          return &c; };
        gf["double&"] =                     (cf_t)+[](dims_t) { static DoubleRefConverter c{};       return &c; };
        gf["const double&"] =               (cf_t)+[](dims_t) { static ConstDoubleRefConverter c{};  return &c; };
        gf["long double"] =                 (cf_t)+[](dims_t) { static LDoubleConverter c{};         return &c; };
        gf["const long double&"] =          (cf_t)+[](dims_t) { static ConstLDoubleRefConverter c{}; return &c; };
        gf["long double&"] =                (cf_t)+[](dims_t) { static LDoubleRefConverter c{};      return &c; };
        gf["std::complex<double>"] =        (cf_t)+[](dims_t) { return new ComplexDConverter{}; };
        gf["complex<double>"] =             (cf_t)+[](dims_t) { return new ComplexDConverter{}; };
        gf["const std::complex<double>&"] = (cf_t)+[](dims_t) { return new ComplexDConverter{}; };
        gf["const complex<double>&"] =      (cf_t)+[](dims_t) { return new ComplexDConverter{}; };
        gf["void"] =                        (cf_t)+[](dims_t) { static VoidConverter c{};            return &c; };

    // pointer/array factories
        gf["bool*"] =                       (cf_t)+[](dims_t d) { return new BoolArrayConverter{d}; };
        gf["bool**"] =                      (cf_t)+[](dims_t d) { return new BoolArrayPtrConverter{d}; };
        gf["const signed char[]"] =         (cf_t)+[](dims_t d) { return new SCharArrayConverter{d}; };
        gf["signed char[]"] =               (cf_t)+[](dims_t d) { return new SCharArrayConverter{d}; };
        gf["signed char**"] =               (cf_t)+[](dims_t d) { return new SCharArrayPtrConverter{d}; };
        gf["const unsigned char*"] =        (cf_t)+[](dims_t d) { return new UCharArrayConverter{d}; };
        gf["unsigned char*"] =              (cf_t)+[](dims_t d) { return new UCharArrayConverter{d}; };
        gf["UCharAsInt*"] =                 (cf_t)+[](dims_t d) { return new UCharArrayConverter{d}; };
        gf["unsigned char**"] =             (cf_t)+[](dims_t d) { return new UCharArrayPtrConverter{d}; };
#if __cplusplus > 201402L
        gf["byte*"] =                       (cf_t)+[](dims_t d) { return new ByteArrayConverter{d}; };
        gf["byte**"] =                      (cf_t)+[](dims_t d) { return new ByteArrayPtrConverter{d}; };
#endif
        gf["short*"] =                      (cf_t)+[](dims_t d) { return new ShortArrayConverter{d}; };
        gf["short**"] =                     (cf_t)+[](dims_t d) { return new ShortArrayPtrConverter{d}; };
        gf["unsigned short*"] =             (cf_t)+[](dims_t d) { return new UShortArrayConverter{d}; };
        gf["unsigned short**"] =            (cf_t)+[](dims_t d) { return new UShortArrayPtrConverter{d}; };
        gf["int*"] =                        (cf_t)+[](dims_t d) { return new IntArrayConverter{d}; };
        gf["int**"] =                       (cf_t)+[](dims_t d) { return new IntArrayPtrConverter{d}; };
        gf["unsigned int*"] =               (cf_t)+[](dims_t d) { return new UIntArrayConverter{d}; };
        gf["unsigned int**"] =              (cf_t)+[](dims_t d) { return new UIntArrayPtrConverter{d}; };
        gf["long*"] =                       (cf_t)+[](dims_t d) { return new LongArrayConverter{d}; };
        gf["long**"] =                      (cf_t)+[](dims_t d) { return new LongArrayPtrConverter{d}; };
        gf["unsigned long*"] =              (cf_t)+[](dims_t d) { return new ULongArrayConverter{d}; };
        gf["unsigned long**"] =             (cf_t)+[](dims_t d) { return new ULongArrayPtrConverter{d}; };
        gf["long long*"] =                  (cf_t)+[](dims_t d) { return new LLongArrayConverter{d}; };
        gf["long long**"] =                 (cf_t)+[](dims_t d) { return new LLongArrayPtrConverter{d}; };
        gf["unsigned long long*"] =         (cf_t)+[](dims_t d) { return new ULLongArrayConverter{d}; };
        gf["unsigned long long**"] =        (cf_t)+[](dims_t d) { return new ULLongArrayPtrConverter{d}; };
        gf["float*"] =                      (cf_t)+[](dims_t d) { return new FloatArrayConverter{d}; };
        gf["float**"] =                     (cf_t)+[](dims_t d) { return new FloatArrayPtrConverter{d}; };
        gf["double*"] =                     (cf_t)+[](dims_t d) { return new DoubleArrayConverter{d}; };
        gf["double**"] =                    (cf_t)+[](dims_t d) { return new DoubleArrayPtrConverter{d}; };
        gf["long double*"] =                (cf_t)+[](dims_t d) { return new LDoubleArrayConverter{d}; };
        gf["long double**"] =               (cf_t)+[](dims_t d) { return new LDoubleArrayPtrConverter{d}; };
        gf["std::complex<double>*"] =       (cf_t)+[](dims_t d) { return new ComplexDArrayConverter{d}; };
        gf["complex<double>*"] =            (cf_t)+[](dims_t d) { return new ComplexDArrayConverter{d}; };
        gf["std::complex<double>**"] =      (cf_t)+[](dims_t d) { return new ComplexDArrayPtrConverter{d}; };
        gf["void*"] =                       (cf_t)+[](dims_t d) { return new VoidArrayConverter{(bool)d}; };

    // aliases
        gf["signed char"] =                 gf["char"];
        gf["const signed char&"] =          gf["const char&"];
#if __cplusplus > 201402L
        gf["byte"] =                        gf["uint8_t"];
        gf["const byte&"] =                 gf["const uint8_t&"];
        gf["byte&"] =                       gf["uint8&"];
#endif
        gf["internal_enum_type_t"] =        gf["int"];
        gf["internal_enum_type_t&"] =       gf["int&"];
        gf["const internal_enum_type_t&"] = gf["const int&"];
        gf["Long64_t"] =                    gf["long long"];
        gf["Long64_t*"] =                   gf["long long*"];
        gf["Long64_t&"] =                   gf["long long&"];
        gf["const Long64_t&"] =             gf["const long long&"];
        gf["ULong64_t"] =                   gf["unsigned long long"];
        gf["ULong64_t*"] =                  gf["unsigned long long*"];
        gf["ULong64_t&"] =                  gf["unsigned long long&"];
        gf["const ULong64_t&"] =            gf["const unsigned long long&"];
        gf["Float16_t"] =                   gf["float"];
        gf["const Float16_t&"] =            gf["const float&"];
        gf["Double32_t"] =                  gf["double"];
        gf["Double32_t&"] =                 gf["double&"];
        gf["const Double32_t&"] =           gf["const double&"];

    // factories for special cases
        gf["TString"] =                     (cf_t)+[](dims_t) { return new TStringConverter{}; };
        gf["TString&"] =                    gf["TString"];
        gf["const TString&"] =              gf["TString"];
        gf["nullptr_t"] =                   (cf_t)+[](dims_t) { static NullptrConverter c{};        return &c;};
        gf["const char*"] =                 (cf_t)+[](dims_t) { return new CStringConverter{}; };
        gf["const signed char*"] =          gf["const char*"];
        gf["const char[]"] =                (cf_t)+[](dims_t) { return new CStringConverter{}; };
        gf["char*"] =                       (cf_t)+[](dims_t) { return new NonConstCStringConverter{}; };
        gf["signed char*"] =                gf["char*"];
        gf["wchar_t*"] =                    (cf_t)+[](dims_t) { return new WCStringConverter{}; };
        gf["char16_t*"] =                   (cf_t)+[](dims_t) { return new CString16Converter{}; };
        gf["char32_t*"] =                   (cf_t)+[](dims_t) { return new CString32Converter{}; };
    // TODO: the following are handled incorrectly upstream (char16_t** where char16_t* intended)?!
        gf["char16_t**"] =                  gf["char16_t*"];
        gf["char32_t**"] =                  gf["char32_t*"];
        gf["const char**"] =                (cf_t)+[](dims_t d) { return new CStringArrayConverter{d}; };
        gf["char**"] =                      gf["const char**"];
        gf["const char*[]"] =               gf["const char**"];
        gf["char*[]"] =                     gf["const char*[]"];
        gf["std::string"] =                 (cf_t)+[](dims_t) { return new STLStringConverter{}; };
        gf["string"] =                      gf["std::string"];
        gf["const std::string&"] =          gf["std::string"];
        gf["const string&"] =               gf["std::string"];
        gf["string&&"] =                    (cf_t)+[](dims_t) { return new STLStringMoveConverter{}; };
        gf["std::string&&"] =               gf["string&&"];
        gf["std::string_view"] =            (cf_t)+[](dims_t) { return new STLStringViewConverter{}; };
        gf["string_view"] =                 gf["std::string_view"];
        gf[STRINGVIEW] =                    gf["std::string_view"];
        gf["experimental::" STRINGVIEW] =   gf["std::string_view"];
        gf["std::string_view&"] =           gf["std::string_view"];
        gf["const string_view&"] =          gf["std::string_view"];
        gf["const " STRINGVIEW "&"] =       gf["std::string_view"];
        gf["std::wstring"] =                (cf_t)+[](dims_t) { return new STLWStringConverter{}; };
        gf[WSTRING] =                       gf["std::wstring"];
        gf["std::" WSTRING] =               gf["std::wstring"];
        gf["const std::wstring&"] =         gf["std::wstring"];
        gf["const std::" WSTRING "&"] =     gf["std::wstring"];
        gf["const " WSTRING "&"] =          gf["std::wstring"];
        gf["void*&"] =                      (cf_t)+[](dims_t) { static VoidPtrRefConverter c{};     return &c; };
        gf["void**"] =                      (cf_t)+[](dims_t d) { return new VoidPtrPtrConverter{size_t((d && d[0] != -1) ? d[1] : -1)}; };
        gf["void*[]"] =                     (cf_t)+[](dims_t d) { return new VoidPtrPtrConverter{size_t((d && d[0] != -1) ? d[1] : -1)}; };
        gf["PyObject*"] =                   (cf_t)+[](dims_t) { static PyObjectConverter c{};       return &c; };
        gf["_object*"] =                    gf["PyObject*"];
        gf["FILE*"] =                       (cf_t)+[](dims_t) { return new VoidArrayConverter{}; };
    }
} initConvFactories_;

} // unnamed namespace
