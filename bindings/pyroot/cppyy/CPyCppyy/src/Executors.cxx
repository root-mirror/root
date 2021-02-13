// Bindings
#include "CPyCppyy.h"
#include "DeclareExecutors.h"
#include "CPPInstance.h"
#include "LowLevelViews.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <cstring>
#include <map>
#include <new>
#include <sstream>
#include <utility>
#include <sys/types.h>


//- data _____________________________________________________________________
namespace CPyCppyy {
    typedef Executor* (*ef_t) ();
    typedef std::map<std::string, ef_t> ExecFactories_t;
    static ExecFactories_t gExecFactories;

    extern PyObject* gNullPtrObject;

    extern std::set<std::string> gIteratorTypes;
}


//- helpers ------------------------------------------------------------------
namespace {

#ifdef WITH_THREAD
    class GILControl {
    public:
        GILControl() : fSave(PyEval_SaveThread()) { }
        ~GILControl() {
            PyEval_RestoreThread(fSave);
        }
    private:
        PyThreadState* fSave;
    };
#endif

} // unnamed namespace

#ifdef WITH_THREAD
#define CPPYY_IMPL_GILCALL(rtype, tcode)                                     \
static inline rtype GILCall##tcode(                                          \
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CPyCppyy::CallContext* ctxt)\
{                                                                            \
    if (!ReleasesGIL(ctxt))                                                  \
        return Cppyy::Call##tcode(method, self, ctxt->GetSize(), ctxt->GetArgs());\
    GILControl gc{};                                                         \
    return Cppyy::Call##tcode(method, self, ctxt->GetSize(), ctxt->GetArgs());\
}
#else
#define CPPYY_IMPL_GILCALL(rtype, tcode)                                     \
static inline rtype GILCall##tcode(                                          \
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CPyCppyy::CallContext* ctxt)\
{                                                                            \
    return Cppyy::Call##tcode(method, self, ctxt->GetSize(), ctxt->GetArgs());\
}
#endif

CPPYY_IMPL_GILCALL(void,          V)
CPPYY_IMPL_GILCALL(unsigned char, B)
CPPYY_IMPL_GILCALL(char,          C)
CPPYY_IMPL_GILCALL(short,         H)
CPPYY_IMPL_GILCALL(Int_t,         I)
CPPYY_IMPL_GILCALL(Long_t,        L)
CPPYY_IMPL_GILCALL(Long64_t,      LL)
CPPYY_IMPL_GILCALL(float,         F)
CPPYY_IMPL_GILCALL(double,        D)
CPPYY_IMPL_GILCALL(LongDouble_t,  LD)
CPPYY_IMPL_GILCALL(void*,         R)

static inline Cppyy::TCppObject_t GILCallO(Cppyy::TCppMethod_t method,
    Cppyy::TCppObject_t self, CPyCppyy::CallContext* ctxt, Cppyy::TCppType_t klass)
{
#ifdef WITH_THREAD
    if (!ReleasesGIL(ctxt))
#endif
        return Cppyy::CallO(method, self, ctxt->GetSize(), ctxt->GetArgs(), klass);
#ifdef WITH_THREAD
    GILControl gc{};
    return Cppyy::CallO(method, self, ctxt->GetSize(), ctxt->GetArgs(), klass);
#endif
}

static inline Cppyy::TCppObject_t GILCallConstructor(
    Cppyy::TCppMethod_t method, Cppyy::TCppType_t klass, CPyCppyy::CallContext* ctxt)
{
#ifdef WITH_THREAD
    if (!ReleasesGIL(ctxt))
#endif
        return Cppyy::CallConstructor(method, klass, ctxt->GetSize(), ctxt->GetArgs());
#ifdef WITH_THREAD
    GILControl gc{};
    return Cppyy::CallConstructor(method, klass, ctxt->GetSize(), ctxt->GetArgs());
#endif
}

static inline PyObject* CPyCppyy_PyText_FromLong(long cl)
{
// python chars are range(256)
    if (cl < -256 || cl > 255) {
        PyErr_SetString(PyExc_ValueError, "char conversion out of range");
        return nullptr;
    }
    int c = (int)cl;
    if (c < 0) return CPyCppyy_PyText_FromFormat("%c", 256 - std::abs(c));
    return CPyCppyy_PyText_FromFormat("%c", c);
}

static inline PyObject* CPyCppyy_PyText_FromULong(unsigned long uc)
{
// TODO: range check here?
    if (255 < uc) {
        PyErr_SetString(PyExc_ValueError, "char conversion out of range");
        return nullptr;
    }
    int c = (int)uc;
    return CPyCppyy_PyText_FromFormat("%c", c);
}

static inline PyObject* CPyCppyy_PyBool_FromLong(long b)
{
    PyObject* result = (bool)b ? Py_True : Py_False;
    Py_INCREF(result);
    return result;
}


//- base executor implementation ---------------------------------------------
CPyCppyy::Executor::~Executor()
{
    /* empty */
}

//- executors for built-ins --------------------------------------------------
PyObject* CPyCppyy::BoolExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python bool return value
    bool retval = GILCallB(method, self, ctxt);
    PyObject* result = retval ? Py_True : Py_False;
    Py_INCREF(result);
    return result;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::BoolConstRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python bool return value
    return CPyCppyy_PyBool_FromLong(*((bool*)GILCallR(method, self, ctxt)));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CharExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method with argument <self, ctxt>, construct python string return value
// with the single char
    return CPyCppyy_PyText_FromLong((int)GILCallC(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CharConstRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python string return value
// with the single char
    return CPyCppyy_PyText_FromLong(*((char*)GILCallR(method, self, ctxt)));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::UCharExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, args>, construct python string return value
// with the single char
    return CPyCppyy_PyText_FromLong((unsigned char)GILCallB(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::UCharConstRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python string return value
// with the single char from the pointer return
    return CPyCppyy_PyText_FromLong(*((unsigned char*)GILCallR(method, self, ctxt)));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::WCharExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, args>, construct python string return value
// with the single wide char
    wchar_t res = (wchar_t)GILCallL(method, self, ctxt);
    return PyUnicode_FromWideChar(&res, 1);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::Char16Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, args>, construct python string return value
// with the single char16
    char16_t res = (char16_t)GILCallL(method, self, ctxt);
    return PyUnicode_DecodeUTF16((const char*)&res, sizeof(char16_t), nullptr, nullptr);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::Char32Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, args>, construct python string return value
// with the single char32
    char32_t res = (char32_t)GILCallL(method, self, ctxt);
    return PyUnicode_DecodeUTF32((const char*)&res, sizeof(char32_t), nullptr, nullptr);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::IntExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python int return value
    return PyInt_FromLong((int)GILCallI(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::Int8Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python int return value
    return PyInt_FromLong((int8_t)GILCallC(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::UInt8Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python int return value
    return PyInt_FromLong((uint8_t)GILCallB(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::ShortExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python int return value
    return PyInt_FromLong((short)GILCallH(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::LongExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python long return value
    return PyLong_FromLong((Long_t)GILCallL(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::ULongExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python unsigned long return value
    return PyLong_FromUnsignedLong((ULong_t)GILCallLL(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::LongLongExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python long long return value
    Long64_t result = GILCallLL(method, self, ctxt);
    return PyLong_FromLongLong(result);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::ULongLongExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python unsigned long long return value
    ULong64_t result = (ULong64_t)GILCallLL(method, self, ctxt);
    return PyLong_FromUnsignedLongLong(result);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::FloatExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python float return value
    return PyFloat_FromDouble((double)GILCallF(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::DoubleExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python float return value
    return PyFloat_FromDouble((double)GILCallD(method, self, ctxt));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::LongDoubleExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python float return value
    return PyFloat_FromDouble((double)GILCallLD(method, self, ctxt));
}

//----------------------------------------------------------------------------
bool CPyCppyy::RefExecutor::SetAssignable(PyObject* pyobject)
{
// prepare "buffer" for by-ref returns, used with __setitem__
    if (pyobject) {
        Py_INCREF(pyobject);
        fAssignable = pyobject;
        return true;
    }

    fAssignable = nullptr;
    return false;
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_REFEXEC(name, type, stype, F1, F2)                        \
PyObject* CPyCppyy::name##RefExecutor::Execute(                              \
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt) \
{                                                                            \
    type* ref = (type*)GILCallR(method, self, ctxt);                         \
    if (!ref) { /* can happen if wrapper compilation fails */                \
        PyErr_SetString(PyExc_ReferenceError, "attempt to access a null-pointer");\
        return nullptr;                                                      \
    }                                                                        \
    if (!fAssignable)                                                        \
        return F1((stype)*ref);                                              \
    else {                                                                   \
        *ref = (type)F2(fAssignable);                                        \
        Py_DECREF(fAssignable);                                              \
        fAssignable = nullptr;                                               \
        if (*ref == (type)-1 && PyErr_Occurred())                            \
            return nullptr;                                                  \
        Py_INCREF(Py_None);                                                  \
        return Py_None;                                                      \
    }                                                                        \
}

CPPYY_IMPL_REFEXEC(Bool,   bool,   Long_t,   CPyCppyy_PyBool_FromLong, PyLong_AsLong)
CPPYY_IMPL_REFEXEC(Char,   char,   Long_t,   CPyCppyy_PyText_FromLong, PyLong_AsLong)
CPPYY_IMPL_REFEXEC(UChar,  unsigned char,  ULong_t,  CPyCppyy_PyText_FromULong, PyLongOrInt_AsULong)
CPPYY_IMPL_REFEXEC(Int8,   int8_t,  Long_t,  PyInt_FromLong, PyLong_AsLong)
CPPYY_IMPL_REFEXEC(UInt8,  uint8_t, ULong_t, PyInt_FromLong, PyLongOrInt_AsULong)
CPPYY_IMPL_REFEXEC(Short,  short,          Long_t,   PyInt_FromLong,     PyLong_AsLong)
CPPYY_IMPL_REFEXEC(UShort, unsigned short, ULong_t,  PyInt_FromLong,     PyLongOrInt_AsULong)
CPPYY_IMPL_REFEXEC(Int,    Int_t,    Long_t,   PyInt_FromLong,     PyLong_AsLong)
CPPYY_IMPL_REFEXEC(UInt,   UInt_t,   ULong_t,  PyLong_FromUnsignedLong, PyLongOrInt_AsULong)
CPPYY_IMPL_REFEXEC(Long,   Long_t,   Long_t,   PyLong_FromLong,    PyLong_AsLong)
CPPYY_IMPL_REFEXEC(ULong,  ULong_t,  ULong_t,  PyLong_FromUnsignedLong, PyLongOrInt_AsULong)
CPPYY_IMPL_REFEXEC(LongLong,  Long64_t,  Long64_t,   PyLong_FromLongLong,         PyLong_AsLongLong)
CPPYY_IMPL_REFEXEC(ULongLong, ULong64_t, ULong64_t,  PyLong_FromUnsignedLongLong, PyLongOrInt_AsULong64)
CPPYY_IMPL_REFEXEC(Float,  float,  double, PyFloat_FromDouble, PyFloat_AsDouble)
CPPYY_IMPL_REFEXEC(Double, double, double, PyFloat_FromDouble, PyFloat_AsDouble)
CPPYY_IMPL_REFEXEC(LongDouble, LongDouble_t, LongDouble_t, PyFloat_FromDouble, PyFloat_AsDouble)

template<typename T>
static inline PyObject* PyComplex_FromComplex(const std::complex<T>& c) {
    return PyComplex_FromDoubles(c.real(), c.imag());
}

template<typename T>
static inline std::complex<T> PyComplex_AsComplex(PyObject* pycplx) {
    Py_complex cplx = PyComplex_AsCComplex(pycplx);
    return std::complex<T>(cplx.real, cplx.imag);
}

CPPYY_IMPL_REFEXEC(ComplexD, std::complex<double>,
    std::complex<double>, PyComplex_FromComplex<double>, PyComplex_AsComplex<double>)


//----------------------------------------------------------------------------
PyObject* CPyCppyy::STLStringRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, return python string return value
    std::string* result = (std::string*)GILCallR(method, self, ctxt);
    if (!fAssignable)
        return CPyCppyy_PyText_FromStringAndSize(result->c_str(), result->size());

    *result = std::string(
        CPyCppyy_PyText_AsString(fAssignable), CPyCppyy_PyText_GET_SIZE(fAssignable));

    Py_DECREF(fAssignable);
    fAssignable = nullptr;

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::VoidExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, return None
    GILCallV(method, self, ctxt);
    Py_RETURN_NONE;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CStringExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python string return value
    char* result = (char*)GILCallR(method, self, ctxt);
    if (!result) {
        Py_INCREF(PyStrings::gEmptyString);
        return PyStrings::gEmptyString;
    }

    return CPyCppyy_PyText_FromString(result);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::WCStringExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python unicode return value
    wchar_t* result = (wchar_t*)GILCallR(method, self, ctxt);
    if (!result) {
        wchar_t w = L'\0';
        return PyUnicode_FromWideChar(&w, 0);
    }

    return PyUnicode_FromWideChar(result, wcslen(result));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CString16Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python unicode return value
    char16_t* result = (char16_t*)GILCallR(method, self, ctxt);
    if (!result) {
        char16_t w = u'\0';
        return PyUnicode_DecodeUTF16((const char*)&w, 0, nullptr, nullptr);
    }

    return PyUnicode_DecodeUTF16((const char*)result,
        std::char_traits<char16_t>::length(result)*sizeof(char16_t), nullptr, nullptr);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CString32Executor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python unicode return value
    char32_t* result = (char32_t*)GILCallR(method, self, ctxt);
    if (!result) {
        char32_t w = U'\0';
        return PyUnicode_DecodeUTF32((const char*)&w, 0, nullptr, nullptr);
    }
 
    return PyUnicode_DecodeUTF32((const char*)result,
        std::char_traits<char32_t>::length(result)*sizeof(char32_t), nullptr, nullptr);
}


//- pointer/array executors --------------------------------------------------
PyObject* CPyCppyy::VoidArrayExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python long return value
    Long_t* result = (Long_t*)GILCallR(method, self, ctxt);
    if (!result) {
        Py_INCREF(gNullPtrObject);
        return gNullPtrObject;
    }
    return CreatePointerView(result);
}

//----------------------------------------------------------------------------
#define CPPYY_IMPL_ARRAY_EXEC(name, type)                                    \
PyObject* CPyCppyy::name##ArrayExecutor::Execute(                            \
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt) \
{                                                                            \
    return CreateLowLevelView((type*)GILCallR(method, self, ctxt));          \
}

CPPYY_IMPL_ARRAY_EXEC(Bool,     bool)
CPPYY_IMPL_ARRAY_EXEC(UChar,    unsigned char)
#if __cplusplus > 201402L
CPPYY_IMPL_ARRAY_EXEC(Byte,     std::byte)
#endif
CPPYY_IMPL_ARRAY_EXEC(Short,    short)
CPPYY_IMPL_ARRAY_EXEC(UShort,   unsigned short)
CPPYY_IMPL_ARRAY_EXEC(Int,      int)
CPPYY_IMPL_ARRAY_EXEC(UInt,     unsigned int)
CPPYY_IMPL_ARRAY_EXEC(Long,     long)
CPPYY_IMPL_ARRAY_EXEC(ULong,    unsigned long)
CPPYY_IMPL_ARRAY_EXEC(LLong,    long long)
CPPYY_IMPL_ARRAY_EXEC(ULLong,   unsigned long long)
CPPYY_IMPL_ARRAY_EXEC(Float,    float)
CPPYY_IMPL_ARRAY_EXEC(Double,   double)
CPPYY_IMPL_ARRAY_EXEC(ComplexF, std::complex<float>)
CPPYY_IMPL_ARRAY_EXEC(ComplexD, std::complex<double>)
CPPYY_IMPL_ARRAY_EXEC(ComplexI, std::complex<int>)
CPPYY_IMPL_ARRAY_EXEC(ComplexL, std::complex<long>)


//- special cases ------------------------------------------------------------
#define CPPYY_COMPLEX_EXEC(code, type)                                       \
PyObject* CPyCppyy::Complex##code##Executor::Execute(                        \
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt) \
{                                                                            \
    static Cppyy::TCppScope_t scopeid = Cppyy::GetScope("std::complex<"#type">");\
    std::complex<type>* result =                                             \
        (std::complex<type>*)GILCallO(method, self, ctxt, scopeid);          \
    if (!result) {                                                           \
        PyErr_SetString(PyExc_ValueError, "NULL result where temporary expected");\
        return nullptr;                                                      \
    }                                                                        \
                                                                             \
    PyObject* pyres = PyComplex_FromDoubles(result->real(), result->imag()); \
    ::operator delete(result); /* Cppyy::CallO calls ::operator new */       \
    return pyres;                                                            \
}

CPPYY_COMPLEX_EXEC(D, double)

//----------------------------------------------------------------------------
PyObject* CPyCppyy::STLStringExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python string return value

// TODO: make use of GILLCallS (?!)
    static Cppyy::TCppScope_t sSTLStringScope = Cppyy::GetScope("std::string");
    std::string* result = (std::string*)GILCallO(method, self, ctxt, sSTLStringScope);
    if (!result) {
        Py_INCREF(PyStrings::gEmptyString);
        return PyStrings::gEmptyString;
    }

    PyObject* pyresult =
        CPyCppyy_PyText_FromStringAndSize(result->c_str(), result->size());
    ::operator delete(result); // calls Cppyy::CallO which calls ::operator new

    return pyresult;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::STLWStringExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python string return value
    static Cppyy::TCppScope_t sSTLWStringScope = Cppyy::GetScope("std::wstring");
    std::wstring* result = (std::wstring*)GILCallO(method, self, ctxt, sSTLWStringScope);
    if (!result) {
        wchar_t w = L'\0';
        return PyUnicode_FromWideChar(&w, 0);
    }

    PyObject* pyresult = PyUnicode_FromWideChar(result->c_str(), result->size());
    ::operator delete(result); // calls Cppyy::CallO which calls ::operator new

    return pyresult;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstancePtrExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python proxy object return value
    return BindCppObject((void*)GILCallR(method, self, ctxt), fClass);
}

//----------------------------------------------------------------------------
CPyCppyy::InstanceExecutor::InstanceExecutor(Cppyy::TCppType_t klass) :
    fClass(klass), fFlags(CPPInstance::kIsValue | CPPInstance::kIsOwner)
{
    /* empty */
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execution will bring a temporary in existence
    Cppyy::TCppObject_t value = GILCallO(method, self, ctxt, fClass);

    if (!value) {
        if (!PyErr_Occurred())         // callee may have set a python error itself
            PyErr_SetString(PyExc_ValueError, "nullptr result where temporary expected");
        return nullptr;
    }

// the result can then be bound
    PyObject* pyobj = BindCppObjectNoCast(value, fClass, fFlags);
    if (!pyobj)
        return nullptr;

// python ref counting will now control this object's life span; it will be
// deleted b/c it is marked as a by-value object (unless C++ ownership is set)
    return pyobj;
}


//----------------------------------------------------------------------------
CPyCppyy::IteratorExecutor::IteratorExecutor(Cppyy::TCppType_t klass) :
    InstanceExecutor(klass)
{
    fFlags = CPPInstance::kNoWrapConv;
}


//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// executor binds the result to the left-hand side, overwriting if an old object
    PyObject* result = BindCppObject((void*)GILCallR(method, self, ctxt), fClass);
    if (!result || !fAssignable)
        return result;
    else {
    // this generic code is quite slow compared to its C++ equivalent ...
        PyObject* assign = PyObject_GetAttr(result, PyStrings::gAssign);
        if (!assign) {
            PyErr_Clear();
            PyObject* descr = PyObject_Str(result);
            if (descr && CPyCppyy_PyText_CheckExact(descr)) {
                PyErr_Format(PyExc_TypeError, "cannot assign to return object (%s)",
                             CPyCppyy_PyText_AsString(descr));
            } else {
                PyErr_SetString(PyExc_TypeError, "cannot assign to result");
            }
            Py_XDECREF(descr);
            Py_DECREF(result);
            Py_DECREF(fAssignable); fAssignable = nullptr;
            return nullptr;
        }

        PyObject* res2 = PyObject_CallFunction(assign, const_cast<char*>("O"), fAssignable);

        Py_DECREF(assign);
        Py_DECREF(result);
        Py_DECREF(fAssignable); fAssignable = nullptr;

        if (res2) {
            Py_DECREF(res2);            // typically, *this from operator=()
            Py_RETURN_NONE;
        }

        return nullptr;
    }
}

//----------------------------------------------------------------------------
static inline PyObject* SetInstanceCheckError(PyObject* pyobj) {
    PyObject* pystr = PyObject_Str(pyobj);
    if (pystr) {
        PyErr_Format(PyExc_TypeError,
           "C++ object expected, got %s", CPyCppyy_PyText_AsString(pystr));
        Py_DECREF(pystr);
    } else
        PyErr_SetString(PyExc_TypeError, "C++ object expected");
    return nullptr;
}

PyObject* CPyCppyy::InstancePtrPtrExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python C++ proxy object
// return ptr value
    if (fAssignable && !CPPInstance_Check(fAssignable))
        return SetInstanceCheckError(fAssignable);

    void** result = (void**)GILCallR(method, self, ctxt);
    if (!fAssignable)
        return BindCppObject((void*)result, fClass,
                             CPPInstance::kIsPtrPtr | CPPInstance::kIsReference);

    CPPInstance* cppinst = (CPPInstance*)fAssignable;
    *result = cppinst->GetObject();

    Py_DECREF(fAssignable);
    fAssignable = nullptr;

    Py_RETURN_NONE;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstancePtrRefExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct python C++ proxy object
// ignoring ref) return ptr value
    if (fAssignable && !CPPInstance_Check(fAssignable))
        return SetInstanceCheckError(fAssignable);

    void** result = (void**)GILCallR(method, self, ctxt);
    if (!fAssignable)
        return BindCppObject(*result, fClass);

    CPPInstance* cppinst = (CPPInstance*)fAssignable;
    *result = cppinst->GetObject();;

    Py_DECREF(fAssignable);
    fAssignable = nullptr;

    Py_RETURN_NONE;
}


//----------------------------------------------------------------------------
PyObject* CPyCppyy::InstanceArrayExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, construct TupleOfInstances from
// return value
    dim_t dims[] = {1, (dim_t)fArraySize};
    return BindCppObjectArray((void*)GILCallR(method, self, ctxt), fClass, dims);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::ConstructorExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t klass, CallContext* ctxt)
{
// package return address in PyObject* for caller to handle appropriately (see
// CPPConstructor for the actual build of the PyObject)
    return (PyObject*)GILCallConstructor(method, (Cppyy::TCppType_t)klass, ctxt);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::PyObjectExecutor::Execute(
    Cppyy::TCppMethod_t method, Cppyy::TCppObject_t self, CallContext* ctxt)
{
// execute <method> with argument <self, ctxt>, return python object
    return (PyObject*)GILCallR(method, self, ctxt);
}


//- factories ----------------------------------------------------------------
CPyCppyy::Executor* CPyCppyy::CreateExecutor(const std::string& fullType)
{
// The matching of the fulltype to an executor factory goes through up to 4 levels:
//   1) full, qualified match
//   2) drop '&' as by ref/full type is often pretty much the same python-wise
//   3) C++ classes, either by ref/ptr or by value
//   4) additional special case for enums
//
// If all fails, void is used, which will cause the return type to be ignored on use

// an exactly matching executor is best
    ExecFactories_t::iterator h = gExecFactories.find(fullType);
    if (h != gExecFactories.end())
        return (h->second)();

// resolve typedefs etc.
    const std::string& resolvedType = Cppyy::ResolveName(fullType);

// a full, qualified matching executor is preferred
    if (resolvedType != fullType) {
         h = gExecFactories.find(resolvedType);
         if (h != gExecFactories.end())
              return (h->second)();
    }

//-- nothing? ok, collect information about the type and possible qualifiers/decorators
    bool isConst = strncmp(resolvedType.c_str(), "const", 5)  == 0;
    const std::string& cpd = Utility::Compound(resolvedType);
    std::string realType = TypeManip::clean_type(resolvedType, false);

// accept unqualified type (as python does not know about qualifiers)
    h = gExecFactories.find(realType + cpd);
    if (h != gExecFactories.end())
        return (h->second)();

// drop const, as that is mostly meaningless to python (with the exception
// of c-strings, but those are specialized in the converter map)
    if (isConst) {
        realType = TypeManip::remove_const(realType);
        h = gExecFactories.find(realType + cpd);
        if (h != gExecFactories.end())
            return (h->second)();
    }

//-- still nothing? try pointer instead of array (for builtins)
    if (cpd == "[]") {
        h = gExecFactories.find(realType + "*");
        if (h != gExecFactories.end())
            return (h->second)();           // TODO: use array size
    }

// C++ classes and special cases
    Executor* result = 0;
    if (Cppyy::TCppType_t klass = Cppyy::GetScope(realType)) {
        if (resolvedType.find("iterator") != std::string::npos || gIteratorTypes.find(fullType) != gIteratorTypes.end()) {
            if (cpd == "")
                return new IteratorExecutor(klass);
        }

        if (cpd == "")
            result = new InstanceExecutor(klass);
        else if (cpd == "&")
            result = new InstanceRefExecutor(klass);
        else if (cpd == "**" || cpd == "*[]" || cpd == "&*")
            result = new InstancePtrPtrExecutor(klass);
        else if (cpd == "*&")
            result = new InstancePtrRefExecutor(klass);
        else if (cpd == "[]") {
            Py_ssize_t asize = Utility::ArraySize(resolvedType);
            if (0 < asize)
                result = new InstanceArrayExecutor(klass, asize);
            else
                result = new InstancePtrRefExecutor(klass);
        } else
            result = new InstancePtrExecutor(klass);
    } else {
    // unknown: void* may work ("user knows best"), void will fail on use of return value
        h = (cpd == "") ? gExecFactories.find("void") : gExecFactories.find("void*");
    }

    if (!result && h != gExecFactories.end())
    // executor factory available, use it to create executor
        result = (h->second)();

   return result;                  // may still be null
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
void CPyCppyy::DestroyExecutor(Executor* p)
{
    if (p && p->HasState())
        delete p;  // state-less executors are always shared
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
bool CPyCppyy::RegisterExecutor(const std::string& name, ef_t fac)
{
// register a custom executor
    auto f = gExecFactories.find(name);
    if (f != gExecFactories.end())
        return false;

    gExecFactories[name] = fac;
    return true;
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
bool CPyCppyy::UnregisterExecutor(const std::string& name)
{
// remove a custom executor
    auto f = gExecFactories.find(name);
    if (f != gExecFactories.end()) {
        gExecFactories.erase(f);
        return true;
    }
    return false;
}

//----------------------------------------------------------------------------
CPYCPPYY_EXPORT
void* CPyCppyy::CallVoidP(Cppyy::TCppMethod_t meth, Cppyy::TCppObject_t obj, CallContext* ctxt)
{
     return GILCallR(meth, obj, ctxt);
}


//----------------------------------------------------------------------------
namespace {

using namespace CPyCppyy;

#define WSTRING "basic_string<wchar_t,char_traits<wchar_t>,allocator<wchar_t> >"

struct InitExecFactories_t {
public:
    InitExecFactories_t() {
    // load all executor factories in the global map 'gExecFactories'
        CPyCppyy::ExecFactories_t& gf = gExecFactories;

    // factories for built-ins
        gf["bool"] =                        (ef_t)+[]() { static BoolExecutor e{};          return &e; };
        gf["bool&"] =                       (ef_t)+[]() { return new BoolRefExecutor{}; };
        gf["const bool&"] =                 (ef_t)+[]() { static BoolConstRefExecutor e{};  return &e; };
        gf["char"] =                        (ef_t)+[]() { static CharExecutor e{};          return &e; };
        gf["signed char"] =                 gf["char"];
        gf["unsigned char"] =               (ef_t)+[]() { static UCharExecutor e{};         return &e; };
        gf["char&"] =                       (ef_t)+[]() { return new CharRefExecutor{}; };
        gf["signed char&"] =                gf["char&"];
        gf["unsigned char&"] =              (ef_t)+[]() { return new UCharRefExecutor{}; };
        gf["const char&"] =                 (ef_t)+[]() { static CharConstRefExecutor e{};  return &e; };
        gf["const signed char&"] =          gf["const char&"];
        gf["const unsigned char&"] =        (ef_t)+[]() { static UCharConstRefExecutor e{}; return &e; };
        gf["wchar_t"] =                     (ef_t)+[]() { static WCharExecutor e{};         return &e; };
        gf["char16_t"] =                    (ef_t)+[]() { static Char16Executor e{};        return &e; };
        gf["char32_t"] =                    (ef_t)+[]() { static Char32Executor e{};        return &e; };
        gf["int8_t"] =                      (ef_t)+[]() { static Int8Executor e{};          return &e; };
        gf["int8_t&"] =                     (ef_t)+[]() { return new Int8RefExecutor{}; };
        gf["const int8_t&"] =               (ef_t)+[]() { static Int8RefExecutor e{};       return &e; };
        gf["uint8_t"] =                     (ef_t)+[]() { static UInt8Executor e{};         return &e; };
        gf["uint8_t&"] =                    (ef_t)+[]() { return new UInt8RefExecutor{}; };
        gf["const uint8_t&"] =              (ef_t)+[]() { static UInt8RefExecutor e{};      return &e; };
        gf["short"] =                       (ef_t)+[]() { static ShortExecutor e{};         return &e; };
        gf["short&"] =                      (ef_t)+[]() { return new ShortRefExecutor{}; };
        gf["int"] =                         (ef_t)+[]() { static IntExecutor e{};           return &e; };
        gf["int&"] =                        (ef_t)+[]() { return new IntRefExecutor{}; };
        gf["unsigned short"] =              gf["int"];
        gf["unsigned short&"] =             (ef_t)+[]() { return new UShortRefExecutor{}; };
        gf["unsigned long"] =               (ef_t)+[]() { static ULongExecutor e{};         return &e; };
        gf["unsigned long&"] =              (ef_t)+[]() { return new ULongRefExecutor{}; };
        gf["unsigned int"] =                gf["unsigned long"];
        gf["unsigned int&"] =               (ef_t)+[]() { return new UIntRefExecutor{}; };
        gf["long"] =                        (ef_t)+[]() { static LongExecutor e{};          return &e; };
        gf["long&"] =                       (ef_t)+[]() { return new LongRefExecutor{}; };
        gf["unsigned long"] =               (ef_t)+[]() { static ULongExecutor e{};         return &e; };
        gf["unsigned long&"] =              (ef_t)+[]() { return new ULongRefExecutor{}; };
        gf["long long"] =                   (ef_t)+[]() { static LongLongExecutor e{};      return &e; };
        gf["long long&"] =                  (ef_t)+[]() { return new LongLongRefExecutor{}; };
        gf["unsigned long long"] =          (ef_t)+[]() { static ULongLongExecutor e{};     return &e; };
        gf["unsigned long long&"] =         (ef_t)+[]() { return new ULongLongRefExecutor{}; };

        gf["float"] =                       (ef_t)+[]() { static FloatExecutor e{};      return &e; };
        gf["float&"] =                      (ef_t)+[]() { return new FloatRefExecutor{}; };
        gf["double"] =                      (ef_t)+[]() { static DoubleExecutor e{};     return &e; };
        gf["double&"] =                     (ef_t)+[]() { return new DoubleRefExecutor{}; };
        gf["long double"] =                 (ef_t)+[]() { static LongDoubleExecutor e{}; return &e; }; // TODO: lost precision
        gf["long double&"] =                (ef_t)+[]() { return new LongDoubleRefExecutor{}; };
        gf["void"] =                        (ef_t)+[]() { static VoidExecutor e{};       return &e; };

    // pointer/array factories
        gf["void*"] =                       (ef_t)+[]() { static VoidArrayExecutor e{};     return &e; };
        gf["bool*"] =                       (ef_t)+[]() { static BoolArrayExecutor e{};     return &e; };
        gf["unsigned char*"] =              (ef_t)+[]() { static UCharArrayExecutor e{};    return &e; };
        gf["const unsigned char*"] =        gf["unsigned char*"];
#if __cplusplus > 201402L
        gf["byte*"] =                       (ef_t)+[]() { static ByteArrayExecutor e{};    return &e; };
        gf["const byte*"] =                 gf["byte*"];
#endif
        gf["short*"] =                      (ef_t)+[]() { static ShortArrayExecutor e{};    return &e; };
        gf["unsigned short*"] =             (ef_t)+[]() { static UShortArrayExecutor e{};   return &e; };
        gf["int*"] =                        (ef_t)+[]() { static IntArrayExecutor e{};      return &e; };
        gf["unsigned int*"] =               (ef_t)+[]() { static UIntArrayExecutor e{};     return &e; };
        gf["long*"] =                       (ef_t)+[]() { static LongArrayExecutor e{};     return &e; };
        gf["unsigned long*"] =              (ef_t)+[]() { static ULongArrayExecutor e{};    return &e; };
        gf["long long*"] =                  (ef_t)+[]() { static LLongArrayExecutor e{};    return &e; };
        gf["unsigned long long*"] =         (ef_t)+[]() { static ULLongArrayExecutor e{};   return &e; };
        gf["float*"] =                      (ef_t)+[]() { static FloatArrayExecutor e{};    return &e; };
        gf["double*"] =                     (ef_t)+[]() { static DoubleArrayExecutor e{};   return &e; };
        gf["complex<float>*"] =             (ef_t)+[]() { static ComplexFArrayExecutor e{}; return &e; };
        gf["complex<double>*"] =            (ef_t)+[]() { static ComplexDArrayExecutor e{}; return &e; };
        gf["complex<int>*"] =               (ef_t)+[]() { static ComplexIArrayExecutor e{}; return &e; };
        gf["complex<long>*"] =              (ef_t)+[]() { static ComplexLArrayExecutor e{}; return &e; };

    // aliases
        gf["internal_enum_type_t"] =        gf["int"];
        gf["internal_enum_type_t&"] =       gf["int&"];
        gf["internal_enum_type_t*"] =       gf["int*"];
#if __cplusplus > 201402L
        gf["byte"] =                        gf["uint8_t"];
        gf["byte&"] =                       gf["uint8_t&"];
        gf["const byte&"] =                 gf["const uint8_t&"];
#endif
        gf["Long64_t"] =                    gf["long long"];
        gf["Long64_t&"] =                   gf["long long&"];
        gf["Long64_t*"] =                   gf["long long*"];
        gf["ULong64_t"] =                   gf["unsigned long long"];
        gf["ULong64_t&"] =                  gf["unsigned long long&"];
        gf["ULong64_t*"] =                  gf["unsigned long long*"];
        gf["Float16_t"] =                   gf["float"];
        gf["Float16_t&"] =                  gf["float&"];
        gf["Double32_t"] =                  gf["double"];
        gf["Double32_t&"] =                 gf["double&"];

    // factories for special cases
        gf["const char*"] =                 (ef_t)+[]() { static CStringExecutor e{};     return &e; };
        gf["char*"] =                       gf["const char*"];
        gf["const signed char*"] =          gf["const char*"];
        gf["signed char*"] =                gf["char*"];
        gf["wchar_t*"] =                    (ef_t)+[]() { static WCStringExecutor e{};    return &e;};
        gf["char16_t*"] =                   (ef_t)+[]() { static CString16Executor e{};   return &e;};
        gf["char32_t*"] =                   (ef_t)+[]() { static CString32Executor e{};   return &e;};
        gf["std::string"] =                 (ef_t)+[]() { static STLStringExecutor e{};   return &e; };
        gf["string"] =                      gf["std::string"];
        gf["std::string&"] =                (ef_t)+[]() { return new STLStringRefExecutor{}; };
        gf["string&"] =                     gf["std::string&"];
        gf["std::wstring"] =                (ef_t)+[]() { static STLWStringExecutor e{};  return &e; };
        gf["std::" WSTRING] =               gf["std::wstring"];
        gf[WSTRING] =                       gf["std::wstring"];
        gf["complex<double>"] =             (ef_t)+[]() { static ComplexDExecutor e{};    return &e; };
        gf["complex<double>&"] =            (ef_t)+[]() { return new ComplexDRefExecutor{}; };
        gf["__init__"] =                    (ef_t)+[]() { static ConstructorExecutor e{}; return &e; };
        gf["PyObject*"] =                   (ef_t)+[]() { static PyObjectExecutor e{};    return &e; };
        gf["_object*"] =                    gf["PyObject*"];
        gf["FILE*"] =                       gf["void*"];
    }
} initExecvFactories_;

} // unnamed namespace
