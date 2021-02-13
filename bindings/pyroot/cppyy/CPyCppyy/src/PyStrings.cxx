// Bindings
#include "CPyCppyy.h"
#include "PyStrings.h"


//- data _____________________________________________________________________
PyObject* CPyCppyy::PyStrings::gAssign           = nullptr;
PyObject* CPyCppyy::PyStrings::gBases            = nullptr;
PyObject* CPyCppyy::PyStrings::gBase             = nullptr;
PyObject* CPyCppyy::PyStrings::gCppName          = nullptr;
PyObject* CPyCppyy::PyStrings::gDeref            = nullptr;
PyObject* CPyCppyy::PyStrings::gPreInc           = nullptr;
PyObject* CPyCppyy::PyStrings::gPostInc          = nullptr;
PyObject* CPyCppyy::PyStrings::gDict             = nullptr;
PyObject* CPyCppyy::PyStrings::gEmptyString      = nullptr;
PyObject* CPyCppyy::PyStrings::gEq               = nullptr;
PyObject* CPyCppyy::PyStrings::gFollow           = nullptr;
PyObject* CPyCppyy::PyStrings::gGetItem          = nullptr;
PyObject* CPyCppyy::PyStrings::gGetNoCheck       = nullptr;
PyObject* CPyCppyy::PyStrings::gInit             = nullptr;
PyObject* CPyCppyy::PyStrings::gIter             = nullptr;
PyObject* CPyCppyy::PyStrings::gLen              = nullptr;
PyObject* CPyCppyy::PyStrings::gLifeLine         = nullptr;
PyObject* CPyCppyy::PyStrings::gModule           = nullptr;
PyObject* CPyCppyy::PyStrings::gMRO              = nullptr;
PyObject* CPyCppyy::PyStrings::gName             = nullptr;
PyObject* CPyCppyy::PyStrings::gNe               = nullptr;
PyObject* CPyCppyy::PyStrings::gTypeCode         = nullptr;
PyObject* CPyCppyy::PyStrings::gCTypesType       = nullptr;

PyObject* CPyCppyy::PyStrings::gUnderlying       = nullptr;

PyObject* CPyCppyy::PyStrings::gAdd              = nullptr;
PyObject* CPyCppyy::PyStrings::gSub              = nullptr;
PyObject* CPyCppyy::PyStrings::gMul              = nullptr;
PyObject* CPyCppyy::PyStrings::gDiv              = nullptr;

PyObject* CPyCppyy::PyStrings::gLShift           = nullptr;
PyObject* CPyCppyy::PyStrings::gLShiftC          = nullptr;

PyObject* CPyCppyy::PyStrings::gAt               = nullptr;
PyObject* CPyCppyy::PyStrings::gBegin            = nullptr;
PyObject* CPyCppyy::PyStrings::gEnd              = nullptr;
PyObject* CPyCppyy::PyStrings::gFirst            = nullptr;
PyObject* CPyCppyy::PyStrings::gSecond           = nullptr;
PyObject* CPyCppyy::PyStrings::gSize             = nullptr;
PyObject* CPyCppyy::PyStrings::gTemplate         = nullptr;
PyObject* CPyCppyy::PyStrings::gVectorAt         = nullptr;

PyObject* CPyCppyy::PyStrings::gCppReal          = nullptr;
PyObject* CPyCppyy::PyStrings::gCppImag          = nullptr;

PyObject* CPyCppyy::PyStrings::gThisModule       = nullptr;

PyObject* CPyCppyy::PyStrings::gNoImplicit       = nullptr;
PyObject* CPyCppyy::PyStrings::gDispInit         = nullptr;

PyObject* CPyCppyy::PyStrings::gExPythonize      = nullptr;
PyObject* CPyCppyy::PyStrings::gPythonize        = nullptr;


//-----------------------------------------------------------------------------
#define CPPYY_INITIALIZE_STRING(var, str)                                     \
    if (!(PyStrings::var = CPyCppyy_PyText_InternFromString((char*)#str)))    \
        return false

bool CPyCppyy::CreatePyStrings() {
// Build cache of commonly used python strings (the cache is python intern, so
// all strings are shared python-wide, not just in cppyy).
    CPPYY_INITIALIZE_STRING(gAssign,         __assign__);
    CPPYY_INITIALIZE_STRING(gBases,          __bases__);
    CPPYY_INITIALIZE_STRING(gBase,           __base__);
    CPPYY_INITIALIZE_STRING(gCppName,        __cpp_name__);
    CPPYY_INITIALIZE_STRING(gDeref,          __deref__);
    CPPYY_INITIALIZE_STRING(gPreInc,         __preinc__);
    CPPYY_INITIALIZE_STRING(gPostInc,        __postinc__);
    CPPYY_INITIALIZE_STRING(gDict,           __dict__);
    if (!(PyStrings::gEmptyString = CPyCppyy_PyText_FromString((char*)"")))
        return false;
    CPPYY_INITIALIZE_STRING(gEq,             __eq__);
    CPPYY_INITIALIZE_STRING(gFollow,         __follow__);
    CPPYY_INITIALIZE_STRING(gGetItem,        __getitem__);
    CPPYY_INITIALIZE_STRING(gGetNoCheck,     _getitem__unchecked);
    CPPYY_INITIALIZE_STRING(gInit,           __init__);
    CPPYY_INITIALIZE_STRING(gIter,           __iter__);
    CPPYY_INITIALIZE_STRING(gLen,            __len__);
    CPPYY_INITIALIZE_STRING(gLifeLine,       __lifeline);
    CPPYY_INITIALIZE_STRING(gModule,         __module__);
    CPPYY_INITIALIZE_STRING(gMRO,            __mro__);
    CPPYY_INITIALIZE_STRING(gName,           __name__);
    CPPYY_INITIALIZE_STRING(gNe,             __ne__);
    CPPYY_INITIALIZE_STRING(gTypeCode,       typecode);
    CPPYY_INITIALIZE_STRING(gCTypesType,     _type_);

    CPPYY_INITIALIZE_STRING(gUnderlying,     __underlying);

    CPPYY_INITIALIZE_STRING(gAdd,            __add__);
    CPPYY_INITIALIZE_STRING(gSub,            __sub__);
    CPPYY_INITIALIZE_STRING(gMul,            __mul__);
    CPPYY_INITIALIZE_STRING(gDiv,            CPPYY__div__);

    CPPYY_INITIALIZE_STRING(gLShift,         __lshift__);
    CPPYY_INITIALIZE_STRING(gLShiftC,        __lshiftc__);

    CPPYY_INITIALIZE_STRING(gAt,             at);
    CPPYY_INITIALIZE_STRING(gBegin,          begin);
    CPPYY_INITIALIZE_STRING(gEnd,            end);
    CPPYY_INITIALIZE_STRING(gFirst,          first);
    CPPYY_INITIALIZE_STRING(gSecond,         second);
    CPPYY_INITIALIZE_STRING(gSize,           size);
    CPPYY_INITIALIZE_STRING(gTemplate,       Template);
    CPPYY_INITIALIZE_STRING(gVectorAt,       _vector__at);

    CPPYY_INITIALIZE_STRING(gCppReal,        __cpp_real);
    CPPYY_INITIALIZE_STRING(gCppImag,        __cpp_imag);

    CPPYY_INITIALIZE_STRING(gThisModule,     cppyy);

    CPPYY_INITIALIZE_STRING(gNoImplicit,     __cppyy_no_implicit);
    CPPYY_INITIALIZE_STRING(gDispInit,       _init_dispatchptr);

    CPPYY_INITIALIZE_STRING(gExPythonize,    __cppyy_explicit_pythonize__);
    CPPYY_INITIALIZE_STRING(gPythonize,      __cppyy_pythonize__);

    return true;
}


//-----------------------------------------------------------------------------
PyObject* CPyCppyy::DestroyPyStrings() {
// Remove all cached python strings.
    Py_DECREF(PyStrings::gBases);       PyStrings::gBases       = nullptr;
    Py_DECREF(PyStrings::gBase);        PyStrings::gBase        = nullptr;
    Py_DECREF(PyStrings::gCppName);     PyStrings::gCppName     = nullptr;
    Py_DECREF(PyStrings::gDeref);       PyStrings::gDeref       = nullptr;
    Py_DECREF(PyStrings::gPreInc);      PyStrings::gPreInc      = nullptr;
    Py_DECREF(PyStrings::gPostInc);     PyStrings::gPostInc     = nullptr;
    Py_DECREF(PyStrings::gDict);        PyStrings::gDict        = nullptr;
    Py_DECREF(PyStrings::gEmptyString); PyStrings::gEmptyString = nullptr;
    Py_DECREF(PyStrings::gEq);          PyStrings::gEq          = nullptr;
    Py_DECREF(PyStrings::gFollow);      PyStrings::gFollow      = nullptr;
    Py_DECREF(PyStrings::gGetItem);     PyStrings::gGetItem     = nullptr;
    Py_DECREF(PyStrings::gGetNoCheck);  PyStrings::gGetNoCheck  = nullptr;
    Py_DECREF(PyStrings::gInit);        PyStrings::gInit        = nullptr;
    Py_DECREF(PyStrings::gIter);        PyStrings::gIter        = nullptr;
    Py_DECREF(PyStrings::gLen);         PyStrings::gLen         = nullptr;
    Py_DECREF(PyStrings::gLifeLine);    PyStrings::gLifeLine    = nullptr;
    Py_DECREF(PyStrings::gModule);      PyStrings::gModule      = nullptr;
    Py_DECREF(PyStrings::gMRO);         PyStrings::gMRO         = nullptr;
    Py_DECREF(PyStrings::gName);        PyStrings::gName        = nullptr;
    Py_DECREF(PyStrings::gNe);          PyStrings::gNe          = nullptr;
    Py_DECREF(PyStrings::gTypeCode);    PyStrings::gTypeCode    = nullptr;
    Py_DECREF(PyStrings::gCTypesType);  PyStrings::gCTypesType  = nullptr;

    Py_DECREF(PyStrings::gUnderlying);  PyStrings::gUnderlying  = nullptr;

    Py_DECREF(PyStrings::gAdd);         PyStrings::gAdd         = nullptr;
    Py_DECREF(PyStrings::gSub);         PyStrings::gSub         = nullptr;
    Py_DECREF(PyStrings::gMul);         PyStrings::gMul         = nullptr;
    Py_DECREF(PyStrings::gDiv);         PyStrings::gDiv         = nullptr;

    Py_DECREF(PyStrings::gLShift);      PyStrings::gLShift      = nullptr;
    Py_DECREF(PyStrings::gLShiftC);     PyStrings::gLShiftC     = nullptr;

    Py_DECREF(PyStrings::gAt);          PyStrings::gAt          = nullptr;
    Py_DECREF(PyStrings::gBegin);       PyStrings::gBegin       = nullptr;
    Py_DECREF(PyStrings::gEnd);         PyStrings::gEnd         = nullptr;
    Py_DECREF(PyStrings::gFirst);       PyStrings::gFirst       = nullptr;
    Py_DECREF(PyStrings::gSecond);      PyStrings::gSecond      = nullptr;
    Py_DECREF(PyStrings::gSize);        PyStrings::gSize        = nullptr;
    Py_DECREF(PyStrings::gTemplate);    PyStrings::gTemplate    = nullptr;
    Py_DECREF(PyStrings::gVectorAt);    PyStrings::gVectorAt    = nullptr;

    Py_DECREF(PyStrings::gCppReal);     PyStrings::gCppReal     = nullptr;
    Py_DECREF(PyStrings::gCppImag);     PyStrings::gCppImag     = nullptr;

    Py_DECREF(PyStrings::gThisModule);  PyStrings::gThisModule  = nullptr;

    Py_DECREF(PyStrings::gNoImplicit);  PyStrings::gNoImplicit  = nullptr;
    Py_DECREF(PyStrings::gDispInit);    PyStrings::gDispInit    = nullptr;

    Py_DECREF(PyStrings::gExPythonize); PyStrings::gExPythonize = nullptr;
    Py_DECREF(PyStrings::gPythonize);   PyStrings::gPythonize   = nullptr;

    Py_RETURN_NONE;
}
