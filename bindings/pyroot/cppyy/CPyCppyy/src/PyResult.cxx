// Bindings
#include "CPyCppyy.h"
#define CPYCPPYY_INTERNAL 1
#include "CPyCppyy/PyResult.h"
#undef CPYCPPYY_INTERNAL

#include "CPPInstance.h"

// Standard
#include <stdexcept>


//______________________________________________________________________________
//                        Python expression eval result
//                        =============================
//
// Transport class for bringing objects from python (dynamically typed) to Cling
// (statically typed). It is best to immediately cast a PyResult to the real
// type, either implicitly (for builtin types) or explicitly (through a void*
// cast for pointers to C++ objects).


//- constructors/destructor --------------------------------------------------
CPyCppyy::PyResult::PyResult()
{
// Construct a PyResult object from Py_None.
    Py_INCREF(Py_None);
    fPyObject = Py_None;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::PyResult(PyObject* pyobject)
{
// Construct a PyResult from a python object. The python object may represent
// a C++ object. Steals reference to given python object.
    if (!pyobject) {
        Py_INCREF(Py_None);
        fPyObject = Py_None;
    } else
        fPyObject = pyobject;           // steals reference
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::PyResult(const PyResult& other)
{
// Copy constructor. Applies python object reference counting.
    Py_INCREF(other.fPyObject);
    fPyObject = other.fPyObject;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult& CPyCppyy::PyResult::operator=(const PyResult& other)
{
// Assignment operator. Applies python object reference counting.
    if (this != &other) {
        Py_INCREF(other.fPyObject);
        Py_DECREF(fPyObject);
        fPyObject = other.fPyObject;
    }

    return *this;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::~PyResult()
{
// Destructor. Reference counting for the held python object is in effect.
    Py_DECREF(fPyObject);
}


//- public members -----------------------------------------------------------
CPyCppyy::PyResult::operator char*() const
{
// Cast python return value to C-style string (may fail).
    return (char*)((const char*)*this);
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator const char*() const
{
// Cast python return value to C-style string (may fail).
    if (fPyObject == Py_None)      // for void returns
        return nullptr;

    const char* s = CPyCppyy_PyText_AsString(fPyObject);
    if (PyErr_Occurred()) {
        PyErr_Print();
        return nullptr;
    }

    return s;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator char() const
{
// Cast python return value to C++ char (may fail).
    std::string s = operator const char*();
    if (s.size())
        return s[0];

    return '\0';
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator long() const
{
// Cast python return value to C++ long (may fail).
    long l = PyLong_AsLong(fPyObject);

    if (PyErr_Occurred())
        PyErr_Print();

    return l;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator unsigned long() const
{
// Cast python return value to C++ unsigned long (may fail).
    unsigned long ul = PyLong_AsUnsignedLong(fPyObject);

    if (PyErr_Occurred())
        PyErr_Print();

    return ul;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator double() const
{
// Cast python return value to C++ double (may fail).
    double d = PyFloat_AsDouble(fPyObject);

    if (PyErr_Occurred())
        PyErr_Print();

    return d;
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator void*() const
{
// Cast python return value to C++ object with dictionary (may fail; note that
// you have to use the void* converter, as CINT will not call any other).
    if (fPyObject == Py_None)
        return nullptr;

    if (CPyCppyy::CPPInstance_Check(fPyObject)) {
        ((CPyCppyy::CPPInstance*)fPyObject)->CppOwns();
        return ((CPyCppyy::CPPInstance*)fPyObject)->GetObject();
    } else
        return fPyObject;              // borrows reference
}

//----------------------------------------------------------------------------
CPyCppyy::PyResult::operator PyObject*() const
{
// Direct return of the held PyObject; note the new reference.
    if (fPyObject == Py_None)
        return nullptr;

    Py_INCREF(fPyObject);
    return fPyObject;
}
