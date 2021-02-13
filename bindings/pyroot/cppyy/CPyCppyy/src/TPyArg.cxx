// Bindings
#include "CPyCppyy.h"
#define CPYCPPYY_INTERNAL 1
#include "CPyCppyy/TPyArg.h"
#undef CPYCPPYY_INTERNAL


//______________________________________________________________________________
//                        Generic wrapper for arguments
//                        =============================
//
// Transport class for bringing C++ values and objects from Cling to Python. It
// provides, from the selected constructor, the proper conversion to a PyObject.
// In principle, there should be no need to use this class directly: it relies
// on implicit conversions.


//- constructor dispatcher ---------------------------------------------------
void TPyArg::CallConstructor(
    PyObject*& pyself, PyObject* pyclass, const std::vector<TPyArg>& args)
{
    int nArgs = (int)args.size();
    PyObject* pyargs = PyTuple_New(nArgs);
    for (int i = 0; i < nArgs; ++i)
        PyTuple_SET_ITEM(pyargs, i, (PyObject*)args[i]);
    pyself = PyObject_Call(pyclass, pyargs, nullptr);
    Py_DECREF(pyargs);
}

//----------------------------------------------------------------------------
void CallConstructor(PyObject*& pyself, PyObject* pyclass)
{
    PyObject* pyargs = PyTuple_New(0);
    pyself = PyObject_Call(pyclass, pyargs, nullptr);
    Py_DECREF(pyargs);
}

//- generic dispatcher -------------------------------------------------------
PyObject* TPyArg::CallMethod(PyObject* pymeth, const std::vector<TPyArg>& args)
{
    int nArgs = (int)args.size();
    PyObject* pyargs = PyTuple_New(nArgs);
    for (int i = 0; i < nArgs; ++i)
        PyTuple_SET_ITEM(pyargs, i, (PyObject*)args[i]);
    PyObject* result = PyObject_Call(pymeth, pyargs, nullptr);
    Py_DECREF(pyargs);
    return result;
}

//- destructor dispatcher ----------------------------------------------------
void TPyArg::CallDestructor(PyObject*& pyself, PyObject*, const std::vector<TPyArg>&)
{
    Py_XDECREF(pyself);       // calls actual dtor if ref-count down to 0
}

//----------------------------------------------------------------------------
void TPyArg::CallDestructor(PyObject*& pyself)
{
    Py_XDECREF(pyself);
}

//- constructors/destructor --------------------------------------------------
TPyArg::TPyArg(PyObject* pyobject)
{
// Construct a TPyArg from a python object.
    Py_XINCREF(pyobject);
    fPyObject = pyobject;
}

//----------------------------------------------------------------------------
TPyArg::TPyArg(int value)
{
// Construct a TPyArg from an integer value.
    fPyObject = PyInt_FromLong(value);
}

//----------------------------------------------------------------------------
TPyArg::TPyArg(long value)
{
// Construct a TPyArg from an integer value.
    fPyObject = PyLong_FromLong(value);
}

//----------------------------------------------------------------------------
TPyArg::TPyArg(double value)
{
// Construct a TPyArg from a double value.
    fPyObject = PyFloat_FromDouble(value);
}

//----------------------------------------------------------------------------
TPyArg::TPyArg(const char* value)
{
// Construct a TPyArg from a C-string.
    fPyObject = CPyCppyy_PyText_FromString(value);
}

//----------------------------------------------------------------------------
TPyArg::TPyArg(const TPyArg& s)
{
// Copy constructor.
    Py_XINCREF(s.fPyObject);
    fPyObject = s.fPyObject;
}

//----------------------------------------------------------------------------
TPyArg& TPyArg::operator=(const TPyArg& s)
{
// Assignment operator.
    if (this != &s) {
        Py_XINCREF(s.fPyObject);
        fPyObject = s.fPyObject;
    }
    return *this;
}

//----------------------------------------------------------------------------
TPyArg::~TPyArg()
{
// Done with held PyObject.
    Py_XDECREF(fPyObject);
    fPyObject = nullptr;
}

//- public members -----------------------------------------------------------
TPyArg::operator PyObject*() const
{
// Extract the python object.
    Py_XINCREF(fPyObject);
    return fPyObject;
}
