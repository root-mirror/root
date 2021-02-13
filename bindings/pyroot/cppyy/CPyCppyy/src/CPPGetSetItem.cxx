// Bindings
#include "CPyCppyy.h"
#include "CPPGetSetItem.h"
#include "Executors.h"


//- protected members ---------------------------------------------------------
bool CPyCppyy::CPPSetItem::InitExecutor_(Executor*& executor, CallContext*)
{
// basic call will do
    if (!CPPMethod::InitExecutor_(executor))
        return false;

// check to make sure we're dealing with a RefExecutor
    if (!dynamic_cast<RefExecutor*>(executor)) {
        PyErr_Format(PyExc_NotImplementedError,
            "no __setitem__ handler for return type (%s)",
            this->GetReturnTypeName().c_str());
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPSetItem::PreProcessArgs(
    CPPInstance*& self, PyObject* args, PyObject* kwds)
{
// Prepare executor with a buffer for the return value.
    Py_ssize_t nArgs = PyTuple_GET_SIZE(args);
    if (nArgs <= 1) {
        PyErr_SetString(PyExc_TypeError, "insufficient arguments to __setitem__");
        return nullptr;
    }

// strip the last element of args to be used on return
    ((RefExecutor*)this->GetExecutor())->SetAssignable(PyTuple_GET_ITEM(args, nArgs-1));
    PyObject* subset = PyTuple_GetSlice(args, 0, nArgs-1);

// see whether any of the arguments is a tuple itself
    Py_ssize_t realsize = 0;
    for (Py_ssize_t i = 0; i < nArgs - 1; ++i) {
        PyObject* item = PyTuple_GET_ITEM(subset, i);
        realsize += PyTuple_Check(item) ? PyTuple_GET_SIZE(item) : 1;
    }

// unroll any tuples, if present in the arguments
    PyObject* unrolled = 0;
    if (realsize != nArgs-1) {
        unrolled = PyTuple_New(realsize);

        int current = 0;
        for (int i = 0; i < nArgs - 1; ++i, ++current) {
            PyObject* item = PyTuple_GET_ITEM(subset, i);
            if (PyTuple_Check(item)) {
                for (int j = 0; j < PyTuple_GET_SIZE(item); ++j, ++current) {
                    PyObject* subitem = PyTuple_GET_ITEM(item, j);
                    Py_INCREF(subitem);
                    PyTuple_SET_ITEM(unrolled, current, subitem);
                }
            } else {
                Py_INCREF(item);
                PyTuple_SET_ITEM(unrolled, current, item);
            }
        }
    }

// continue normal method processing
    PyObject* result = CPPMethod::PreProcessArgs(self, unrolled ? unrolled : subset, kwds);

    Py_XDECREF(unrolled);
    Py_DECREF(subset);
    return result;
}


//-----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPGetItem::PreProcessArgs(
    CPPInstance*& self, PyObject* args, PyObject* kwds)
{
// Unroll tuples for call, otherwise just like CPPMethod (this is very similar
// to the code in CPPSetItem above, but subtly different in the details, hence
// not factored out).
    Py_ssize_t nArgs = PyTuple_GET_SIZE(args);

// see whether any of the arguments is a tuple itself
    Py_ssize_t realsize = 0;
    for (Py_ssize_t i = 0; i < nArgs; ++i) {
        PyObject* item = PyTuple_GET_ITEM(args, i);
        realsize += PyTuple_Check(item) ? PyTuple_GET_SIZE(item) : 1;
    }

// unroll any tuples, if present in the arguments
    PyObject* unrolled = 0;
    if (realsize != nArgs-1) {
        unrolled = PyTuple_New(realsize);

        int current = 0;
        for (int i = 0; i < nArgs; ++i, ++current) {
            PyObject* item = PyTuple_GET_ITEM(args, i);
            if (PyTuple_Check(item)) {
                for (int j = 0; j < PyTuple_GET_SIZE(item); ++j, ++current) {
                    PyObject* subitem = PyTuple_GET_ITEM(item, j);
                    Py_INCREF(subitem);
                    PyTuple_SET_ITEM(unrolled, current, subitem);
                }
            } else {
                Py_INCREF(item);
                PyTuple_SET_ITEM(unrolled, current, item);
            }
        }
    }

// continue normal method processing
    PyObject* result = CPPMethod::PreProcessArgs(self, unrolled ? unrolled : args, kwds);

    Py_XDECREF(unrolled);
    return result;
}
