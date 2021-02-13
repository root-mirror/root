#ifndef CPYCPPYY_DISPATCHPTR_H
#define CPYCPPYY_DISPATCHPTR_H

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// DispatchPtr                                                              //
//                                                                          //
// Smart pointer for reference management and C++ instance tracking when    //
// cross-inheriting. The carried pointer is always expected to be derived   //
// from CPPInstance, and the DispatchPtr to be embedded in the C++ instance //
// derived dispatcher to which it points (ownership is two-way; life-times  //
// are equal). The C++ dispatcher then uses the DispatchPtr to call Python  //
// functions for virtual methods.                                           //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

// Bindings
#include "CPyCppyy/CommonDefs.h"


namespace CPyCppyy {

class CPYCPPYY_CLASS_EXTERN DispatchPtr {
public:
// Default constructor: only ever created from C++, as by definition, creation
// from the Python side makes the relevant Python instance available. Calls to
// the default ctor happen eg. in STL containers. It is expected that the
// pointer to the Python object is filled in later, eg. through assign().
    DispatchPtr() : fPyHardRef(nullptr), fPyWeakRef(nullptr) {}

// Conversion constructor: called with C++ object construction when the PyObject
// is known (eg. when instantiating from Python), with pyobj the Python-side
// representation of the C++ object.
    explicit DispatchPtr(PyObject* pyobj);

// Copy constructor: only ever called from C++. The Python object needs to be
// copied, in case it has added state, and rebound to the new C++ instance.
    DispatchPtr(const DispatchPtr& other, void* cppinst);

// Assignment: only ever called from C++. Similarly to the copy constructor, the
// Pythonb object needs to be copied and rebound.
    DispatchPtr& assign(const DispatchPtr& other, void* cppinst);

// Do not otherwise allow straight copies/assignment.
    DispatchPtr(DispatchPtr&)  = delete;
    DispatchPtr(DispatchPtr&&) = delete;
    DispatchPtr& operator=(const DispatchPtr& other) = delete;

// lifetime is directly bound to the lifetime of the dispatcher object
    ~DispatchPtr() {
        Py_XDECREF(fPyWeakRef);
        Py_XDECREF(fPyHardRef);
    }

// either C++ owns the Python object through a reference count (on fPyHardRef) or
// Python owns the C++ object and we only have a weak reference (through fPyWeakRef)
    void PythonOwns();
    void CppOwns();

// access to underlying object: cast and dereferencing
    operator PyObject*() const {
        return Get();
    }

    PyObject* operator->() const {
        return Get();
    }

private:
    PyObject* Get() const;

private:
    PyObject* fPyHardRef;
    PyObject* fPyWeakRef;
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_DISPATCHPTR_H
