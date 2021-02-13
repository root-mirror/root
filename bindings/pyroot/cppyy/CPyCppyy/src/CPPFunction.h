#ifndef CPYCPPYY_CPPFUNCTION_H
#define CPYCPPYY_CPPFUNCTION_H

// Bindings
#include "CPPMethod.h"


namespace CPyCppyy {

// Wrapper for global free/static C++ functions
class CPPFunction : public CPPMethod {
public:
    using CPPMethod::CPPMethod;

    virtual PyCallable* Clone() { return new CPPFunction(*this); }

    virtual PyObject* Call(
        CPPInstance*&, PyObject* args, PyObject* kwds, CallContext* ctx = nullptr);

protected:
    virtual PyObject* PreProcessArgs(CPPInstance*& self, PyObject* args, PyObject* kwds);
};

// Wrapper for global binary operators that swap arguments
class CPPReverseBinary : public CPPFunction {
public:
    using CPPFunction::CPPFunction;

    virtual PyCallable* Clone() { return new CPPFunction(*this); }

    virtual PyObject* Call(
        CPPInstance*&, PyObject* args, PyObject* kwds, CallContext* ctx = nullptr);

protected:
    virtual PyObject* PreProcessArgs(CPPInstance*& self, PyObject* args, PyObject* kwds);
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_CPPFUNCTION_H
