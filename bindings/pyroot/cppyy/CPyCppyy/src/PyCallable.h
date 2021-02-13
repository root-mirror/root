#ifndef CPYCPPYY_PYCALLABLE_H
#define CPYCPPYY_PYCALLABLE_H

// Bindings
#include "CallContext.h"


namespace CPyCppyy {

class CPPInstance;

class PyCallable {
public:
    virtual ~PyCallable() {}

public:
    virtual PyObject* GetSignature(bool show_formalargs = true) = 0;
    virtual PyObject* GetPrototype(bool show_formalargs = true) = 0;
    virtual PyObject* GetDocString() { return GetPrototype(); }

    virtual int GetPriority() = 0;
    virtual bool IsGreedy() = 0;

    virtual int GetMaxArgs() = 0;
    virtual PyObject* GetCoVarNames() = 0;
    virtual PyObject* GetArgDefault(int /* iarg */) = 0;

    virtual PyObject* GetScopeProxy() = 0;
    virtual Cppyy::TCppFuncAddr_t GetFunctionAddress() = 0;

    virtual PyCallable* Clone() = 0;

public:
    virtual PyObject* Call(
        CPPInstance*& self, PyObject* args, PyObject* kwds, CallContext* ctxt = nullptr) = 0;
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_PYCALLABLE_H
