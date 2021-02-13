#ifndef CPYCPPYY_CPPMETHOD_H
#define CPYCPPYY_CPPMETHOD_H

// Bindings
#include "PyCallable.h"

// Standard
#include <map>
#include <string>
#include <vector>


namespace CPyCppyy {

class Executor;
class Converter;

class CPPMethod : public PyCallable {
public:
    CPPMethod(Cppyy::TCppScope_t scope, Cppyy::TCppMethod_t method);
    CPPMethod(const CPPMethod&);
    CPPMethod& operator=(const CPPMethod&);
    virtual ~CPPMethod();

public:
    virtual PyObject* GetSignature(bool show_formalargs = true);
    virtual PyObject* GetPrototype(bool show_formalargs = true);
    virtual int       GetPriority();
    virtual bool IsGreedy();

    virtual int       GetMaxArgs();
    virtual PyObject* GetCoVarNames();
    virtual PyObject* GetArgDefault(int iarg);
    virtual PyObject* GetScopeProxy();
    virtual Cppyy::TCppFuncAddr_t GetFunctionAddress();

    virtual PyCallable* Clone() { return new CPPMethod(*this); }

public:
    virtual PyObject* Call(
        CPPInstance*& self, PyObject* args, PyObject* kwds, CallContext* ctxt = nullptr);

protected:
    virtual PyObject* PreProcessArgs(CPPInstance*& self, PyObject* args, PyObject* kwds);

    bool      Initialize(CallContext* ctxt = nullptr);
    PyObject* ProcessKeywords(PyObject* self, PyObject* args, PyObject* kwds);
    bool      ConvertAndSetArgs(PyObject* args, CallContext* ctxt = nullptr);
    PyObject* Execute(void* self, ptrdiff_t offset, CallContext* ctxt = nullptr);

    Cppyy::TCppMethod_t GetMethod()   { return fMethod; }
// TODO: the following is a special case to allow shimming of the
// constructor; there's probably a better way ...
    void SetMethod(Cppyy::TCppMethod_t m) { fMethod = m; }
    Cppyy::TCppScope_t  GetScope()    { return fScope; }
    Executor*           GetExecutor() { return fExecutor; }
    std::string         GetSignatureString(bool show_formalargs = true);
    std::string         GetReturnTypeName();

    virtual bool InitExecutor_(Executor*&, CallContext* ctxt = nullptr);

private:
    void Copy_(const CPPMethod&);
    void Destroy_() const;

    PyObject* ExecuteFast(void*, ptrdiff_t, CallContext*);
    PyObject* ExecuteProtected(void*, ptrdiff_t, CallContext*);

    bool InitConverters_();

    void SetPyError_(PyObject* msg);

private:
// representation
    Cppyy::TCppMethod_t fMethod;
    Cppyy::TCppScope_t  fScope;
    Executor*           fExecutor;

// call dispatch buffers
    std::vector<Converter*>     fConverters;
    std::map<std::string, int>* fArgIndices;

protected:
// cached value that doubles as initialized flag (uninitialized if -1)
    int fArgsRequired;
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_CPPMETHOD_H
