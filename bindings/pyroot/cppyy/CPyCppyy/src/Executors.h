#ifndef CPYCPPYY_EXECUTORS_H
#define CPYCPPYY_EXECUTORS_H

// Standard
#include <string>


namespace CPyCppyy {

struct CallContext;

class CPYCPPYY_CLASS_EXPORT Executor {
public:
    virtual ~Executor();
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*) = 0;
    virtual bool HasState() { return false; }
};

// special case needed for CPPSetItem
class RefExecutor : public Executor {
public:
    RefExecutor() : fAssignable(nullptr) {}
    virtual bool SetAssignable(PyObject*);
    virtual bool HasState() { return true; }

protected:
    PyObject* fAssignable;
};

// create/destroy executor from fully qualified type (public API)
CPYCPPYY_EXPORT Executor* CreateExecutor(const std::string& fullType);
CPYCPPYY_EXPORT void DestroyExecutor(Executor* p);
typedef Executor* (*ef_t) ();
CPYCPPYY_EXPORT bool RegisterExecutor(const std::string& name, ef_t fac);
CPYCPPYY_EXPORT bool UnregisterExecutor(const std::string& name);

// helper for the actual call
CPYCPPYY_EXPORT void* CallVoidP(Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);

} // namespace CPyCppyy

#endif // !CPYCPPYY_EXECUTORS_H
