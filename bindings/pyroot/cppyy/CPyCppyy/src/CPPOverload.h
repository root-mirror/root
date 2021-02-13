#ifndef CPYCPPYY_CPPOVERLOAD_H
#define CPYCPPYY_CPPOVERLOAD_H

// Bindings
#include "PyCallable.h"

// Standard
#include <string>
#include <utility>
#include <vector>


namespace CPyCppyy {

// signature hashes are also used by TemplateProxy
inline uint64_t HashSignature(PyObject* args)
{
// Build a hash from the types of the given python function arguments.
    uint64_t hash = 0;

    int nargs = (int)PyTuple_GET_SIZE(args);
    for (int i = 0; i < nargs; ++i) {
    // TODO: hashing in the ref-count is for moves; resolve this together with the
    // improved overloads for implicit conversions
        PyObject* pyobj = PyTuple_GET_ITEM(args, i);
        hash += (uint64_t)Py_TYPE(pyobj);
        hash += (uint64_t)(pyobj->ob_refcnt == 1 ? 1 : 0);
        hash += (hash << 10); hash ^= (hash >> 6);
    }

    hash += (hash << 3); hash ^= (hash >> 11); hash += (hash << 15);

    return hash;
}

class CPPOverload {
public:
    typedef std::vector<std::pair<uint64_t, PyCallable*>> DispatchMap_t;
    typedef std::vector<PyCallable*> Methods_t;

    struct MethodInfo_t {
        MethodInfo_t() : fFlags(CallContext::kNone) { fRefCount = new int(1); }
        ~MethodInfo_t();

        std::string                 fName;
        CPPOverload::DispatchMap_t  fDispatchMap;
        CPPOverload::Methods_t      fMethods;
        uint64_t                    fFlags;

        int* fRefCount;

    private:
        MethodInfo_t(const MethodInfo_t&) = delete;
        MethodInfo_t& operator=(const MethodInfo_t&) = delete;
    };

public:
    void Set(const std::string& name, std::vector<PyCallable*>& methods);
    void AdoptMethod(PyCallable* pc);
    void MergeOverload(CPPOverload* meth);

    const std::string& GetName() const { return fMethodInfo->fName; }
    bool HasMethods() const { return !fMethodInfo->fMethods.empty(); }

public:                 // public, as the python C-API works with C structs
    PyObject_HEAD
    CPPInstance*   fSelf;         // must be first (same layout as TemplateProxy)
    MethodInfo_t*  fMethodInfo;

private:
    CPPOverload() = delete;
};


//- method proxy type and type verification ----------------------------------
CPYCPPYY_IMPORT PyTypeObject CPPOverload_Type;

template<typename T>
inline bool CPPOverload_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &CPPOverload_Type);
}

template<typename T>
inline bool CPPOverload_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &CPPOverload_Type;
}

//- creation -----------------------------------------------------------------
inline CPPOverload* CPPOverload_New(
    const std::string& name, std::vector<PyCallable*>& methods)
{
// Create and initialize a new method proxy from the overloads.
    CPPOverload* pymeth = (CPPOverload*)CPPOverload_Type.tp_new(&CPPOverload_Type, nullptr, nullptr);
    pymeth->Set(name, methods);
    return pymeth;
}

inline CPPOverload* CPPOverload_New(const std::string& name, PyCallable* method)
{
// Create and initialize a new method proxy from the method.
    std::vector<PyCallable*> p;
    p.push_back(method);
    return CPPOverload_New(name, p);
}

} // namespace CPyCppyy

#endif // !CPYCPPYY_CPPOVERLOAD_H
