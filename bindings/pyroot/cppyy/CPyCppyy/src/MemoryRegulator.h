#ifndef CPYCPPYY_MEMORYREGULATOR_H
#define CPYCPPYY_MEMORYREGULATOR_H

#include <functional>
#include <utility>

namespace CPyCppyy {

class CPPInstance;

typedef std::function<std::pair<bool, bool>(Cppyy::TCppObject_t, Cppyy::TCppType_t)> MemHook_t;

class MemoryRegulator {
private:
    static MemHook_t registerHook, unregisterHook;

public:
    MemoryRegulator();

// callback from C++-side frameworks
    static bool RecursiveRemove(Cppyy::TCppObject_t cppobj, Cppyy::TCppType_t klass);

// called when a new python proxy object is created
    static bool RegisterPyObject(CPPInstance* pyobj, void* cppobj);

// called when a the python proxy object is about to be garbage collected or when it is
// about to delete the proxied C++ object, if owned
    static bool UnregisterPyObject(CPPInstance* pyobj, PyObject* pyclass);

// new reference to python object matching cppobj, or 0 on failure
    static PyObject* RetrievePyObject(Cppyy::TCppObject_t cppobj, PyObject* pyclass);

// set hooks for custom memory regulation
    static void SetRegisterHook(MemHook_t h);
    static void SetUnregisterHook(MemHook_t h);
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_MEMORYREGULATOR_H
