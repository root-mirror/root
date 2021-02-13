#ifndef CPYCPPYY_CONVERTERS_H
#define CPYCPPYY_CONVERTERS_H

// Standard
#include <string>


namespace CPyCppyy {

struct Parameter;
struct CallContext;

class CPYCPPYY_CLASS_EXPORT Converter {
public:
    virtual ~Converter();

public:
    virtual bool SetArg(PyObject*, Parameter&, CallContext* = nullptr) = 0;
    virtual PyObject* FromMemory(void* address);
    virtual bool ToMemory(PyObject* value, void* address);
    virtual bool HasState() { return false; }
};

// create/destroy converter from fully qualified type (public API)
CPYCPPYY_EXPORT Converter* CreateConverter(const std::string& fullType, dims_t dims = nullptr);
CPYCPPYY_EXPORT void DestroyConverter(Converter* p);
typedef Converter* (*cf_t)(dims_t d);
CPYCPPYY_EXPORT bool RegisterConverter(const std::string& name, cf_t fac);
CPYCPPYY_EXPORT bool UnregisterConverter(const std::string& name);


// converters for special cases (only here b/c of external use of StrictInstancePtrConverter)
class VoidArrayConverter : public Converter {
public:
    VoidArrayConverter(bool keepControl = true) { fKeepControl = keepControl; }

public:
    virtual bool SetArg(PyObject*, Parameter&, CallContext* = nullptr);
    virtual PyObject* FromMemory(void* address);
    virtual bool ToMemory(PyObject* value, void* address);
    virtual bool HasState() { return true; }

protected:
    virtual bool GetAddressSpecialCase(PyObject* pyobject, void*& address);
    bool KeepControl() { return fKeepControl; }

private:
    bool fKeepControl;
};

class InstancePtrConverter : public VoidArrayConverter {
public:
    InstancePtrConverter(Cppyy::TCppType_t klass, bool keepControl = false) :
        VoidArrayConverter(keepControl), fClass(klass) {}

public:
    virtual bool SetArg(PyObject*, Parameter&, CallContext* = nullptr);
    virtual PyObject* FromMemory(void* address);
    virtual bool ToMemory(PyObject* value, void* address);

protected:
    Cppyy::TCppType_t fClass;
};

class StrictInstancePtrConverter : public InstancePtrConverter {
public:
    using InstancePtrConverter::InstancePtrConverter;

protected:
    virtual bool GetAddressSpecialCase(PyObject*, void*&) { return false; }
};

} // namespace CPyCppyy

#endif // !CPYCPPYY_CONVERTERS_H
