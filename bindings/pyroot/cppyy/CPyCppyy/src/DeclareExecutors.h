#ifndef CPYCPPYY_DECLAREEXECUTORS_H
#define CPYCPPYY_DECLAREEXECUTORS_H

// Bindings
#include "Executors.h"
#include "CallContext.h"

#if __cplusplus > 201402L
#include <cstddef>
#endif


namespace CPyCppyy {

namespace {

#define CPPYY_DECL_EXEC(name)                                                \
class name##Executor : public Executor {                                     \
public:                                                                      \
    virtual PyObject* Execute(                                               \
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);             \
}

// executors for built-ins
CPPYY_DECL_EXEC(Bool);
CPPYY_DECL_EXEC(BoolConstRef);
CPPYY_DECL_EXEC(Char);
CPPYY_DECL_EXEC(CharConstRef);
CPPYY_DECL_EXEC(UChar);
CPPYY_DECL_EXEC(UCharConstRef);
CPPYY_DECL_EXEC(WChar);
CPPYY_DECL_EXEC(Char16);
CPPYY_DECL_EXEC(Char32);
CPPYY_DECL_EXEC(Int8);
CPPYY_DECL_EXEC(UInt8);
CPPYY_DECL_EXEC(Short);
CPPYY_DECL_EXEC(Int);
CPPYY_DECL_EXEC(Long);
CPPYY_DECL_EXEC(ULong);
CPPYY_DECL_EXEC(LongLong);
CPPYY_DECL_EXEC(ULongLong);
CPPYY_DECL_EXEC(Float);
CPPYY_DECL_EXEC(Double);
CPPYY_DECL_EXEC(LongDouble);
CPPYY_DECL_EXEC(Void);
CPPYY_DECL_EXEC(CString);
CPPYY_DECL_EXEC(WCString);
CPPYY_DECL_EXEC(CString16);
CPPYY_DECL_EXEC(CString32);

// pointer/array executors
CPPYY_DECL_EXEC(VoidArray);
CPPYY_DECL_EXEC(BoolArray);
CPPYY_DECL_EXEC(UCharArray);
#if __cplusplus > 201402L
CPPYY_DECL_EXEC(ByteArray);
#endif
CPPYY_DECL_EXEC(ShortArray);
CPPYY_DECL_EXEC(UShortArray);
CPPYY_DECL_EXEC(IntArray);
CPPYY_DECL_EXEC(UIntArray);
CPPYY_DECL_EXEC(LongArray);
CPPYY_DECL_EXEC(ULongArray);
CPPYY_DECL_EXEC(LLongArray);
CPPYY_DECL_EXEC(ULLongArray);
CPPYY_DECL_EXEC(FloatArray);
CPPYY_DECL_EXEC(DoubleArray);
CPPYY_DECL_EXEC(ComplexFArray);
CPPYY_DECL_EXEC(ComplexDArray);
CPPYY_DECL_EXEC(ComplexIArray);
CPPYY_DECL_EXEC(ComplexLArray);

// special cases
CPPYY_DECL_EXEC(ComplexD);
CPPYY_DECL_EXEC(STLString);
CPPYY_DECL_EXEC(STLWString);

class InstancePtrExecutor : public Executor {
public:
    InstancePtrExecutor(Cppyy::TCppType_t klass) : fClass(klass) {}
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);
    virtual bool HasState() { return true; }

protected:
    Cppyy::TCppType_t fClass;
};

class InstanceExecutor : public Executor {
public:
    InstanceExecutor(Cppyy::TCppType_t klass);
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);
    virtual bool HasState() { return true; }

protected:
    Cppyy::TCppType_t fClass;
    unsigned int      fFlags;
};

class IteratorExecutor : public InstanceExecutor {
public:
    IteratorExecutor(Cppyy::TCppType_t klass);
};

CPPYY_DECL_EXEC(Constructor);
CPPYY_DECL_EXEC(PyObject);

#define CPPYY_DECL_REFEXEC(name)                                             \
class name##RefExecutor : public RefExecutor {                               \
public:                                                                      \
    virtual PyObject* Execute(                                               \
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);             \
}

CPPYY_DECL_REFEXEC(Bool);
CPPYY_DECL_REFEXEC(Char);
CPPYY_DECL_REFEXEC(UChar);
CPPYY_DECL_REFEXEC(Int8);
CPPYY_DECL_REFEXEC(UInt8);
CPPYY_DECL_REFEXEC(Short);
CPPYY_DECL_REFEXEC(UShort);
CPPYY_DECL_REFEXEC(Int);
CPPYY_DECL_REFEXEC(UInt);
CPPYY_DECL_REFEXEC(Long);
CPPYY_DECL_REFEXEC(ULong);
CPPYY_DECL_REFEXEC(LongLong);
CPPYY_DECL_REFEXEC(ULongLong);
CPPYY_DECL_REFEXEC(Float);
CPPYY_DECL_REFEXEC(Double);
CPPYY_DECL_REFEXEC(LongDouble);
CPPYY_DECL_REFEXEC(ComplexD);
CPPYY_DECL_REFEXEC(STLString);

// special cases
class InstanceRefExecutor : public RefExecutor {
public:
    InstanceRefExecutor(Cppyy::TCppType_t klass) : fClass(klass) {}
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);

protected:
    Cppyy::TCppType_t fClass;
};

class InstancePtrPtrExecutor : public InstanceRefExecutor {
public:
    using InstanceRefExecutor::InstanceRefExecutor;
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);
};

class InstancePtrRefExecutor : public InstanceRefExecutor {
public:
    using InstanceRefExecutor::InstanceRefExecutor;
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);
};

class InstanceArrayExecutor : public InstancePtrExecutor {
public:
    InstanceArrayExecutor(Cppyy::TCppType_t klass, Py_ssize_t array_size)
        : InstancePtrExecutor(klass), fArraySize(array_size) {}
    virtual PyObject* Execute(
        Cppyy::TCppMethod_t, Cppyy::TCppObject_t, CallContext*);

protected:
    Py_ssize_t fArraySize;
};

} // unnamed namespace

} // namespace CPyCppyy

#endif // !CPYCPPYY_DECLAREEXECUTORS_H
