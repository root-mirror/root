#ifndef CPYCPPYY_CPPDATAMEMBER_H
#define CPYCPPYY_CPPDATAMEMBER_H

// Bindings
#include "Converters.h"

// Standard
#include <string>


namespace CPyCppyy {

class CPPInstance;

class CPPDataMember {
public:
    void Set(Cppyy::TCppScope_t scope, Cppyy::TCppIndex_t idata);
    void Set(Cppyy::TCppScope_t scope, const std::string& name, void* address);

    std::string GetName() { return CPyCppyy_PyText_AsString(fName); }
    void* GetAddress(CPPInstance* pyobj /* owner */);

public:                 // public, as the python C-API works with C structs
    PyObject_HEAD
    intptr_t           fOffset;
    Long_t             fFlags;
    Converter*         fConverter;
    Cppyy::TCppScope_t fEnclosingScope;
    PyObject*          fName;

private:                // private, as the python C-API will handle creation
    CPPDataMember() = delete;
};


//- property proxy for C++ data members, type and type verification ----------
extern PyTypeObject CPPDataMember_Type;

template<typename T>
inline bool CPPDataMember_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &CPPDataMember_Type);
}

template<typename T>
inline bool CPPDataMember_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &CPPDataMember_Type;
}

//- creation -----------------------------------------------------------------
inline CPPDataMember* CPPDataMember_New(
    Cppyy::TCppScope_t scope, Cppyy::TCppIndex_t idata)
{
// Create an initialize a new property descriptor, given the C++ datum.
    CPPDataMember* pyprop =
        (CPPDataMember*)CPPDataMember_Type.tp_new(&CPPDataMember_Type, nullptr, nullptr);
    pyprop->Set(scope, idata);
    return pyprop;
}

inline CPPDataMember* CPPDataMember_NewConstant(
    Cppyy::TCppScope_t scope, const std::string& name, void* address)
{
// Create an initialize a new property descriptor, given the C++ datum.
    CPPDataMember* pyprop =
        (CPPDataMember*)CPPDataMember_Type.tp_new(&CPPDataMember_Type, nullptr, nullptr);
    pyprop->Set(scope, name, address);
    return pyprop;
}

} // namespace CPyCppyy

#endif // !CPYCPPYY_CPPDATAMEMBER_H
