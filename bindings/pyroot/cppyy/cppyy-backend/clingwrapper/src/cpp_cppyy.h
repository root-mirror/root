#ifndef CPYCPPYY_CPPYY_H
#define CPYCPPYY_CPPYY_H

// Standard
#include <set>
#include <string>
#include <vector>
#include <stddef.h>
#include <stdint.h>


// ROOT types
typedef long long            Long64_t;
typedef unsigned long long   ULong64_t;
typedef long double          LongDouble_t;

namespace Cppyy {
    typedef size_t      TCppScope_t;
    typedef TCppScope_t TCppType_t;
    typedef void*       TCppEnum_t;
    typedef void*       TCppObject_t;
    typedef intptr_t    TCppMethod_t;

    typedef size_t      TCppIndex_t;
    typedef void*       TCppFuncAddr_t;

// direct interpreter access -------------------------------------------------
    RPY_EXPORTED
    bool Compile(const std::string& code);

// name to opaque C++ scope representation -----------------------------------
    RPY_EXPORTED
    std::string ResolveName(const std::string& cppitem_name);
    RPY_EXPORTED
    std::string ResolveEnum(const std::string& enum_type);
    RPY_EXPORTED
    TCppScope_t GetScope(const std::string& scope_name);
    RPY_EXPORTED
    TCppType_t  GetActualClass(TCppType_t klass, TCppObject_t obj);
    RPY_EXPORTED
    size_t      SizeOf(TCppType_t klass);
    RPY_EXPORTED
    size_t      SizeOf(const std::string& type_name);

    RPY_EXPORTED
    bool        IsBuiltin(const std::string& type_name);
    RPY_EXPORTED
    bool        IsComplete(const std::string& type_name);

    RPY_EXPORTED
    TCppScope_t gGlobalScope;      // for fast access

// memory management ---------------------------------------------------------
    RPY_EXPORTED
    TCppObject_t Allocate(TCppType_t type);
    RPY_EXPORTED
    void         Deallocate(TCppType_t type, TCppObject_t instance);
    RPY_EXPORTED
    TCppObject_t Construct(TCppType_t type);
    RPY_EXPORTED
    void         Destruct(TCppType_t type, TCppObject_t instance);

// method/function dispatching -----------------------------------------------
    RPY_EXPORTED
    void          CallV(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    unsigned char CallB(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    char          CallC(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    short         CallH(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    int           CallI(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    long          CallL(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    Long64_t      CallLL(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    float         CallF(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    double        CallD(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    LongDouble_t  CallLD(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    void*         CallR(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    RPY_EXPORTED
    char*         CallS(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args, size_t* length);
    RPY_EXPORTED
    TCppObject_t  CallConstructor(TCppMethod_t method, TCppType_t type, size_t nargs, void* args);
    RPY_EXPORTED
    void          CallDestructor(TCppType_t type, TCppObject_t self);
    RPY_EXPORTED
    TCppObject_t  CallO(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args, TCppType_t result_type);

    RPY_EXPORTED
    TCppFuncAddr_t GetFunctionAddress(TCppMethod_t method, bool check_enabled=true);

// handling of function argument buffer --------------------------------------
    RPY_EXPORTED
    void*  AllocateFunctionArgs(size_t nargs);
    RPY_EXPORTED
    void   DeallocateFunctionArgs(void* args);
    RPY_EXPORTED
    size_t GetFunctionArgSizeof();
    RPY_EXPORTED
    size_t GetFunctionArgTypeoffset();

// scope reflection information ----------------------------------------------
    RPY_EXPORTED
    bool IsNamespace(TCppScope_t scope);
    RPY_EXPORTED
    bool IsTemplate(const std::string& template_name);
    RPY_EXPORTED
    bool IsAbstract(TCppType_t type);
    RPY_EXPORTED
    bool IsEnum(const std::string& type_name);

    RPY_EXPORTED
    void GetAllCppNames(TCppScope_t scope, std::set<std::string>& cppnames);

// namespace reflection information ------------------------------------------
    RPY_EXPORTED
    std::vector<TCppScope_t> GetUsingNamespaces(TCppScope_t);

// class reflection information ----------------------------------------------
    RPY_EXPORTED
    std::string GetFinalName(TCppType_t type);
    RPY_EXPORTED
    std::string GetScopedFinalName(TCppType_t type);
    RPY_EXPORTED
    bool        HasVirtualDestructor(TCppType_t type);
    RPY_EXPORTED
    bool        HasComplexHierarchy(TCppType_t type);
    RPY_EXPORTED
    TCppIndex_t GetNumBases(TCppType_t type);
    RPY_EXPORTED
    std::string GetBaseName(TCppType_t type, TCppIndex_t ibase);
    RPY_EXPORTED
    bool        IsSubtype(TCppType_t derived, TCppType_t base);
    RPY_EXPORTED
    bool        IsSmartPtr(TCppType_t type);
    RPY_EXPORTED
    bool        GetSmartPtrInfo(const std::string&, TCppType_t* raw, TCppMethod_t* deref);
    RPY_EXPORTED
    void        AddSmartPtrType(const std::string&);

// calculate offsets between declared and actual type, up-cast: direction > 0; down-cast: direction < 0
    RPY_EXPORTED
    ptrdiff_t GetBaseOffset(
        TCppType_t derived, TCppType_t base, TCppObject_t address, int direction, bool rerror = false);

// method/function reflection information ------------------------------------
    RPY_EXPORTED
    TCppIndex_t GetNumMethods(TCppScope_t scope);
    RPY_EXPORTED
    std::vector<TCppIndex_t> GetMethodIndicesFromName(TCppScope_t scope, const std::string& name);

    RPY_EXPORTED
    TCppMethod_t GetMethod(TCppScope_t scope, TCppIndex_t imeth);

    RPY_EXPORTED
    std::string GetMethodName(TCppMethod_t);
    RPY_EXPORTED
    std::string GetMethodFullName(TCppMethod_t);
    RPY_EXPORTED
    std::string GetMethodMangledName(TCppMethod_t);
    RPY_EXPORTED
    std::string GetMethodResultType(TCppMethod_t);
    RPY_EXPORTED
    TCppIndex_t GetMethodNumArgs(TCppMethod_t);
    RPY_EXPORTED
    TCppIndex_t GetMethodReqArgs(TCppMethod_t);
    RPY_EXPORTED
    std::string GetMethodArgName(TCppMethod_t, TCppIndex_t iarg);
    RPY_EXPORTED
    std::string GetMethodArgType(TCppMethod_t, TCppIndex_t iarg);
    RPY_EXPORTED
    std::string GetMethodArgDefault(TCppMethod_t, TCppIndex_t iarg);
    RPY_EXPORTED
    std::string GetMethodSignature(TCppMethod_t, bool show_formalargs, TCppIndex_t maxargs = (TCppIndex_t)-1);
    RPY_EXPORTED
    std::string GetMethodPrototype(TCppScope_t scope, TCppMethod_t, bool show_formalargs);
    RPY_EXPORTED
    bool        IsConstMethod(TCppMethod_t);

    RPY_EXPORTED
    TCppIndex_t GetNumTemplatedMethods(TCppScope_t scope);
    RPY_EXPORTED
    std::string GetTemplatedMethodName(TCppScope_t scope, TCppIndex_t imeth);
    RPY_EXPORTED
    bool        IsTemplatedConstructor(TCppScope_t scope, TCppIndex_t imeth);
    RPY_EXPORTED
    bool        ExistsMethodTemplate(TCppScope_t scope, const std::string& name);
    RPY_EXPORTED
    bool        IsMethodTemplate(TCppScope_t scope, TCppIndex_t imeth);
    RPY_EXPORTED
    TCppMethod_t GetMethodTemplate(
        TCppScope_t scope, const std::string& name, const std::string& proto);

    RPY_EXPORTED
    TCppIndex_t  GetGlobalOperator(
        TCppType_t scope, const std::string& lc, const std::string& rc, const std::string& op);

// method properties ---------------------------------------------------------
    RPY_EXPORTED
    bool IsPublicMethod(TCppMethod_t method);
    RPY_EXPORTED
    bool IsProtectedMethod(TCppMethod_t method);
    RPY_EXPORTED
    bool IsConstructor(TCppMethod_t method);
    RPY_EXPORTED
    bool IsDestructor(TCppMethod_t method);
    RPY_EXPORTED
    bool IsStaticMethod(TCppMethod_t method);

// data member reflection information ----------------------------------------
    RPY_EXPORTED
    TCppIndex_t GetNumDatamembers(TCppScope_t scope);
    RPY_EXPORTED
    std::string GetDatamemberName(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    std::string GetDatamemberType(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    intptr_t    GetDatamemberOffset(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    TCppIndex_t GetDatamemberIndex(TCppScope_t scope, const std::string& name);

// data member properties ----------------------------------------------------
    RPY_EXPORTED
    bool IsPublicData(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    bool IsProtectedData(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    bool IsStaticData(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    bool IsConstData(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    bool IsEnumData(TCppScope_t scope, TCppIndex_t idata);
    RPY_EXPORTED
    int  GetDimensionSize(TCppScope_t scope, TCppIndex_t idata, int dimension);

// enum properties -----------------------------------------------------------
    RPY_EXPORTED
    TCppEnum_t  GetEnum(TCppScope_t scope, const std::string& enum_name);
    RPY_EXPORTED
    TCppIndex_t GetNumEnumData(TCppEnum_t);
    RPY_EXPORTED
    std::string GetEnumDataName(TCppEnum_t, TCppIndex_t idata);
    RPY_EXPORTED
    long long   GetEnumDataValue(TCppEnum_t, TCppIndex_t idata);

} // namespace Cppyy

#endif // !CPYCPPYY_CPPYY_H
