#ifndef CPYCPPYY_CPPYY_H
#define CPYCPPYY_CPPYY_H

// Standard
#include <set>
#include <string>
#include <vector>
#include <stddef.h>
#include <stdint.h>

// import/export (after precommondefs.h from PyPy)
#ifdef _MSC_VER
#define CPPYY_IMPORT extern __declspec(dllimport)
#else
#define CPPYY_IMPORT extern
#endif


namespace Cppyy {

    typedef size_t      TCppScope_t;
    typedef TCppScope_t TCppType_t;
    typedef void*       TCppEnum_t;
    typedef void*       TCppObject_t;
    typedef intptr_t    TCppMethod_t;

    typedef size_t      TCppIndex_t;
    typedef void*       TCppFuncAddr_t;

// direct interpreter access -------------------------------------------------
    CPPYY_IMPORT
    bool Compile(const std::string& code);

// name to opaque C++ scope representation -----------------------------------
    CPPYY_IMPORT
    std::string ResolveName(const std::string& cppitem_name);
    CPPYY_IMPORT
    std::string ResolveEnum(const std::string& enum_type);
    CPPYY_IMPORT
    TCppScope_t GetScope(const std::string& scope_name);
    CPPYY_IMPORT
    TCppType_t  GetActualClass(TCppType_t klass, TCppObject_t obj);
    CPPYY_IMPORT
    size_t      SizeOf(TCppType_t klass);
    CPPYY_IMPORT
    size_t      SizeOf(const std::string& type_name);

    CPPYY_IMPORT
    bool        IsBuiltin(const std::string& type_name);
    CPPYY_IMPORT
    bool        IsComplete(const std::string& type_name);

    CPPYY_IMPORT
    TCppScope_t gGlobalScope;      // for fast access

// memory management ---------------------------------------------------------
    CPPYY_IMPORT
    TCppObject_t Allocate(TCppType_t type);
    CPPYY_IMPORT
    void         Deallocate(TCppType_t type, TCppObject_t instance);
    CPPYY_IMPORT
    TCppObject_t Construct(TCppType_t type);
    CPPYY_IMPORT
    void         Destruct(TCppType_t type, TCppObject_t instance);

// method/function dispatching -----------------------------------------------
    CPPYY_IMPORT
    void          CallV(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    unsigned char CallB(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    char          CallC(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    short         CallH(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    int           CallI(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    long          CallL(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    Long64_t      CallLL(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    float         CallF(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    double        CallD(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    LongDouble_t  CallLD(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    void*         CallR(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args);
    CPPYY_IMPORT
    char*         CallS(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args, size_t* length);
    CPPYY_IMPORT
    TCppObject_t  CallConstructor(TCppMethod_t method, TCppType_t type, size_t nargs, void* args);
    CPPYY_IMPORT
    void          CallDestructor(TCppType_t type, TCppObject_t self);
    CPPYY_IMPORT
    TCppObject_t  CallO(TCppMethod_t method, TCppObject_t self, size_t nargs, void* args, TCppType_t result_type);

    CPPYY_IMPORT
    TCppFuncAddr_t GetFunctionAddress(TCppMethod_t method, bool check_enabled=true);

// handling of function argument buffer --------------------------------------
    CPPYY_IMPORT
    void*  AllocateFunctionArgs(size_t nargs);
    CPPYY_IMPORT
    void   DeallocateFunctionArgs(void* args);
    CPPYY_IMPORT
    size_t GetFunctionArgSizeof();
    CPPYY_IMPORT
    size_t GetFunctionArgTypeoffset();

// scope reflection information ----------------------------------------------
    CPPYY_IMPORT
    bool IsNamespace(TCppScope_t scope);
    CPPYY_IMPORT
    bool IsTemplate(const std::string& template_name);
    CPPYY_IMPORT
    bool IsAbstract(TCppType_t type);
    CPPYY_IMPORT
    bool IsEnum(const std::string& type_name);

    CPPYY_IMPORT
    void GetAllCppNames(TCppScope_t scope, std::set<std::string>& cppnames);

// namespace reflection information ------------------------------------------
    CPPYY_IMPORT
    std::vector<TCppScope_t> GetUsingNamespaces(TCppScope_t);

// class reflection information ----------------------------------------------
    CPPYY_IMPORT
    std::string GetFinalName(TCppType_t type);
    CPPYY_IMPORT
    std::string GetScopedFinalName(TCppType_t type);
    CPPYY_IMPORT
    bool        HasVirtualDestructor(TCppType_t type);
    CPPYY_IMPORT
    bool        HasComplexHierarchy(TCppType_t type);
    CPPYY_IMPORT
    TCppIndex_t GetNumBases(TCppType_t type);
    CPPYY_IMPORT
    std::string GetBaseName(TCppType_t type, TCppIndex_t ibase);
    CPPYY_IMPORT
    bool        IsSubtype(TCppType_t derived, TCppType_t base);
    CPPYY_IMPORT
    bool        IsSmartPtr(TCppType_t type);
    CPPYY_IMPORT
    bool        GetSmartPtrInfo(const std::string&, TCppType_t* raw, TCppMethod_t* deref);
    CPPYY_IMPORT
    void        AddSmartPtrType(const std::string&);

// calculate offsets between declared and actual type, up-cast: direction > 0; down-cast: direction < 0
    CPPYY_IMPORT
    ptrdiff_t GetBaseOffset(
        TCppType_t derived, TCppType_t base, TCppObject_t address, int direction, bool rerror = false);

// method/function reflection information ------------------------------------
    CPPYY_IMPORT
    TCppIndex_t GetNumMethods(TCppScope_t scope);
    CPPYY_IMPORT
    std::vector<TCppIndex_t> GetMethodIndicesFromName(TCppScope_t scope, const std::string& name);

    CPPYY_IMPORT
    TCppMethod_t GetMethod(TCppScope_t scope, TCppIndex_t imeth);

    CPPYY_IMPORT
    std::string GetMethodName(TCppMethod_t);
    CPPYY_IMPORT
    std::string GetMethodFullName(TCppMethod_t);
    CPPYY_IMPORT
    std::string GetMethodMangledName(TCppMethod_t);
    CPPYY_IMPORT
    std::string GetMethodResultType(TCppMethod_t);
    CPPYY_IMPORT
    TCppIndex_t GetMethodNumArgs(TCppMethod_t);
    CPPYY_IMPORT
    TCppIndex_t GetMethodReqArgs(TCppMethod_t);
    CPPYY_IMPORT
    std::string GetMethodArgName(TCppMethod_t, TCppIndex_t iarg);
    CPPYY_IMPORT
    std::string GetMethodArgType(TCppMethod_t, TCppIndex_t iarg);
    CPPYY_IMPORT
    std::string GetMethodArgDefault(TCppMethod_t, TCppIndex_t iarg);
    CPPYY_IMPORT
    std::string GetMethodSignature(TCppMethod_t, bool show_formalargs, TCppIndex_t maxargs = (TCppIndex_t)-1);
    CPPYY_IMPORT
    std::string GetMethodPrototype(TCppScope_t scope, TCppMethod_t, bool show_formalargs);
    CPPYY_IMPORT
    bool        IsConstMethod(TCppMethod_t);

    CPPYY_IMPORT
    TCppIndex_t GetNumTemplatedMethods(TCppScope_t scope);
    CPPYY_IMPORT
    std::string GetTemplatedMethodName(TCppScope_t scope, TCppIndex_t imeth);
    CPPYY_IMPORT
    bool        IsTemplatedConstructor(TCppScope_t scope, TCppIndex_t imeth);
    CPPYY_IMPORT
    bool        ExistsMethodTemplate(TCppScope_t scope, const std::string& name);
    CPPYY_IMPORT
    bool        IsMethodTemplate(TCppScope_t scope, TCppIndex_t imeth);
    CPPYY_IMPORT
    TCppMethod_t GetMethodTemplate(
        TCppScope_t scope, const std::string& name, const std::string& proto);

    CPPYY_IMPORT
    TCppIndex_t  GetGlobalOperator(
        TCppType_t scope, const std::string& lc, const std::string& rc, const std::string& op);

// method properties ---------------------------------------------------------
    CPPYY_IMPORT
    bool IsPublicMethod(TCppMethod_t method);
    CPPYY_IMPORT
    bool IsProtectedMethod(TCppMethod_t method);
    CPPYY_IMPORT
    bool IsConstructor(TCppMethod_t method);
    CPPYY_IMPORT
    bool IsDestructor(TCppMethod_t method);
    CPPYY_IMPORT
    bool IsStaticMethod(TCppMethod_t method);

// data member reflection information ----------------------------------------
    CPPYY_IMPORT
    TCppIndex_t GetNumDatamembers(TCppScope_t scope);
    CPPYY_IMPORT
    std::string GetDatamemberName(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    std::string GetDatamemberType(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    intptr_t    GetDatamemberOffset(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    TCppIndex_t GetDatamemberIndex(TCppScope_t scope, const std::string& name);

// data member properties ----------------------------------------------------
    CPPYY_IMPORT
    bool IsPublicData(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    bool IsProtectedData(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    bool IsStaticData(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    bool IsConstData(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    bool IsEnumData(TCppScope_t scope, TCppIndex_t idata);
    CPPYY_IMPORT
    int  GetDimensionSize(TCppScope_t scope, TCppIndex_t idata, int dimension);

// enum properties -----------------------------------------------------------
    CPPYY_IMPORT
    TCppEnum_t  GetEnum(TCppScope_t scope, const std::string& enum_name);
    CPPYY_IMPORT
    TCppIndex_t GetNumEnumData(TCppEnum_t);
    CPPYY_IMPORT
    std::string GetEnumDataName(TCppEnum_t, TCppIndex_t idata);
    CPPYY_IMPORT
    long long   GetEnumDataValue(TCppEnum_t, TCppIndex_t idata);

} // namespace Cppyy

#endif // !CPYCPPYY_CPPYY_H
