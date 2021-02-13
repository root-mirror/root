#ifndef CPYCPPYY_TEMPLATEPROXY_H
#define CPYCPPYY_TEMPLATEPROXY_H

// Bindings
#include "CPPScope.h"
#include "Utility.h"

// Standard
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>


namespace CPyCppyy {

class PyCallable;
class CPPOverload;

/** Template proxy object to return functions and methods
 */

typedef std::pair<uint64_t, CPPOverload*> TP_DispatchEntry_t;
typedef std::map<std::string, std::vector<TP_DispatchEntry_t>> TP_DispatchMap_t;

class TemplateInfo {
public:
    TemplateInfo();
    TemplateInfo(const TemplateInfo&) = delete;
    TemplateInfo& operator=(const TemplateInfo&) = delete;
    ~TemplateInfo();

public:
    PyObject* fCppName;
    PyObject* fPyName;
    PyObject* fPyClass;
    CPPOverload* fNonTemplated;   // holder for non-template overloads
    CPPOverload* fTemplated;      // holder for templated overloads
    CPPOverload* fLowPriority;    // low priority overloads such as void*/void**
    PyObject* fWeakrefList;

    TP_DispatchMap_t fDispatchMap;

    uint64_t fFlags;              // collective for all methods
};

typedef std::shared_ptr<TemplateInfo> TP_TInfo_t;

class TemplateProxy {
private:
    friend TemplateProxy* TemplateProxy_New(
        const std::string& cppname, const std::string& pyname, PyObject* pyclass);
    void Set(const std::string& cppname, const std::string& pyname, PyObject* pyclass);

public:                 // public, as the python C-API works with C structs
    PyObject_HEAD
    PyObject* fSelf;              // must be first (same layout as CPPOverload)
    PyObject* fTemplateArgs;
    PyObject* fWeakrefList;
    TP_TInfo_t fTI;

public:
    void MergeOverload(CPPOverload* mp);
    void AdoptMethod(PyCallable* pc);
    void AdoptTemplate(PyCallable* pc);
    PyObject* Instantiate(const std::string& fname,
        PyObject* tmplArgs, Utility::ArgPreference, int* pcnt = nullptr);

private:                // private, as the python C-API will handle creation
    TemplateProxy() = delete;
};


//- template proxy type and type verification --------------------------------
extern PyTypeObject TemplateProxy_Type;

template<typename T>
inline bool TemplateProxy_Check(T* object)
{
    return object && PyObject_TypeCheck(object, &TemplateProxy_Type);
}

template<typename T>
inline bool TemplateProxy_CheckExact(T* object)
{
    return object && Py_TYPE(object) == &TemplateProxy_Type;
}

//- creation -----------------------------------------------------------------
inline TemplateProxy* TemplateProxy_New(
    const std::string& cppname, const std::string& pyname, PyObject* pyclass)
{
// Create and initialize a new template method proxy for the class.
    if (!CPPScope_Check(pyclass)) return nullptr;

    TemplateProxy* pytmpl =
        (TemplateProxy*)TemplateProxy_Type.tp_new(&TemplateProxy_Type, nullptr, nullptr);
    pytmpl->Set(cppname, pyname, pyclass);
    return pytmpl;
}

} // namespace CPyCppyy

#endif // !CPYCPPYY_TEMPLATEPROXY_H
