// Bindings
#include "CPyCppyy.h"
#include "ProxyWrappers.h"
#include "CPPClassMethod.h"
#include "CPPConstructor.h"
#include "CPPDataMember.h"
#include "CPPExcInstance.h"
#include "CPPFunction.h"
#include "CPPGetSetItem.h"
#include "CPPInstance.h"
#include "CPPMethod.h"
#include "CPPOverload.h"
#include "CPPScope.h"
#include "MemoryRegulator.h"
#include "PyStrings.h"
#include "Pythonize.h"
#include "TemplateProxy.h"
#include "TupleOfInstances.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <algorithm>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <vector>


//- data _______________________________________________________________________
namespace CPyCppyy {
    extern PyObject* gThisModule;
    extern PyObject* gPyTypeMap;
    extern std::set<Cppyy::TCppType_t> gPinnedTypes;
}

// to prevent having to walk scopes, track python classes by C++ class
typedef std::map<Cppyy::TCppScope_t, PyObject*> PyClassMap_t;
static PyClassMap_t gPyClasses;


//- helpers --------------------------------------------------------------------

namespace CPyCppyy {

typedef struct {
    PyObject_HEAD
    PyObject *dict;
} proxyobject;

// helper for creating new C++ proxy python types
static PyObject* CreateNewCppProxyClass(Cppyy::TCppScope_t klass, PyObject* pybases)
{
// Create a new python shadow class with the required hierarchy and meta-classes.
    PyObject* pymetabases = PyTuple_New(PyTuple_GET_SIZE(pybases));
    for (int i = 0; i < PyTuple_GET_SIZE(pybases); ++i) {
        PyObject* btype = (PyObject*)Py_TYPE(PyTuple_GetItem(pybases, i));
        Py_INCREF(btype);
        PyTuple_SET_ITEM(pymetabases, i, btype);
    }

    std::string name = Cppyy::GetFinalName(klass);

// create meta-class, add a dummy __module__ to pre-empt the default setting
    PyObject* args = Py_BuildValue((char*)"sO{}", (name+"_meta").c_str(), pymetabases);
    PyDict_SetItem(PyTuple_GET_ITEM(args, 2), PyStrings::gModule, Py_True);
    Py_DECREF(pymetabases);

    PyObject* pymeta = (PyObject*)CPPScopeMeta_New(klass, args);
    Py_DECREF(args);
    if (!pymeta) {
        PyErr_Print();
        return nullptr;
    }

// alright, and now we really badly want to get rid of the dummy ...
    PyObject* dictproxy = PyObject_GetAttr(pymeta, PyStrings::gDict);
    PyDict_DelItem(((proxyobject*)dictproxy)->dict, PyStrings::gModule);

// create actual class
    args = Py_BuildValue((char*)"sO{}", name.c_str(), pybases);
    PyObject* pyclass =
        ((PyTypeObject*)pymeta)->tp_new((PyTypeObject*)pymeta, args, nullptr);

    Py_DECREF(args);
    Py_DECREF(pymeta);

    return pyclass;
}

static inline
void AddPropertyToClass(PyObject* pyclass,
    Cppyy::TCppScope_t scope, Cppyy::TCppIndex_t idata)
{
    CPyCppyy::CPPDataMember* property = CPyCppyy::CPPDataMember_New(scope, idata);
    PyObject* pname = CPyCppyy_PyText_InternFromString(const_cast<char*>(property->GetName().c_str()));

// allow access at the instance level
    PyType_Type.tp_setattro(pyclass, pname, (PyObject*)property);

// allow access at the class level (always add after setting instance level)
    if (Cppyy::IsStaticData(scope, idata))
        PyType_Type.tp_setattro((PyObject*)Py_TYPE(pyclass), pname, (PyObject*)property);

// cleanup
    Py_DECREF(pname);
    Py_DECREF(property);
}

static inline
void AddScopeToParent(PyObject* parent, const std::string& name, PyObject* newscope)
{
    PyObject* pyname = CPyCppyy_PyText_InternFromString((char*)name.c_str());
    if (CPPScope_Check(parent)) PyType_Type.tp_setattro(parent, pyname, newscope);
    else PyObject_SetAttr(parent, pyname, newscope);
    Py_DECREF(pyname);
}

} // namespace CPyCppyy


//- public functions ---------------------------------------------------------
namespace CPyCppyy {

static inline void sync_templates(
    PyObject* pyclass, const std::string& mtCppName, const std::string& mtName)
{
    PyObject* dct = PyObject_GetAttr(pyclass, PyStrings::gDict);
    PyObject* pyname = CPyCppyy_PyText_InternFromString(const_cast<char*>(mtName.c_str()));
    PyObject* attr = PyObject_GetItem(dct, pyname);
    if (!attr) PyErr_Clear();
    Py_DECREF(dct);
    if (!TemplateProxy_Check(attr)) {
        TemplateProxy* pytmpl = TemplateProxy_New(mtCppName, mtName, pyclass);
        if (CPPOverload_Check(attr)) pytmpl->MergeOverload((CPPOverload*)attr);
        PyType_Type.tp_setattro(pyclass, pyname, (PyObject*)pytmpl);
        Py_DECREF(pytmpl);
    }
    Py_XDECREF(attr);
    Py_DECREF(pyname);
}

static int BuildScopeProxyDict(Cppyy::TCppScope_t scope, PyObject* pyclass)
{
// Collect methods and data for the given scope, and add them to the given python
// proxy object.

// some properties that'll affect building the dictionary
    bool isNamespace = Cppyy::IsNamespace(scope);
    bool isAbstract  = Cppyy::IsAbstract(scope);
    bool hasConstructor = false;
    Cppyy::TCppMethod_t potGetItem = (Cppyy::TCppMethod_t)0;

// load all public methods and data members
    typedef std::vector<PyCallable*> Callables_t;
    typedef std::map<std::string, Callables_t> CallableCache_t;
    CallableCache_t cache;

// bypass custom __getattr__ for efficiency
    getattrofunc oldgetattro = Py_TYPE(pyclass)->tp_getattro;
    Py_TYPE(pyclass)->tp_getattro = PyType_Type.tp_getattro;

// functions in namespaces are properly found through lazy lookup, so do not
// create them until needed (the same is not true for data members)
    const Cppyy::TCppIndex_t nMethods = isNamespace ? 0 : Cppyy::GetNumMethods(scope);
    for (Cppyy::TCppIndex_t imeth = 0; imeth < nMethods; ++imeth) {
        Cppyy::TCppMethod_t method = Cppyy::GetMethod(scope, imeth);

    // process the method based on its name
        std::string mtCppName = Cppyy::GetMethodName(method);

    // special case trackers
        bool setupSetItem = false;
        bool isConstructor = Cppyy::IsConstructor(method);
        bool isTemplate = isConstructor ? false : Cppyy::IsMethodTemplate(scope, imeth);

    // filter empty names (happens for namespaces, is bug?)
        if (mtCppName == "")
            continue;

    // filter C++ destructors
        if (mtCppName[0] == '~')
            continue;

    // translate operators
        std::string mtName = Utility::MapOperatorName(mtCppName, Cppyy::GetMethodNumArgs(method));
        if (mtName.empty())
            continue;

    // operator[]/() returning a reference type will be used for __setitem__
        bool isCall = mtName == "__call__";
        if (isCall || mtName == "__getitem__") {
            const std::string& qual_return = Cppyy::ResolveName(Cppyy::GetMethodResultType(method));
            const std::string& cpd = Utility::Compound(qual_return);
            if (!cpd.empty() && cpd[cpd.size()- 1] == '&' && \
                    qual_return.find("const", 0, 5) == std::string::npos) {
                if (isCall && !potGetItem) potGetItem = method;
                setupSetItem = true;     // will add methods as overloads
            } else if (isCall) {
            // not a non-const by-ref return, thus better __getitem__ candidate
                potGetItem = method;
            }
        }

    // do not expose private methods as the Cling wrappers for them won't compile
        if (!Cppyy::IsPublicMethod(method))
            continue;

    // template members; handled by adding a dispatcher to the class
        bool storeOnTemplate =
            isTemplate ? true : (!isConstructor && Cppyy::ExistsMethodTemplate(scope, mtCppName));
        if (storeOnTemplate) {
            sync_templates(pyclass, mtCppName, mtName);
        // continue processing to actually add the method so that the proxy can find
        // it on the class when called explicitly
        }

    // construct the holder
        PyCallable* pycall = nullptr;
        if (Cppyy::IsStaticMethod(method))  // class method
            pycall = new CPPClassMethod(scope, method);
        else if (isNamespace)               // free function
            pycall = new CPPFunction(scope, method);
        else if (isConstructor) {           // ctor
            mtName = "__init__";
            hasConstructor = true;
            if (!isAbstract)
                pycall = new CPPConstructor(scope, method);
            else
                pycall = new CPPAbstractClassConstructor(scope, method);
        } else                               // member function
            pycall = new CPPMethod(scope, method);

        if (storeOnTemplate) {
        // template proxy was already created in sync_templates call above, so
        // add only here, not to the cache of collected methods
            PyObject* attr = PyObject_GetAttrString(pyclass, const_cast<char*>(mtName.c_str()));
            if (isTemplate) ((TemplateProxy*)attr)->AdoptTemplate(pycall);
            else ((TemplateProxy*)attr)->AdoptMethod(pycall);
            Py_DECREF(attr);

        // for operator[]/() that returns by ref, also add __setitem__
            if (setupSetItem) {
                TemplateProxy* pysi = (TemplateProxy*)PyObject_GetAttrString(pyclass, const_cast<char*>("__setitem__"));
                if (!pysi) {
                     pysi = TemplateProxy_New(mtCppName, "__setitem__", pyclass);
                     PyObject_SetAttrString(pyclass, const_cast<char*>("__setitem__"), (PyObject*)pysi);
                }
                if (isTemplate) pysi->AdoptTemplate(new CPPSetItem(scope, method));
                else pysi->AdoptMethod(new CPPSetItem(scope, method));
                Py_XDECREF(pysi);
            }

        } else {
        // lookup method dispatcher and store method
            Callables_t& md = (*(cache.insert(
                std::make_pair(mtName, Callables_t())).first)).second;
            md.push_back(pycall);

        // special case for operator[]/() that returns by ref, use for getitem/call and setitem
            if (setupSetItem) {
                Callables_t& setitem = (*(cache.insert(
                    std::make_pair(std::string("__setitem__"), Callables_t())).first)).second;
                setitem.push_back(new CPPSetItem(scope, method));
            }
        }
    }

// add proxies for un-instantiated/non-overloaded templated methods
    const Cppyy::TCppIndex_t nTemplMethods = isNamespace ? 0 : Cppyy::GetNumTemplatedMethods(scope);
    for (Cppyy::TCppIndex_t imeth = 0; imeth < nTemplMethods; ++imeth) {
        const std::string mtCppName = Cppyy::GetTemplatedMethodName(scope, imeth);
    // the number of arguments isn't known until instantiation and as far as C++ is concerned, all
    // same-named operators are simply overloads; so will pre-emptively add both names if with and
    // without arguments differ, letting the normal overload mechanism resolve on call
        bool isConstructor = Cppyy::IsTemplatedConstructor(scope, imeth);

    // first add with no arguments
        std::string mtName0 = isConstructor ? "__init__" : Utility::MapOperatorName(mtCppName, false);
        sync_templates(pyclass, mtCppName, mtName0);

    // then add when taking arguments, if this method is different
        if (!isConstructor) {
            std::string mtName1 = Utility::MapOperatorName(mtCppName, true);
            if (mtName0 != mtName1)
                sync_templates(pyclass, mtCppName, mtName1);
        }
    }

// add a pseudo-default ctor, if none defined
    if (!hasConstructor) {
        PyCallable* defctor = nullptr;
        if (isAbstract)
            defctor = new CPPAbstractClassConstructor(scope, (Cppyy::TCppMethod_t)0);
        else if (isNamespace)
            defctor = new CPPNamespaceConstructor(scope, (Cppyy::TCppMethod_t)0);
        else if (!Cppyy::IsComplete(Cppyy::GetScopedFinalName(scope))) {
            ((CPPScope*)pyclass)->fFlags |= CPPScope::kIsInComplete;
            defctor = new CPPIncompleteClassConstructor(scope, (Cppyy::TCppMethod_t)0);
        } else
            defctor = new CPPConstructor(scope, (Cppyy::TCppMethod_t)0);
        cache["__init__"].push_back(defctor);
    }

// map __call__ to __getitem__ if also mapped to __setitem__
    if (potGetItem) {
        Callables_t& getitem = (*(cache.insert(
           std::make_pair(std::string("__getitem__"), Callables_t())).first)).second;
        getitem.push_back(new CPPGetItem(scope, potGetItem));
    }

// add the methods to the class dictionary
    PyObject* dct = PyObject_GetAttr(pyclass, PyStrings::gDict);
    for (CallableCache_t::iterator imd = cache.begin(); imd != cache.end(); ++imd) {
    // in order to prevent removing templated editions of this method (which were set earlier,
    // above, as a different proxy object), we'll check and add this method flagged as a generic
    // one (to be picked up by the templated one as appropriate) if a template exists
        PyObject* pyname = CPyCppyy_PyText_FromString(const_cast<char*>(imd->first.c_str()));
        PyObject* attr = PyObject_GetItem(dct, pyname);
        Py_DECREF(pyname);
        if (TemplateProxy_Check(attr)) {
        // template exists, supply it with the non-templated method overloads
            for (auto cit : imd->second)
                ((TemplateProxy*)attr)->AdoptMethod(cit);
        } else {
            if (!attr) PyErr_Clear();
        // normal case, add a new method
            CPPOverload* method = CPPOverload_New(imd->first, imd->second);
            PyObject* pymname = CPyCppyy_PyText_InternFromString(const_cast<char*>(method->GetName().c_str()));
            PyType_Type.tp_setattro(pyclass, pymname, (PyObject*)method);
            Py_DECREF(pymname);
            Py_DECREF(method);
        }

        Py_XDECREF(attr);         // could have been found in base class or non-existent
    }
    Py_DECREF(dct);

 // collect data members (including enums)
    const Cppyy::TCppIndex_t nDataMembers = Cppyy::GetNumDatamembers(scope);
    for (Cppyy::TCppIndex_t idata = 0; idata < nDataMembers; ++idata) {
    // allow only public members
        if (!Cppyy::IsPublicData(scope, idata))
            continue;

    // enum datamembers (this in conjunction with previously collected enums above)
        if (Cppyy::IsEnumData(scope, idata) && Cppyy::IsStaticData(scope, idata)) {
        // some implementation-specific data members have no address: ignore them
            if (!Cppyy::GetDatamemberOffset(scope, idata))
                continue;

        // two options: this is a static variable, or it is the enum value, the latter
        // already exists, so check for it and move on if set
            PyObject* eset = PyObject_GetAttrString(pyclass,
                const_cast<char*>(Cppyy::GetDatamemberName(scope, idata).c_str()));
            if (eset) {
                Py_DECREF(eset);
                continue;
            }

            PyErr_Clear();

        // it could still be that this is an anonymous enum, which is not in the list
        // provided by the class
            if (strstr(Cppyy::GetDatamemberType(scope, idata).c_str(), "(anonymous)") != 0) {
                AddPropertyToClass(pyclass, scope, idata);
                continue;
            }
        }

    // properties (aka public (static) data members)
        AddPropertyToClass(pyclass, scope, idata);
    }

// restore custom __getattr__
    Py_TYPE(pyclass)->tp_getattro = oldgetattro;

// all ok, done
    return 0;
}

//----------------------------------------------------------------------------
static void CollectUniqueBases(Cppyy::TCppType_t klass, std::deque<std::string>& uqb)
{
// collect bases in acceptable mro order, while removing duplicates (this may
// break the overload resolution in esoteric cases, but otherwise the class can
// not be used at all, as CPython will refuse the mro).
    size_t nbases = Cppyy::GetNumBases(klass);

    std::deque<Cppyy::TCppType_t> bids;
    for (size_t ibase = 0; ibase < nbases; ++ibase) {
        const std::string& name = Cppyy::GetBaseName(klass, ibase);
        int decision = 2;
        Cppyy::TCppType_t tp = Cppyy::GetScope(name);
        if (!tp) continue;   // means this base with not be available Python-side
        for (size_t ibase2 = 0; ibase2 < uqb.size(); ++ibase2) {
            if (uqb[ibase2] == name) {         // not unique ... skip
                decision = 0;
                break;
            }

            if (Cppyy::IsSubtype(tp, bids[ibase2])) {
            // mro requirement: sub-type has to follow base
                decision = 1;
                break;
            }
        }

        if (decision == 1) {
            uqb.push_front(name);
            bids.push_front(tp);
        } else if (decision == 2) {
            uqb.push_back(name);
            bids.push_back(tp);
        }
    // skipped if decision == 0 (not unique)
    }
}

static PyObject* BuildCppClassBases(Cppyy::TCppType_t klass)
{
// Build a tuple of python proxy classes of all the bases of the given 'klass'.
    std::deque<std::string> uqb;
    CollectUniqueBases(klass, uqb);

// allocate a tuple for the base classes, special case for first base
    size_t nbases = uqb.size();

    PyObject* pybases = PyTuple_New(nbases ? nbases : 1);
    if (!pybases)
        return nullptr;

// build all the bases
    if (nbases == 0) {
        Py_INCREF((PyObject*)(void*)&CPPInstance_Type);
        PyTuple_SET_ITEM(pybases, 0, (PyObject*)(void*)&CPPInstance_Type);
    } else {
        for (std::deque<std::string>::size_type ibase = 0; ibase < nbases; ++ibase) {
            PyObject* pyclass = CreateScopeProxy(uqb[ibase]);
            if (!pyclass) {
                Py_DECREF(pybases);
                return nullptr;
            }

            PyTuple_SET_ITEM(pybases, ibase, pyclass);
        }

    // special case, if true python types enter the hierarchy, make sure that
    // the first base seen is still the CPPInstance_Type
        if (!PyObject_IsSubclass(PyTuple_GET_ITEM(pybases, 0), (PyObject*)&CPPInstance_Type)) {
            PyObject* newpybases = PyTuple_New(nbases+1);
            Py_INCREF((PyObject*)(void*)&CPPInstance_Type);
            PyTuple_SET_ITEM(newpybases, 0, (PyObject*)(void*)&CPPInstance_Type);
            for (int ibase = 0; ibase < (int)nbases; ++ibase) {
                PyObject* pyclass = PyTuple_GET_ITEM(pybases, ibase);
                Py_INCREF(pyclass);
                PyTuple_SET_ITEM(newpybases, ibase+1, pyclass);
            }
            Py_DECREF(pybases);
            pybases = newpybases;
        }
    }

    return pybases;
}

} // namespace CPyCppyy

//----------------------------------------------------------------------------
PyObject* CPyCppyy::GetScopeProxy(Cppyy::TCppScope_t scope)
{
// Retrieve scope proxy from the known ones.
    PyClassMap_t::iterator pci = gPyClasses.find(scope);
    if (pci != gPyClasses.end()) {
        PyObject* pyclass = PyWeakref_GetObject(pci->second);
        if (pyclass != Py_None) {
            Py_INCREF(pyclass);
            return pyclass;
        }
    }

    return nullptr;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CreateScopeProxy(Cppyy::TCppScope_t scope)
{
// Convenience function with a lookup first through the known existing proxies.
    PyObject* pyclass = GetScopeProxy(scope);
    if (pyclass)
        return pyclass;

    return CreateScopeProxy(Cppyy::GetScopedFinalName(scope));
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CreateScopeProxy(PyObject*, PyObject* args)
{
// Build a python shadow class for the named C++ class.
    std::string cname = CPyCppyy_PyText_AsString(PyTuple_GetItem(args, 0));
    if (PyErr_Occurred())
        return nullptr;

   return CreateScopeProxy(cname);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CreateScopeProxy(const std::string& name, PyObject* parent)
{
// Build a python shadow class for the named C++ class or namespace.

// determine complete scope name, if a python parent has been given
    std::string scName = "";
    if (parent) {
        if (CPPScope_Check(parent))
            scName = Cppyy::GetScopedFinalName(((CPPScope*)parent)->fCppType);
        else {
            PyObject* parname = PyObject_GetAttr(parent, PyStrings::gName);
            if (!parname) {
                PyErr_Format(PyExc_SystemError, "given scope has no name for %s", name.c_str());
                return nullptr;
            }

        // should be a string
            scName = CPyCppyy_PyText_AsString(parname);
            Py_DECREF(parname);
            if (PyErr_Occurred())
                return nullptr;
        }

    // accept this parent scope and use it's name for prefixing
        Py_INCREF(parent);
    }

// retrieve C++ class (this verifies name, and is therefore done first)
    const std::string& lookup = scName.empty() ? name : (scName+"::"+name);
    Cppyy::TCppScope_t klass = Cppyy::GetScope(lookup);

    if (!(bool)klass && Cppyy::IsTemplate(lookup)) {
    // a "naked" templated class is requested: return callable proxy for instantiations
        PyObject* pytcl = PyObject_GetAttr(gThisModule, PyStrings::gTemplate);
        PyObject* pytemplate = PyObject_CallFunction(
            pytcl, const_cast<char*>("s"), const_cast<char*>(lookup.c_str()));
        Py_DECREF(pytcl);

    // cache the result
        AddScopeToParent(parent ? parent : gThisModule, name, pytemplate);

    // done, next step should be a call into this template
        Py_XDECREF(parent);
        return pytemplate;
    }

    if (!(bool)klass) {
    // could be an enum, which are treated seperately in CPPScope (TODO: maybe they
    // should be handled here instead anyway??)
        if (Cppyy::IsEnum(lookup))
            return nullptr;

    // final possibility is a typedef of a builtin; these are mapped on the python side
        std::string resolved = Cppyy::ResolveName(lookup);
        if (gPyTypeMap) {
            PyObject* tc = PyDict_GetItemString(gPyTypeMap, resolved.c_str()); // borrowed
            if (tc && PyCallable_Check(tc)) {
                PyObject* nt = PyObject_CallFunction(tc, (char*)"ss", name.c_str(), scName.c_str());
                if (nt) {
                    if (parent) {
                        AddScopeToParent(parent, name, nt);
                        Py_DECREF(parent);
                    }
                    return nt;
                }
                PyErr_Clear();
            }
        }

    // all options have been exhausted: it doesn't exist as such
        PyErr_Format(PyExc_TypeError, "\'%s\' is not a known C++ class", lookup.c_str());
        Py_XDECREF(parent);
        return nullptr;
    }

// locate class by ID, if possible, to prevent parsing scopes/templates anew
    PyObject* pyscope = GetScopeProxy(klass);
    if (pyscope) {
        if (parent) {
            AddScopeToParent(parent, name, pyscope);
            Py_DECREF(parent);
        }
        return pyscope;
    }

// now have a class ... get the actual, fully scoped class name, so that typedef'ed
// classes are created in the right place
    const std::string& actual = Cppyy::GetScopedFinalName(klass);
    if (actual != lookup) {
        pyscope = CreateScopeProxy(actual);
        if (!pyscope) PyErr_Clear();
    }

// locate the parent, if necessary, for memoizing the class if not specified
    std::string::size_type last = 0;
    if (!parent) {
    // TODO: move this to TypeManip, which already has something similar in
    // the form of 'extract_namespace'
    // need to deal with template parameters that can have scopes themselves
        int tpl_open = 0;
        for (std::string::size_type pos = 0; pos < name.size(); ++pos) {
            std::string::value_type c = name[pos];

        // count '<' and '>' to be able to skip template contents
            if (c == '<')
                ++tpl_open;
            else if (c == '>')
                --tpl_open;

      // by only checking for "::" the last part (class name) is dropped
            else if (tpl_open == 0 && \
                c == ':' && pos+1 < name.size() && name[ pos+1 ] == ':') {
            // found a new scope part
                const std::string& part = name.substr(last, pos-last);

                PyObject* next = PyObject_GetAttrString(
                    parent ? parent : gThisModule, const_cast<char*>(part.c_str()));

                if (!next) {            // lookup failed, try to create it
                    PyErr_Clear();
                    next = CreateScopeProxy(part, parent);
                }
                Py_XDECREF(parent);

                if (!next)              // create failed, give up
                    return nullptr;

            // found scope part
                parent = next;

            // done with part (note that pos is moved one ahead here)
                last = pos+2; ++pos;
            }
        }

        if (parent && !CPPScope_Check(parent)) {
        // Special case: parent found is not one of ours (it's e.g. a pure Python module), so
        // continuing would fail badly. One final lookup, then out of here ...
            std::string unscoped = name.substr(last, std::string::npos);
            PyObject* ret = PyObject_GetAttrString(parent, unscoped.c_str());
            Py_DECREF(parent);
            return ret;
        }
    }

// use the module as a fake scope if no outer scope found
    if (!parent) {
        Py_INCREF(gThisModule);
        parent = gThisModule;
    }

// if the scope was earlier found as actual, then we're done already, otherwise
// build a new scope proxy
    if (!pyscope) {
    // construct the base classes
        PyObject* pybases = BuildCppClassBases(klass);
        if (pybases != 0) {
        // create a fresh Python class, given bases, name, and empty dictionary
            pyscope = CreateNewCppProxyClass(klass, pybases);
            Py_DECREF(pybases);
        }

    // fill the dictionary, if successful
        if (pyscope) {
            if (BuildScopeProxyDict(klass, pyscope)) {
            // something failed in building the dictionary
                Py_DECREF(pyscope);
                pyscope = nullptr;
            }
        }

    // store a ref from cppyy scope id to new python class
        if (pyscope && !(((CPPScope*)pyscope)->fFlags & CPPScope::kIsInComplete)) {
            gPyClasses[klass] = PyWeakref_NewRef(pyscope, nullptr);

            if (!(((CPPScope*)pyscope)->fFlags & CPPScope::kIsNamespace)) {
            // add python-style features to classes only
                if (!Pythonize(pyscope, Cppyy::GetScopedFinalName(klass))) {
                    Py_DECREF(pyscope);
                    pyscope = nullptr;
                }
            } else {
            // add to sys.modules to allow importing from this namespace
                PyObject* pyfullname = PyObject_GetAttr(pyscope, PyStrings::gModule);
                CPyCppyy_PyText_AppendAndDel(&pyfullname, CPyCppyy_PyText_FromString("."));
                CPyCppyy_PyText_AppendAndDel(&pyfullname, PyObject_GetAttr(pyscope, PyStrings::gName));
                PyObject* modules = PySys_GetObject(const_cast<char*>("modules"));
                if (modules && PyDict_Check(modules))
                    PyDict_SetItem(modules, pyfullname, pyscope);
                Py_DECREF(pyfullname);
            }
        }
    }

// store on parent if found/created and complete
    if (pyscope && !(((CPPScope*)pyscope)->fFlags & CPPScope::kIsInComplete))
        AddScopeToParent(parent, name, pyscope);
    Py_DECREF(parent);

// all done
    return pyscope;
}


//----------------------------------------------------------------------------
PyObject* CPyCppyy::CreateExcScopeProxy(PyObject* pyscope, PyObject* pyname, PyObject* parent)
{
// To allow use of C++ exceptions in lieue of Python exceptions, they need to
// derive from BaseException, which can not mix with the normal CPPInstance and
// use of the meta-class. Instead, encapsulate them in a forwarding class that
// derives from Pythons Exception class

// start with creation of CPPExcInstance type base classes
    std::deque<std::string> uqb;
    CollectUniqueBases(((CPPScope*)pyscope)->fCppType, uqb);
    size_t nbases = uqb.size();

// Support for multiple bases actually can not actually work as-is: the reason
// for deriving from BaseException is to guarantee the layout needed for storing
// traces. If some other base is std::exception (as e.g. boost::bad_any_cast) or
// also derives from std::exception, then there are two trace locations. OTOH,
// if the other class is a non-exception type, then the exception class does not
// need to derive from it because it can never be caught as that type forwarding
// to the proxy will work as expected, through, which is good enough).
//
// The code below restricts the hierarchy to a single base class, picking the
// "best" by filtering std::exception and non-exception bases.

    PyObject* pybases = PyTuple_New(1);
    if (nbases == 0) {
        Py_INCREF((PyObject*)(void*)&CPPExcInstance_Type);
        PyTuple_SET_ITEM(pybases, 0, (PyObject*)(void*)&CPPExcInstance_Type);
    } else {
        PyObject* best_base = nullptr;

        for (std::deque<std::string>::size_type ibase = 0; ibase < nbases; ++ibase) {
        // retrieve bases through their enclosing scope to guarantee treatment as
        // exception classes and proper caching
            const std::string& finalname = Cppyy::GetScopedFinalName(Cppyy::GetScope(uqb[ibase]));
            const std::string& parentname = TypeManip::extract_namespace(finalname);
            PyObject* base_parent = CreateScopeProxy(parentname);
            if (!base_parent) {
                Py_DECREF(pybases);
                return nullptr;
            }

            PyObject* excbase = PyObject_GetAttrString(base_parent,
                parentname.empty() ? finalname.c_str() : finalname.substr(parentname.size()+2, std::string::npos).c_str());
            Py_DECREF(base_parent);
            if (!excbase) {
                Py_DECREF(pybases);
                return nullptr;
            }

            if (PyType_IsSubtype((PyTypeObject*)excbase, &CPPExcInstance_Type)) {
                Py_XDECREF(best_base);
                best_base = excbase;
                if (finalname != "std::exception")
                    break;
            } else {
            // just skip: there will be at least one exception derived base class
                Py_DECREF(excbase);
            }
        }

        PyTuple_SET_ITEM(pybases, 0, best_base);
    }

    PyObject* args = Py_BuildValue((char*)"OO{}", pyname, pybases);

// meta-class attributes (__cpp_name__, etc.) can not be resolved lazily so add
// them directly instead in case they are needed
    PyObject* dct = PyTuple_GET_ITEM(args, 2);
    PyDict_SetItem(dct, PyStrings::gUnderlying, pyscope);
    PyDict_SetItem(dct, PyStrings::gName,    PyObject_GetAttr(pyscope, PyStrings::gName));
    PyDict_SetItem(dct, PyStrings::gCppName, PyObject_GetAttr(pyscope, PyStrings::gCppName));
    PyDict_SetItem(dct, PyStrings::gModule,  PyObject_GetAttr(pyscope, PyStrings::gModule));

// create the actual exception class
    PyObject* exc_pyscope = PyType_Type.tp_new(&PyType_Type, args, nullptr);
    Py_DECREF(args);
    Py_DECREF(pybases);

// cache the result for future lookups and return
    PyType_Type.tp_setattro(parent, pyname, exc_pyscope);
    return exc_pyscope;
}


//----------------------------------------------------------------------------
PyObject* CPyCppyy::BindCppObjectNoCast(Cppyy::TCppObject_t address,
        Cppyy::TCppType_t klass, const unsigned flags)
{
// only known or knowable objects will be bound (null object is ok)
    if (!klass) {
        PyErr_SetString(PyExc_TypeError, "attempt to bind C++ object w/o class");
        return nullptr;
    }

// retrieve python class
    PyObject* pyclass = CreateScopeProxy(klass);
    if (!pyclass)
        return nullptr;                 // error has been set in CreateScopeProxy

    bool isRef   = flags & CPPInstance::kIsReference;
    bool isValue = flags & CPPInstance::kIsValue;

// TODO: make sure that a consistent address is used (may have to be done in BindCppObject)
    if (address && !isValue /* always fresh */ && !(flags & (CPPInstance::kNoWrapConv|CPPInstance::kNoMemReg))) {
        PyObject* oldPyObject = MemoryRegulator::RetrievePyObject(
            isRef ? *(void**)address : address, pyclass);

    // ptr-ptr requires old object to be a reference to enable re-use
        if (oldPyObject && (!(flags & CPPInstance::kIsPtrPtr) ||
                ((CPPInstance*)oldPyObject)->fFlags & CPPInstance::kIsReference)) {
            return oldPyObject;
        }
    }

// if smart, instantiate a Python-side object of the underlying type, carrying the smartptr
    PyObject* smart_type = (flags != CPPInstance::kNoWrapConv && (((CPPClass*)pyclass)->fFlags & CPPScope::kIsSmart)) ? pyclass : nullptr;
    if (smart_type) {
        pyclass = CreateScopeProxy(((CPPSmartClass*)smart_type)->fUnderlyingType);
        if (!pyclass) {
        // simply restore and expose as the actual smart pointer class
            pyclass = smart_type;
            smart_type = nullptr;
        }
    }

// instantiate an object of this class
    PyObject* args = PyTuple_New(0);
    CPPInstance* pyobj =
        (CPPInstance*)((PyTypeObject*)pyclass)->tp_new((PyTypeObject*)pyclass, args, nullptr);
    Py_DECREF(args);

// bind, register and return if successful
    if (pyobj != 0) { // fill proxy value?
        unsigned objflags =
            (isRef ? CPPInstance::kIsReference : 0) | (isValue ? CPPInstance::kIsValue : 0) | (flags & CPPInstance::kIsOwner);
        pyobj->Set(address, (CPPInstance::EFlags)objflags);

        if (smart_type)
            pyobj->SetSmart(smart_type);

    // do not register null pointers, references (?), or direct usage of smart pointers or iterators
        if (address && !isRef && !(flags & (CPPInstance::kNoWrapConv|CPPInstance::kNoMemReg)))
            MemoryRegulator::RegisterPyObject(pyobj, pyobj->GetObject());
    }

// successful completion; wrap exception options to make them raiseable, normal return otherwise
    if (((CPPClass*)pyclass)->fFlags & CPPScope::kIsException) {
        PyObject* exc_obj = CPPExcInstance_Type.tp_new(&CPPExcInstance_Type, nullptr, nullptr);
        ((CPPExcInstance*)exc_obj)->fCppInstance = (PyObject*)pyobj;
        Py_DECREF(pyclass);
        return exc_obj;
    }

    Py_DECREF(pyclass);

    return (PyObject*)pyobj;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::BindCppObject(Cppyy::TCppObject_t address,
        Cppyy::TCppType_t klass, const unsigned flags)
{
// if the object is a null pointer, return a typed one (as needed for overloading)
    if (!address)
        return BindCppObjectNoCast(address, klass, flags);

// only known or knowable objects will be bound
    if (!klass) {
        PyErr_SetString(PyExc_TypeError, "attempt to bind C++ object w/o class");
        return nullptr;
    }

    bool isRef = flags & CPPInstance::kIsReference;

// get actual class for recycling checking and/or downcasting
    Cppyy::TCppType_t clActual = isRef ? 0 : Cppyy::GetActualClass(klass, address);

// downcast to real class for object returns, unless pinned
    if (clActual && klass != clActual) {
        auto pci = gPinnedTypes.find(klass);
        if (pci == gPinnedTypes.end()) {
            intptr_t offset = Cppyy::GetBaseOffset(
                clActual, klass, address, -1 /* down-cast */, true /* report errors */);
            if (offset != -1) {   // may fail if clActual not fully defined
                address = (void*)((intptr_t)address + offset);
                klass = clActual;
            }
        }
    }

// actual binding (returned object may be zero w/ a python exception set)
    return BindCppObjectNoCast(address, klass, flags);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::BindCppObjectArray(
    Cppyy::TCppObject_t address, Cppyy::TCppType_t klass, dims_t dims)
{
// TODO: this function exists for symmetry; need to figure out if it's useful
    return TupleOfInstances_New(address, klass, dims[0], dims+1);
}
