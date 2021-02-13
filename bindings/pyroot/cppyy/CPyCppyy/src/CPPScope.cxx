// Bindings
#include "CPyCppyy.h"
#include "CPPScope.h"
#include "CPPDataMember.h"
#include "CPPFunction.h"
#include "CPPOverload.h"
#include "CustomPyTypes.h"
#include "Dispatcher.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"
#include "TemplateProxy.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <string.h>
#include <algorithm>         // for for_each
#include <set>
#include <string>
#include <vector>


namespace CPyCppyy {

extern PyTypeObject CPPInstance_Type;

//- helpers ------------------------------------------------------------------
static inline PyObject* add_template(PyObject* pyclass,
    const std::string& name, std::vector<PyCallable*>* overloads = nullptr)
{
// If templated, the user-facing function must be the template proxy, but the
// specific lookup must be the current overload, if already found.
    TemplateProxy* pytmpl = nullptr;

    const std::string& ncl = TypeManip::clean_type(name);
    if (ncl != name) {
        PyObject* pyncl = CPyCppyy_PyText_InternFromString(ncl.c_str());
        pytmpl = (TemplateProxy*)PyType_Type.tp_getattro((PyObject*)Py_TYPE(pyclass), pyncl);
        if (!pytmpl) {
            PyErr_Clear();
            pytmpl = TemplateProxy_New(ncl, ncl, pyclass);
        // cache the template on its clean name
            PyType_Type.tp_setattro((PyObject*)Py_TYPE(pyclass), pyncl, (PyObject*)pytmpl);
        }
        Py_DECREF(pyncl);
    }

    if (pytmpl) {
        if (!TemplateProxy_CheckExact((PyObject*)pytmpl)) {
            Py_DECREF(pytmpl);
            return nullptr;
        }
    } else
        pytmpl = TemplateProxy_New(ncl, ncl, pyclass);

    if (overloads) {
    // adopt the new overloads
        if (ncl != name)
            for (auto clb : *overloads) pytmpl->AdoptTemplate(clb);
        else
            for (auto clb : *overloads) pytmpl->AdoptMethod(clb);
    }

    if (ncl == name)
        return (PyObject*)pytmpl;

    Py_DECREF(pytmpl);
    return nullptr;       // so that caller caches the method on full name
}

//----------------------------------------------------------------------------
static int enum_setattro(PyObject* /* pyclass */, PyObject* /* pyname */, PyObject* /* pyval */)
{
// Helper to make enums read-only.
    PyErr_SetString(PyExc_TypeError, "enum values are read-only");
    return -1;
}


//= CPyCppyy type proxy construction/destruction =============================
static PyObject* meta_alloc(PyTypeObject* meta, Py_ssize_t nitems)
{
// pure memory allocation; object initialization is in pt_new
    return PyType_Type.tp_alloc(meta, nitems);
}

//----------------------------------------------------------------------------
static void meta_dealloc(CPPScope* scope)
{
    if (scope->fFlags & CPPScope::kIsNamespace) {
        if (scope->fImp.fUsing) {
            for (auto pyobj : *scope->fImp.fUsing) Py_DECREF(pyobj);
            delete scope->fImp.fUsing; scope->fImp.fUsing = nullptr;
        }
    } else if (!(scope->fFlags & CPPScope::kIsPython)) {
        delete scope->fImp.fCppObjects; scope->fImp.fCppObjects = nullptr;
    }
    delete scope->fOperators;
    free(scope->fModuleName);
    return PyType_Type.tp_dealloc((PyObject*)scope);
}

//-----------------------------------------------------------------------------
static PyObject* meta_getcppname(CPPScope* scope, void*)
{
    if ((void*)scope == (void*)&CPPInstance_Type)
        return CPyCppyy_PyText_FromString("CPPInstance_Type");
    return CPyCppyy_PyText_FromString(Cppyy::GetScopedFinalName(scope->fCppType).c_str());
}

//-----------------------------------------------------------------------------
static PyObject* meta_getmodule(CPPScope* scope, void*)
{
    if ((void*)scope == (void*)&CPPInstance_Type)
        return CPyCppyy_PyText_FromString("cppyy.gbl");

    if (scope->fModuleName)
        return CPyCppyy_PyText_FromString(scope->fModuleName);

// get C++ representation of outer scope
    std::string modname =
        TypeManip::extract_namespace(Cppyy::GetScopedFinalName(scope->fCppType));
    if (modname.empty())
        return CPyCppyy_PyText_FromString(const_cast<char*>("cppyy.gbl"));

// now peel scopes one by one, pulling in the python naming (which will
// simply recurse if not overridden in python)
    PyObject* pymodule = nullptr;
    PyObject* pyscope = CPyCppyy::GetScopeProxy(Cppyy::GetScope(modname));
    if (pyscope) {
    // get the module of our module
        pymodule = PyObject_GetAttr(pyscope, PyStrings::gModule);
        if (pymodule) {
        // append name of our module
            PyObject* pymodname = PyObject_GetAttr(pyscope, PyStrings::gName);
            if (pymodname) {
                CPyCppyy_PyText_AppendAndDel(&pymodule, CPyCppyy_PyText_FromString("."));
                CPyCppyy_PyText_AppendAndDel(&pymodule, pymodname);
            }
        }
        Py_DECREF(pyscope);
    }

    if (pymodule)
        return pymodule;
    PyErr_Clear();

// lookup through python failed, so simply cook up a '::' -> '.' replacement
    TypeManip::cppscope_to_pyscope(modname);
    return CPyCppyy_PyText_FromString(("cppyy.gbl."+modname).c_str());
}

//-----------------------------------------------------------------------------
static int meta_setmodule(CPPScope* scope, PyObject* value, void*)
{
    if ((void*)scope == (void*)&CPPInstance_Type) {
        PyErr_SetString(PyExc_AttributeError,
            "attribute \'__module__\' of 'cppyy.CPPScope\' objects is not writable");
        return -1;
    }

    const char* newname = CPyCppyy_PyText_AsStringChecked(value);
    if (!value)
        return -1;

    free(scope->fModuleName);
    Py_ssize_t sz = CPyCppyy_PyText_GET_SIZE(value);
    scope->fModuleName = (char*)malloc(sz+1);
    memcpy(scope->fModuleName, newname, sz+1);

    return 0;
}

//----------------------------------------------------------------------------
static PyObject* meta_repr(CPPScope* scope)
{
// Specialized b/c type_repr expects __module__ to live in the dictionary,
// whereas it is a property (to save memory).
    if ((void*)scope == (void*)&CPPInstance_Type)
        return CPyCppyy_PyText_FromFormat(
            const_cast<char*>("<class cppyy.CPPInstance at %p>"), scope);

    if (scope->fFlags & (CPPScope::kIsMeta | CPPScope::kIsPython)) {
    // either meta type or Python-side derived class: use default type printing
        return PyType_Type.tp_repr((PyObject*)scope);
    }

// printing of C++ classes
    PyObject* modname = meta_getmodule(scope, nullptr);
    std::string clName = Cppyy::GetFinalName(scope->fCppType);
    const char* kind = (scope->fFlags & CPPScope::kIsNamespace) ? "namespace" : "class";

    PyObject* repr = CPyCppyy_PyText_FromFormat("<%s %s.%s at %p>",
        kind, CPyCppyy_PyText_AsString(modname), clName.c_str(), scope);

    Py_DECREF(modname);
    return repr;
}


//= CPyCppyy type metaclass behavior =========================================
static PyObject* pt_new(PyTypeObject* subtype, PyObject* args, PyObject* kwds)
{
// Called when CPPScope acts as a metaclass; since type_new always resets
// tp_alloc, and since it does not call tp_init on types, the metaclass is
// being fixed up here, and the class is initialized here as well.

// fixup of metaclass (left permanent, and in principle only called once b/c
// cppyy caches python classes)
    subtype->tp_alloc   = (allocfunc)meta_alloc;
    subtype->tp_dealloc = (destructor)meta_dealloc;

// creation of the python-side class; extend the size if this is a smart ptr
    Cppyy::TCppType_t raw{0}; Cppyy::TCppMethod_t deref{0};
    if (CPPScope_CheckExact(subtype)) {
        if (Cppyy::GetSmartPtrInfo(Cppyy::GetScopedFinalName(((CPPScope*)subtype)->fCppType), &raw, &deref))
            subtype->tp_basicsize = sizeof(CPPSmartClass);
    }
    CPPScope* result = (CPPScope*)PyType_Type.tp_new(subtype, args, kwds);
    if (!result)
        return nullptr;

    result->fFlags      = CPPScope::kNone;
    result->fOperators  = nullptr;
    result->fModuleName = nullptr;

    if (raw && deref) {
        result->fFlags |= CPPScope::kIsSmart;
        ((CPPSmartClass*)result)->fUnderlyingType = raw;
        ((CPPSmartClass*)result)->fDereferencer   = deref;
    }

// initialization of class (based on metatype)
    const char* mp = strstr(subtype->tp_name, "_meta");
    if (!mp || !CPPScope_CheckExact(subtype)) {
    // there has been a user meta class override in a derived class, so do
    // the consistent thing, thus allowing user control over naming
        result->fCppType = Cppyy::GetScope(
            CPyCppyy_PyText_AsString(PyTuple_GET_ITEM(args, 0)));
    } else {
    // coming here from cppyy or from sub-classing in python; take the
    // C++ type from the meta class to make sure that the latter category
    // has fCppType properly set (it inherits the meta class, but has an
    // otherwise unknown (or wrong) C++ type)
        result->fCppType = ((CPPScope*)subtype)->fCppType;

    // the following is not robust, but by design, C++ classes get their
    // dictionaries filled after creation (chicken & egg problem as they
    // can return themselves in methods), whereas a derived Python class
    // with method overrides will have a non-empty dictionary (even if it
    // has no methods, it will at least have a module name)
        if (3 <= PyTuple_GET_SIZE(args)) {
            PyObject* dct = PyTuple_GET_ITEM(args, 2);
            Py_ssize_t sz = PyDict_Size(dct);
            if (0 < sz && !Cppyy::IsNamespace(result->fCppType)) {
                result->fFlags |= CPPScope::kIsPython;
                std::ostringstream errmsg;
                if (!InsertDispatcher(result, PyTuple_GET_ITEM(args, 1), dct, errmsg)) {
                    PyErr_Format(PyExc_TypeError, "no python-side overrides supported (%s)", errmsg.str().c_str());
                    return nullptr;
                } else {
                // the direct base can be useful for some templates, such as shared_ptrs,
                // so make it accessible (the __cpp_cross__ data member also signals that
                // this is a cross-inheritance class)
                    PyObject* bname = CPyCppyy_PyText_FromString(Cppyy::GetBaseName(result->fCppType, 0).c_str());
                    if (PyObject_SetAttrString((PyObject*)result, "__cpp_cross__", bname) == -1)
                        PyErr_Clear();
                    Py_DECREF(bname);
                }
            } else if (sz == (Py_ssize_t)-1)
                PyErr_Clear();
        }
    }

// maps for using namespaces and tracking objects
    if (!Cppyy::IsNamespace(result->fCppType)) {
        static Cppyy::TCppType_t exc_type = (Cppyy::TCppType_t)Cppyy::GetScope("std::exception");
        if (Cppyy::IsSubtype(result->fCppType, exc_type))
            result->fFlags |= CPPScope::kIsException;
        if (!(result->fFlags & CPPScope::kIsPython))
            result->fImp.fCppObjects = new CppToPyMap_t;
        else {
        // special case: the C++ objects should be stored with the associated C++, not Python, type
            CPPClass* kls = (CPPClass*)GetScopeProxy(result->fCppType);
            if (kls) {
                result->fImp.fCppObjects = kls->fImp.fCppObjects;
                Py_DECREF(kls);
            } else
                result->fImp.fCppObjects = nullptr;
        }
    } else {
        result->fImp.fUsing = nullptr;
        result->fFlags |= CPPScope::kIsNamespace;
    }

    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        return nullptr;
    }
    return (PyObject*)result;
}


//----------------------------------------------------------------------------
static PyObject* meta_getattro(PyObject* pyclass, PyObject* pyname)
{
// normal type-based lookup
    PyObject* attr = PyType_Type.tp_getattro(pyclass, pyname);
    if (attr || pyclass == (PyObject*)&CPPInstance_Type)
        return attr;

    if (!CPyCppyy_PyText_CheckExact(pyname) || !CPPScope_Check(pyclass))
        return nullptr;

// filter for python specials
    std::string name = CPyCppyy_PyText_AsString(pyname);
    if (name.size() >= 2 && name.compare(0, 2, "__") == 0 &&
            name.compare(name.size()-2, name.size(), "__") == 0)
        return nullptr;

// more elaborate search in case of failure (eg. for inner classes on demand)
    std::vector<Utility::PyError_t> errors;
    Utility::FetchError(errors);
    attr = CreateScopeProxy(name, pyclass);
    if (CPPScope_Check(attr) && (((CPPScope*)attr)->fFlags & CPPScope::kIsException)) {
    // Instead of the CPPScope, return a fresh exception class derived from CPPExcInstance.
        return CreateExcScopeProxy(attr, pyname, pyclass);
    }

    CPPScope* klass = ((CPPScope*)pyclass);
    if (!attr) {
        Utility::FetchError(errors);
        Cppyy::TCppScope_t scope = klass->fCppType;

    // namespaces may have seen updates in their list of global functions, which
    // are available as "methods" even though they're not really that
        if (klass->fFlags & CPPScope::kIsNamespace) {
        // tickle lazy lookup of functions
            const std::vector<Cppyy::TCppIndex_t> methods =
                Cppyy::GetMethodIndicesFromName(scope, name);
            if (!methods.empty()) {
            // function exists, now collect overloads
                std::vector<PyCallable*> overloads;
                for (auto idx : methods) {
                    overloads.push_back(
                        new CPPFunction(scope, Cppyy::GetMethod(scope, idx)));
                }

            // Note: can't re-use Utility::AddClass here, as there's the risk of
            // a recursive call. Simply add method directly, as we're guaranteed
            // that it doesn't exist yet.
                if (Cppyy::ExistsMethodTemplate(scope, name))
                    attr = add_template(pyclass, name, &overloads);

                if (!attr)    // add_template can fail if the method can not be added
                    attr = (PyObject*)CPPOverload_New(name, overloads);
            }

        // tickle lazy lookup of data members
            if (!attr) {
                Cppyy::TCppIndex_t dmi = Cppyy::GetDatamemberIndex(scope, name);
                if (dmi != (Cppyy::TCppIndex_t)-1) attr = (PyObject*)CPPDataMember_New(scope, dmi);
            }
        }

    // this may be a typedef that resolves to a sugared type
        if (!attr) {
            const std::string& lookup = Cppyy::GetScopedFinalName(klass->fCppType) + "::" + name;
            const std::string& resolved = Cppyy::ResolveName(lookup);
            if (resolved != lookup) {
                const std::string& cpd = Utility::Compound(resolved);
                if (cpd == "*") {
                    const std::string& clean = TypeManip::clean_type(resolved, false, true);
                    Cppyy::TCppType_t tcl = Cppyy::GetScope(clean);
                    if (tcl) {
                        typedefpointertoclassobject* tpc =
                            PyObject_GC_New(typedefpointertoclassobject, &TypedefPointerToClass_Type);
                        tpc->fType = tcl;
                        attr = (PyObject*)tpc;
                    }
                }
            }
        }

    // function templates that have not been instantiated
        if (!attr) {
            if (Cppyy::ExistsMethodTemplate(scope, name))
                attr = add_template(pyclass, name);
            else {
            // for completeness in error reporting
                PyErr_Format(PyExc_TypeError, "\'%s\' is not a known C++ template", name.c_str());
                Utility::FetchError(errors);
            }
        }

    // enums types requested as type (rather than the constants)
        if (!attr) {
        // TODO: IsEnum should deal with the scope, using klass->GetListOfEnums()->FindObject()
            if (Cppyy::IsEnum(scope == Cppyy::gGlobalScope ? name : Cppyy::GetScopedFinalName(scope)+"::"+name)) {
            // enum types (incl. named and class enums)
                Cppyy::TCppEnum_t etype = Cppyy::GetEnum(scope, name);
                if (etype) {
                // create new enum type with labeled values in place, with a meta-class
                // to make sure the enum values are read-only
                    PyObject* pymetabases = PyTuple_New(1);
                    PyObject* btype = (PyObject*)Py_TYPE(&PyInt_Type);
                    Py_INCREF(btype);
                    PyTuple_SET_ITEM(pymetabases, 0, btype);

                    PyObject* args = Py_BuildValue((char*)"sO{}", (name+"_meta").c_str(), pymetabases);
                    Py_DECREF(pymetabases);
                    PyObject* pymeta = PyType_Type.tp_new(Py_TYPE(&PyInt_Type), args, nullptr);
                    ((PyTypeObject*)pymeta)->tp_setattro = enum_setattro;
                    Py_DECREF(args);

                // prepare the base class
                    PyObject* pybases = PyTuple_New(1);
                    Py_INCREF(&PyInt_Type);
                    PyTuple_SET_ITEM(pybases, 0, (PyObject*)&PyInt_Type);

                // collect the enum values
                    Cppyy::TCppIndex_t ndata = Cppyy::GetNumEnumData(etype);
                    PyObject* dct = PyDict_New();
                    for (Cppyy::TCppIndex_t idata = 0; idata < ndata; ++idata) {
                        PyObject* val = PyLong_FromLongLong(Cppyy::GetEnumDataValue(etype, idata));
                        PyDict_SetItemString(dct, Cppyy::GetEnumDataName(etype, idata).c_str(), val);
                        Py_DECREF(val);
                    }

                // add the __cpp_name__ for templates
                    PyObject* cppname = nullptr;
                    if (scope == Cppyy::gGlobalScope) {
                        Py_INCREF(pyname);
                        cppname = pyname;
                    } else
                        cppname = CPyCppyy_PyText_FromString((Cppyy::GetScopedFinalName(scope)+"::"+name).c_str());
                    PyDict_SetItem(dct, PyStrings::gCppName, cppname);
                    Py_DECREF(cppname);

                // create the actual enum class
                    args = Py_BuildValue((char*)"sOO", name.c_str(), pybases, dct);
                    Py_DECREF(pybases);
                    Py_DECREF(dct);
                    attr = ((PyTypeObject*)pymeta)->tp_new((PyTypeObject*)pymeta, args, nullptr);

                // final cleanup
                    Py_DECREF(args);
                    Py_DECREF(pymeta);

                } else {
                // presumably not a class enum; simply pretend int
                    Py_INCREF(&PyInt_Type);
                    attr = (PyObject*)&PyInt_Type;
                }
            } else {
            // for completeness in error reporting
                PyErr_Format(PyExc_TypeError, "\'%s\' is not a known C++ enum", name.c_str());
                Utility::FetchError(errors);
            }
        }

        if (attr) {
        // cache the result
            if (CPPDataMember_Check(attr)) {
                PyType_Type.tp_setattro((PyObject*)Py_TYPE(pyclass), pyname, attr);
                Py_DECREF(attr);
                attr = PyType_Type.tp_getattro(pyclass, pyname);
            } else
                PyType_Type.tp_setattro(pyclass, pyname, attr);

        } else {
            Utility::FetchError(errors);
        }
    }

    if (!attr && (klass->fFlags & CPPScope::kIsNamespace)) {
    // refresh using list as necessary
        const std::vector<Cppyy::TCppScope_t>& uv = Cppyy::GetUsingNamespaces(klass->fCppType);
        if (!klass->fImp.fUsing || uv.size() != klass->fImp.fUsing->size()) {
            if (klass->fImp.fUsing) {
                for (auto pyref : *klass->fImp.fUsing) Py_DECREF(pyref);
                klass->fImp.fUsing->clear();
            } else
                klass->fImp.fUsing = new std::vector<PyObject*>;

        // reload and reset weak refs
            for (auto uid : uv) {
                std::string uname = Cppyy::GetScopedFinalName(uid);
                PyObject* pyuscope = CreateScopeProxy(uname);
                if (pyuscope) {
                    klass->fImp.fUsing->push_back(PyWeakref_NewRef(pyuscope, nullptr));
                // the namespace may not otherwise be held, so tie the lifetimes
                    PyObject* llname = CPyCppyy_PyText_FromString(("__lifeline_"+uname).c_str());
                    PyType_Type.tp_setattro(pyclass, llname, pyuscope);
                    Py_DECREF(llname);
                    Py_DECREF(pyuscope);
                }
            }
        }

    // try all outstanding using namespaces in turn to find the attribute (will cache
    // locally later; TODO: doing so may cause pathological cases)
        for (auto pyref : *klass->fImp.fUsing) {
            PyObject* pyuscope = PyWeakref_GetObject(pyref);
            if (pyuscope) {
                attr = PyObject_GetAttr(pyuscope, pyname);
                if (attr) break;
                PyErr_Clear();
            }
        }
    }

    if (attr) {
        std::for_each(errors.begin(), errors.end(), Utility::PyError_t::Clear);
        PyErr_Clear();
    } else {
    // not found: prepare a full error report
        PyObject* topmsg = nullptr;
        PyObject* sklass = PyObject_Str(pyclass);
        if (sklass) {
            topmsg = CPyCppyy_PyText_FromFormat("%s has no attribute \'%s\'. Full details:",
                CPyCppyy_PyText_AsString(sklass), CPyCppyy_PyText_AsString(pyname));
            Py_DECREF(sklass);
        } else {
            topmsg = CPyCppyy_PyText_FromFormat("no such attribute \'%s\'. Full details:",
                CPyCppyy_PyText_AsString(pyname));
        }
        SetDetailedException(errors, topmsg /* steals */, PyExc_AttributeError /* default error */);
    }

    return attr;
}

//----------------------------------------------------------------------------
static int meta_setattro(PyObject* pyclass, PyObject* pyname, PyObject* pyval)
{
// Global data and static data in namespaces is found lazily, thus if the first
// use is setting of the global data by the user, it will not be reflected on
// the C++ side, b/c there is no descriptor yet. This triggers the creation for
// for such data as necessary. The many checks to narrow down the specific case
// are needed to prevent unnecessary lookups and recursion.
    if (((CPPScope*)pyclass)->fFlags & CPPScope::kIsNamespace) {
    // skip if the given pyval is a descriptor already, or an unassignable class
        if (!CPyCppyy::CPPDataMember_Check(pyval) && !CPyCppyy::CPPScope_Check(pyval)) {
            std::string name = CPyCppyy_PyText_AsString(pyname);
            Cppyy::TCppIndex_t dmi = Cppyy::GetDatamemberIndex(((CPPScope*)pyclass)->fCppType, name);
            if (dmi != (Cppyy::TCppIndex_t)-1)
                meta_getattro(pyclass, pyname);       // triggers creation
        }
    }

    return PyType_Type.tp_setattro(pyclass, pyname, pyval);
}


//----------------------------------------------------------------------------
// p2.7 does not have a __dir__ in object, and object.__dir__ in p3 does not
// quite what I'd expected of it, so the following pulls in the internal code
#include "PyObjectDir27.inc"

static PyObject* meta_dir(CPPScope* klass)
{
// Collect a list of everything (currently) available in the namespace.
// The backend can filter by returning empty strings. Special care is
// taken for functions, which need not be unique (overloading).
    using namespace Cppyy;

    if ((void*)klass == (void*)&CPPInstance_Type)
        return PyList_New(0);

    if (!CPyCppyy::CPPScope_Check((PyObject*)klass)) {
        PyErr_SetString(PyExc_TypeError, "C++ proxy scope expected");
        return nullptr;
    }

    PyObject* dirlist = _generic_dir((PyObject*)klass);
    if (!(klass->fFlags & CPPScope::kIsNamespace))
        return dirlist;

    std::set<std::string> cppnames;
    Cppyy::GetAllCppNames(klass->fCppType, cppnames);

// get rid of duplicates
    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(dirlist); ++i)
        cppnames.insert(CPyCppyy_PyText_AsString(PyList_GET_ITEM(dirlist, i)));

    Py_DECREF(dirlist);
    dirlist = PyList_New(cppnames.size());

// copy total onto python list
    Py_ssize_t i = 0;
    for (const auto& name : cppnames) {
        PyList_SET_ITEM(dirlist, i++, CPyCppyy_PyText_FromString(name.c_str()));
    }
    return dirlist;
}

//-----------------------------------------------------------------------------
static PyMethodDef meta_methods[] = {
    {(char*)"__dir__", (PyCFunction)meta_dir, METH_NOARGS, nullptr},
    {(char*)nullptr, nullptr, 0, nullptr}
};


//-----------------------------------------------------------------------------
static PyGetSetDef meta_getset[] = {
    {(char*)"__cpp_name__", (getter)meta_getcppname, nullptr, nullptr, nullptr},
    {(char*)"__module__",   (getter)meta_getmodule,  (setter)meta_setmodule, nullptr, nullptr},
    {(char*)nullptr, nullptr, nullptr, nullptr, nullptr}
};


//= CPyCppyy object proxy type type ==========================================
PyTypeObject CPPScope_Type = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    (char*)"cppyy.CPPScope",       // tp_name
    sizeof(CPyCppyy::CPPScope),    // tp_basicsize
    0,                             // tp_itemsize
    0,                             // tp_dealloc
    0,                             // tp_print
    0,                             // tp_getattr
    0,                             // tp_setattr
    0,                             // tp_compare
    (reprfunc)meta_repr,           // tp_repr
    0,                             // tp_as_number
    0,                             // tp_as_sequence
    0,                             // tp_as_mapping
    0,                             // tp_hash
    0,                             // tp_call
    0,                             // tp_str
    (getattrofunc)meta_getattro,   // tp_getattro
    (setattrofunc)meta_setattro,   // tp_setattro
    0,                             // tp_as_buffer
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE
#if PY_VERSION_HEX >= 0x03040000
        | Py_TPFLAGS_TYPE_SUBCLASS
#endif
        ,                          // tp_flags
    (char*)"CPyCppyy metatype (internal)",        // tp_doc
    0,                             // tp_traverse
    0,                             // tp_clear
    0,                             // tp_richcompare
    0,                             // tp_weaklistoffset
    0,                             // tp_iter
    0,                             // tp_iternext
    meta_methods,                  // tp_methods
    0,                             // tp_members
    meta_getset,                   // tp_getset
    &PyType_Type,                  // tp_base
    0,                             // tp_dict
    0,                             // tp_descr_get
    0,                             // tp_descr_set
    0,                             // tp_dictoffset
    0,                             // tp_init
    0,                             // tp_alloc
    (newfunc)pt_new,               // tp_new
    0,                             // tp_free
    0,                             // tp_is_gc
    0,                             // tp_bases
    0,                             // tp_mro
    0,                             // tp_cache
    0,                             // tp_subclasses
    0                              // tp_weaklist
#if PY_VERSION_HEX >= 0x02030000
    , 0                            // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
    , 0                            // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x03040000
    , 0                            // tp_finalize
#endif
};

} // namespace CPyCppyy
