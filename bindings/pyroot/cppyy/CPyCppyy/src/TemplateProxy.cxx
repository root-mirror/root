// Bindings
#include "CPyCppyy.h"
#include "TemplateProxy.h"
#include "CPPClassMethod.h"
#include "CPPConstructor.h"
#include "CPPFunction.h"
#include "CPPMethod.h"
#include "CPPOverload.h"
#include "PyCallable.h"
#include "PyStrings.h"
#include "Utility.h"

// Standard
#include <algorithm>


namespace CPyCppyy {

//- helper for ctypes conversions --------------------------------------------
static PyObject* TC2CppName(PyObject* pytc, const char* cpd, bool allow_voidp)
{
    const char* name = nullptr;
    if (CPyCppyy_PyText_Check(pytc)) {
        char tc = ((char*)CPyCppyy_PyText_AsString(pytc))[0];
        switch (tc) {
            case '?': name = "bool";               break;
            case 'c': name = "char";               break;
            case 'b': name = "char";               break;
            case 'B': name = "unsigned char";      break;
            case 'h': name = "short";              break;
            case 'H': name = "unsigned short";     break;
            case 'i': name = "int";                break;
            case 'I': name = "unsigned int";       break;
            case 'l': name = "long";               break;
            case 'L': name = "unsigned long";      break;
            case 'q': name = "long long";          break;
            case 'Q': name = "unsigned long long"; break;
            case 'f': name = "float";              break;
            case 'd': name = "double";             break;
            case 'g': name = "long double";        break;
            default:  name = (allow_voidp ? "void*" : nullptr); break;
        }
    }

    if (name)
        return CPyCppyy_PyText_FromString((std::string{name}+cpd).c_str());
    return nullptr;
}

//----------------------------------------------------------------------------
TemplateInfo::TemplateInfo() : fCppName(nullptr), fPyName(nullptr), fPyClass(nullptr),
    fNonTemplated(nullptr), fTemplated(nullptr), fLowPriority(nullptr)
{
    /* empty */
}

//----------------------------------------------------------------------------
TemplateInfo::~TemplateInfo()
{
    Py_XDECREF(fCppName);
    Py_XDECREF(fPyName);
    Py_XDECREF(fPyClass);

    Py_DECREF(fNonTemplated);
    Py_DECREF(fTemplated);
    Py_DECREF(fLowPriority);

    for (const auto& p : fDispatchMap) {
        for (const auto& c : p.second) {
            Py_DECREF(c.second);
        }
    }
}

//----------------------------------------------------------------------------
void TemplateProxy::Set(const std::string& cppname, const std::string& pyname, PyObject* pyclass)
{
// Initialize the proxy for the given 'pyclass.'
    fSelf         = nullptr;
    fTemplateArgs = nullptr;

    fTI->fCppName = CPyCppyy_PyText_FromString(const_cast<char*>(cppname.c_str()));
    fTI->fPyName  = CPyCppyy_PyText_FromString(const_cast<char*>(pyname.c_str()));
    Py_XINCREF(pyclass);
    fTI->fPyClass = pyclass;

    std::vector<PyCallable*> dummy;
    fTI->fNonTemplated = CPPOverload_New(pyname, dummy);
    fTI->fTemplated    = CPPOverload_New(pyname, dummy);
    fTI->fLowPriority  = CPPOverload_New(pyname, dummy);
}

//----------------------------------------------------------------------------
void TemplateProxy::MergeOverload(CPPOverload* mp) {
// Store overloads of this templated method.
    bool isGreedy = false;
    for (auto pc : mp->fMethodInfo->fMethods) {
        if (pc->IsGreedy()) {
            isGreedy = true;
            break;
        }
    }

    CPPOverload* cppol = isGreedy ? fTI->fLowPriority : fTI->fNonTemplated;
    cppol->MergeOverload(mp);
}

void TemplateProxy::AdoptMethod(PyCallable* pc) {
// Store overload of this templated method.
    CPPOverload* cppol = pc->IsGreedy() ? fTI->fLowPriority : fTI->fNonTemplated;
    cppol->AdoptMethod(pc);
}

void TemplateProxy::AdoptTemplate(PyCallable* pc)
{
// Store known template methods.
    fTI->fTemplated->AdoptMethod(pc);
}

//----------------------------------------------------------------------------
PyObject* TemplateProxy::Instantiate(const std::string& fname,
    PyObject* args, Utility::ArgPreference pref, int* pcnt)
{
// Instantiate (and cache) templated methods, return method if any
    std::string proto = "";

    Py_ssize_t nArgs = PyTuple_GET_SIZE(args);
    if (nArgs != 0) {
        PyObject* tpArgs = PyTuple_New(nArgs);
        for (int i = 0; i < nArgs; ++i) {
            PyObject* itemi = PyTuple_GET_ITEM(args, i);

            bool bArgSet = false;

        // special case for arrays
            PyObject* pytc = PyObject_GetAttr(itemi, PyStrings::gTypeCode);
            if (pytc) {
                PyObject* pyptrname = TC2CppName(pytc, "*", true);
                if (pyptrname) {
                    PyTuple_SET_ITEM(tpArgs, i, pyptrname);
                    bArgSet = true;
                // string added, but not counted towards nStrings
                }
                Py_DECREF(pytc); pytc = nullptr;
            } else
                PyErr_Clear();

        // if not arg set, try special case for ctypes
            if (!bArgSet) pytc = PyObject_GetAttr(itemi, PyStrings::gCTypesType);

            if (!bArgSet && pytc) {
                PyObject* pyactname = TC2CppName(pytc, "&", false);
                if (!pyactname) {
                // _type_ of a pointer to c_type is that type, which will have a type
                    PyObject* newpytc = PyObject_GetAttr(pytc, PyStrings::gCTypesType);
                    Py_DECREF(pytc);
                    pytc = newpytc;
                    if (pytc) {
                        pyactname = TC2CppName(pytc, "*", false);
                    } else
                        PyErr_Clear();
                }
                Py_DECREF(pytc); pytc = nullptr;
                if (pyactname) {
                    PyTuple_SET_ITEM(tpArgs, i, pyactname);
                    bArgSet = true;
                // string added, but not counted towards nStrings
                }
            } else
                PyErr_Clear();

            if (!bArgSet) {
            // normal case (may well fail)
                PyErr_Clear();
                PyObject* tp = (PyObject*)Py_TYPE(itemi);
                Py_INCREF(tp);
                PyTuple_SET_ITEM(tpArgs, i, tp);
            }
        }

        const std::string& name_v1 = Utility::ConstructTemplateArgs(nullptr, tpArgs, args, pref, 0, pcnt);
        Py_DECREF(tpArgs);
        if (name_v1.size())
            proto = name_v1.substr(1, name_v1.size()-2);
    }

// the following causes instantiation as necessary
    Cppyy::TCppScope_t scope = ((CPPClass*)fTI->fPyClass)->fCppType;
    Cppyy::TCppMethod_t cppmeth = Cppyy::GetMethodTemplate(scope, fname, proto);
    if (cppmeth) {    // overload stops here
    // A successful instantiation needs to be cached to pre-empt future instantiations. There
    // are two names involved, the original asked (which may be partial) and the received.
    //
    // Caching scheme: if the match is exact, simply add the overload to the pre-existing
    // one, or create a new overload for later lookups. If the match is not exact, do the
    // same, but also create an alias. Only add exact matches to the set of known template
    // instantiations, to prevent piling on from different partial instantiations.
    //
    // TODO: this caches the lookup method before the call, meaning that failing overloads
    // can add already existing overloads to the set of methods.

        std::string resname = Cppyy::GetMethodFullName(cppmeth);

    // An initializer_list is preferred for the argument types, but should not leak into
    // the argument types. If it did, replace with vector and lookup anew.
        if (resname.find("initializer_list") != std::string::npos) {
            auto pos = proto.find("initializer_list");
	    while (pos != std::string::npos) {
		proto.replace(pos, 16, "vector");
		pos = proto.find("initializer_list", pos + 6);
	    }

            Cppyy::TCppMethod_t m2 = Cppyy::GetMethodTemplate(scope, fname, proto);
            if (m2 && m2 != cppmeth) {
            // replace if the new method with vector was found; otherwise just continue
            // with the previously found method with initializer_list.
                cppmeth = m2;
                resname = Cppyy::GetMethodFullName(cppmeth);
            }
        }

        bool bExactMatch = fname == resname;

    // lookup on existing name in case this was an overload, not a caching, failure
        PyObject* dct = PyObject_GetAttr(fTI->fPyClass, PyStrings::gDict);
        PyObject* pycachename = CPyCppyy_PyText_InternFromString(fname.c_str());
        PyObject* pyol = PyObject_GetItem(dct, pycachename);
        if (!pyol) PyErr_Clear();
        bool bIsCppOL = CPPOverload_Check(pyol);

        if (pyol && !bIsCppOL && !TemplateProxy_Check(pyol)) {
        // unknown object ... leave well alone
            Py_DECREF(pyol);
            Py_DECREF(pycachename);
            Py_DECREF(dct);
            return nullptr;
        }

    // find the full name if the requested one was partial
        PyObject* exact = nullptr;
        PyObject* pyresname = CPyCppyy_PyText_FromString(resname.c_str());
        if (!bExactMatch) {
            exact = PyObject_GetItem(dct, pyresname);
            if (!exact) PyErr_Clear();
        }
        Py_DECREF(dct);

        bool bIsConstructor = false, bNeedsRebind = true;

        PyCallable* meth = nullptr;
        if (Cppyy::IsNamespace(scope)) {
            meth = new CPPFunction(scope, cppmeth);
            bNeedsRebind = false;
        } else if (Cppyy::IsStaticMethod(cppmeth)) {
            meth = new CPPClassMethod(scope, cppmeth);
            bNeedsRebind = false;
        } else if (Cppyy::IsConstructor(cppmeth)) {
            bIsConstructor = true;
            meth = new CPPConstructor(scope, cppmeth);
        } else
            meth = new CPPMethod(scope, cppmeth);

    // Case 1/2: method simply did not exist before
        if (!pyol) {
        // actual overload to use (now owns meth)
            pyol = (PyObject*)CPPOverload_New(fname, meth);
            if (bIsConstructor) {
            // TODO: this is an ugly hack :(
                ((CPPOverload*)pyol)->fMethodInfo->fFlags |= \
                    CallContext::kIsCreator | CallContext::kIsConstructor;
            }

        // add to class dictionary
            PyType_Type.tp_setattro(fTI->fPyClass, pycachename, pyol);
        }

    // Case 3/4: pre-existing method that was either not found b/c the full
    // templated name was constructed in this call or it failed as overload
        else if (bIsCppOL) {
        // TODO: see above, since the call hasn't happened yet, this overload may
        // already exist and fail again.
            ((CPPOverload*)pyol)->AdoptMethod(meth);   // takes ownership
        }

    // Case 5: must be a template proxy, meaning that current template name is not
    // a template overload
        else {
            ((TemplateProxy*)pyol)->AdoptTemplate(meth->Clone());
            Py_DECREF(pyol);
            pyol = (PyObject*)CPPOverload_New(fname, meth);      // takes ownership
        }

    // Special Case if name was aliased (e.g. typedef in template instantiation)
        if (!exact && !bExactMatch) {
            PyType_Type.tp_setattro(fTI->fPyClass, pyresname, pyol);
        }

    // cleanup
        Py_DECREF(pyresname);
        Py_DECREF(pycachename);

    // retrieve fresh (for boundedness) and call
        PyObject* pymeth =
            CPPOverload_Type.tp_descr_get(pyol, bNeedsRebind ? fSelf : nullptr, (PyObject*)&CPPOverload_Type);
        Py_DECREF(pyol);
        return pymeth;
    }

    PyErr_Format(PyExc_TypeError, "Failed to instantiate \"%s(%s)\"", fname.c_str(), proto.c_str());
    return nullptr;
}


//= CPyCppyy template proxy construction/destruction =========================
static TemplateProxy* tpp_new(PyTypeObject*, PyObject*, PyObject*)
{
// Create a new empty template method proxy.
    TemplateProxy* pytmpl = PyObject_GC_New(TemplateProxy, &TemplateProxy_Type);
    pytmpl->fSelf         = nullptr;
    pytmpl->fTemplateArgs = nullptr;
    pytmpl->fWeakrefList  = nullptr;
    new (&pytmpl->fTI) TP_TInfo_t{};
    pytmpl->fTI = std::make_shared<TemplateInfo>();

    PyObject_GC_Track(pytmpl);
    return pytmpl;
}

//----------------------------------------------------------------------------
static int tpp_clear(TemplateProxy* pytmpl)
{
// Garbage collector clear of held python member objects.
    Py_CLEAR(pytmpl->fSelf);
    Py_CLEAR(pytmpl->fTemplateArgs);

    return 0;
}

//----------------------------------------------------------------------------
static void tpp_dealloc(TemplateProxy* pytmpl)
{
// Destroy the given template method proxy.
    if (pytmpl->fWeakrefList)
        PyObject_ClearWeakRefs((PyObject*)pytmpl);
    PyObject_GC_UnTrack(pytmpl);
    tpp_clear(pytmpl);
    pytmpl->fTI.~TP_TInfo_t();
    PyObject_GC_Del(pytmpl);
}

//----------------------------------------------------------------------------
static int tpp_traverse(TemplateProxy* pytmpl, visitproc visit, void* arg)
{
// Garbage collector traverse of held python member objects.
    Py_VISIT(pytmpl->fSelf);
    Py_VISIT(pytmpl->fTemplateArgs);

    return 0;
}

//----------------------------------------------------------------------------
static PyObject* tpp_doc(TemplateProxy* pytmpl, void*)
{
// Forward to method proxies to doc all overloads
    PyObject* doc = nullptr;
    if (pytmpl->fTI->fNonTemplated->HasMethods())
        doc = PyObject_GetAttrString((PyObject*)pytmpl->fTI->fNonTemplated, "__doc__");
    if (pytmpl->fTI->fTemplated->HasMethods()) {
        PyObject* doc2 = PyObject_GetAttrString((PyObject*)pytmpl->fTI->fTemplated, "__doc__");
        if (doc && doc2) {
            CPyCppyy_PyText_AppendAndDel(&doc, CPyCppyy_PyText_FromString("\n"));
            CPyCppyy_PyText_AppendAndDel(&doc, doc2);
        } else if (!doc && doc2) {
            doc = doc2;
        }
    }
    if (pytmpl->fTI->fLowPriority->HasMethods()) {
        PyObject* doc2 = PyObject_GetAttrString((PyObject*)pytmpl->fTI->fLowPriority, "__doc__");
        if (doc && doc2) {
            CPyCppyy_PyText_AppendAndDel(&doc, CPyCppyy_PyText_FromString("\n"));
            CPyCppyy_PyText_AppendAndDel(&doc, doc2);
        } else if (!doc && doc2) {
            doc = doc2;
        }
    }

    if (doc)
        return doc;

    return CPyCppyy_PyText_FromString(TemplateProxy_Type.tp_doc);
}

//----------------------------------------------------------------------------
static PyObject* tpp_repr(TemplateProxy* pytmpl)
{
// Simply return the doc string as that's the most useful info (this will appear
// on clsses on calling help()).
     return tpp_doc(pytmpl, nullptr);
}


//= CPyCppyy template proxy callable behavior ================================
static inline std::string targs2str(TemplateProxy* pytmpl)
{
    if (!pytmpl || !pytmpl->fTemplateArgs) return "";
    return CPyCppyy_PyText_AsString(pytmpl->fTemplateArgs);
}

static inline void UpdateDispatchMap(TemplateProxy* pytmpl, bool use_targs, uint64_t sighash, CPPOverload* pymeth)
{
// memoize a method in the dispatch map after successful call; replace old if need be (may be
// with the same CPPOverload, just with more methods)
    bool bInserted = false;
    auto& v = pytmpl->fTI->fDispatchMap[use_targs ? targs2str(pytmpl) : ""];

    Py_INCREF(pymeth);
    for (auto& p : v) {
        if (p.first == sighash) {
            Py_DECREF(p.second);
            p.second = pymeth;
            bInserted = true;
        }
    }
    if (!bInserted) v.push_back(std::make_pair(sighash, pymeth));
}

static inline PyObject* CallMethodImp(TemplateProxy* pytmpl, PyObject*& pymeth,
    PyObject* args, PyObject* kwds, bool impOK, uint64_t sighash)
{
// Actual call of a given overload: takes care of handlign of "self" and
// dereferences the overloaded method after use.
    PyObject* result;
    if (!impOK) PyDict_SetItem(kwds, PyStrings::gNoImplicit, Py_True);
    bool isNS = (((CPPScope*)pytmpl->fTI->fPyClass)->fFlags & CPPScope::kIsNamespace);
    if (isNS && pytmpl->fSelf) {
    // this is a global method added a posteriori to the class
        Py_ssize_t sz = PyTuple_GET_SIZE(args);
        PyObject* newArgs = PyTuple_New(sz+1);
        for (int i = 0; i < sz; ++i) {
            PyObject* item = PyTuple_GET_ITEM(args, i);
            Py_INCREF(item);
            PyTuple_SET_ITEM(newArgs, i+1, item);
        }
        Py_INCREF((PyObject*)pytmpl->fSelf);
        PyTuple_SET_ITEM(newArgs, 0, (PyObject*)pytmpl->fSelf);
        result = CPPOverload_Type.tp_call(pymeth, newArgs, kwds);
        Py_DECREF(newArgs);
    } else
        result = CPPOverload_Type.tp_call(pymeth, args, kwds);

    if (result) {
        Py_XDECREF(((CPPOverload*)pymeth)->fSelf); ((CPPOverload*)pymeth)->fSelf = nullptr;    // unbind
        UpdateDispatchMap(pytmpl, true, sighash, (CPPOverload*)pymeth);
    }

    Py_DECREF(pymeth); pymeth = nullptr;
    return result;
}

#define TPPCALL_RETURN                                                       \
{ if (!errors.empty())                                                       \
      std::for_each(errors.begin(), errors.end(), Utility::PyError_t::Clear);\
  Py_DECREF(kwds);                                                           \
  return result; }

static PyObject* tpp_call(TemplateProxy* pytmpl, PyObject* args, PyObject* kwds)
{
// Dispatcher to the actual member method, several uses possible; in order:
//
// case 1: explicit template previously selected through subscript
//
// case 2: select known non-template overload
//
//    obj.method(a0, a1, ...)
//       => obj->method(a0, a1, ...)        // non-template
//
// case 3: select known template overload
//
//    obj.method(a0, a1, ...)
//       => obj->method(a0, a1, ...)        // all known templates
//
// case 4: auto-instantiation from types of arguments
//
//    obj.method(a0, a1, ...)
//       => obj->method<type(a0), type(a1), ...>(a0, a1, ...)
//
// Note: explicit instantiation needs to use [] syntax:
//
//    obj.method[type<a0>, type<a1>, ...](a0, a1, ...)
//
// case 5: low priority methods, such as ones that take void* arguments
//

// TODO: should previously instantiated templates be considered first?

// container for collecting errors
    std::vector<Utility::PyError_t> errors;

    PyObject* pymeth = nullptr, *result = nullptr;

// short-cut through memoization map
    uint64_t sighash = HashSignature(args);

    if (!pytmpl->fTemplateArgs) {
    // look for known signatures ...
        CPPOverload* ol = nullptr;
        auto& v = pytmpl->fTI->fDispatchMap[targs2str(pytmpl)];
        for (const auto& p : v) {
            if (p.first == sighash) {
                ol = p.second;
                break;
            }
        }

        if (ol != nullptr) {
            if (!pytmpl->fSelf) {
                result = CPPOverload_Type.tp_call((PyObject*)ol, args, kwds);
            } else {
                pymeth = CPPOverload_Type.tp_descr_get(
                    (PyObject*)ol, pytmpl->fSelf, (PyObject*)&CPPOverload_Type);
                result = CPPOverload_Type.tp_call(pymeth, args, kwds);
                Py_DECREF(pymeth); pymeth = nullptr;
            }
            if (result)
                return result;
            Utility::FetchError(errors);
        }
    }

// do not mix template instantiations with implicit conversions
    if (!kwds) kwds = PyDict_New();
    else {
        Py_INCREF(kwds);
    }

// case 1: explicit template previously selected through subscript
    if (pytmpl->fTemplateArgs) {
    // instantiate explicitly
        PyObject* pyfullname = CPyCppyy_PyText_FromString(
            CPyCppyy_PyText_AsString(pytmpl->fTI->fCppName));
        CPyCppyy_PyText_Append(&pyfullname, pytmpl->fTemplateArgs);

    // first, lookup by full name, if previously stored
        bool isNS = (((CPPScope*)pytmpl->fTI->fPyClass)->fFlags & CPPScope::kIsNamespace);
        pymeth = PyObject_GetAttr((pytmpl->fSelf && !isNS) ? pytmpl->fSelf : pytmpl->fTI->fPyClass, pyfullname);

    // attempt call if found (this may fail if there are specializations)
        if (CPPOverload_Check(pymeth)) {
        // since the template args are fully explicit, allow implicit conversion of arguments
            result = CallMethodImp(pytmpl, pymeth, args, kwds, true, sighash);
            if (result) {
                Py_DECREF(pyfullname);
                TPPCALL_RETURN;
            }
            Utility::FetchError(errors);
        } else if (pymeth && PyCallable_Check(pymeth)) {
        // something different (user provided?)
            result = PyObject_CallObject(pymeth, args);
            Py_DECREF(pymeth);
            if (result) {
                Py_DECREF(pyfullname);
                TPPCALL_RETURN;
            }
            Utility::FetchError(errors);
        } else if (!pymeth)
            PyErr_Clear();

    // not cached or failed call; try instantiation
        pymeth = pytmpl->Instantiate(
            CPyCppyy_PyText_AsString(pyfullname), args, Utility::kNone);
        if (pymeth) {
        // attempt actuall call; same as above, allow implicit conversion of arguments
            result = CallMethodImp(pytmpl, pymeth, args, kwds, true, sighash);
            if (result) {
                Py_DECREF(pyfullname);
                TPPCALL_RETURN;
            }
        }

    // no drop through if failed (if implicit was desired, don't provide template args)
        Utility::FetchError(errors);
        PyObject* topmsg = CPyCppyy_PyText_FromFormat("Could not instantiate %s:", CPyCppyy_PyText_AsString(pyfullname));
        Py_DECREF(pyfullname);
        Utility::SetDetailedException(errors, topmsg /* steals */, PyExc_TypeError /* default error */);

        Py_DECREF(kwds);
        return nullptr;
    }

// case 2: select known non-template overload
    if (pytmpl->fTI->fNonTemplated->HasMethods()) {
    // simply forward the call: all non-templated methods are defined on class definition
    // and thus already available
        pymeth = CPPOverload_Type.tp_descr_get(
           (PyObject*)pytmpl->fTI->fNonTemplated, pytmpl->fSelf, (PyObject*)&CPPOverload_Type);
    // now call the method with the arguments (loops internally and implicit is okay as
    // these are not templated methods that should match exactly)
        result = CPPOverload_Type.tp_call(pymeth, args, kwds);
        Py_DECREF(pymeth); pymeth = nullptr;
        if (result) {
            UpdateDispatchMap(pytmpl, false, sighash, pytmpl->fTI->fNonTemplated);
            TPPCALL_RETURN;
        }
        Utility::FetchError(errors);
    }

// case 3: select known template overload
    if (pytmpl->fTI->fTemplated->HasMethods()) {
    // simply forward the call
        pymeth = CPPOverload_Type.tp_descr_get(
            (PyObject*)pytmpl->fTI->fTemplated, pytmpl->fSelf, (PyObject*)&CPPOverload_Type);
    // now call the method with the arguments (loops internally)
        PyDict_SetItem(kwds, PyStrings::gNoImplicit, Py_True);
        result = CPPOverload_Type.tp_call(pymeth, args, kwds);
        Py_DECREF(pymeth); pymeth = nullptr;
        if (result) {
            UpdateDispatchMap(pytmpl, true, sighash, pytmpl->fTI->fTemplated);
            TPPCALL_RETURN;
        }
        Utility::FetchError(errors);
    }

// case 4: auto-instantiation from types of arguments
    for (auto pref : {Utility::kReference, Utility::kPointer, Utility::kValue}) {
        // TODO: no need to loop if there are no non-instance arguments; also, should any
        // failed lookup be removed?
        int pcnt = 0;
        pymeth = pytmpl->Instantiate(
            CPyCppyy_PyText_AsString(pytmpl->fTI->fCppName), args, pref, &pcnt);
        if (pymeth) {
        // attempt actuall call; argument based, so do not allow implicit conversions
            result = CallMethodImp(pytmpl, pymeth, args, kwds, false, sighash);
            if (result) TPPCALL_RETURN;
        }
        Utility::FetchError(errors);
        if (!pcnt) break;         // preference never used; no point trying others
    }

// case 5: low priority methods, such as ones that take void* arguments
    if (pytmpl->fTI->fLowPriority->HasMethods()) {
    // simply forward the call
        pymeth = CPPOverload_Type.tp_descr_get(
            (PyObject*)pytmpl->fTI->fLowPriority, pytmpl->fSelf, (PyObject*)&CPPOverload_Type);
    // now call the method with the arguments (loops internally)
        PyDict_SetItem(kwds, PyStrings::gNoImplicit, Py_True);
        result = CPPOverload_Type.tp_call(pymeth, args, kwds);
        Py_DECREF(pymeth); pymeth = nullptr;
        if (result) {
            UpdateDispatchMap(pytmpl, false, sighash, pytmpl->fTI->fLowPriority);
            TPPCALL_RETURN;
        }
        Utility::FetchError(errors);
    }

// error reporting is fraud, given the numerous steps taken, but more details seems better
    if (!errors.empty()) {
        PyObject* topmsg = CPyCppyy_PyText_FromString("Template method resolution failed:");
        Utility::SetDetailedException(errors, topmsg /* steals */, PyExc_TypeError /* default error */);
    } else {
        PyErr_Format(PyExc_TypeError, "cannot resolve method template call for \'%s\'",
            CPyCppyy_PyText_AsString(pytmpl->fTI->fPyName));
    }

    Py_DECREF(kwds);
    return nullptr;
}

//----------------------------------------------------------------------------
static TemplateProxy* tpp_descrget(TemplateProxy* pytmpl, PyObject* pyobj, PyObject*)
{
// create and use a new template proxy (language requirement)
    TemplateProxy* newPyTmpl = (TemplateProxy*)TemplateProxy_Type.tp_alloc(&TemplateProxy_Type, 0);

// new method is to be bound to current object (may be nullptr)
    Py_XINCREF(pyobj);
    newPyTmpl->fSelf = pyobj;

    Py_XINCREF(pytmpl->fTemplateArgs);
    newPyTmpl->fTemplateArgs = pytmpl->fTemplateArgs;

// copy name, class, etc. pointers
    new (&newPyTmpl->fTI) std::shared_ptr<TemplateInfo>{pytmpl->fTI};

    return newPyTmpl;
}

//----------------------------------------------------------------------------
static PyObject* tpp_subscript(TemplateProxy* pytmpl, PyObject* args)
{
// Explicit template member lookup/instantiation; works by re-bounding. This method can
// not cache overloads as instantiations need not be unique for the argument types due
// to template specializations.
    TemplateProxy* typeBoundMethod = tpp_descrget(pytmpl, pytmpl->fSelf, nullptr);
    Py_XDECREF(typeBoundMethod->fTemplateArgs);
    typeBoundMethod->fTemplateArgs = CPyCppyy_PyText_FromString(
        Utility::ConstructTemplateArgs(nullptr, args).c_str());
    return (PyObject*)typeBoundMethod;
}

//-----------------------------------------------------------------------------
static PyObject* tpp_getuseffi(CPPOverload*, void*)
{   
    return PyInt_FromLong(0); // dummy (__useffi__ unused)
}   
    
//-----------------------------------------------------------------------------
static int tpp_setuseffi(CPPOverload*, PyObject*, void*)
{   
    return 0;                 // dummy (__useffi__ unused)
}   


//----------------------------------------------------------------------------
static PyMappingMethods tpp_as_mapping = {
    nullptr, (binaryfunc)tpp_subscript, nullptr
};

static PyGetSetDef tpp_getset[] = {
    {(char*)"__doc__", (getter)tpp_doc, nullptr, nullptr, nullptr},
    {(char*)"__useffi__", (getter)tpp_getuseffi, (setter)tpp_setuseffi,
      (char*)"unused", nullptr},
    {(char*)nullptr,   nullptr,         nullptr, nullptr, nullptr}
};


//= CPyCppyy template proxy type =============================================
PyTypeObject TemplateProxy_Type = {
   PyVarObject_HEAD_INIT(&PyType_Type, 0)
   (char*)"cppyy.TemplateProxy", // tp_name
   sizeof(TemplateProxy),     // tp_basicsize
   0,                         // tp_itemsize
   (destructor)tpp_dealloc,   // tp_dealloc
   0,                         // tp_print
   0,                         // tp_getattr
   0,                         // tp_setattr
   0,                         // tp_compare
   (reprfunc)tpp_repr,        // tp_repr
   0,                         // tp_as_number
   0,                         // tp_as_sequence
   &tpp_as_mapping,           // tp_as_mapping
   0,                         // tp_hash
   (ternaryfunc)tpp_call,     // tp_call
   0,                         // tp_str
   0,                         // tp_getattro
   0,                         // tp_setattro
   0,                         // tp_as_buffer
   Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_GC,      // tp_flags
   (char*)"cppyy template proxy (internal)",     // tp_doc
   (traverseproc)tpp_traverse,// tp_traverse
   (inquiry)tpp_clear,        // tp_clear
   0,                         // tp_richcompare
   offsetof(TemplateProxy, fWeakrefList),        // tp_weaklistoffset
   0,                         // tp_iter
   0,                         // tp_iternext
   0,                         // tp_methods
   0,                         // tp_members
   tpp_getset,                // tp_getset
   0,                         // tp_base
   0,                         // tp_dict
   (descrgetfunc)tpp_descrget,// tp_descr_get
   0,                         // tp_descr_set
   0,                         // tp_dictoffset
   0,                         // tp_init
   0,                         // tp_alloc
   (newfunc)tpp_new,          // tp_new
   0,                         // tp_free
   0,                         // tp_is_gc
   0,                         // tp_bases
   0,                         // tp_mro
   0,                         // tp_cache
   0,                         // tp_subclasses
   0                          // tp_weaklist
#if PY_VERSION_HEX >= 0x02030000
   , 0                        // tp_del
#endif
#if PY_VERSION_HEX >= 0x02060000
   , 0                        // tp_version_tag
#endif
#if PY_VERSION_HEX >= 0x03040000
   , 0                        // tp_finalize
#endif
};

} // namespace CPyCppyy
