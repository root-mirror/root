// Bindings
#include "CPyCppyy.h"
#include "CPPMethod.h"
#include "CPPExcInstance.h"
#include "CPPInstance.h"
#include "Converters.h"
#include "Executors.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"
#include "TypeManip.h"
#include "SignalTryCatch.h"
#include "Utility.h"

#include "CPyCppyy/PyException.h"

// Standard
#include <assert.h>
#include <string.h>
#include <exception>
#include <sstream>
#include <string>
#include <typeinfo>
#include <memory>


//- data and local helpers ---------------------------------------------------
namespace CPyCppyy {
    extern PyObject* gThisModule;
    extern PyObject* gBusException;
    extern PyObject* gSegvException;
    extern PyObject* gIllException;
    extern PyObject* gAbrtException;
}


//- private helpers ----------------------------------------------------------
inline void CPyCppyy::CPPMethod::Copy_(const CPPMethod& /* other */)
{
// fScope and fMethod handled separately

// do not copy caches
    fExecutor     = nullptr;
    fArgIndices   = nullptr;
    fArgsRequired = -1;
}

//----------------------------------------------------------------------------
inline void CPyCppyy::CPPMethod::Destroy_() const
{
// destroy executor and argument converters
    if (fExecutor && fExecutor->HasState()) delete fExecutor;

    for (auto p : fConverters) {
        if (p && p->HasState()) delete p;
    }

    delete fArgIndices;
}

//----------------------------------------------------------------------------
inline PyObject* CPyCppyy::CPPMethod::ExecuteFast(
    void* self, ptrdiff_t offset, CallContext* ctxt)
{
// call into C++ through fExecutor; abstracted out from Execute() to prevent some
// code duplication with ProtectedCall()
    PyObject* result = nullptr;

    try {       // C++ try block
        result = fExecutor->Execute(fMethod, (Cppyy::TCppObject_t)((intptr_t)self+offset), ctxt);
    } catch (PyException&) {
        result = nullptr;           // error already set
    } catch (std::exception& e) {
    // attempt to set the exception to the actual type, to allow catching with the Python C++ type
        static Cppyy::TCppType_t exc_type = (Cppyy::TCppType_t)Cppyy::GetScope("std::exception");

        PyObject* pyexc_type = nullptr;
        PyObject* pyexc_obj  = nullptr;

    // TODO: factor this code with the same in ProxyWrappers (and cache it there to be able to
    // look up based on TCppType_t):
        Cppyy::TCppType_t actual = Cppyy::GetActualClass(exc_type, &e);
        const std::string& finalname = Cppyy::GetScopedFinalName(actual);
        const std::string& parentname = TypeManip::extract_namespace(finalname);
        PyObject* parent = CreateScopeProxy(parentname);
        if (parent) {
            pyexc_type = PyObject_GetAttrString(parent,
                parentname.empty() ? finalname.c_str() : finalname.substr(parentname.size()+2, std::string::npos).c_str());
            Py_DECREF(parent);
        }

        if (pyexc_type) {
        // create a copy of the exception (TODO: factor this code with the same in ProxyWrappers)
            PyObject* pyclass = CPyCppyy::GetScopeProxy(actual);
            PyObject* source = BindCppObjectNoCast(&e, actual);
            PyObject* pyexc_copy = PyObject_CallFunctionObjArgs(pyclass, source, nullptr);
            Py_DECREF(source);
            Py_DECREF(pyclass);
            if (pyexc_copy) {
                pyexc_obj = CPPExcInstance_Type.tp_new((PyTypeObject*)pyexc_type, nullptr, nullptr);
                ((CPPExcInstance*)pyexc_obj)->fCppInstance = (PyObject*)pyexc_copy;
            } else
                PyErr_Clear();
        } else
            PyErr_Clear();

        if (pyexc_type && pyexc_obj) {
            PyErr_SetObject(pyexc_type, pyexc_obj);
            Py_DECREF(pyexc_obj);
            Py_DECREF(pyexc_type);
        } else {
            PyErr_Format(PyExc_Exception, "%s (C++ exception)", e.what());
            Py_XDECREF(pyexc_obj);
            Py_XDECREF(pyexc_type);
        }

        result = nullptr;
    } catch (...) {
        PyErr_SetString(PyExc_Exception, "unhandled, unknown C++ exception");
        result = nullptr;
    }

// TODO: covers the PyException throw case, which does not seem to work on Windows, so
// instead leaves the error be
#ifdef _WIN32
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = nullptr;
    }
#endif

    return result;
}

//----------------------------------------------------------------------------
inline PyObject* CPyCppyy::CPPMethod::ExecuteProtected(
    void* self, ptrdiff_t offset, CallContext* ctxt)
{
// helper code to prevent some code duplication; this code embeds a "try/catch"
// block that saves the call environment for restoration in case of an otherwise
// fatal signal
    PyObject* result = 0;

    TRY {     // copy call environment to be able to jump back on signal
        result = ExecuteFast(self, offset, ctxt);
    } CATCH(excode) {
    // Unfortunately, the excodes are not the ones from signal.h, but enums from TSysEvtHandler.h
        if (excode == 0)
            PyErr_SetString(gBusException, "bus error in C++; program state was reset");
        else if (excode == 1)
            PyErr_SetString(gSegvException, "segfault in C++; program state was reset");
        else if (excode == 4)
            PyErr_SetString(gIllException, "illegal instruction in C++; program state was reset");
        else if (excode == 5)
            PyErr_SetString(gAbrtException, "abort from C++; program state was reset");
        else if (excode == 12)
            PyErr_SetString(PyExc_FloatingPointError, "floating point exception in C++; program state was reset");
        else
            PyErr_SetString(PyExc_SystemError, "problem in C++; program state was reset");
        result = 0;
    } ENDTRY;

    return result;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CPPMethod::InitConverters_()
{
// build buffers for argument dispatching
    const size_t nArgs = Cppyy::GetMethodNumArgs(fMethod);
    fConverters.resize(nArgs);

// setup the dispatch cache
    for (int iarg = 0; iarg < (int)nArgs; ++iarg) {
        const std::string& fullType = Cppyy::GetMethodArgType(fMethod, iarg);
        Converter* conv = CreateConverter(fullType);
        if (!conv) {
            PyErr_Format(PyExc_TypeError, "argument type %s not handled", fullType.c_str());
            return false;
        }

        fConverters[iarg] = conv;
    }

    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CPPMethod::InitExecutor_(Executor*& executor, CallContext* /* ctxt */)
{
// install executor conform to the return type
    executor = CreateExecutor(
        (bool)fMethod == true ? Cppyy::GetMethodResultType(fMethod) \
                              : Cppyy::GetScopedFinalName(fScope));

    if (!executor)
        return false;

    return true;
}

//----------------------------------------------------------------------------
std::string CPyCppyy::CPPMethod::GetSignatureString(bool fa)
{
// built a signature representation (used for doc strings)
    std::stringstream sig; sig << "(";
    int count = 0;
    const size_t nArgs = Cppyy::GetMethodNumArgs(fMethod);
    for (int iarg = 0; iarg < (int)nArgs; ++iarg) {
        if (count) sig << (fa ? ", " : ",");

        sig << Cppyy::GetMethodArgType(fMethod, iarg);

        if (fa) {
            const std::string& parname = Cppyy::GetMethodArgName(fMethod, iarg);
            if (!parname.empty())
                sig << " " << parname;

            const std::string& defvalue = Cppyy::GetMethodArgDefault(fMethod, iarg);
            if (!defvalue.empty())
                sig << " = " << defvalue;
        }
        count++;
    }
    sig << ")";
    return sig.str();
}

//----------------------------------------------------------------------------
void CPyCppyy::CPPMethod::SetPyError_(PyObject* msg)
{
// helper to report errors in a consistent format (derefs msg)
    std::string details{};

    PyObject *etype = nullptr, *evalue = nullptr;
    if (PyErr_Occurred()) {
        PyObject* etrace = nullptr;

        PyErr_Fetch(&etype, &evalue, &etrace);

        if (evalue) {
            PyObject* descr = PyObject_Str(evalue);
            if (descr) {
                details = CPyCppyy_PyText_AsString(descr);
                Py_DECREF(descr);
            }
        }

        Py_XDECREF(etrace);
    }

    PyObject* doc = GetDocString();
    PyObject* errtype = etype;
    if (!errtype)
        errtype = PyExc_TypeError;
    PyObject* pyname = PyObject_GetAttr(errtype, PyStrings::gName);
    const char* cname = pyname ? CPyCppyy_PyText_AsString(pyname) : "Exception";

    if (!PyType_IsSubtype((PyTypeObject*)errtype, &CPPExcInstance_Type)) {
        if (details.empty()) {
            PyErr_Format(errtype, "%s =>\n    %s: %s", CPyCppyy_PyText_AsString(doc),
                 cname, msg ? CPyCppyy_PyText_AsString(msg) : "");
        } else if (msg) {
            PyErr_Format(errtype, "%s =>\n    %s: %s (%s)",
                CPyCppyy_PyText_AsString(doc), cname, CPyCppyy_PyText_AsString(msg),
                details.c_str());
        } else {
            PyErr_Format(errtype, "%s =>\n    %s: %s",
                CPyCppyy_PyText_AsString(doc), cname, details.c_str());
        }
    } else {
        Py_XDECREF(((CPPExcInstance*)evalue)->fTopMessage);
        if (msg) {
            ((CPPExcInstance*)evalue)->fTopMessage = CPyCppyy_PyText_FromFormat(\
                "%s =>\n    %s: %s | ", CPyCppyy_PyText_AsString(doc), cname, CPyCppyy_PyText_AsString(msg));
        } else {
            ((CPPExcInstance*)evalue)->fTopMessage = CPyCppyy_PyText_FromFormat(\
                 "%s =>\n    %s: ", CPyCppyy_PyText_AsString(doc), cname);
        }
        PyErr_SetObject(errtype, evalue);
    }

    Py_XDECREF(pyname);
    Py_XDECREF(evalue);
    Py_XDECREF(etype);
    Py_DECREF(doc);
    Py_XDECREF(msg);
}

//- constructors and destructor ----------------------------------------------
CPyCppyy::CPPMethod::CPPMethod(
        Cppyy::TCppScope_t scope, Cppyy::TCppMethod_t method) :
    fMethod(method), fScope(scope), fExecutor(nullptr), fArgIndices(nullptr),
    fArgsRequired(-1)
{
   // empty
}

//----------------------------------------------------------------------------
CPyCppyy::CPPMethod::CPPMethod(const CPPMethod& other) :
    PyCallable(other), fMethod(other.fMethod), fScope(other.fScope)
{
    Copy_(other);
}

//----------------------------------------------------------------------------
CPyCppyy::CPPMethod& CPyCppyy::CPPMethod::operator=(const CPPMethod& other)
{
    if (this != &other) {
        Destroy_();
        Copy_(other);
        fScope  = other.fScope;
        fMethod = other.fMethod;
    }

    return *this;
}

//----------------------------------------------------------------------------
CPyCppyy::CPPMethod::~CPPMethod()
{
    Destroy_();
}


//- public members -----------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::GetPrototype(bool fa)
{
// construct python string from the method's prototype
    return CPyCppyy_PyText_FromFormat("%s%s %s::%s%s",
        (Cppyy::IsStaticMethod(fMethod) ? "static " : ""),
        Cppyy::GetMethodResultType(fMethod).c_str(),
        Cppyy::GetScopedFinalName(fScope).c_str(), Cppyy::GetMethodName(fMethod).c_str(),
        GetSignatureString(fa).c_str());
}

//----------------------------------------------------------------------------
int CPyCppyy::CPPMethod::GetPriority()
{
// To help with overload selection, methods are given a priority based on the
// affinity of Python and C++ types. Priority only matters for methods that have
// an equal number of arguments and types that are possible substitutes (the
// normal selection mechanisms would simply distinguish them otherwise).

// The following types are ordered, in favor (variants implicit):
//
//  bool >> long >> int >> short
//  double >> long double >> float
//  const char* >> char
//
// Further, all integer types are preferred over floating point b/c int to float
// is allowed implicitly, float to int is not.
//
// Special cases that are disliked include void* and unknown/incomplete types.
// Also, moves are preferred over references. std::initializer_list is not a nice
// conversion candidate either, but needs to be higher priority to mix well with
// implicit conversions.
// TODO: extend this to favour classes that are not bases.
// TODO: profile this method (it's expensive, but should be called too often)

    int priority = 0;

    const size_t nArgs = Cppyy::GetMethodNumArgs(fMethod);
    for (int iarg = 0; iarg < (int)nArgs; ++iarg) {
        const std::string aname = Cppyy::GetMethodArgType(fMethod, iarg);

        if (Cppyy::IsBuiltin(aname)) {
        // integer types
            if (strstr(aname.c_str(), "bool"))
                priority +=    1;      // bool over int (does accept 1 and 0)
            else if (strstr(aname.c_str(), "long long"))
                priority +=   -5;      // will very likely fit
            else if (strstr(aname.c_str(), "long"))
                priority +=  -10;      // most affine integer type
            // no need to compare with int; leave at zero
            else if (strstr(aname.c_str(), "short"))
                priority +=  -50;      // not really relevant as a type

        // floating point types (note all numbers lower than integer types)
            else if (strstr(aname.c_str(), "float"))
                priority += -100;      // not really relevant as a type
            else if (strstr(aname.c_str(), "long double"))
                priority +=  -90;      // fits double with least loss of precision
            else if (strstr(aname.c_str(), "double"))
                priority +=  -80;      // most affine floating point type

        // string/char types
            else if (strstr(aname.c_str(), "char") && aname[aname.size()-1] != '*')
                priority += -60;       // prefer (const) char* over char

        // oddball
            else if (strstr(aname.c_str(), "void*"))
                priority -= 1000;      // void*/void** shouldn't be too greedy

        } else {
        // This is a user-defined type (class, struct, enum, etc.).

        // There's a bit of hysteresis here for templates: once GetScope() is called, their
        // IsComplete() succeeds, the other way around it does not. Since GetPriority() is
        // likely called several times in a sort, the GetScope() _must_ come first, or
        // different GetPriority() calls may return different results (since the 2nd time,
        // GetScope() will have been called from the first), killing the stable_sort.

        // prefer more derived classes
            const std::string& clean_name = TypeManip::clean_type(aname, false);
            Cppyy::TCppScope_t scope = Cppyy::GetScope(clean_name);
            if (scope)
                priority += (int)Cppyy::GetNumBases(scope);

            if (Cppyy::IsEnum(clean_name))
                priority -= 100;

        // a couple of special cases as explained above
            if (aname.find("initializer_list") != std::string::npos) {
                priority +=   150;     // needed for proper implicit conversion rules
            } else if (aname.rfind("&&", aname.size()-2) != std::string::npos) {
                priority +=   100;     // prefer moves over other ref/ptr
            } else if (!aname.empty() && !Cppyy::IsComplete(aname)) {
            // class is known, but no dictionary available, 2 more cases: * and &
                if (aname[aname.size() - 1] == '&')
                    priority += -5000;
                else
                    priority += -2000; // prefer pointer passing over reference
            }
        }
    }

// prefer methods w/o optional arguments b/c ones with optional arguments are easier to
// select by providing the optional arguments explicitly
    priority += ((int)Cppyy::GetMethodReqArgs(fMethod) - (int)nArgs);

// add a small penalty to prefer non-const methods over const ones for get/setitem
    if (Cppyy::IsConstMethod(fMethod) && Cppyy::GetMethodName(fMethod) == "operator[]")
        priority += -10;

    return priority;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CPPMethod::IsGreedy()
{
// Methods will all void*-like arguments should be sorted after template
// instanstations, so that they don't greedily take over pointers to object.
// GetPriority() is too heavy-handed, as it will pull in all the argument
// types, so use this cheaper check.
    const size_t nArgs = Cppyy::GetMethodReqArgs(fMethod);
    if (!nArgs) return false;

    for (int iarg = 0; iarg < (int)nArgs; ++iarg) {
        const std::string aname = Cppyy::GetMethodArgType(fMethod, iarg);
        if (aname.find("void*") != 0)
            return false;
    }
    return true;
}


//----------------------------------------------------------------------------
int CPyCppyy::CPPMethod::GetMaxArgs()
{
    return (int)Cppyy::GetMethodNumArgs(fMethod);
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::GetCoVarNames()
{
// Build a tuple of the argument types/names.
    int co_argcount = (int)GetMaxArgs() /* +1 for self */;

// TODO: static methods need no 'self' (but is harmless otherwise)

    PyObject* co_varnames = PyTuple_New(co_argcount+1 /* self */);
    PyTuple_SET_ITEM(co_varnames, 0, CPyCppyy_PyText_FromString("self"));
    for (int iarg = 0; iarg < co_argcount; ++iarg) {
        std::string argrep = Cppyy::GetMethodArgType(fMethod, iarg);
        const std::string& parname = Cppyy::GetMethodArgName(fMethod, iarg);
        if (!parname.empty()) {
            argrep += " ";
            argrep += parname;
        }

        PyObject* pyspec = CPyCppyy_PyText_FromString(argrep.c_str());
        PyTuple_SET_ITEM(co_varnames, iarg+1, pyspec);
    }

    return co_varnames;
}

PyObject* CPyCppyy::CPPMethod::GetArgDefault(int iarg)
{
// get the default value (if any) of argument iarg of this method
    if (iarg >= (int)GetMaxArgs())
        return nullptr;

    const std::string& defvalue = Cppyy::GetMethodArgDefault(fMethod, iarg);
    if (!defvalue.empty()) {

    // attempt to evaluate the string representation (will work for all builtin types)
        PyObject* pyval = (PyObject*)PyRun_String(
            (char*)defvalue.c_str(), Py_eval_input, gThisModule, gThisModule);
        if (!pyval && PyErr_Occurred()) {
            PyErr_Clear();
            return CPyCppyy_PyText_FromString(defvalue.c_str());
        }

        return pyval;
    }

    return nullptr;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::GetScopeProxy()
{
// Get or build the scope of this method.
    return CreateScopeProxy(fScope);
}


//----------------------------------------------------------------------------
Cppyy::TCppFuncAddr_t CPyCppyy::CPPMethod::GetFunctionAddress()
{
// Return the C++ pointer of this function
    return Cppyy::GetFunctionAddress(fMethod, false /* don't check fast path envar */);
}


//----------------------------------------------------------------------------
bool CPyCppyy::CPPMethod::Initialize(CallContext* ctxt)
{
// done if cache is already setup
    if (fArgsRequired != -1)
        return true;

    if (!InitConverters_())
        return false;

    if (!InitExecutor_(fExecutor, ctxt))
        return false;

// minimum number of arguments when calling
    fArgsRequired = (int)((bool)fMethod == true ? Cppyy::GetMethodReqArgs(fMethod) : 0);

    return true;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::ProcessKeywords(PyObject* self, PyObject* args, PyObject* kwds)
{
    if (!PyDict_CheckExact(kwds)) {
        SetPyError_(CPyCppyy_PyText_FromString("received unknown keyword arguments object"));
        return nullptr;
    }

    if (PyDict_Size(kwds) == 0 && !self) {
        Py_INCREF(args);
        return args;
    }

    if (!fArgIndices) {
        fArgIndices = new std::map<std::string, int>{};
        for (int iarg = 0; iarg < (int)Cppyy::GetMethodNumArgs(fMethod); ++iarg)
            (*fArgIndices)[Cppyy::GetMethodArgName(fMethod, iarg)] = iarg;
    }

    Py_ssize_t nKeys = PyDict_Size(kwds);
    Py_ssize_t nArgs = PyTuple_GET_SIZE(args) + (self ? 1 : 0);
    if (nKeys+nArgs < fArgsRequired) {
        SetPyError_(CPyCppyy_PyText_FromFormat(
            "takes at least %d arguments (%zd given)", fArgsRequired, nKeys+nArgs));
        return nullptr;
    }

    PyObject* newArgs = PyTuple_New(nArgs+nKeys);

// set all values to zero to be able to check them later (this also guarantees normal
// cleanup by the tuple deallocation)
    for (Py_ssize_t i = 0; i < nArgs+nKeys; ++i)
        PyTuple_SET_ITEM(newArgs, i, nullptr);

// next, insert the keyword values
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(kwds, &pos, &key, &value)) {
        const char* ckey = CPyCppyy_PyText_AsStringChecked(key);
        if (!ckey) {
            Py_DECREF(newArgs);
            return nullptr;
        }
        auto p = fArgIndices->find(ckey);
        if (p == fArgIndices->end()) {
            SetPyError_(CPyCppyy_PyText_FromFormat("%s::%s got an unexpected keyword argument \'%s\'",
                Cppyy::GetFinalName(fScope).c_str(), Cppyy::GetMethodName(fMethod).c_str(), ckey));
            Py_DECREF(newArgs);
            return nullptr;
        }
        Py_INCREF(value);
        PyTuple_SetItem(newArgs, (*fArgIndices)[ckey], value);
    }

// fill out the rest of the arguments
    Py_ssize_t start = 0;
    if (self) {
        Py_INCREF(self);
        PyTuple_SET_ITEM(newArgs, 0, self);
        start = 1;
    }

    for (Py_ssize_t i = start; i < nArgs; ++i) {
        if (PyTuple_GET_ITEM(newArgs, i)) {
            SetPyError_(CPyCppyy_PyText_FromFormat("%s::%s got multiple values for argument %d",
                Cppyy::GetFinalName(fScope).c_str(), Cppyy::GetMethodName(fMethod).c_str(), (int)i+1));
            Py_DECREF(newArgs);
            return nullptr;
        }

        PyObject* item = PyTuple_GET_ITEM(args, i);
        Py_INCREF(item);
        PyTuple_SET_ITEM(newArgs, i, item);
    }

    return newArgs;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::PreProcessArgs(
    CPPInstance*& self, PyObject* args, PyObject* kwds)
{
// verify existence of self, return if ok
    if (self) {
        if (kwds) return ProcessKeywords(nullptr, args, kwds);
        Py_INCREF(args);
        return args;
    }

// otherwise, check for a suitable 'self' in args and update accordingly
    if (PyTuple_GET_SIZE(args) != 0) {
        CPPInstance* pyobj = (CPPInstance*)PyTuple_GET_ITEM(args, 0);

    // demand CPyCppyy object, and an argument that may match down the road
        if (CPPInstance_Check(pyobj) &&
             (fScope == Cppyy::gGlobalScope ||                  // free global
             (pyobj->ObjectIsA() == 0)     ||                   // null pointer or ctor call
             (Cppyy::IsSubtype(pyobj->ObjectIsA(), fScope)))) { // matching types

        // reset self
            Py_INCREF(pyobj);      // corresponding Py_DECREF is in CPPOverload
            self = pyobj;

        // offset args by 1 (new ref)
            PyObject* newArgs = PyTuple_GetSlice(args, 1, PyTuple_GET_SIZE(args));

        // put the keywords, if any, in their places in the arguments array
            if (kwds) {
                args = ProcessKeywords(nullptr, newArgs, kwds);
                Py_DECREF(newArgs);
                newArgs = args;
            }

            return newArgs;  // may be nullptr if kwds insertion failed
        }
    }

// no self, set error and lament
    SetPyError_(CPyCppyy_PyText_FromFormat(
        "unbound method %s::%s must be called with a %s instance as first argument",
        Cppyy::GetFinalName(fScope).c_str(), Cppyy::GetMethodName(fMethod).c_str(),
        Cppyy::GetFinalName(fScope).c_str()));
    return nullptr;
}

//----------------------------------------------------------------------------
bool CPyCppyy::CPPMethod::ConvertAndSetArgs(PyObject* args, CallContext* ctxt)
{
    Py_ssize_t argc = PyTuple_GET_SIZE(args);
    Py_ssize_t argMax = (Py_ssize_t)fConverters.size();

    if (argMax != argc) {
    // argc must be between min and max number of arguments
        if (argc < (Py_ssize_t)fArgsRequired) {
            SetPyError_(CPyCppyy_PyText_FromFormat(
                "takes at least %d arguments (%zd given)", fArgsRequired, argc));
            return false;
        } else if (argMax < argc) {
            SetPyError_(CPyCppyy_PyText_FromFormat(
                "takes at most %zd arguments (%zd given)", argMax, argc));
            return false;
        }
    }

    if (argc == 0)
        return true;

// pass current scope for which the call is made
    ctxt->fCurScope = fScope;

// convert the arguments to the method call array
    bool isOK = true;
    Parameter* cppArgs = ctxt->GetArgs(argc);
    for (int i = 0; i < (int)argc; ++i) {
        if (!fConverters[i]->SetArg(PyTuple_GET_ITEM(args, i), cppArgs[i], ctxt)) {
            SetPyError_(CPyCppyy_PyText_FromFormat("could not convert argument %d", i+1));
            isOK = false;
            break;
        }
    }

    return isOK;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::Execute(void* self, ptrdiff_t offset, CallContext* ctxt)
{
// call the interface method
    PyObject* result = 0;

    if (CallContext::sSignalPolicy != CallContext::kProtected && \
        !(ctxt->fFlags & CallContext::kProtected)) {
    // bypasses try block (i.e. segfaults will abort)
        result = ExecuteFast(self, offset, ctxt);
    } else {
    // at the cost of ~10% performance, don't abort the interpreter on any signal
        result = ExecuteProtected(self, offset, ctxt);
    }

// TODO: the following is dreadfully slow and dead-locks on Apache: revisit
// raising exceptions through callbacks by using magic returns
//    if (result && Utility::PyErr_Occurred_WithGIL()) {
//    // can happen in the case of a CINT error: trigger exception processing
//        Py_DECREF(result);
//        result = 0;
//    } else if (!result && PyErr_Occurred())
    if (!result && PyErr_Occurred())
        SetPyError_(0);

    return result;
}

//----------------------------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::Call(
    CPPInstance*& self, PyObject* args, PyObject* kwds, CallContext* ctxt)
{
// setup as necessary
    if (fArgsRequired == -1 && !Initialize(ctxt))
        return nullptr;

// fetch self, verify, and put the arguments in usable order
    if (!(args = PreProcessArgs(self, args, kwds)))
        return nullptr;

// translate the arguments
    if (fArgsRequired || PyTuple_GET_SIZE(args)) {
        if (!ConvertAndSetArgs(args, ctxt)) {
            Py_DECREF(args);
            return nullptr;
        }
    }

// get the C++ object that this object proxy is a handle for
    void* object = self->GetObject();

// validity check that should not fail
    if (!object) {
        PyErr_SetString(PyExc_ReferenceError, "attempt to access a null-pointer");
        Py_DECREF(args);
        return nullptr;
    }

// get its class
    Cppyy::TCppType_t derived = self->ObjectIsA();

// calculate offset (the method expects 'this' to be an object of fScope)
    ptrdiff_t offset = 0;
    if (derived && derived != fScope)
        offset = Cppyy::GetBaseOffset(derived, fScope, object, 1 /* up-cast */);

// actual call; recycle self instead of returning new object for same address objects
    CPPInstance* pyobj = (CPPInstance*)Execute(object, offset, ctxt);
    Py_DECREF(args);

    if (CPPInstance_Check(pyobj) &&
            derived && pyobj->ObjectIsA() == derived &&
            pyobj->GetObject() == object) {
        Py_INCREF((PyObject*)self);
        Py_DECREF(pyobj);
        return (PyObject*)self;
    }

    return (PyObject*)pyobj;
}

//- protected members --------------------------------------------------------
PyObject* CPyCppyy::CPPMethod::GetSignature(bool fa)
{
// construct python string from the method's signature
    return CPyCppyy_PyText_FromString(GetSignatureString(fa).c_str());
}

//----------------------------------------------------------------------------
std::string CPyCppyy::CPPMethod::GetReturnTypeName()
{
    return Cppyy::GetMethodResultType(fMethod);
}
