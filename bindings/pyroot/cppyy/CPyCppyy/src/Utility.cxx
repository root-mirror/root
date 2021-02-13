// Bindings
#include "CPyCppyy.h"
#include "Utility.h"
#include "CPPFunction.h"
#include "CPPInstance.h"
#include "CPPOverload.h"
#include "ProxyWrappers.h"
#include "PyCallable.h"
#include "PyStrings.h"
#include "CustomPyTypes.h"
#include "TemplateProxy.h"
#include "TypeManip.h"

// Standard
#include <limits.h>
#include <string.h>
#include <algorithm>
#include <list>
#include <mutex>
#include <set>
#include <sstream>
#include <utility>


//- data _____________________________________________________________________
dict_lookup_func CPyCppyy::gDictLookupOrg = 0;
bool CPyCppyy::gDictLookupActive = false;

typedef std::map<std::string, std::string> TC2POperatorMapping_t;
static TC2POperatorMapping_t gC2POperatorMapping;
static std::set<std::string> gOpSkip;
static std::set<std::string> gOpRemove;

namespace {

    using namespace CPyCppyy::Utility;

    struct InitOperatorMapping_t {
    public:
        InitOperatorMapping_t() {
        // Initialize the global map of operator names C++ -> python.

            gOpSkip.insert("[]");      // __s/getitem__, depends on return type
            gOpSkip.insert("+");       // __add__, depends on # of args (see __pos__)
            gOpSkip.insert("-");       // __sub__, id. (eq. __neg__)
            gOpSkip.insert("*");       // __mul__ or __deref__
            gOpSkip.insert("++");      // __postinc__ or __preinc__
            gOpSkip.insert("--");      // __postdec__ or __predec__

            gOpRemove.insert("new");   // this and the following not handled at all
            gOpRemove.insert("new[]");
            gOpRemove.insert("delete");
            gOpRemove.insert("delete[]");

            gC2POperatorMapping["[]"]  = "__getitem__";
            gC2POperatorMapping["()"]  = "__call__";
            gC2POperatorMapping["/"]   = CPPYY__div__;
            gC2POperatorMapping["%"]   = "__mod__";
            gC2POperatorMapping["**"]  = "__pow__";
            gC2POperatorMapping["<<"]  = "__lshift__";
            gC2POperatorMapping[">>"]  = "__rshift__";
            gC2POperatorMapping["&"]   = "__and__";
            gC2POperatorMapping["&&"]  = "__dand__";
            gC2POperatorMapping["|"]   = "__or__";
            gC2POperatorMapping["||"]  = "__dor__";
            gC2POperatorMapping["^"]   = "__xor__";
            gC2POperatorMapping["~"]   = "__invert__";
            gC2POperatorMapping[","]   = "__comma__";
            gC2POperatorMapping["+="]  = "__iadd__";
            gC2POperatorMapping["-="]  = "__isub__";
            gC2POperatorMapping["*="]  = "__imul__";
            gC2POperatorMapping["/="]  = CPPYY__idiv__;
            gC2POperatorMapping["%="]  = "__imod__";
            gC2POperatorMapping["**="] = "__ipow__";
            gC2POperatorMapping["<<="] = "__ilshift__";
            gC2POperatorMapping[">>="] = "__irshift__";
            gC2POperatorMapping["&="]  = "__iand__";
            gC2POperatorMapping["|="]  = "__ior__";
            gC2POperatorMapping["^="]  = "__ixor__";
            gC2POperatorMapping["=="]  = "__eq__";
            gC2POperatorMapping["!="]  = "__ne__";
            gC2POperatorMapping[">"]   = "__gt__";
            gC2POperatorMapping["<"]   = "__lt__";
            gC2POperatorMapping[">="]  = "__ge__";
            gC2POperatorMapping["<="]  = "__le__";

        // the following type mappings are "exact"
            gC2POperatorMapping["const char*"]  = "__str__";
            gC2POperatorMapping["char*"]        = "__str__";
            gC2POperatorMapping["const char *"] = gC2POperatorMapping["const char*"];
            gC2POperatorMapping["char *"]       = gC2POperatorMapping["char*"];
            gC2POperatorMapping["int"]          = "__int__";
            gC2POperatorMapping["long"]         = CPPYY__long__;
            gC2POperatorMapping["double"]       = "__float__";

        // the following type mappings are "okay"; the assumption is that they
        // are not mixed up with the ones above or between themselves (and if
        // they are, that it is done consistently)
            gC2POperatorMapping["short"]              = "__int__";
            gC2POperatorMapping["unsigned short"]     = "__int__";
            gC2POperatorMapping["unsigned int"]       = CPPYY__long__;
            gC2POperatorMapping["unsigned long"]      = CPPYY__long__;
            gC2POperatorMapping["long long"]          = CPPYY__long__;
            gC2POperatorMapping["unsigned long long"] = CPPYY__long__;
            gC2POperatorMapping["float"]              = "__float__";

            gC2POperatorMapping["->"]  = "__follow__";      // not an actual python operator
            gC2POperatorMapping["="]   = "__assign__";      // id.

#if PY_VERSION_HEX < 0x03000000
            gC2POperatorMapping["bool"] = "__nonzero__";
#else
            gC2POperatorMapping["bool"] = "__bool__";
#endif
        }
    } initOperatorMapping_;

// TODO: this should live with Helpers
    inline void RemoveConst(std::string& cleanName) {
        std::string::size_type spos = std::string::npos;
        while ((spos = cleanName.find("const")) != std::string::npos) {
            cleanName.swap(cleanName.erase(spos, 5));
        }
    }

} // unnamed namespace


//- public functions ---------------------------------------------------------
unsigned long CPyCppyy::PyLongOrInt_AsULong(PyObject* pyobject)
{
// Convert <pybject> to C++ unsigned long, with bounds checking, allow int -> ulong.
    unsigned long ul = PyLong_AsUnsignedLong(pyobject);
    if (PyErr_Occurred() && PyInt_Check(pyobject)) {
        PyErr_Clear();
        long i = PyInt_AS_LONG(pyobject);
        if (0 <= i) {
            ul = (unsigned long)i;
        } else {
            PyErr_SetString(PyExc_ValueError,
                "can\'t convert negative value to unsigned long");
            return (unsigned long)-1;
        }
    }

    return ul;
}

//----------------------------------------------------------------------------
ULong64_t CPyCppyy::PyLongOrInt_AsULong64(PyObject* pyobject)
{
// Convert <pyobject> to C++ unsigned long long, with bounds checking.
    ULong64_t ull = PyLong_AsUnsignedLongLong(pyobject);
    if (PyErr_Occurred() && PyInt_Check(pyobject)) {
        PyErr_Clear();
        long i = PyInt_AS_LONG(pyobject);
        if (0 <= i) {
            ull = (ULong64_t)i;
        } else {
            PyErr_SetString(PyExc_ValueError,
                "can\'t convert negative value to unsigned long long");
        }
    }

    return ull;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Utility::AddToClass(
    PyObject* pyclass, const char* label, PyCFunction cfunc, int flags)
{
// Add the given function to the class under name 'label'.

// use list for clean-up (.so's are unloaded only at interpreter shutdown)
    static std::list<PyMethodDef> s_pymeths;

    s_pymeths.push_back(PyMethodDef());
    PyMethodDef* pdef = &s_pymeths.back();
    pdef->ml_name  = const_cast<char*>(label);
    pdef->ml_meth  = cfunc;
    pdef->ml_flags = flags;
    pdef->ml_doc   = nullptr;

    PyObject* func = PyCFunction_New(pdef, nullptr);
    PyObject* name = CPyCppyy_PyText_InternFromString(pdef->ml_name);
    PyObject* method = CustomInstanceMethod_New(func, nullptr, pyclass);
    bool isOk = PyType_Type.tp_setattro(pyclass, name, method) == 0;
    Py_DECREF(method);
    Py_DECREF(name);
    Py_DECREF(func);

    if (PyErr_Occurred())
        return false;

    if (!isOk) {
        PyErr_Format(PyExc_TypeError, "could not add method %s", label);
        return false;
    }

    return true;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Utility::AddToClass(PyObject* pyclass, const char* label, const char* func)
{
// Add the given function to the class under name 'label'.
    PyObject* pyfunc = PyObject_GetAttrString(pyclass, const_cast<char*>(func));
    if (!pyfunc)
        return false;

    PyObject* pylabel = CPyCppyy_PyText_InternFromString(const_cast<char*>(label));
    bool isOk = PyType_Type.tp_setattro(pyclass, pylabel, pyfunc) == 0;
    Py_DECREF(pylabel);

    Py_DECREF(pyfunc);
    return isOk;
}

//----------------------------------------------------------------------------
bool CPyCppyy::Utility::AddToClass(PyObject* pyclass, const char* label, PyCallable* pyfunc)
{
// Add the given function to the class under name 'label'.
    CPPOverload* method =
        (CPPOverload*)PyObject_GetAttrString(pyclass, const_cast<char*>(label));

    if (!method || !CPPOverload_Check(method)) {
    // not adding to existing CPPOverload; add callable directly to the class
        if (PyErr_Occurred())
            PyErr_Clear();
        Py_XDECREF((PyObject*)method);
        method = CPPOverload_New(label, pyfunc);
        PyObject* pylabel = CPyCppyy_PyText_InternFromString(const_cast<char*>(label));
        bool isOk = PyType_Type.tp_setattro(pyclass, pylabel, (PyObject*)method) == 0;
        Py_DECREF(pylabel);
        Py_DECREF(method);
        return isOk;
    }

    method->AdoptMethod(pyfunc);

    Py_DECREF(method);
    return true;
}


//----------------------------------------------------------------------------
static inline
CPyCppyy::PyCallable* BuildOperator(const std::string& lcname, const std::string& rcname,
    const char* op, Cppyy::TCppScope_t scope, bool reverse=false)
{
// Helper to find a function with matching signature in 'funcs'.
    std::string opname = "operator";
    opname += op;

    Cppyy::TCppIndex_t idx = Cppyy::GetGlobalOperator(scope, lcname, rcname, opname);
    if (idx == (Cppyy::TCppIndex_t)-1)
        return nullptr;

    Cppyy::TCppMethod_t meth = Cppyy::GetMethod(scope, idx);
    if (!reverse)
        return new CPyCppyy::CPPFunction(scope, meth);
    return new CPyCppyy::CPPReverseBinary(scope, meth);
}

//----------------------------------------------------------------------------
CPyCppyy::PyCallable* CPyCppyy::Utility::FindUnaryOperator(PyObject* pyclass, const char* op)
{
// Find a callable matching named operator (op) and klass arguments in the global
// namespace or the klass' namespace.
    if (!CPPScope_Check(pyclass))
        return nullptr;

    CPPClass* klass = (CPPClass*)pyclass;
    const std::string& lcname = Cppyy::GetScopedFinalName(klass->fCppType);
    Cppyy::TCppScope_t scope = Cppyy::GetScope(TypeManip::extract_namespace(lcname));
    return FindBinaryOperator(lcname, "", op, scope, false);
}

//----------------------------------------------------------------------------
CPyCppyy::PyCallable* CPyCppyy::Utility::FindBinaryOperator(PyObject* left, PyObject* right,
    const char* op, Cppyy::TCppScope_t scope)
{
// Find a callable matching the named operator (op) and the (left, right)
// arguments in the global or these objects' namespaces.

    bool reverse = false;
    if (!CPPInstance_Check(left)) {
        if (CPPInstance_Check(right))
           reverse = true;
        else
           return nullptr;
    }

// retrieve the class names to match the signature of any found global functions
    const std::string& lcname = ClassName(left);
    const std::string& rcname = ClassName(right);
    return FindBinaryOperator(lcname, rcname, op, scope, reverse);
}

//----------------------------------------------------------------------------
CPyCppyy::PyCallable* CPyCppyy::Utility::FindBinaryOperator(
    const std::string& lcname, const std::string& rcname,
    const char* op, Cppyy::TCppScope_t scope, bool reverse)
{
// Find a global function with a matching signature; search __gnu_cxx, std::__1,
// and __cppyy_internal pro-actively (as there's AFAICS no way to unearth 'using'
// information).

    if (rcname == "<unknown>" || lcname == "<unknown>")
        return nullptr;

    PyCallable* pyfunc = 0;

    const std::string& lnsname = TypeManip::extract_namespace(lcname);
    if (!scope) scope = Cppyy::GetScope(lnsname);
    if (scope)
        pyfunc = BuildOperator(lcname, rcname, op, scope, reverse);

    if (!pyfunc && scope != Cppyy::gGlobalScope)      // search in global scope anyway
        pyfunc = BuildOperator(lcname, rcname, op, Cppyy::gGlobalScope, reverse);

    if (!pyfunc) {
    // For GNU on clang, search the internal __gnu_cxx namespace for binary operators (is
    // typically the case for STL iterators operator==/!=.
    // TODO: only look in __gnu_cxx for iterators (and more generally: do lookups in the
    //       namespace where the class is defined
        static Cppyy::TCppScope_t gnucxx = Cppyy::GetScope("__gnu_cxx");
        if (gnucxx)
            pyfunc = BuildOperator(lcname, rcname, op, gnucxx, reverse);
    }

    if (!pyfunc) {
    // Same for clang (on Mac only?). TODO: find proper pre-processor magic to only use those
    // specific namespaces that are actually around; although to be sure, this isn't expensive.
        static Cppyy::TCppScope_t std__1 = Cppyy::GetScope("std::__1");

        if (std__1
#ifdef __APPLE__
 && lcname.find("__wrap_iter") == std::string::npos   // wrapper call does not compile
#endif
        ) {
            pyfunc = BuildOperator(lcname, rcname, op, std__1, reverse);
        }
    }

    if (!pyfunc) {
    // One more, mostly for Mac, but again not sure whether this is not a general issue. Some
    // operators are declared as friends only in classes, so then they're not found in the
    // global namespace, so this helper let's the compiler resolve the operator.
        static Cppyy::TCppScope_t s_intern = Cppyy::GetScope("__cppyy_internal");
        if (s_intern) {
            std::stringstream fname, proto;
            if (strncmp(op, "==", 2) == 0) { fname << "is_equal<"; }
            else if (strncmp(op, "!=", 2) == 0) { fname << "is_not_equal<"; }
            else { fname << "not_implemented<"; }
            fname << lcname << ", " << rcname << ">";
            proto << "const " << lcname << "&, const " << rcname;
            Cppyy::TCppMethod_t method = Cppyy::GetMethodTemplate(s_intern, fname.str(), proto.str());
            if (method) pyfunc = new CPPFunction(s_intern, method);
        }
    }

    return pyfunc;
}

//----------------------------------------------------------------------------
static bool AddTypeName(std::string& tmpl_name, PyObject* tn, PyObject* arg,
    CPyCppyy::Utility::ArgPreference pref, int* pcnt = nullptr)
{
// Determine the appropriate C++ type for a given Python type; this is a helper because
// it can recurse if the type is list or tuple and needs matching on std::vector.
    using namespace CPyCppyy;
    using namespace CPyCppyy::Utility;

    if (tn == (PyObject*)&PyInt_Type) {
        if (arg) {
#if PY_VERSION_HEX < 0x03000000
            long l = PyInt_AS_LONG(arg);
            tmpl_name.append((l < INT_MIN || INT_MAX < l) ? "long" : "int");
#else
             Long64_t ll = PyLong_AsLongLong(arg);
             if (ll == (Long64_t)-1 && PyErr_Occurred()) {
                 PyErr_Clear();
                 ULong64_t ull = PyLong_AsUnsignedLongLong(arg);
                 if (ull == (ULong64_t)-1 && PyErr_Occurred()) {
                     PyErr_Clear();
                     tmpl_name.append("int");    // still out of range, will fail later
                 } else
                     tmpl_name.append("ULong64_t");   // since already failed long long
             } else
                 tmpl_name.append((ll < INT_MIN || INT_MAX < ll) ? \
                     ((ll < LONG_MIN || LONG_MAX < ll) ? "Long64_t" : "long") : "int");
#endif
        } else
            tmpl_name.append("int");
#if PY_VERSION_HEX < 0x03000000
    } else if (tn == (PyObject*)&PyLong_Type) {
        if (arg) {
             Long64_t ll = PyLong_AsLongLong(arg);
             if (ll == (Long64_t)-1 && PyErr_Occurred()) {
                 PyErr_Clear();
                 ULong64_t ull = PyLong_AsUnsignedLongLong(arg);
                 if (ull == (ULong64_t)-1 && PyErr_Occurred()) {
                     PyErr_Clear();
                     tmpl_name.append("long");   // still out of range, will fail later
                 } else
                     tmpl_name.append("ULong64_t");   // since already failed long long
             } else
                 tmpl_name.append((ll < LONG_MIN || LONG_MAX < ll) ? "Long64_t" : "long");
        } else
            tmpl_name.append("long");
#endif
    } else if (tn == (PyObject*)&PyFloat_Type) {
    // special case for floats (Python-speak for double) if from argument (only)
        tmpl_name.append(arg ? "double" : "float");
#if PY_VERSION_HEX < 0x03000000
    } else if (tn == (PyObject*)&PyString_Type) {
#else
    } else if (tn == (PyObject*)&PyUnicode_Type) {
#endif
        tmpl_name.append("std::string");
    } else if (tn == (PyObject*)&PyList_Type || tn == (PyObject*)&PyTuple_Type) {
        if (arg && PySequence_Size(arg)) {
            std::string subtype{"std::initializer_list<"};
            PyObject* item = PySequence_GetItem(arg, 0);
            ArgPreference subpref = pref == kValue ? kValue : kPointer;
            if (AddTypeName(subtype, (PyObject*)Py_TYPE(item), item, subpref)) {
                tmpl_name.append(subtype);
                tmpl_name.append(">");
            }
            Py_DECREF(item);
        }

    } else if (CPPScope_Check(tn)) {
        tmpl_name.append(Cppyy::GetScopedFinalName(((CPPClass*)tn)->fCppType));
        if (arg) {
        // try to specialize the type match for the given object
            CPPInstance* pyobj = (CPPInstance*)arg;
            if (CPPInstance_Check(pyobj)) {
                if (pyobj->fFlags & CPPInstance::kIsRValue)
                    tmpl_name.append("&&");
                else {
                    if (pcnt) *pcnt += 1;
                    if ((pyobj->fFlags & CPPInstance::kIsReference) || pref == kPointer)
                        tmpl_name.push_back('*');
                    else if (pref != kValue)
                        tmpl_name.push_back('&');
                }
            }
        }
    } else if (PyObject_HasAttr(tn, PyStrings::gCppName)) {
        PyObject* tpName = PyObject_GetAttr(tn, PyStrings::gCppName);
        tmpl_name.append(CPyCppyy_PyText_AsString(tpName));
        Py_DECREF(tpName);
    } else if (PyObject_HasAttr(tn, PyStrings::gName)) {
        PyObject* tpName = PyObject_GetAttr(tn, PyStrings::gName);
        tmpl_name.append(CPyCppyy_PyText_AsString(tpName));
        Py_DECREF(tpName);
    } else if (PyInt_Check(tn) || PyLong_Check(tn) || PyFloat_Check(tn)) {
    // last ditch attempt, works for things like int values; since this is a
    // source of errors otherwise, it is limited to specific types and not
    // generally used (str(obj) can print anything ...)
        PyObject* pystr = PyObject_Str(tn);
        tmpl_name.append(CPyCppyy_PyText_AsString(pystr));
        Py_DECREF(pystr);
    } else {
        return false;
    }

    return true;
}

std::string CPyCppyy::Utility::ConstructTemplateArgs(
    PyObject* pyname, PyObject* tpArgs, PyObject* args, ArgPreference pref, int argoff, int* pcnt)
{
// Helper to construct the "<type, type, ...>" part of a templated name (either
// for a class or method lookup
    bool justOne = !PyTuple_CheckExact(tpArgs);

// Note: directly appending to string is a lot faster than stringstream
    std::string tmpl_name;
    tmpl_name.reserve(128);
    if (pyname)
        tmpl_name.append(CPyCppyy_PyText_AsString(pyname));
    tmpl_name.push_back('<');

    if (pcnt) *pcnt = 0;     // count number of times 'pref' is used

    Py_ssize_t nArgs = justOne ? 1 : PyTuple_GET_SIZE(tpArgs);
    for (int i = argoff; i < nArgs; ++i) {
    // add type as string to name
        PyObject* tn = justOne ? tpArgs : PyTuple_GET_ITEM(tpArgs, i);
        if (CPyCppyy_PyText_Check(tn)) {
            tmpl_name.append(CPyCppyy_PyText_AsString(tn));
    // some commmon numeric types (separated out for performance: checking for
    // __cpp_name__ and/or __name__ is rather expensive)
        } else {
            if (!AddTypeName(tmpl_name, tn, (args ? PyTuple_GET_ITEM(args, i) : nullptr), pref, pcnt)) {
                PyErr_SetString(PyExc_SyntaxError,
                    "could not construct C++ name from provided template argument.");
                return "";
            }
        }

    // add a comma, as needed (no space as internally, final names don't have them)
        if (i != nArgs-1)
            tmpl_name.push_back(',');
    }

// close template name
    tmpl_name.push_back('>');

    return tmpl_name;
}

//----------------------------------------------------------------------------
void CPyCppyy::Utility::ConstructCallbackPreamble(const std::string& retType,
    const std::vector<std::string>& argtypes, std::ostringstream& code)
{
// Generate function setup to be used in callbacks (wrappers and overrides).
    int nArgs = (int)argtypes.size();

// return value and argument type converters
    bool isVoid = retType == "void";
    if (!isVoid)
        code << "    CPYCPPYY_STATIC std::unique_ptr<CPyCppyy::Converter, std::function<void(CPyCppyy::Converter*)>> "
                     "retconv{CPyCppyy::CreateConverter(\""
             << retType << "\"), CPyCppyy::DestroyConverter};\n";
    if (nArgs) {
        code << "    CPYCPPYY_STATIC std::vector<std::unique_ptr<CPyCppyy::Converter, std::function<void(CPyCppyy::Converter*)>>> argcvs;\n"
             << "    if (argcvs.empty()) {\n"
             << "      argcvs.reserve(" << nArgs << ");\n";
        for (int i = 0; i < nArgs; ++i)
            code << "      argcvs.emplace_back(CPyCppyy::CreateConverter(\"" << argtypes[i] << "\"), CPyCppyy::DestroyConverter);\n";
        code << "    }\n";
    }

// declare return value (TODO: this does not work for most non-builtin values)
    if (!isVoid)
        code << "    " << retType << " ret{};\n";

// acquire GIL
    code << "    PyGILState_STATE state = PyGILState_Ensure();\n";

// build argument tuple if needed
    if (nArgs) {
        code << "    std::vector<PyObject*> pyargs;\n";
        code << "    pyargs.reserve(" << nArgs << ");\n"
             << "    try {\n";
        for (int i = 0; i < nArgs; ++i) {
            code << "      pyargs.emplace_back(argcvs[" << i << "]->FromMemory((void*)&arg" << i << "));\n"
                 << "      if (!pyargs.back()) throw " << i << ";\n";
        }
        code << "    } catch(int) {\n"
             << "      for (auto pyarg : pyargs) Py_XDECREF(pyarg);\n"
             << "      PyGILState_Release(state); throw CPyCppyy::PyException{};\n"
             << "    }\n";
    }
}

void CPyCppyy::Utility::ConstructCallbackReturn(const std::string& retType, int nArgs, std::ostringstream& code)
{
// Generate code for return value conversion and error handling.
    bool isVoid = retType == "void";
    bool isPtr  = Cppyy::ResolveName(retType).back() == '*';
    if (nArgs)
        code << "    for (auto pyarg : pyargs) Py_DECREF(pyarg);\n";
    code << "    bool cOk = (bool)pyresult;\n"
            "    if (pyresult) {\n";
    if (isPtr) {
    // If the return type is a CPPInstance, owned by Python, and the ref-count down
    // to 1, the return will hold a dangling pointer, so set it to nullptr instead.
        code << "      if (!CPyCppyy::Instance_IsLively(pyresult))\n"
                "        ret = nullptr;\n"
                "      else {\n";
    }
    code << (isVoid ? "" : "        cOk = retconv->ToMemory(pyresult, &ret);\n")
         <<                "        Py_DECREF(pyresult);\n    }\n";
    if (isPtr) code << "  }\n";
    code << "    if (!cOk) {"     // assume error set when converter failed
// TODO: On Windows, throwing a C++ exception here makes the code hang; leave
// the error be which allows at least one layer of propagation
#ifdef _WIN32
            " /* do nothing */ }\n"
#else
            " PyGILState_Release(state); throw CPyCppyy::PyException{}; }\n"
#endif
            "    PyGILState_Release(state);\n"
            "    return";
    code << (isVoid ? ";\n  }\n" : " ret;\n  }\n");
}

//----------------------------------------------------------------------------
bool CPyCppyy::Utility::InitProxy(PyObject* module, PyTypeObject* pytype, const char* name)
{
// Initialize a proxy class for use by python, and add it to the module.

// finalize proxy type
    if (PyType_Ready(pytype) < 0)
        return false;

// add proxy type to the given module
    Py_INCREF(pytype);       // PyModule_AddObject steals reference
    if (PyModule_AddObject(module, (char*)name, (PyObject*)pytype) < 0) {
        Py_DECREF(pytype);
        return false;
    }

// declare success
    return true;
}

//----------------------------------------------------------------------------
Py_ssize_t CPyCppyy::Utility::GetBuffer(PyObject* pyobject, char tc, int size, void*& buf, bool check)
{
// Retrieve a linear buffer pointer from the given pyobject.

// special case: don't handle character strings here (yes, they're buffers, but not quite)
    if (PyBytes_Check(pyobject))
        return 0;

// special case: bytes array
    if ((!check || tc == '*' || tc == 'B') && PyByteArray_CheckExact(pyobject)) {
        buf = PyByteArray_AS_STRING(pyobject);
        return PyByteArray_GET_SIZE(pyobject);
    }

// new-style buffer interface
    if (PyObject_CheckBuffer(pyobject)) {
        Py_buffer bufinfo;
        memset(&bufinfo, 0, sizeof(Py_buffer));
        if (PyObject_GetBuffer(pyobject, &bufinfo, PyBUF_FORMAT) == 0) {
            if (tc == '*' || strchr(bufinfo.format, tc)
#ifdef _WIN32
            // ctypes is inconsistent in format on Windows; either way these types are the same size
                || (tc == 'I' && strchr(bufinfo.format, 'L')) || (tc == 'i' && strchr(bufinfo.format, 'l'))
#endif
            // allow 'signed char' ('b') from array to pass through '?' (bool as from struct)
                || (tc == '?' && strchr(bufinfo.format, 'b'))
                    ) {
                buf = bufinfo.buf;
                if (buf && bufinfo.ndim == 0) {
                    PyBuffer_Release(&bufinfo);
                    return bufinfo.len/bufinfo.itemsize;
                } else if (buf && bufinfo.ndim == 1) {
                    Py_ssize_t size1d = bufinfo.shape ? bufinfo.shape[0] : bufinfo.len/bufinfo.itemsize;
                    PyBuffer_Release(&bufinfo);
                    return size1d;
                }
            } else {
            // have buf, but format mismatch: bail out now, otherwise the old
            // code will return based on itemsize match
                PyBuffer_Release(&bufinfo);
                return 0;
            }
        }
        PyErr_Clear();
    }

// attempt to retrieve pointer through old-style buffer interface
    PyBufferProcs* bufprocs = Py_TYPE(pyobject)->tp_as_buffer;

    PySequenceMethods* seqmeths = Py_TYPE(pyobject)->tp_as_sequence;
    if (seqmeths != 0 && bufprocs != 0
#if  PY_VERSION_HEX < 0x03000000
         && bufprocs->bf_getwritebuffer != 0
         && (*(bufprocs->bf_getsegcount))(pyobject, 0) == 1
#else
         && bufprocs->bf_getbuffer != 0
#endif
        ) {

   // get the buffer
#if PY_VERSION_HEX < 0x03000000
        Py_ssize_t buflen = (*(bufprocs->bf_getwritebuffer))(pyobject, 0, &buf);
#else
        Py_buffer bufinfo;
        (*(bufprocs->bf_getbuffer))(pyobject, &bufinfo, PyBUF_WRITABLE);
        buf = (char*)bufinfo.buf;
        Py_ssize_t buflen = bufinfo.len;
#if PY_VERSION_HEX < 0x03010000
        PyBuffer_Release(pyobject, &bufinfo);
#else
        PyBuffer_Release(&bufinfo);
#endif
#endif

        if (buf && check == true) {
        // determine buffer compatibility (use "buf" as a status flag)
            PyObject* pytc = PyObject_GetAttr(pyobject, PyStrings::gTypeCode);
            if (pytc != 0) {      // for array objects
                char cpytc = CPyCppyy_PyText_AsString(pytc)[0];
                if (!(cpytc == tc || (tc == '?' && cpytc == 'b')))
                    buf = 0;      // no match
                Py_DECREF(pytc);
            } else if (seqmeths->sq_length &&
                       (int)(buflen/(*(seqmeths->sq_length))(pyobject)) == size) {
            // this is a gamble ... may or may not be ok, but that's for the user
                PyErr_Clear();
            } else if (buflen == size) {
            // also a gamble, but at least 1 item will fit into the buffer, so very likely ok ...
                PyErr_Clear();
            } else {
                buf = 0;                      // not compatible

            // clarify error message
                PyObject* pytype = 0, *pyvalue = 0, *pytrace = 0;
                PyErr_Fetch(&pytype, &pyvalue, &pytrace);
                PyObject* pyvalue2 = CPyCppyy_PyText_FromFormat(
                    (char*)"%s and given element size (%ld) do not match needed (%d)",
                    CPyCppyy_PyText_AsString(pyvalue),
                    seqmeths->sq_length ? (Long_t)(buflen/(*(seqmeths->sq_length))(pyobject)) : (Long_t)buflen,
                    size);
                Py_DECREF(pyvalue);
                PyErr_Restore(pytype, pyvalue2, pytrace);
            }
        }

        if (!buf) return 0;
        return buflen/(size ? size : 1);
    }

    return 0;
}

//----------------------------------------------------------------------------
std::string CPyCppyy::Utility::MapOperatorName(const std::string& name, bool bTakesParams)
{
// Map the given C++ operator name on the python equivalent.
    if (8 < name.size() && name.substr(0, 8) == "operator") {
        std::string op = name.substr(8, std::string::npos);

    // stripping ...
        std::string::size_type start = 0, end = op.size();
        while (start < end && isspace(op[start])) ++start;
        while (start < end && isspace(op[end-1])) --end;
        op = op.substr(start, end - start);

    // certain operators should be removed completely (e.g. operator delete & friends)
        if (gOpRemove.find(op) != gOpRemove.end())
            return "";

    // check first if none, to prevent spurious deserializing downstream
        TC2POperatorMapping_t::iterator pop = gC2POperatorMapping.find(op);
        if (pop == gC2POperatorMapping.end() && gOpSkip.find(op) == gOpSkip.end()) {
            op = Cppyy::ResolveName(op);
            pop = gC2POperatorMapping.find(op);
        }

    // map C++ operator to python equivalent, or made up name if no equivalent exists
        if (pop != gC2POperatorMapping.end()) {
            return pop->second;

        } else if (op == "*") {
        // dereference v.s. multiplication of two instances
            return bTakesParams ? "__mul__" : "__deref__";

        } else if (op == "+") {
        // unary positive v.s. addition of two instances
            return bTakesParams ? "__add__" : "__pos__";

        } else if (op == "-") {
        // unary negative v.s. subtraction of two instances
            return bTakesParams ? "__sub__" : "__neg__";

        } else if (op == "++") {
        // prefix v.s. postfix increment
            return bTakesParams ? "__postinc__" : "__preinc__";

        } else if (op == "--") {
        // prefix v.s. postfix decrement
            return bTakesParams ? "__postdec__" : "__predec__";
        }

    }

// might get here, as not all operator methods are handled (new, delete, etc.)
    return name;
}

//----------------------------------------------------------------------------
const std::string CPyCppyy::Utility::Compound(const std::string& name)
{
// TODO: consolidate with other string manipulations in TypeManip.cxx
// Break down the compound of a fully qualified type name.
    std::string cleanName = name;
    RemoveConst(cleanName);

    std::string compound = "";
    for (int ipos = (int)cleanName.size()-1; 0 <= ipos; --ipos) {
        char c = cleanName[ipos];
        if (isspace(c)) continue;
        if (isalnum(c) || c == '_' || c == '>' || c == ')') break;

        compound = c + compound;
    }

// for arrays (TODO: deal with the actual size)
    if (compound == "]")
        return "[]";

    return compound;
}

//----------------------------------------------------------------------------
Py_ssize_t CPyCppyy::Utility::ArraySize(const std::string& name)
{
// TODO: consolidate with other string manipulations in Helpers.cxx
// Extract size from an array type, if available.
    std::string cleanName = name;
    RemoveConst(cleanName);

    if (cleanName[cleanName.size()-1] == ']') {
        std::string::size_type idx = cleanName.rfind('[');
        if (idx != std::string::npos) {
            const std::string asize = cleanName.substr(idx+1, cleanName.size()-2);
            return strtoul(asize.c_str(), nullptr, 0);
        }
    }

    return -1;
}

//----------------------------------------------------------------------------
std::string CPyCppyy::Utility::ClassName(PyObject* pyobj)
{
// Retrieve the class name from the given Python instance.
    std::string clname = "<unknown>";
    PyObject* pyclass = (PyObject*)Py_TYPE(pyobj);
    PyObject* pyname = PyObject_GetAttr(pyclass, PyStrings::gCppName);
    if (!pyname) {
        PyErr_Clear();
        pyname = PyObject_GetAttr(pyclass, PyStrings::gName);
    }

    if (pyname) {
        clname = CPyCppyy_PyText_AsString(pyname);
        Py_DECREF(pyname);
    } else
        PyErr_Clear();
    return clname;
}


//----------------------------------------------------------------------------
CPyCppyy::Utility::PyOperators::~PyOperators()
{
    Py_XDECREF(fEq);
    Py_XDECREF(fNe);
    Py_XDECREF(fLAdd); Py_XDECREF(fRAdd);
    Py_XDECREF(fSub);
    Py_XDECREF(fLMul); Py_XDECREF(fRMul);
    Py_XDECREF(fDiv);
    Py_XDECREF(fHash);
}


//----------------------------------------------------------------------------
PyObject* CPyCppyy::Utility::PyErr_Occurred_WithGIL()
{
// Re-acquire the GIL before calling PyErr_Occurred() in case it has been
// released; note that the p2.2 code assumes that there are no callbacks in
// C++ to python (or at least none returning errors).
#if PY_VERSION_HEX >= 0x02030000
    PyGILState_STATE gstate = PyGILState_Ensure();
    PyObject* e = PyErr_Occurred();
    PyGILState_Release(gstate);
#else
    if (PyThreadState_GET())
        return PyErr_Occurred();
    PyObject* e = 0;
#endif

    return e;
}


//----------------------------------------------------------------------------
size_t CPyCppyy::Utility::FetchError(std::vector<PyError_t>& errors)
{
// Fetch the current python error, if any, and store it for future use.
    if (PyErr_Occurred()) {
        PyError_t e;
        PyErr_Fetch(&e.fType, &e.fValue, &e.fTrace);
        errors.push_back(e);
    }
    return errors.size();
}

//----------------------------------------------------------------------------
void CPyCppyy::Utility::SetDetailedException(std::vector<PyError_t>& errors, PyObject* topmsg, PyObject* defexc)
{
// Use the collected exceptions to build up a detailed error log.
    if (errors.empty()) {
    // should not happen ...
        PyErr_SetString(defexc, CPyCppyy_PyText_AsString(topmsg));
        Py_DECREF(topmsg);
        return;
    }

// add the details to the topmsg
    PyObject* separator = CPyCppyy_PyText_FromString("\n  ");

    PyObject* exc_type = nullptr;
    for (auto& e : errors) {
        if (!exc_type) exc_type = e.fType;
        else if (exc_type != e.fType) exc_type = defexc;
        CPyCppyy_PyText_Append(&topmsg, separator);
        if (CPyCppyy_PyText_Check(e.fValue)) {
            CPyCppyy_PyText_Append(&topmsg, e.fValue);
        } else if (e.fValue) {
            PyObject* excstr = PyObject_Str(e.fValue);
            if (!excstr) {
                PyErr_Clear();
                excstr = PyObject_Str((PyObject*)Py_TYPE(e.fValue));
            }
            CPyCppyy_PyText_AppendAndDel(&topmsg, excstr);
        } else {
            CPyCppyy_PyText_AppendAndDel(&topmsg,
                CPyCppyy_PyText_FromString("unknown exception"));
        }
    }

    Py_DECREF(separator);
    std::for_each(errors.begin(), errors.end(), PyError_t::Clear);

// set the python exception
    PyErr_SetString(exc_type, CPyCppyy_PyText_AsString(topmsg));
    Py_DECREF(topmsg);
}


//----------------------------------------------------------------------------
static bool includesDone = false;
bool CPyCppyy::Utility::IncludePython()
{
// setup Python API for callbacks
    if (!includesDone) {
        bool okay = Cppyy::Compile(
        // basic API (converters etc.)
            "#include \"CPyCppyy/API.h\"\n"

        // utilities from the CPyCppyy public API
            "#include \"CPyCppyy/DispatchPtr.h\"\n"
            "#include \"CPyCppyy/PyException.h\"\n"
        );
        includesDone = okay;
    }

    return includesDone;
}
