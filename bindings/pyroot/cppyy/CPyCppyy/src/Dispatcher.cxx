// Bindings
#include "CPyCppyy.h"
#include "Dispatcher.h"
#include "CPPScope.h"
#include "PyStrings.h"
#include "ProxyWrappers.h"
#include "TypeManip.h"
#include "Utility.h"

// Standard
#include <set>
#include <sstream>


//----------------------------------------------------------------------------
static inline void InjectMethod(Cppyy::TCppMethod_t method, const std::string& mtCppName, std::ostringstream& code)
{
// inject implementation for an overridden method
    using namespace CPyCppyy;

// method declaration
    std::string retType = Cppyy::GetMethodResultType(method);
    code << "  " << retType << " " << mtCppName << "(";

// build out the signature with predictable formal names
    Cppyy::TCppIndex_t nArgs = Cppyy::GetMethodNumArgs(method);
    std::vector<std::string> argtypes; argtypes.reserve(nArgs);
    for (Cppyy::TCppIndex_t i = 0; i < nArgs; ++i) {
        argtypes.push_back(Cppyy::GetMethodArgType(method, i));
        if (i != 0) code << ", ";
        code << argtypes.back() << " arg" << i;
    }
    code << ") ";
    if (Cppyy::IsConstMethod(method)) code << "const ";
    code << "{\n";

// start function body
    Utility::ConstructCallbackPreamble(retType, argtypes, code);

// perform actual method call
#if PY_VERSION_HEX < 0x03000000
    code << "    PyObject* mtPyName = PyString_FromString(\"" << mtCppName << "\");\n" // TODO: intern?
#else
    code << "    PyObject* mtPyName = PyUnicode_FromString(\"" << mtCppName << "\");\n"
#endif
            "    PyObject* pyresult = PyObject_CallMethodObjArgs((PyObject*)_internal_self, mtPyName";

    for (Cppyy::TCppIndex_t i = 0; i < nArgs; ++i)
        code << ", pyargs[" << i << "]";
    code << ", NULL);\n    Py_DECREF(mtPyName);\n";

// close
    Utility::ConstructCallbackReturn(retType, nArgs, code);
}

//----------------------------------------------------------------------------
namespace {
    struct BaseInfo {
        BaseInfo(Cppyy::TCppType_t t, std::string&& bn, std::string&& bns) :
            btype(t), bname(bn), bname_scoped(bns) {}
        Cppyy::TCppType_t btype;
        std::string bname;
        std::string bname_scoped;
    };

    typedef std::vector<BaseInfo> BaseInfos_t;
} // unnamed namespace

//----------------------------------------------------------------------------
bool CPyCppyy::InsertDispatcher(CPPScope* klass, PyObject* bases, PyObject* dct, std::ostringstream& err)
{
// Scan all methods in dct and where it overloads base methods in klass, create
// dispatchers on the C++ side. Then interject the dispatcher class.
    if (!PyTuple_Check(bases) || !PyTuple_GET_SIZE(bases) || !dct || !PyDict_Check(dct)) {
        err << "internal error: expected tuple of bases and proper dictionary";
        return false;
    }

    if (!Utility::IncludePython()) {
        err << "failed to include Python.h";
        return false;
    }

    const Py_ssize_t nBases = PyTuple_GET_SIZE(bases);
    BaseInfos_t base_infos; base_infos.reserve(nBases);
    for (Py_ssize_t ibase = 0; ibase < nBases; ++ibase) {
	auto currentBase = PyTuple_GET_ITEM(bases, ibase);
        if (!CPPScope_Check(currentBase))
        // Python base class
            continue;

        Cppyy::TCppType_t basetype = ((CPPScope*)currentBase)->fCppType;
        if (!basetype) {
            err << "base class is incomplete";
            break;
        }

        if (Cppyy::IsNamespace(basetype)) {
            err << Cppyy::GetScopedFinalName(basetype) << " is a namespace";
            break;
        }

        if (!Cppyy::HasVirtualDestructor(basetype)) {
            const std::string& bname = Cppyy::GetScopedFinalName(basetype);
            PyErr_Warn(PyExc_RuntimeWarning, (char*)("class \""+bname+"\" has no virtual destructor").c_str());
        }

        base_infos.emplace_back(
            basetype, TypeManip::template_base(Cppyy::GetFinalName(basetype)), Cppyy::GetScopedFinalName(basetype));
    }

    if ((Py_ssize_t)base_infos.size() == 0)
        return false;
    else if ((Py_ssize_t)base_infos.size() > 1) {
        err << "multi cross-inheritance not supported";
        return false;
    }

    const auto& binfo = base_infos[0];
    const Cppyy::TCppType_t baseType  = binfo.btype;
    const std::string& baseName       = binfo.bname;
    const std::string& baseNameScoped = binfo.bname_scoped;

    bool isDeepHierarchy = klass->fCppType && baseType != klass->fCppType;

// once classes can be extended, should consider re-use; for now, since derived
// python classes can differ in what they override, simply use different shims
    static int counter = 0;
    std::ostringstream osname;
    osname << "Dispatcher" << ++counter;
    const std::string& derivedName = osname.str();

// generate proxy class with the relevant method dispatchers
    std::ostringstream code;

// start class declaration
    code << "namespace __cppyy_internal {\n"
         << "class " << derivedName << " : public ::" << baseNameScoped << " {\n";
    if (!isDeepHierarchy)
        code << "protected:\n  CPyCppyy::DispatchPtr _internal_self;\n";
    code << "public:\n";

// add a virtual destructor for good measure
    code << "  virtual ~" << derivedName << "() {}\n";

// methods: first collect all callables, then get overrides from base class, for
// those that are still missing, search the hierarchy
    PyObject* clbs = PyDict_New();
    PyObject* items = PyDict_Items(dct);
    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(items); ++i) {
        PyObject* value = PyTuple_GET_ITEM(PyList_GET_ITEM(items, i), 1);
        if (PyCallable_Check(value))
            PyDict_SetItem(clbs, PyTuple_GET_ITEM(PyList_GET_ITEM(items, i), 0), value);
    }
    Py_DECREF(items);
    if (PyDict_DelItem(clbs, PyStrings::gInit) != 0)
        PyErr_Clear();

// protected methods and data need their access changed in the C++ trampoline and then
// exposed on the Python side; so, collect their names as we go along
    std::set<std::string> protected_names;

// simple case: methods from current class
    bool has_default = false;
    bool has_cctor = false;
    bool has_constructors = false;
    const Cppyy::TCppIndex_t nMethods = Cppyy::GetNumMethods(baseType);
    for (Cppyy::TCppIndex_t imeth = 0; imeth < nMethods; ++imeth) {
        Cppyy::TCppMethod_t method = Cppyy::GetMethod(baseType, imeth);

        if (Cppyy::IsConstructor(method) && (Cppyy::IsPublicMethod(method) || Cppyy::IsProtectedMethod(method))) {
            has_constructors = true;
            Cppyy::TCppIndex_t nreq = Cppyy::GetMethodReqArgs(method);
            if (nreq == 0)
                has_default = true;
            else if (!has_cctor && nreq == 1) {
                const std::string& argtype = Cppyy::GetMethodArgType(method, 0);
                if (Utility::Compound(argtype) == "&" && TypeManip::clean_type(argtype, false) == baseNameScoped)
                    has_cctor = true;
            }
            continue;
        }

        std::string mtCppName = Cppyy::GetMethodName(method);
        PyObject* key = CPyCppyy_PyText_FromString(mtCppName.c_str());
        int contains = PyDict_Contains(dct, key);
        if (contains == -1) PyErr_Clear();
        if (contains != 1) {
            Py_DECREF(key);

        // if the method is protected, we expose it through re-declaration and forwarding (using
        // does not work here b/c there may be private overloads)
            if (Cppyy::IsProtectedMethod(method)) {
                protected_names.insert(mtCppName);

                code << "  " << Cppyy::GetMethodResultType(method) << " " << mtCppName << "(";
                Cppyy::TCppIndex_t nArgs = Cppyy::GetMethodNumArgs(method);
                for (Cppyy::TCppIndex_t i = 0; i < nArgs; ++i) {
                    if (i != 0) code << ", ";
                    code << Cppyy::GetMethodArgType(method, i) << " arg" << i;
                }
                code << ") ";
                if (Cppyy::IsConstMethod(method)) code << "const ";
                code << "{\n    return " << baseName << "::" << mtCppName << "(";
                for (Cppyy::TCppIndex_t i = 0; i < nArgs; ++i) {
                    if (i != 0) code << ", ";
                    code << "arg" << i;
                }
                code << ");\n  }\n";
            }

            continue;
        }

        InjectMethod(method, mtCppName, code);

        if (PyDict_DelItem(clbs, key) != 0)
            PyErr_Clear();        // happens for overloads
        Py_DECREF(key);
    }

// try to locate left-overs in base classes
    if (PyDict_Size(clbs)) {
        size_t nbases = Cppyy::GetNumBases(baseType);
        for (size_t ibase = 0; ibase < nbases; ++ibase) {
            Cppyy::TCppScope_t tbase = (Cppyy::TCppScope_t)Cppyy::GetScope( \
                Cppyy::GetBaseName(baseType, ibase));

            PyObject* keys = PyDict_Keys(clbs);
            for (Py_ssize_t i = 0; i < PyList_GET_SIZE(keys); ++i) {
            // TODO: should probably invert this looping; but that makes handling overloads clunky
                PyObject* key = PyList_GET_ITEM(keys, i);
                std::string mtCppName = CPyCppyy_PyText_AsString(key);
                const auto& v = Cppyy::GetMethodIndicesFromName(tbase, mtCppName);
                for (auto idx : v)
                    InjectMethod(Cppyy::GetMethod(tbase, idx), mtCppName, code);
                if (!v.empty()) {
                    if (PyDict_DelItem(clbs, key) != 0) PyErr_Clear();
                }
             }
             Py_DECREF(keys);
        }
    }
    Py_DECREF(clbs);

// constructors: most are simply inherited, for use by the Python derived class
    code << "  using " << baseName << "::" << baseName << ";\n";

// for working with C++ templates, additional constructors are needed to make
// sure the python object is properly carried, but they can only be generated
// if the base class supports them
    if (has_default || !has_constructors)
        code << "  " << derivedName << "() {}\n";
    if (has_default || has_cctor || !has_constructors) {
        code << "  " << derivedName << "(const " << derivedName << "& other)";
        if (has_cctor)
            code << " : " << baseName << "(other)";
        if (!isDeepHierarchy) {
            if (has_cctor) code << ", ";
            else code << " : ";
            code << "_internal_self(other._internal_self, this)";
        }
        code << " {}\n";
    }

// destructor: default is fine

// pull in data members that are protected
    Cppyy::TCppIndex_t nData = Cppyy::GetNumDatamembers(baseType);
    if (nData) code << "public:\n";
    for (Cppyy::TCppIndex_t idata = 0; idata < nData; ++idata) {
        if (Cppyy::IsProtectedData(baseType, idata)) {
            const std::string dm_name = Cppyy::GetDatamemberName(baseType, idata);
            if (dm_name != "_internal_self") {
                protected_names.insert(dm_name);
                code << "  using " << baseName << "::" << dm_name << ";\n";
            }
        }
    }

// initialize the dispatch pointer for base
    code << "public:\n  static void _init_dispatchptr(" << derivedName << "* inst, PyObject* self) {\n"
            "    new ((void*)&inst->_internal_self) CPyCppyy::DispatchPtr{self};\n"
            "  }\n";

// finish class declaration
    code << "};\n}";

// finally, compile the code
    if (!Cppyy::Compile(code.str())) {
        err << "failed to compile the dispatcher code";
        return false;
    }

// keep track internally of the actual C++ type (this is used in
// CPPConstructor to call the dispatcher's one instead of the base)
    Cppyy::TCppScope_t disp = Cppyy::GetScope("__cppyy_internal::"+derivedName);
    if (!disp) {
        err << "failed to retrieve the internal dispatcher";
        return false;
    }
    klass->fCppType = disp;

// at this point, the dispatcher only lives in C++, as opposed to regular classes
// that are part of the hierarchy in Python, so create it, which will cache it for
// later use by e.g. the MemoryRegulator
    PyObject* disp_proxy = CPyCppyy::CreateScopeProxy(disp);

// finally, to expose protected members, copy them over from the C++ dispatcher base
// to the Python dictionary (the C++ dispatcher's Python proxy is not a base of the
// Python class to keep the inheritance tree intact)
    for (const auto& name : protected_names) {
         PyObject* disp_dct = PyObject_GetAttr(disp_proxy, PyStrings::gDict);
         PyObject* pyf = PyMapping_GetItemString(disp_dct, (char*)name.c_str());
         if (pyf) {
             PyObject_SetAttrString((PyObject*)klass, (char*)name.c_str(), pyf);
             Py_DECREF(pyf);
         }
         Py_DECREF(disp_dct);
    }

    Py_XDECREF(disp_proxy);

    return true;
}
