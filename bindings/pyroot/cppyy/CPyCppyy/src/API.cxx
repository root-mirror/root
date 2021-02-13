// Bindings
#include "CPyCppyy.h"
#define CPYCPPYY_INTERNAL 1
#include "CPyCppyy/API.h"
#undef CPYCPPYY_INTERNAL

#include "CPPInstance.h"
#include "CPPOverload.h"
#include "CPPScope.h"
#include "ProxyWrappers.h"
#include "PyStrings.h"

// Standard
#include <stdio.h>
#include <iostream>
#include <string>

//______________________________________________________________________________
//                  CPyCppyy API: Interpreter and Proxy Access
//                  ==========================================
//
// Access to cppyy Python objects from Cling and C++: allows conversion for
// instances and type checking for scopes, instances, etc.
// Adds a few convenience functions to call Python from Cling and expose Python
// classes to Cling for use in inheritance etc.


//- data ---------------------------------------------------------------------
static PyObject* gMainDict = nullptr;

namespace CPyCppyy {
    extern PyObject* gThisModule;
}


//- private helpers ----------------------------------------------------------
namespace {

static bool Initialize()
{
// Private initialization method: setup the python interpreter and load the
// cppyy module.
    static bool isInitialized = false;
    if (isInitialized)
        return true;

    if (!Py_IsInitialized()) {
    // this happens if Cling comes in first
#if PY_VERSION_HEX < 0x03020000
        PyEval_InitThreads();
#endif
        Py_Initialize();
#if PY_VERSION_HEX >= 0x03020000
        PyEval_InitThreads();
#endif

    // try again to see if the interpreter is initialized
        if (!Py_IsInitialized()) {
        // give up ...
            std::cerr << "Error: python has not been intialized; returning." << std::endl;
            return false;
        }

   // set the command line arguments on python's sys.argv
#if PY_VERSION_HEX < 0x03000000
        char* argv[] = {const_cast<char*>("cppyy")};
#else
        wchar_t* argv[] = {const_cast<wchar_t*>(L"cppyy")};
#endif
        PySys_SetArgv(sizeof(argv)/sizeof(argv[0]), argv);

    // force loading of the cppyy module
        PyRun_SimpleString(const_cast<char*>("import cppyy"));
    }

    if (!gMainDict) {
    // retrieve the main dictionary
        gMainDict = PyModule_GetDict(
            PyImport_AddModule(const_cast<char*>("__main__")));
        Py_INCREF(gMainDict);
    }

// declare success ...
    isInitialized = true;
    return true;
}

} // unnamed namespace


//- C++ access to cppyy objects ---------------------------------------------
void* CPyCppyy::Instance_AsVoidPtr(PyObject* pyobject)
{
// Extract the object pointer held by the CPPInstance pyobject.
    if (!Initialize())
        return nullptr;

// check validity of cast
    if (!CPPInstance_Check(pyobject))
        return nullptr;

// get held object (may be null)
    return ((CPPInstance*)pyobject)->GetObject();
}

//-----------------------------------------------------------------------------
PyObject* CPyCppyy::Instance_FromVoidPtr(
    void* addr, const std::string& classname, bool python_owns)
{
// Bind the addr to a python object of class defined by classname.
    if (!Initialize())
        return nullptr;

// perform cast (the call will check TClass and addr, and set python errors)
    PyObject* pyobject = BindCppObjectNoCast(addr, Cppyy::GetScope(classname), false);

// give ownership, for ref-counting, to the python side, if so requested
    if (python_owns && CPPInstance_Check(pyobject))
        ((CPPInstance*)pyobject)->PythonOwns();

    return pyobject;
}


//-----------------------------------------------------------------------------
bool CPyCppyy::Scope_Check(PyObject* pyobject)
{
// Test if the given object is of a CPPScope derived type.
    if (!Initialize())
        return false;

    return CPPScope_Check(pyobject);
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Scope_CheckExact(PyObject* pyobject)
{
// Test if the given object is of a CPPScope type.
    if (!Initialize())
        return false;

    return CPPScope_CheckExact(pyobject);
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Instance_Check(PyObject* pyobject)
{
// Test if the given pyobject is of CPPInstance derived type.
    if (!Initialize())
        return false;

// detailed walk through inheritance hierarchy
    return CPPInstance_Check(pyobject);
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Instance_CheckExact(PyObject* pyobject)
{
// Test if the given pyobject is of CPPInstance type.
    if (!Initialize())
        return false;

// direct pointer comparison of type member
    return CPPInstance_CheckExact(pyobject);
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Instance_IsLively(PyObject* pyobject)
{
// Test whether the given instance can safely return to C++, or whether
    if (!CPPInstance_Check(pyobject))
        return true;    // simply don't know

// the instance fails the lively test if it owns the C++ object while having a
// reference count of 1 (meaning: it could delete the C++ instance any moment)
    if (pyobject->ob_refcnt <= 1 && (((CPPInstance*)pyobject)->fFlags & CPPInstance::kIsOwner))
        return false;

    return true;
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Overload_Check(PyObject* pyobject)
{
// Test if the given pyobject is of CPPOverload derived type.
    if (!Initialize())
        return false;

// detailed walk through inheritance hierarchy
    return CPPOverload_Check(pyobject);
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Overload_CheckExact(PyObject* pyobject)
{
// Test if the given pyobject is of CPPOverload type.
    if (!Initialize())
        return false;

// direct pointer comparison of type member
    return CPPOverload_CheckExact(pyobject);
}


//- access to the python interpreter ----------------------------------------
bool CPyCppyy::Import(const std::string& mod_name)
{
// Import the named python module and create Cling equivalents for its classes.
    if (!Initialize())
        return false;

    PyObject* mod = PyImport_ImportModule(mod_name.c_str());
    if (!mod) {
        PyErr_Print();
        return false;
    }

// allow finding to prevent creation of a python proxy for the C++ proxy
    Py_INCREF(mod);
    PyModule_AddObject(gThisModule, mod_name.c_str(), mod);

// force creation of the module as a namespace
// TODO: the following is broken (and should live in Cppyy.cxx)
//   TClass::GetClass(mod_name, true);

    PyObject* dct = PyModule_GetDict(mod);

// create Cling classes for all new python classes
    PyObject* values = PyDict_Values(dct);
    for (int i = 0; i < PyList_GET_SIZE(values); ++i) {
        PyObject* value = PyList_GET_ITEM(values, i);
        Py_INCREF(value);

    // collect classes
        if (PyClass_Check(value) || PyObject_HasAttr(value, PyStrings::gBases)) {
        // get full class name (including module)
            PyObject* pyClName = PyObject_GetAttr(value, PyStrings::gName);
            if (PyErr_Occurred())
                PyErr_Clear();

        // build full, qualified name
            std::string fullname = mod_name;
            fullname += ".";
            fullname += CPyCppyy_PyText_AsString(pyClName);

      // force class creation (this will eventually call TPyClassGenerator)
      // TODO: the following is broken (and should live in Cppyy.cxx) to
      //         TClass::GetClass(fullname.c_str(), true);

            Py_XDECREF(pyClName);
        }

        Py_DECREF(value);
    }

    Py_DECREF(values);

// TODO: mod "leaks" here
    if (PyErr_Occurred())
        return false;
    return true;
}

//-----------------------------------------------------------------------------
void CPyCppyy::ExecScript(const std::string& name, const std::vector<std::string>& args)
{
// Execute a python stand-alone script, with argv CLI arguments.
//
// example of use:
//    CPyCppyy::ExecScript("test.py", {"1", "2", "3"});

    if (!Initialize())
        return;

// verify arguments
    if (name.empty()) {
        std::cerr << "Error: no file name specified." << std::endl;
        return;
    }

    FILE* fp = fopen(name.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: could not open file \"" << name << "\"." << std::endl;
        return;
    }

// store a copy of the old cli for restoration
    PyObject* oldargv = PySys_GetObject(const_cast<char*>("argv"));   // borrowed
    if (!oldargv)                                 // e.g. apache
        PyErr_Clear();
    else {
        PyObject* l = PyList_New(PyList_GET_SIZE(oldargv));
        for (int i = 0; i < PyList_GET_SIZE(oldargv); ++i) {
            PyObject* item = PyList_GET_ITEM(oldargv, i);
            Py_INCREF(item);
            PyList_SET_ITEM(l, i, item);          // steals ref
        }
        oldargv = l;
    }

// create and set (add progam name) the new command line
#if PY_VERSION_HEX < 0x03000000
    int argc = args.size() + 1;
    const char** argv = new const char*[argc];
    for (int i = 1; i < argc; ++i) argv[i] = args[i-1].c_str();
    argv[0] = Py_GetProgramName();
    PySys_SetArgv(argc, const_cast<char**>(argv));
    delete [] argv;
#else
// TODO: fix this to work like above ...
    (void)args;
#endif

// actual script execution
    PyObject* gbl = PyDict_Copy(gMainDict);
    PyObject* result =       // PyRun_FileEx closes fp (b/c of last argument "1")
        PyRun_FileEx(fp, const_cast<char*>(name.c_str()), Py_file_input, gbl, gbl, 1);
    if (!result)
        PyErr_Print();
    Py_XDECREF(result);
    Py_DECREF(gbl);

// restore original command line
    if (oldargv) {
        PySys_SetObject(const_cast<char*>("argv"), oldargv);
        Py_DECREF(oldargv);
    }
}

//-----------------------------------------------------------------------------
bool CPyCppyy::Exec(const std::string& cmd)
{
// Execute a python statement (e.g. "import noddy").
    if (!Initialize())
        return false;

// execute the command
    PyObject* result =
        PyRun_String(const_cast<char*>(cmd.c_str()), Py_file_input, gMainDict, gMainDict);

// test for error
    if (result) {
        Py_DECREF(result);
        return true;
    }

    PyErr_Print();
    return false;
}

//-----------------------------------------------------------------------------
const CPyCppyy::PyResult CPyCppyy::Eval(const std::string& expr)
{
// Evaluate a python expression.
//
// Caution: do not hold on to the return value: either store it in a builtin
// type (implicit casting will work), or in a pointer to a cppyy object (explicit
// casting to a void* is required).
    if (!Initialize())
        return PyResult();

// evaluate the expression
    PyObject* result =
        PyRun_String(const_cast<char*>(expr.c_str()), Py_eval_input, gMainDict, gMainDict);

// report errors as appropriate; return void
    if (!result) {
        PyErr_Print();
        return PyResult();
    }

// results that require no convserion
    if (result == Py_None || CPPInstance_Check(result) ||
            PyBytes_Check(result) ||
            PyFloat_Check(result) || PyLong_Check(result) || PyInt_Check(result))
        return PyResult(result);

// explicit conversion for python type required
    PyObject* pyclass = (PyObject*)Py_TYPE(result);

// retrieve class name and the module in which it resides
    PyObject* name = PyObject_GetAttr(pyclass, PyStrings::gName);
    PyObject* module = PyObject_GetAttr(pyclass, PyStrings::gModule);

 // concat name
    std::string qname =
        std::string(CPyCppyy_PyText_AsString(module)) + \
                    '.' + CPyCppyy_PyText_AsString(name);
    Py_DECREF(module);
    Py_DECREF(name);

// locate cppyy style class with this name
    // TODO: use Cppyy.cxx ...
    //TClass* klass = TClass::GetClass(qname.c_str());
    void* klass = nullptr;

// construct general cppyy python object that pretends to be of class 'klass'
    if (klass)
        return PyResult(result);

// no conversion, return null pointer object
    Py_DECREF(result);
    return PyResult();
}

//-----------------------------------------------------------------------------
void CPyCppyy::Prompt() {
// Enter an interactive python session (exit with ^D). State is preserved
// between successive calls.
    if (!Initialize())
        return;

// enter i/o interactive mode
    PyRun_InteractiveLoop(stdin, const_cast<char*>("\0"));
}
