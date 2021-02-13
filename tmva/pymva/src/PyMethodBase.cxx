// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015, Stefan Wunsch 2017

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : PyMethodBase                                                          *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method based on python                     *
 *                                                                                *
 **********************************************************************************/

#include <Python.h> // Needs to be included first to avoid redefinition of _POSIX_C_SOURCE
#include <TMVA/PyMethodBase.h>

#include "TMVA/DataSet.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/MsgLogger.h"
#include "TMVA/Results.h"
#include "TMVA/Timer.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <cwchar>

using namespace TMVA;

namespace TMVA {
namespace Internal {
class PyGILRAII {
   PyGILState_STATE m_GILState;

public:
   PyGILRAII() : m_GILState(PyGILState_Ensure()) {}
   ~PyGILRAII() { PyGILState_Release(m_GILState); }
};
} // namespace Internal
} // namespace TMVA

ClassImp(PyMethodBase);

// NOTE: Introduce here nothing that breaks if multiple instances
// of the same method share these objects, e.g., the local namespace.
PyObject *PyMethodBase::fModuleBuiltin = NULL;
PyObject *PyMethodBase::fEval = NULL;
PyObject *PyMethodBase::fOpen = NULL;

PyObject *PyMethodBase::fModulePickle = NULL;
PyObject *PyMethodBase::fPickleDumps = NULL;
PyObject *PyMethodBase::fPickleLoads = NULL;

PyObject *PyMethodBase::fMain = NULL;
PyObject *PyMethodBase::fGlobalNS = NULL;

///////////////////////////////////////////////////////////////////////////////

PyMethodBase::PyMethodBase(const TString &jobName, Types::EMVA methodType, const TString &methodTitle, DataSetInfo &dsi,
                           const TString &theOption)
   : MethodBase(jobName, methodType, methodTitle, dsi, theOption),
      fClassifier(NULL)
{
   if (!PyIsInitialized()) {
      PyInitialize();
   }

   // Set up private local namespace for each method instance
   fLocalNS = PyDict_New();
   if (!fLocalNS) {
      Log() << kFATAL << "Can't init local namespace" << Endl;
   }
}

///////////////////////////////////////////////////////////////////////////////

PyMethodBase::PyMethodBase(Types::EMVA methodType,
                           DataSetInfo &dsi,
                           const TString &weightFile): MethodBase(methodType, dsi, weightFile),
   fClassifier(NULL)
{
   if (!PyIsInitialized()) {
      PyInitialize();
   }

   // Set up private local namespace for each method instance
   fLocalNS = PyDict_New();
   if (!fLocalNS) {
      Log() << kFATAL << "Can't init local namespace" << Endl;
   }
}

///////////////////////////////////////////////////////////////////////////////

PyMethodBase::~PyMethodBase()
{
}

///////////////////////////////////////////////////////////////////////////////
/// Evaluate Python code
///
/// \param[in] code Python code as string
/// \return Python object from evaluation of code line
///
/// Take a Python code as input and evaluate it in the local namespace. Then,
/// return the result as Python object.

PyObject *PyMethodBase::Eval(TString code)
{
   if(!PyIsInitialized()) PyInitialize();
   PyObject *pycode = Py_BuildValue("(sOO)", code.Data(), fGlobalNS, fLocalNS);
   PyObject *result = PyObject_CallObject(fEval, pycode);
   Py_DECREF(pycode);
   return result;
}

///////////////////////////////////////////////////////////////////////////////
/// Initialize Python interpreter
///
/// NOTE: We introduce a shared global namespace `fGlobalNS`, but using
/// a private local namespace `fLocalNS`. This prohibits the interference
/// of instances of the same method with the same factory, e.g., by overriding
/// variables in the same local namespace.

void PyMethodBase::PyInitialize()
{
   TMVA::MsgLogger Log;

   bool pyIsInitialized = PyIsInitialized();
   if (!pyIsInitialized) {
      Py_Initialize();
   }

   TMVA::Internal::PyGILRAII raii;
   if (!pyIsInitialized) {
      _import_array();
   }

   fMain = PyImport_AddModule("__main__");
   if (!fMain) {
      Log << kFATAL << "Can't import __main__" << Endl;
      Log << Endl;
   }

   fGlobalNS = PyModule_GetDict(fMain);
   if (!fGlobalNS) {
      Log << kFATAL << "Can't init global namespace" << Endl;
      Log << Endl;
   }

   #if PY_MAJOR_VERSION < 3
   //preparing objects for eval
   PyObject *bName =  PyUnicode_FromString("__builtin__");
   // Import the file as a Python module.
   fModuleBuiltin = PyImport_Import(bName);
   if (!fModuleBuiltin) {
      Log << kFATAL << "Can't import __builtin__" << Endl;
      Log << Endl;
   }
   #else
   //preparing objects for eval
   PyObject *bName =  PyUnicode_FromString("builtins");
   // Import the file as a Python module.
   fModuleBuiltin = PyImport_Import(bName);
   if (!fModuleBuiltin) {
      Log << kFATAL << "Can't import builtins" << Endl;
      Log << Endl;
   }
   #endif

   PyObject *mDict = PyModule_GetDict(fModuleBuiltin);
   fEval = PyDict_GetItemString(mDict, "eval");
   fOpen = PyDict_GetItemString(mDict, "open");

   Py_DECREF(bName);
   Py_DECREF(mDict);
   //preparing objects for pickle
   PyObject *pName = PyUnicode_FromString("pickle");
   // Import the file as a Python module.
   fModulePickle = PyImport_Import(pName);
   if (!fModulePickle) {
      Log << kFATAL << "Can't import pickle" << Endl;
      Log << Endl;
   }
   PyObject *pDict = PyModule_GetDict(fModulePickle);
   fPickleDumps = PyDict_GetItemString(pDict, "dump");
   fPickleLoads = PyDict_GetItemString(pDict, "load");

   Py_DECREF(pName);
   Py_DECREF(pDict);
}

///////////////////////////////////////////////////////////////////////////////
// Finalize Python interpreter

void PyMethodBase::PyFinalize()
{
   Py_Finalize();
   if (fEval) Py_DECREF(fEval);
   if (fModuleBuiltin) Py_DECREF(fModuleBuiltin);
   if (fPickleDumps) Py_DECREF(fPickleDumps);
   if (fPickleLoads) Py_DECREF(fPickleLoads);
   if(fMain) Py_DECREF(fMain);//objects fGlobalNS and fLocalNS will be free here
}

///////////////////////////////////////////////////////////////////////////////
/// Set program name for Python interpeter
///
/// \param[in] name Program name

void PyMethodBase::PySetProgramName(TString name)
{
   #if PY_MAJOR_VERSION < 3
   Py_SetProgramName(const_cast<char*>(name.Data()));
   #else
   Py_SetProgramName((wchar_t *)name.Data());
   #endif
}


///////////////////////////////////////////////////////////////////////////////

size_t mystrlen(const char* s) { return strlen(s); }

///////////////////////////////////////////////////////////////////////////////

size_t mystrlen(const wchar_t* s) { return wcslen(s); }

///////////////////////////////////////////////////////////////////////////////
/// Get program name from Python interpreter
///
/// \return Program name

TString PyMethodBase::Py_GetProgramName()
{
   auto progName = ::Py_GetProgramName();
   return std::string(progName, progName + mystrlen(progName));
}

///////////////////////////////////////////////////////////////////////////////
/// Check Python interpreter initialization status
///
/// \return Boolean whether interpreter is initialized

int PyMethodBase::PyIsInitialized()
{
   if (!Py_IsInitialized()) return kFALSE;
   if (!fEval) return kFALSE;
   if (!fModuleBuiltin) return kFALSE;
   if (!fPickleDumps) return kFALSE;
   if (!fPickleLoads) return kFALSE;
   return kTRUE;
}

///////////////////////////////////////////////////////////////////////////////
/// Serialize Python object
///
/// \param[in] path Path where object is written to file
/// \param[in] obj Python object
///
/// The input Python object is serialized and written to a file. The Python
/// module `pickle` is used to do so.

void PyMethodBase::Serialize(TString path, PyObject *obj)
{
   if(!PyIsInitialized()) PyInitialize();

   PyObject *file_arg = Py_BuildValue("(ss)", path.Data(),"wb");
   PyObject *file = PyObject_CallObject(fOpen,file_arg);
   PyObject *model_arg = Py_BuildValue("(OO)", obj,file);
   PyObject *model_data = PyObject_CallObject(fPickleDumps , model_arg);

   Py_DECREF(file_arg);
   Py_DECREF(file);
   Py_DECREF(model_arg);
   Py_DECREF(model_data);
}

///////////////////////////////////////////////////////////////////////////////
/// Unserialize Python object
///
/// \param[in] path Path to serialized Python object
/// \param[in] obj Python object where the unserialized Python object is loaded
///  \return Error code

Int_t PyMethodBase::UnSerialize(TString path, PyObject **obj)
{
   // Load file
   PyObject *file_arg = Py_BuildValue("(ss)", path.Data(),"rb");
   PyObject *file = PyObject_CallObject(fOpen,file_arg);
   if(!file) return 1;

   // Load object from file using pickle
   PyObject *model_arg = Py_BuildValue("(O)", file);
   *obj = PyObject_CallObject(fPickleLoads , model_arg);
   if(!obj) return 2;

   Py_DECREF(file_arg);
   Py_DECREF(file);
   Py_DECREF(model_arg);

   return 0;
}

///////////////////////////////////////////////////////////////////////////////
/// Execute Python code from string
///
/// \param[in] code Python code as string
/// \param[in] errorMessage Error message which shall be shown if the execution fails
/// \param[in] start Start symbol
///
/// Helper function to run python code from string in local namespace with
/// error handling
/// `start` defines the start symbol defined in PyRun_String (Py_eval_input,
/// Py_single_input, Py_file_input)

void PyMethodBase::PyRunString(TString code, TString errorMessage, int start) {
   fPyReturn = PyRun_String(code, start, fGlobalNS, fLocalNS);
   if (!fPyReturn) {
      Log() << kWARNING << "Failed to run python code: " << code << Endl;
      Log() << kWARNING << "Python error message:" << Endl;
      PyErr_Print();
      Log() << kFATAL << errorMessage << Endl;
   }
}
