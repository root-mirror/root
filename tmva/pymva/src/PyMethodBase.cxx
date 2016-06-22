// @(#)root/tmva/pymva $Id$
// Authors: Omar Zapata, Lorenzo Moneta, Sergei Gleyzer 2015

/**********************************************************************************
 * Project: TMVA - a Root-integrated toolkit for multivariate data analysis       *
 * Package: TMVA                                                                  *
 * Class  : PyMethodBase                                                          *
 *                                                                                *
 * Description:                                                                   *
 *      Virtual base class for all MVA method based on python                     *
 *                                                                                *
 **********************************************************************************/

#include<TMVA/PyMethodBase.h>
#include<TApplication.h>

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-function"
#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#if PY_VERSION_HEX >= 0x03000000
#define PyString_FromString PyBytes_FromString
#endif


using namespace TMVA;

ClassImp(PyMethodBase)

PyObject *PyMethodBase::fModuleBuiltin = NULL;
PyObject *PyMethodBase::fEval = NULL;
PyObject *PyMethodBase::fModulePickle = NULL;
PyObject *PyMethodBase::fPickleDumps = NULL;
PyObject *PyMethodBase::fPickleLoads = NULL;

//_______________________________________________________________________
PyMethodBase::PyMethodBase(const TString &jobName,
                           Types::EMVA methodType,
                           const TString &methodTitle,
                           DataSetInfo &dsi,
                           const TString &theOption ,
                           TDirectory *theBaseDir): MethodBase(jobName, methodType, methodTitle, dsi, theOption, theBaseDir),
   fClassifier(NULL)
{
   if (!PyIsInitialized()) {
      PyInitialize();
   }
}

//_______________________________________________________________________
PyMethodBase::PyMethodBase(Types::EMVA methodType,
                           DataSetInfo &dsi,
                           const TString &weightFile,
                           TDirectory *theBaseDir): MethodBase(methodType, dsi, weightFile, theBaseDir),
   fClassifier(NULL)
{
   if (!PyIsInitialized()) {
      PyInitialize();
   }
}

//_______________________________________________________________________
PyMethodBase::~PyMethodBase()
{
   if (PyIsInitialized()) {
//     PyFinalize();
   }
}

//_______________________________________________________________________
PyObject *PyMethodBase::Eval(TString code)
{
   PyObject *main = PyImport_AddModule("__main__");
   PyObject *global = PyModule_GetDict(main);
   PyObject *local = PyDict_New();

   PyObject *pycode = Py_BuildValue("(sOO)", code.Data(), global, local);
   PyObject *result = PyObject_CallObject(fEval, pycode);
//     Py_DECREF(main);
//     Py_DECREF(global);
//     Py_DECREF(local);
   Py_DECREF(pycode);
   return result;
}

//_______________________________________________________________________
void PyMethodBase::PyInitialize()
{
   TMVA::MsgLogger Log;
   if (!PyIsInitialized()) {
      Py_Initialize();
      _import_array();
      import_array();
   }
   //preparing objects for eval
   PyObject *bName = PyString_FromString("__builtin__");
   // Import the file as a Python module.
   fModuleBuiltin = PyImport_Import(bName);
   if (!fModuleBuiltin) {
      Log << kFATAL << "Can't import __builtin__" << Endl;
      Log << Endl;
   }
   PyObject *mDict = PyModule_GetDict(fModuleBuiltin);
   fEval = PyDict_GetItemString(mDict, "eval");

   Py_DECREF(bName);
   Py_DECREF(mDict);
   //preparing objects for pickle
   PyObject *pName = PyString_FromString("pickle");
   // Import the file as a Python module.
   fModulePickle = PyImport_Import(pName);
   if (!fModulePickle) {
      Log << kFATAL << "Can't import pickle" << Endl;
      Log << Endl;
   }
   PyObject *pDict = PyModule_GetDict(fModulePickle);
   fPickleDumps = PyDict_GetItemString(pDict, "dumps");
   fPickleLoads = PyDict_GetItemString(pDict, "loads");

   Py_DECREF(pName);
   Py_DECREF(pDict);


}

//_______________________________________________________________________
void PyMethodBase::PyFinalize()
{
   Py_Finalize();
   if (fEval) delete fEval;
   if (fModuleBuiltin) delete fModuleBuiltin;
   if (fPickleDumps) delete fPickleDumps;
   if (fPickleLoads) delete fPickleLoads;

}
void PyMethodBase::PySetProgramName(TString name)
{
   #if PY_VERSION_HEX >= 0x03000000
   const char* char_name = const_cast<char *>(name.Data());
   const size_t size = strlen(char_name)+1;
   wchar_t wchar_name[size];
   mbstowcs(wchar_name, char_name, size);
   Py_SetProgramName(wchar_name);
   #else
   Py_SetProgramName(const_cast<char *>(name.Data()));
   #endif
}
//_______________________________________________________________________
TString PyMethodBase::Py_GetProgramName()
{
   return Py_GetProgramName();
}
//_______________________________________________________________________
int  PyMethodBase::PyIsInitialized()
{
   if (!Py_IsInitialized()) return kFALSE;
   if (!fEval) return kFALSE;
   if (!fModuleBuiltin) return kFALSE;
   if (!fPickleDumps) return kFALSE;
   if (!fPickleLoads) return kFALSE;
   return kTRUE;
}


