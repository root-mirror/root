// @(#)root/tmva/sofie $Id$
// Author: Sanjiban Sengupta, 2021


#ifndef TMVA_SOFIE_PYINITIALIZE
#define TMVA_SOFIE_PYINITIALIZE

#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include "TMVA/RTensor.hxx"
#include "TString.h"


#include "TMVA/RModel.hxx"
#include "TMVA/SOFIE_common.hxx"


namespace TMVA{

// Declaring Global variables
PyObject *fModuleBuiltin = NULL;
PyObject *fEval = NULL;
PyObject *fOpen = NULL;


PyObject *fMain = NULL;
PyObject *fGlobalNS = NULL;
PyObject *fPyReturn = NULL;


///////////////////////////////////////////////////////////////////////////////
/// Check Python interpreter initialization status
///
/// \return Boolean whether interpreter is initialized
bool PyIsInitialized()
{
   if (!Py_IsInitialized()) return false;
   if (!fEval) return false;
   if (!fModuleBuiltin) return false;
   return true;
}


void PyInitialize()
{
   bool pyIsInitialized = PyIsInitialized();
   if (!pyIsInitialized) {
      Py_Initialize();
   }

   if (!pyIsInitialized) {
      _import_array();
   }

   fMain = PyImport_AddModule("__main__");
   if (!fMain) {
      std::cout<<"Python Error: Cannot import __main__\n";

   }

   fGlobalNS = PyModule_GetDict(fMain);
   if (!fGlobalNS) {
      std::cout<<"Python Error: Cannot init global namespace\n";
   }

   #if PY_MAJOR_VERSION < 3
   //preparing objects for eval
   PyObject *bName =  PyUnicode_FromString("__builtin__");
   // Import the file as a Python module.
   fModuleBuiltin = PyImport_Import(bName);
   if (!fModuleBuiltin) {
      std::cout<<"Python Error: Cannot import __builtin__\n";
   }
   #else
   //preparing objects for eval
   PyObject *bName =  PyUnicode_FromString("builtins");
   // Import the file as a Python module.
   fModuleBuiltin = PyImport_Import(bName);
   if (!fModuleBuiltin) {
      std::cout<<"Python Error: Cannot import __builtin__\n";
   }
   #endif

   PyObject *mDict = PyModule_GetDict(fModuleBuiltin);
   fEval = PyDict_GetItemString(mDict, "eval");
   fOpen = PyDict_GetItemString(mDict, "open");

   Py_DECREF(bName);
   Py_DECREF(mDict);
}




///////////////////////////////////////////////////////////////////////////////
// Finalize Python interpreter

void PyFinalize()
{
   Py_Finalize();
   if (fEval) Py_DECREF(fEval);
   if (fModuleBuiltin) Py_DECREF(fModuleBuiltin);
   if(fMain) Py_DECREF(fMain);//objects fGlobalNS and fLocalNS will be free here
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

void PyRunString(TString code, PyObject *fLocalNS, TString errorMessage="Failed to run python code", int start=Py_single_input) {
   fPyReturn = PyRun_String(code, start, fGlobalNS, fLocalNS);
   if (!fPyReturn) {
      std::cout<<"Failed to run python code: " << code <<"\n";
      std::cout<< "Python error message:\n";
      PyErr_Print();
      std::cout<<errorMessage;
   }
 }


namespace Experimental{
   namespace SOFIE{
      RTensor<float> getArray(PyObject* value){
         //Check and modify the function signature
         PyArrayObject* weightArray = (PyArrayObject*)value;
         std::vector<std::size_t>shape;
         std::vector<std::size_t>strides;

         //Preparing the shape vector
         for(npy_intp* j=PyArray_SHAPE(weightArray); j<PyArray_SHAPE(weightArray)+PyArray_NDIM(weightArray); ++j){
            shape.push_back((std::size_t)(*j));
            }

         //Preparing the strides vector
         for(npy_intp* k=PyArray_STRIDES(weightArray); k<PyArray_STRIDES(weightArray)+PyArray_NDIM(weightArray); ++k){
            strides.push_back((std::size_t)(*k));
            }

         //Declaring the RTensor object for storing weights values.
         RTensor<float>x((float*)PyArray_DATA(weightArray),shape,strides);
         return x;
         }
         }
         }

const char* PyString_AsString(PyObject* str){
   #if PY_MAJOR_VERSION < 3   // for Python2
      const char *stra_name = PyBytes_AsString(str);
      // need to add string delimiter for Python2
      TString sname = TString::Format("'%s'",stra_name);
      const char * name = sname.Data();
#else   // for Python3
      PyObject* repr = PyObject_Repr(str);
      PyObject* stra = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
      const char *name = PyBytes_AsString(stra);
#endif
return name;
}
}

#endif //TMVA_SOFIE_PYINITIALIZE
